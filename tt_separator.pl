#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my @files          = @ARGV;
my $NUM_OF_SAMPLES = 288;
my $K = 5;
my @CANCERS = ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca");

my %FILEHANDLES = ();
foreach my $i (0..$K) {
    my $dir ="set" . $i . "/";
    mkdir $dir;
    foreach my $c (@CANCERS) {
        my $test_f  = $dir . "test_$c.tbl";
        my $train_f = $dir . "train_$c.tbl";

        open $FILEHANDLES{"set$i"}->{$c}->{"test"}, ">", $test_f
            or die "Can't open $test_f : $!\n";

        open $FILEHANDLES{"set$i"}->{$c}->{"train"}, ">", $train_f
            or die "Can't open $train_f : $!\n";

    }

}

foreach my $filename (@files) {
    my $file = master_key($filename);
	my $train_filename = $filename;
    my $test_filename  = $filename;
    my $cancer = $filename;
    $cancer =~ s/.+\/(.+)\.tbl/$1/g;
	$train_filename    =~ s/.+\/(.*?)$/\/train_$1/g;
    $test_filename     =~ s/.+\/(.*?)$/\/test_$1/g;

   	my $first_line = <$file>;
    chomp($first_line);
    my @samples = split /\t/, $first_line;
    @samples    = @samples[0..$NUM_OF_SAMPLES - 1];
    my @shuff_samples = @samples;
    array_shuffle(\@shuff_samples);
    my %prev_order = ();
    my %new_order  = ();

    foreach my $i (0..$#samples) {
        $prev_order{$samples[$i]} = $i;
    }

    foreach my $i (0..$#shuff_samples) {
        $new_order{$shuff_samples[$i]} = $i;
    }

    my %new_to_prev = ();
    foreach my $samples (@samples) {
        my $prev_idx = $prev_order{$samples};
        my $new_idx  = $new_order{$samples};
        $new_to_prev{$new_idx} = $prev_idx ;
    }



    my %header_done = ();
    while (<$file>){
	   chomp;
	   my ($gene, @expr) = split /\t/;

	   my @new_expr = map {
		   if     ($_ eq "NA")  { $_ }
		   elsif  ($_ eq "Inf") { $_ }
		   elsif  ($_ > 2)  { "up" }
		   elsif  ($_ < -2) { "down" }
		   elsif  ($_ > -2 and $_ < 2) { "nochange" }
		   else { print STDERR "Error...\n"; "ERROR" }
	   } @expr;

       @new_expr = @new_expr[0..$NUM_OF_SAMPLES -1 ];


       foreach my $k (0..$K) {
           my $step     = $k * 48;
           my %testing  = ();
           my %training = ();
           my $initial = 0 + $step;
           my $max_to_test = 95;

           # THESE ARE FOR THE TESTING SET
           my $train_f = $FILEHANDLES{"set$k"}->{$cancer}->{"train"};
           my $test_f  = $FILEHANDLES{"set$k"}->{$cancer}->{"test"};
           foreach my $i ($initial..$#shuff_samples) {
               if (scalar(keys %testing) >= $max_to_test) {
                   last;
               }
               my $prev = $new_to_prev{$i};
               $testing{$i} = $new_expr[$prev];
           }

           # THESE ARE FOR THE TRAINING SET
           foreach my $i (0..$#shuff_samples) {
               next if exists $testing{$i};
               my $prev = $new_to_prev{$i};
               $training{$i} = $new_expr[$prev];
           }

           # Print testing and testing sample names
           # only once for CANCER
           if (not exists $header_done{$filename}->{$k}) {
               $header_done{$filename}->{$k} = 1;
               foreach my $i (0..$#shuff_samples) {
                   if (exists $testing{$i}) {
                       print $test_f "$shuff_samples[$i]\t";
                   } else {
                       print $train_f "$shuff_samples[$i]\t";
                   }
               }
           }

           print $test_f "\n";
           print $train_f "\n";

           print $test_f "$gene\t";
           print $train_f "$gene\t";

           foreach my $i (sort keys %testing) {
               print $test_f "$testing{$i}\t";
           }

           foreach my $i (sort keys %training) {
               print $train_f "$training{$i}\t";
           }

       }

	}

}


############################### FUNCTIONS #############################
sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}

sub array_shuffle { # F-Y shuffle
    my $array = shift;
    my $i     = @$array;
    while ( --$i ) {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
    return;
}
