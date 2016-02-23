#/usr/bin/perl

=head 1 NAME

tt_separator.pl

=head 1 VERSION

v.0.1.0

=head 1 Usage

=cut

#===============================================================================
# VARIABLES AND OPTIONS
#===============================================================================


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use diagnostics;

my %options;
GetOptions (
    \%options    ,
    "help|?"     ,
    "directory=s",
);

if ($options{"help"}) {
    print "sh perl tt_separator folder_data/\n";
    exit(0);
} elsif (not $options{"directory"}) {
    die "You have to introduce the directory where the files are.\n";
}

#===============================================================================
# MAIN
#===============================================================================

# Dealing with imput

my $direct          = $options{"directory"};
my @files;

foreach my $file (glob("$direct/*")) {
    push (@files, $file);
}


# Predefined parameters

my $NUM_OF_SAMPLES = 288; # number of samples in the shortest cancer file
my $K = 4; # number of k-folds to test
my @CANCERS = ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca"); # cancer types

# Splitting the dataset

my $FILEHANDLES = filehandler(\@CANCERS,$K);


foreach my $filename (@files) {
    # Creating the set files
    my($file,$train_filename,$test_filename,$cancer) = fileOmatic($filename);

    my $first_line = <$file>;
    chomp($first_line);
    my @samples = split /\t/, $first_line;
    @samples    = @samples[0..$NUM_OF_SAMPLES - 1];
   
    # Shuffling samples to obtain a random set
    my ($new_to_prev,$shuff_samples) = shuffler(\@samples);

    my %header_done = ();
    while (<$file>){
     chomp;
     my ($gene, @expr) = split /\t/;

     # Discretization of expression values and dealing with NA and Inf
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
           my $step     = $k * 56;
           my %testing  = ();
           my %training = ();
           my $max_to_test = $step + 56;

           my $train_f = $FILEHANDLES->{"set$k"}->{$cancer}->{"train"};
           my $test_f  = $FILEHANDLES->{"set$k"}->{$cancer}->{"test"};

           # Creating the test subsample
           foreach my $i ($step..$max_to_test) {
               my $prev = $new_to_prev->{$i};
               $testing{$i} = $new_expr[$prev];
           }

           # Creating the test subsample
           foreach my $i (0..$#{ $shuff_samples }) {
               next if exists $testing{$i};
               my $prev = $new_to_prev->{$i};
               $training{$i} = $new_expr[$prev];
           }

           # Print testing and testing sample names
           # only once for CANCER
           if (not exists $header_done{$filename}->{$k}) {
               $header_done{$filename}->{$k} = 1;
               foreach my $i (0..$#{ $shuff_samples }) {
                   if (exists $testing{$i}) {
                       print $test_f "@{ $shuff_samples }[$i]\t";
                   } else {
                       print $train_f "@{ $shuff_samples }[$i]\t";
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


#--------------------------------------------------------------------------------
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

sub fill_hash{
    my $limit = shift;
    my %hash;
    foreach my $i (0..$#{ $limit }) {
        $hash{ @{ $limit }[$i] } = $i;
    }
    return (\%hash);
}

sub filehandler{
    my $CANCERS = shift;
    my $K = shift;
    my %FILEHANDLES = ();
    foreach my $i (0..$K) {
        my $dir ="set" . $i . "/";
        mkdir $dir;
        foreach my $c (@{ $CANCERS }) {
            my $test_f  = $dir . "test_$c.tbl";
            my $train_f = $dir . "train_$c.tbl";

            open $FILEHANDLES{"set$i"}->{$c}->{"test"}, ">", $test_f
                or die "Can't open $test_f : $!\n";

            open $FILEHANDLES{"set$i"}->{$c}->{"train"}, ">", $train_f
                or die "Can't open $train_f : $!\n";
        }
    }
    return (\%FILEHANDLES);
}

sub fileOmatic{
    my $filename = shift;
    my $file = master_key($filename);
    my $train_filename = $filename;
    my $test_filename  = $filename;
    my $cancer = $filename;
    $cancer =~ s/.+\/(.+)\.tbl/$1/g;
    $train_filename    =~ s/.+\/(.*?)$/\/train_$1/g;
    $test_filename     =~ s/.+\/(.*?)$/\/test_$1/g;
    return($file,$train_filename,$test_filename,$cancer);
}

sub shuffler{
    # Shuffle the samples in order to have a randomised sampling
    my $samples = shift;
    my @shuff_samples = @{ $samples };
    array_shuffle(\@shuff_samples);
    
    # Using index relationship between samples and shuffled samples
    # to get a table of correspondence: position of the samples:order of sampling
    my $prev_order = fill_hash($samples);
    my $new_order  = fill_hash(\@shuff_samples);

    my %new_to_prev = ();
    foreach my $sample (@{ $samples }) {
        my $prev_idx = $prev_order->{$sample};
        my $new_idx  = $new_order->{$sample};
        $new_to_prev{$new_idx} = $prev_idx;
    }
    return (\%new_to_prev,\@shuff_samples);
}