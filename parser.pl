#!/usr/bin/perl
use strict;
use warnings;

my @files          = @ARGV;
my $NUM_OF_SAMPLES = 288;

mkdir "set1";
mkdir "set2";

foreach my $filename (@files) {
    my $file = master_key($filename);
	my $train_filename = $filename;
    my $test_filename  = $filename;
	$train_filename    =~ s/.+\/(.*?)$/\/train_$1/g;
    $test_filename     =~ s/.+\/(.*?)$/\/test_$1/g;

	open my $train_fh1, ">", "set1/" . "$train_filename"
		or die "Can't write to $train_filename : $!\n";

    open my $test_fh1, ">", "set1/" . "$test_filename"
    	or die "Can't write to $test_filename : $!\n";

    open my $train_fh2, ">", "set2/" . "$train_filename"
    	or die "Can't write to $train_filename : $!\n";

    open my $test_fh2, ">", "set2/" . "$test_filename"
    	or die "Can't write to $test_filename : $!\n";

   	my $first_line = <$file>;
    my @samples = split /\t/, $first_line;
    @samples = @samples[0..$NUM_OF_SAMPLES - 1];
    foreach my $i (0..$#samples) {
        if ($i % 2 == 0) {
            print $train_fh1 "$samples[$i]\t";
            print $test_fh2  "$samples[$i]\t";
        } else {
            print $test_fh1  "$samples[$i]\t";
            print $train_fh2 "$samples[$i]\t";
        }
    }
    print $train_fh1 "\n";
    print $test_fh1  "\n";
    print $train_fh2 "\n";
    print $test_fh2  "\n";

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

       print $train_fh1 "$gene\t";
       print $test_fh1  "$gene\t";
       print $train_fh2 "$gene\t";
       print $test_fh2  "$gene\t";

       foreach my $i (0..$#new_expr) {
           if ($i % 2 == 0) {
               print $train_fh1 "$new_expr[$i]\t";
               print $test_fh2 "$new_expr[$i]\t";
           } else {
               print $test_fh1 "$new_expr[$i]\t";
               print $train_fh2 "$new_expr[$i]\t";
           }
       }
       print $train_fh1 "\n";
       print $test_fh1  "\n";
       print $train_fh2 "\n";
       print $test_fh2  "\n";
	}

}


############################### FUNCTIONS #############################
sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
