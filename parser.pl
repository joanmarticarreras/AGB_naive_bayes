#!/usr/bin/perl
use strict;
use warnings;

my @files = @ARGV;

foreach my $filename (@files) {
    my $file = master_key($filename);
	my $outfilename = $filename;
	$outfilename    =~ s/\/(.*?)/\/f_$1/g;
	open my $out_fh, ">", "$outfilename"
		or die "Can't write to $outfilename : $!\n";

   	my $first_line = <$file>;
   	print $out_fh "$first_line";
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
	   print $out_fh "$gene\t", join("\t", @new_expr), "\n";
	}

}









############################### FUNCTIONS #############################
sub master_key {+
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
