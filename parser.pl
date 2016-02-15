#!/usr/bin/perl
use strict;
use warnings;

my $work_file = shift (@ARGV);

my $file = master_key($work_file);

my $first_line = <$file>;
print "$first_line";
while (<$file>){
	chomp;
	my ($gene, @expr) = split /\t/;
    my @new_expr = map {
        if    ($_ > 2)  { "up" }
		elsif ($_ < -2) { "down" }
		elsif  ($_ > -2 and $_ < 2) { "nochange" }
		elsif  ($_ eq "Inf") { $_ }
		else { "NA" }
    } @expr;

	print "$gene\t", join("\t", @new_expr), "\n";

}








############################### FUNCTIONS #############################
sub master_key {+
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
