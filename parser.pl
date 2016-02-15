#!/usr/bin/perl
use strict;
use warnings;

my $work_file = shift (@ARGV);

my $brca = master_key($work_file);


print (<$brca>);
while (< $brca >){
	chomp;
	my @line = split /\t/,$brca;
	
}







############################### FUNCTIONS #############################
sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}