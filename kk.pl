#/usr/bin/perl

use warnings;
use strict;

my $file = shift @ARGV;

my $value = shift @ARGV;

my $fh = master_key($file);

while (<$fh>){
	my @cols = split /\t/;
	if ($cols[3] > $value){
		print $cols[0]."\t".$cols[1]."\t".$cols[2]."\t".$cols[3]."\t".$cols[4];
	}
}























sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
