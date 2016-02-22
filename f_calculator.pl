#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;


my $file = shift @ARGV;

my @CANCERS = ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca");

open my $fh, "<", $file
    or die "Can't open $file :$!\n";

my %results = ();
# Initialize cancers:

foreach my $c (@CANCERS) {
    $results{$c}->{TP} = 0;
    $results{$c}->{FP} = 0;
    $results{$c}->{TN} = 0;
    $results{$c}->{FN} = 0;

}

while (<$fh>) {
    chomp;
    my ($s, $real, $prediction) = split /\t/;
    foreach my $cancer (@CANCERS) {
        if ($real eq $cancer) {
            # TP or FN
            if ($real eq $prediction) {
                $results{$cancer}->{TP}++;
            } else {
                $results{$cancer}->{FN}++;
            }
        } else {
            # FP or TN
            if ($prediction eq $cancer) {
                $results{$cancer}->{FP}++;
            } else {
                $results{$cancer}->{TN}++;
            }

        }
    }
}

foreach my $cancer (keys %results) {
    my $precision = $results{$cancer}->{TP} / ($results{$cancer}->{TP} + $results{$cancer}->{FP});
    my $recall    = $results{$cancer}->{TP} / ($results{$cancer}->{TP} + $results{$cancer}->{FN});
    my $f_measure = 2 / ((1/$precision) + (1/$recall));

    printf ("%s\t%.4f\t%.4f\t%.4f\n", $cancer, $precision, $recall, $f_measure);
}
