#!/usr/bin/perl

=head1 NAME

f_calculator.pl

=head1 VERSION

v.0.1.0

=head1 DESCRIPTION

This program takes a list of results files from NotSoNaive.pl and computes Precision/Recall/F for each tumor type

=head1 USAGE

perl tt_separator.pl files

=head1 OPTIONS

=over 8

=item B<-h>, B<-help>

Shows this help.

=back

=head1 AUTHORS

Joan Marti i Carreras, Sergio Castillo Lara

=cut


use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;

pod2usage( -verbose => 0,
           -output  => \*STDOUT   ) unless @ARGV;

my %options = ();

GetOptions (
    \%options    ,
    "help|?"     ,
);

pod2usage( -verbose => 1,
           -output  => \*STDERR   ) if $options{help};

my @files = @ARGV;

print "CANCER\tMEASURE\tVALUE\n";
foreach my $file (@files) {

    my @CANCERS = ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca");

    open my $fh, "<", $file
        or die "Can't open $file :$!\n";
    print"$file\n";

    # Initialize cancers:
    my %results = ();
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
        my ($precision, $recall, $f_measure);
        if (($results{$cancer}->{TP} + $results{$cancer}->{FP}) == 0) {
            $precision = 0;
        } else {
            $precision = $results{$cancer}->{TP} / ($results{$cancer}->{TP} + $results{$cancer}->{FP});
        }
        if (($results{$cancer}->{TP} + $results{$cancer}->{FN}) == 0) {
            $recall = 0;
        } else {
            $recall = $results{$cancer}->{TP} / ($results{$cancer}->{TP} + $results{$cancer}->{FN});

        }

        if ($precision == 0 or $recall == 0) {
            print STDERR "$file\n";
            $f_measure = 0;
        } else {
            $f_measure = 2 / ((1/$precision) + (1/$recall));
        }
        printf ("%s\t%s\t%.4f\n", $cancer, "PRECISION", $precision);
        printf ("%s\t%s\t%.4f\n", $cancer, "RECALL", $recall);
        printf ("%s\t%s\t%.4f\n", $cancer, "F", $f_measure);
    }

        print("###################################################\n");

}
