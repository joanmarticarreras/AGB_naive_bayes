#/usr/bin/perl

=head 1 NAME

NotSoNaive.pl

=head 1 VERSION

v.0.1.0

=head 1 Usage

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

my %options;
GetOptions (
    \%options   ,
    "help|?"    ,
    "train=s"   ,
    "test=s"    ,
);

if ($options{"help"}) {
    print "BLA";
    exit(0);
} elsif (not $options{"train"}) {
    die "You have to introduce training files separated by commas\n";
} elsif (not $options{"test"}) {
    die "You have to introduce one test file to make the predictions\n";
}

my @train_files = split /,/, $options{"train"};
my %model = ();

foreach my $t_file (@train_files) {
    model_charger(\%model, $t_file);
    print Dumper(\%model);
}

sub model_charger {
    my $model  = shift;
    my $t_file = shift;

    my $fh = master_key($t_file);

    my $first_line = <$fh>;
    while (<$fh>) {
        chomp;
        my ($gene, @expr) = split /\t/;

        $model->{$t_file} = ()
            unless exists $model->{$t_file};
        $model->{$t_file}->{$gene} = ()
            unless exists $model->{$t_file}->{$gene};
        foreach my $ex_val (@expr) {
            $model->{$t_file}->{$gene}->{$ex_val}++;
            $model->{$t_file}->{$gene}->{'total'}++;
        }

    }


}

sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
