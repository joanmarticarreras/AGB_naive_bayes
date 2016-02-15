#/usr/bin/perl

=head 1 NAME

NotSoNaive.pl

=head 1 VERSION

v.0.1.0

=head 1 Usage

=cut


#===============================================================================
# VARIABLES AND OPTIONS
#===============================================================================

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;

my %options;
GetOptions (
    \%options    ,
    "help|?"     ,
    "train=s"    ,
    "test=s"     ,
    "notvalid=s" ,
);

if ($options{"help"}) {
    print "BLA";
    exit(0);
} elsif (not $options{"train"}) {
    die "You have to introduce training files separated by commas\n";
} elsif (not $options{"test"}) {
    die "You have to introduce one test file to make the predictions\n";
}


#===============================================================================
# MAIN
#===============================================================================

my @train_files = split /,/, $options{"train"};
my $not_valid   = read_not_valid($options{"notvalid"});
my %model       = ();

foreach my $t_file (@train_files) {
    model_charger(\%model, $t_file);
}

mutual_information(\%model, \@train_files);

print Dumper(\%model);

#===============================================================================
# FUNCTIONS
#===============================================================================
#--------------------------------------------------------------------------------
sub read_not_valid {
    my $file = shift;
    my $fh   = master_key($file);
    my %not_valid = ();

    while (<$fh>) {
        chomp;
        $not_valid{$_} = 1;
    }

    return \%not_valid;
}

#--------------------------------------------------------------------------------
sub model_charger {
    my $model     = shift;
    my $t_file    = shift;
    my $not_valid = shift;

    my $fh = master_key($t_file);

    my $first_line = <$fh>;
    while (<$fh>) {
        chomp;
        my ($gene, @expr) = split /\t/;
        next if exists $not_valid->{$gene};

        $model->{$t_file} = ()
            unless exists $model->{$t_file};
        $model->{$t_file}->{$gene} = ()
            unless exists $model->{$t_file}->{$gene};
        foreach my $ex_val (@expr) {
            $model->{$t_file}->{$gene}->{$ex_val}++;
            $model->{$t_file}->{$gene}->{'total'}++;
        }
        add_pseudocounts(\%model, $t_file, $gene);

    }
    return;
}

#--------------------------------------------------------------------------------
sub add_pseudocounts {
    my $model  = shift;
    my $cancer = shift;
    my $gene   = shift;

    $model->{$cancer}->{$gene}->{"up"}++;
    $model->{$cancer}->{$gene}->{"down"}++;
    $model->{$cancer}->{$gene}->{"nochange"}++;
    $model->{$cancer}->{$gene}->{"total"} += 3;

    return;
}

#--------------------------------------------------------------------------------
sub compute_probabilities {
    my $model  = shift;
    my $cancer = shift;
    my $gene   = shift;

    $model->{$cancer}->{$gene}->{"up"}       = $model->{$cancer}->{$gene}->{"up"} / $model->{$cancer}->{$gene}->{"total"};
    $model->{$cancer}->{$gene}->{"down"}     = $model->{$cancer}->{$gene}->{"down"} / $model->{$cancer}->{$gene}->{"total"};
    $model->{$cancer}->{$gene}->{"nochange"} = $model->{$cancer}->{$gene}->{"nochange"} / $model->{$cancer}->{$gene}->{"total"};
    return;
}

#--------------------------------------------------------------------------------
sub mutual_information {
    my $model       = shift;
    my $train_files = shift;
    my $class_entropy = entropy($train_files);
    my $condi_entropy = conditional_entropy($model);

}

#--------------------------------------------------------------------------------
sub entropy {
    my $classes = shift;
    my $total   = scalar(@{ $classes });
    my $result  = log2($total);

    return $result;
}

#--------------------------------------------------------------------------------
sub log2 {
    my $n = shift;
    return log($n) / log(2);
}

#--------------------------------------------------------------------------------
sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
