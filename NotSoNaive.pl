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
    \%options   ,
    "help|?"    ,
    "train=s"   ,
    "test=s"    ,
    "notvalid"  ,
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
    print Dumper(\%model);
}



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
        
        my %check_pseudo  = ("up"      => 0, "down"  => 0,
                            "nochange" => 0, "total" => 0);
        $model->{$t_file} = ()
            unless exists $model->{$t_file};
        $model->{$t_file}->{$gene} = ()
            unless exists $model->{$t_file}->{$gene};
        foreach my $ex_val (@expr) {
            $model->{$t_file}->{$gene}->{$ex_val}++;
            $check_pseudo{$ex_val} = 1;
            $model->{$t_file}->{$gene}->{'total'}++;
        }
        if (keys  %check_pseudo != 4) {
            check_pseudocounts(\%model);
        }

    }
    return;
}

#--------------------------------------------------------------------------------
sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
