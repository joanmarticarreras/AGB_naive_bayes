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
my %likelihoods = ();
my @cancer_probs = (
    1/8,
    1/8,
    1/8,
    1/8,
    1/8,
    1/8,
    1/8,
    1/8
);

foreach my $t_file (@train_files) {
    model_charger(\%likelihoods, $t_file, $not_valid);
}

mutual_information(\%likelihoods, \@cancer_probs);
print Dumper(\%likelihoods);


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
    my $likelihoods = shift;
    my $t_file      = shift;
    my $not_valid   = shift;
    my $cancer      = $t_file;
    $cancer =~ s/.+\_(.+)\.tbl/$1/g;
    my $fh = master_key($t_file);

    my $first_line = <$fh>;
    while (<$fh>) {
        chomp;
        next unless /[^\s]/;
        my ($gene, @expr) = split /\t/;
        next if exists $not_valid->{$gene};

        $likelihoods->{$cancer} = ()
            unless exists $likelihoods->{$cancer};
        $likelihoods->{$cancer}->{$gene} = ()
            unless exists $likelihoods->{$cancer}->{$gene};
        foreach my $ex_val (@expr) {
            $likelihoods->{$cancer}->{$gene}->{$ex_val}++;
            $likelihoods->{$cancer}->{$gene}->{'total'}++;
        }
        add_pseudocounts($likelihoods, $cancer, $gene);

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
    my $likelihoods     = shift;
    my $c_probabilities = shift;
    my $class_entropy   = entropy($c_probabilities);

    foreach my $gene (keys %{ $likelihoods->{"brca"} }) {
        # we use brca to get the genes, but we could use whatever cancer
        # we want.
        my ($condi_entropy, $total_expr) = conditional_entropy($gene, $likelihoods);
        my $up_entropy       = compute_entropy("up", $gene, $likelihoods, $total_expr);
        my $down_entropy     = compute_entropy("down", $gene, $likelihoods, $total_expr);
        my $nochange_entropy = compute_entropy("nochange", $gene, $likelihoods, $total_expr);

        my $information_gain = $class_entropy - ($condi_entropy->{"up"} * $up_entropy
                                                + $condi_entropy->{"down"} * $down_entropy
                                                + $condi_entropy->{"nochange"} * $nochange_entropy);

        # print "\n$gene:\n";
        # print "CLASS: $class_entropy\n";
        # print "Condi_up: $condi_entropy->{up}\n";
        # print "Condi_down: $condi_entropy->{down}\n";
        # print "Condi_nochange: $condi_entropy->{nochange}\n";
        # print "UP_entropu: $up_entropy\n";
        # print "DOWN_entropu: $down_entropy\n";
        # print "NOCHANGE_entropu: $nochange_entropy\n";
        # print "IG: $information_gain\n";
        # print "\n----\n";
    }

    return;
}

#--------------------------------------------------------------------------------
sub compute_entropy {
    my $expr        = shift;
    my $gene        = shift;
    my $likelihoods = shift;
    my $total_expr  = shift;

    my @probabilities = ();
    foreach my $cancer ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca") {
        my $prob = $likelihoods->{"$cancer"}->{"$gene"}->{$expr} / $total_expr->{$expr};
        push @probabilities, $prob;
    }
    my $entropy = entropy(\@probabilities);

    return($entropy)
}


#--------------------------------------------------------------------------------
sub conditional_entropy {
    my $gene        = shift;
    my $likelihoods = shift;

    my $total_people = 0;
    my %result       = ();
    foreach my $expr ("up", "down", "nochange") {
        my $total_expr = 0;
        foreach my $cancer ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca") {
            $total_expr   += $likelihoods->{"$cancer"}->{"$gene"}->{"$expr"};
            $total_people += $likelihoods->{"$cancer"}->{"$gene"}->{"$expr"};
        }
        $result{$expr} = $total_expr;

    }
    my %total_expr = ();
    $total_expr{"up"}       = $result{"up"};
    $total_expr{"down"}     = $result{"down"};
    $total_expr{"nochange"} = $result{"nochange"};

    $result{"up"}       = $result{"up"} /$total_people;
    $result{"down"}     = $result{"down"} /$total_people;
    $result{"nochange"} = $result{"nochange"} /$total_people;

    # This returns something like:
    # {up}       => number
    # {down}     => number
    # {nochange} => number
    return \%result, \%total_expr;
}

#--------------------------------------------------------------------------------
sub entropy {
    my $probabilities = shift;
    my $result = 0;
    foreach my $prob (@{$probabilities}) {
        $result -= $prob * log2($prob);
    }

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
