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
    "directory=s",
    "notvalid=s" ,
    "igthreshold=s",
);

if (not $options{"igthreshold"}) {
    # This magic number is the 3rd percentile of our computed IGs
    $options{"igthreshold"} = "0.241287";
}

if ($options{"help"}) {
    print "BLA";
    exit(0);
} elsif (not $options{"directory"}) {
    die "You have to introduce the directory where the files are.\n";
}

#===============================================================================
# MAIN
#===============================================================================

my @CANCERS = ("brca", "coad", "hnsc", "kric", "luad", "lusc", "prad", "thca");
my $directory = $options{"directory"};
my $not_valid   = read_not_valid($options{"notvalid"});
my %likelihoods = ();
my %priors = (
    "brca" => 1/8,
    "coad" => 1/8,
    "hnsc" => 1/8,
    "kric" => 1/8,
    "luad" => 1/8,
    "lusc" => 1/8,
    "prad" => 1/8,
    "thca" => 1/8
);
my @priors_arr = values %priors;

foreach my $train_files (glob("$directory/train_*")) {
    model_charger(\%likelihoods, $train_files, $not_valid);
}

my $rel_likelihoods = mutual_information(\%likelihoods, \@priors_arr);
undef %likelihoods;

my @test_files = glob("$directory/test_*");

predict_cancer($rel_likelihoods, \@test_files, \%priors);

#print Dumper($rel_likelihoods);


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
        next unless /[^\s\t]/;
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
    my $likelihoods          = shift;
    my $c_probabilities      = shift;
    my $relevant_likelihoods = ();
    my $class_entropy        = entropy($c_probabilities);

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

        if ($information_gain >= $options{"igthreshold"}) {
            foreach my $cancer (@CANCERS) {
                $relevant_likelihoods->{$cancer}->{$gene}->{"up"}
                    = log($likelihoods->{$cancer}->{$gene}->{"up"} / $likelihoods->{$cancer}->{$gene}->{"total"});
                $relevant_likelihoods->{$cancer}->{$gene}->{"down"}
                    = log($likelihoods->{$cancer}->{$gene}->{"down"} / $likelihoods->{$cancer}->{$gene}->{"total"});
                $relevant_likelihoods->{$cancer}->{$gene}->{"nochange"}
                    = log($likelihoods->{$cancer}->{$gene}->{"nochange"} / $likelihoods->{$cancer}->{$gene}->{"total"});
            }
        }

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

    return $relevant_likelihoods;
}

#--------------------------------------------------------------------------------
sub compute_entropy {
    my $expr        = shift;
    my $gene        = shift;
    my $likelihoods = shift;
    my $total_expr  = shift;

    my @probabilities = ();
    foreach my $cancer (@CANCERS) {
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
        foreach my $cancer (@CANCERS) {
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
sub predict_cancer {
    my $likelihoods  = shift;
    my $test_files   = shift;
    my $priors       = shift;
    my %cancer_probs = ();

    foreach my $t_file (@{ $test_files }) {
        my $cancer = $t_file;
        $cancer =~ s/.+\_(.+)\.tbl/$1/g;
        my $fh = master_key($t_file);
        my $first = <$fh>;
        chomp($first);
        my @samples = split /\t/, $first;

        SAMPLE:
        foreach my $s (@samples) {
            foreach my $c (@CANCERS) {
                next SAMPLE unless $s =~ /[^\s\t]/; # avoid whitespace
                $cancer_probs{$s}->{$c} = log($priors->{$c});
            }
        }

        while (<$fh>) {
            chomp;
            my ($gene, @expr) = split /\t/;
            next unless /[^\s\t]/; # avoid whitespace
            next unless exists $likelihoods->{"brca"}->{$gene};

            foreach my $i (0..$#expr) {
                foreach my $c (@CANCERS) {
                    next unless defined $samples[$i];
                    $cancer_probs{$samples[$i]}->{$c} += $likelihoods->{$c}->{$gene}->{$expr[$i]};
                }
            }
        }

        foreach my $samp (@samples) {
            my @sorted_cancers = sort {$cancer_probs{$samp}->{$b} <=> $cancer_probs{$samp}->{$a}} keys %{ $cancer_probs{$samp} };
            my $best_prediction = $sorted_cancers[0];
            my $prob = compute_prob($best_prediction, $samp, \%cancer_probs, $cancer);
            print "$samp\t$cancer\t$best_prediction\t$cancer_probs{$samp}->{$cancer}\t$prob\n";

        }

    }

}

sub compute_prob {
    my $best          = shift;
    my $sample        = shift;
    my $probabilities = shift;
    my $cancer        = shift;
    my $final_p       = 0;
    my $denominator   = 0;

    foreach my $cancer (@CANCERS) {
        $denominator += 10**($probabilities->{$sample}->{$cancer});
    }
    if ($denominator == 0) {
        $denominator = 1;
        # We did this to avoid floating point overflow that
        # produces an ilegal division by zero. All the probabilities will be wrong
        # for this sample.
        print STDERR "$sample probability will be wrong because of floating point problems: limit of Perl precision.\n";
    }

    $final_p = 10**($probabilities->{$sample}->{$best}) / $denominator;
    return $final_p;
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
