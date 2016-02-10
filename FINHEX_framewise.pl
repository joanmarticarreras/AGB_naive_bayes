#!/usr/bin/perl
use strict;
use warnings;

############### DOCUMENTATION ###############

=pod

=head1

This is the document information concerning the script FINHEX_framewise.pl (FINHEX comes from FINd HEXamers).

=head2

This scripts allowes to find coding genes in several sequences using hexamer frequencies algorithm as a guide. This algorithm:

- Uses hexamer as a unit, due to the probability of a given hexamer is lower than for a given codon (a trimer). Lowering the probability to be found in a non-coding context, increases the chance to find a real coding gene.

- It suposes indepency between hexamers.

- Uses a log-likelihood ratio between a coding model and a non-coding model (intronic model, which it has been shown to be better than completly random model).

- "Training files" must be from the same species or at closely related to the query sequences, as codon bias changes between species.

=head3

- FINHEX_framewise.pl searches in every frame of the query sequences and uses the best score from the six frames to classify the sequences between coding and non-coding.

=head4

For further information on this topic, consult the README file.

In order to use it properly, introduce into the terminal in strict order:

$ perl FINHEX_framewise.pl <coding fasta file> <non-coding fasta file> <query fasta file>

All files must be in fasta format. In addition, the "training files" for coding and non-coding background, must provide the frame of lecture.

All fields at the fasta title must be separated by spaces. At natural position 10th (9th if starting from 0), must be placed the number of frame (0,1 or 2).


Tests performed with "training sets" of 2000 sequences (variable length) each and a quert test of 20 (variable length).

Statistical performance performance:

H0: Sequence is non-coding     HA: Sequences is coding

- By sequences

                        BLAST
                        CODING  NON-CODING
PREDICTED   CODING      9       0
            NON-CODING  5       6

- By probabilities

                        BLAST
                        CODING  NON-CODING
PREDICTED   CODING      0.45    0
            NON-CODING  0.25    0.3

STAT POWER = 0.45
TYPE I ERROR = 0.25
TYPE II ERROR = 0
CONFIDENCE = 0.3

CPU time and resources: 1 wallclock secs ( 0.37 usr +  0.00 sys =  0.37 CPU)

NOTE: BLAST algorithm was used on the query sequences in order to compare the efficiency of the script.

This script is an intellectual property which is owned by Joan Marti i Carreras. Even if my intention as its author is to make it open source, please acknowledge the authorship when being used.

Kind regards,

JMC

=cut

###############    MAIN    ###############
# Defining the 3 files with which the .pl will work
my ($cod_file, $ncod_file, $input_file) = @ARGV;

# Open files
my $cod_open = master_key($cod_file);
my $ncod_open = master_key($ncod_file);
my $input_open = master_key($input_file);

# Training sets by coding and non-coding background
my ($train_coding) = count_hex($cod_open);
my ($train_ncoding) = count_hex($ncod_open);

# Compute HEXAMER SCORE foreach sequences using coding and non-coding context to predict the outcome
$input_open = master_key($input_file);
local $/ = ">";
<$input_open>; # As using > as end of line (above command), by default, the first line is empty. Skip it.
my ($id, $seq, @seq, @aux, $inv_seq, $length, @scores, $p_scores, $n_scores);
while (<$input_open>){
    chomp $_;
    @aux = split (/\n/,$_);
    ($id, @seq) = @aux;
    $seq = join ("", @seq);
    $seq =~ tr/atgc/ATGC/;
    $inv_seq = reverse ($seq);
    $inv_seq =~ tr/ATGC/TACG/;
    $length = length($seq);
    if ($length > 5){ # Checks no sequences is shorter than 6 nt
        $p_scores = framer($seq, $length);
        $n_scores = framer($inv_seq, $length);
        @scores = (@{ $p_scores }, @{ $n_scores });
        @scores = sort {$b <=> $a} (@scores);
        if ($scores[0] >0){
            print $id."_CODING\n$seq\n";
        }else{
            print $id."_NON-CODING\n$seq\n";
        }
    }
}

############### FUNCTIONS ###############
sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}

sub count_hex{
    my $fh = shift;
    local $/ = ">";
    <$fh>; # As using > as end of line (above command), by default, the first line is empty. Skip it.
    my ($id, @seq,$seq, $length, %count_hex, $counter);
    while (<$fh>){
        chomp $_;
        my @aux = split (/\n/,$_);
        ($id, @seq) = @aux;
        $seq = join ("", @seq);
        $seq =~ tr/atgc/ATGC/;
        $length = length($seq);
        if ($length > 5){ # Checks no sequences is shorter than 6 nt
            for (my $i = 0; $i < $length; $i+=3) {
                my $sub_seq = substr($seq,$i,6);
                if (length ($sub_seq) > 5){ # Check no mer < 6-mer is taken into account
                $count_hex{$sub_seq}++; 
                $counter++;     
                }
            }       
        }
    }
    foreach my $keys1 (keys %count_hex){
        $count_hex{$keys1} = $count_hex{$keys1}/$counter;
    }
    return (\%count_hex);
}

sub framer{
    my $seq = shift;
    my $length = shift;
    my ($scr, $c, @scrs);
    for (my $frame = 0; $frame < 3; $frame++){
        for (my $i = $frame; $i < $length; $i+=3) {
          	my $sub_seq = substr($seq,$i,6);
            if (length($sub_seq) > 5){
            	$c++;
                $scr += score ($train_coding->{$sub_seq}, $train_ncoding->{$sub_seq});
            }
        }   
        $scr = $scr/$c;
        push (@scrs, $scr);
        $c = 0;
    }
    return(\@scrs);
}
    
sub score{
    my $cod = shift;
    my $ncod = shift;
    my $n;
    if (not defined $cod){ # Needs a defined, as not all hexamers found in the query sequences can be found in each frame of the training
        $n = -4; # log of 0 in base 2 is -inf, then is better to give an alternative negative return is recomended
    }
    elsif(not defined $ncod){
        $n = 4; # 1/0 does not exists, but implies that the hexamer is only from coding origin, an alternative positive return is recommended
    }
    else{
        my $ratio = $cod/$ncod;
        $n = log ($ratio)/log (2); 
    }
    return $n;  
}