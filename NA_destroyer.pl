use warnings;
use strict;
use Data::Dumper;

my @files = @ARGV;

foreach my $filename (@files) {
    my $fh = master_key($filename);
	my $outfilename = $filename;
	$outfilename    =~ s/\/(.*?)/\/${1}_NA/g;
	open my $out_fh, ">", "$outfilename"
		or die "Can't write to $outfilename : $!\n";

    # READ ALL THE SAMPLES
    my %data = ();
    my $first_line = <$fh>;
    my @samples = split /\t/, $first_line;
    foreach my $i (0..$#samples) {
        $data{$i}  = ();
    }

    # ADD DATA FOR EACH SAMPLE
    while (<$fh>) {
        chomp;
        my ($gene, @expr) = split /\t/;
        for (my $i = 0; $i < $#expr; $i++) {
            if (exists $data{$i}) {
                if ($expr[$i] eq "NA" or $expr[$i] eq "Inf") {
                    # OJO CON LOS INFINITOS
                    delete $data{$i};
                    next;
                } else {
                    $data{$i}->{$gene} = $expr[$i];
                }
            }

        }

    }
    print Dumper(%data);

}

sub master_key {
    my $file = shift or die "Usage: $0 FILE\n";
    open my $file_fh, '<', $file or die "Could not open '$file' $!";
    return $file_fh;
}
