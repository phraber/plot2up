#!/opt/local/bin/perl -w

use vars qw ( $perllibdir );

BEGIN {
  $perllibdir="$ENV{BPDIR}" if (defined $ENV{BPDIR});
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;

use strict;

use Bio::Seq;
use Bio::SeqIO;

my $in = Bio::SeqIO->newFh(-fh => \*ARGV, '-format' => 'Fasta');
#$in = Bio::SeqIO->newFh('-format' => 'Fasta') if ($#ARGV<0);

my $out = Bio::SeqIO->newFh('-format' => 'Fasta');

my $currseq;
my $firstseq;

while ( eval { $currseq = <$in> } ) {
    my $temp = uc($currseq->seq());
    $currseq->seq($temp);
#    $temp =~ s/-//g;
#    $out->write_seq($currseq) if (length($temp) > $minLength);
#    print $currseq->id(), "\t";
#    print $currseq->length(), "\n";

    if (! defined $firstseq) {
	print $out $currseq;
	$firstseq = $currseq;
    }
    else {

	my $newseq = "";
	for (my $site = 0; $site < length($firstseq->seq); $site++) {
	    if (substr($currseq->seq, $site, 1) ne substr($firstseq->seq, $site, 1)) {
		$newseq .= substr($currseq->seq, $site, 1);
	    }
	    else {
		$newseq .= ".";
	    }
	}
	my $s = Bio::Seq->new ( -seq => $newseq, -id => $currseq->id);
	print $out $s;
    }
} # end of while loop
