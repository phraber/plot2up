#!/usr/local/bin/perl -w

# just translate a fasta file of nt sequences to aa sequences

# derived from codon-consensus.pl
# this routine attempts to report the most frequent (within-patient) codon
# so that the consensus nt sequence matches the consensus translation

#   INPUTS: a fasta-formatted file

#  OUTPUTS: a fasta-formatted file with one more sequence than the
#  input, the consensus sequence, at the beginning of the file

# USAGE: translate.pl [infile] [outfile]

# DEFAULTS: without arguments reads from and writes to standard
# output (like a unix filter), uses "consensus" to identify
# consensus sequence, and +1 frame for translation

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $DEBUG = 0;

my $in = Bio::SeqIO->newFh(-fh => \*ARGV, '-format' => 'Fasta');

my $unknown_aa = "X";
my $stop_aa = "Z";
my $myframe = 0;
my $out = Bio::SeqIO->newFh('-format' => 'Fasta');
my $nt;
while ( eval { $nt = <$in> } ) {

    my $aa = $nt->translate(-terminator => $stop_aa,
			    -unknown => $unknown_aa,
			    -codontable_id => 1,
			    -frame => $myframe);

    print $out $aa;
}

