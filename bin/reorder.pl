#!/usr/local/bin/perl -w

#    NAME: sigtreed.pl
# PURPOSE: report phylogenetic distribution of signatures
#  INPUTS: tree file & fasta sequence file
# OUTPUTS: tree file annotated with signature counts
#  AUTHOR: Peter Hraber - phraber@lanl.gov
# VERSION: 1.3, 02/26/2008
#
# HISTORY: 
# 1.3 - 02/26/08 - Enabled option to reroot trees.
# 1.2 - 02/25/08 - Implemented empirical signature function.
# 1.1 - 02/25/08 - Fix to include all taxa below focal node, not just lineage.
# 1.0 - 02/22/08 - Initial prototype. Report depth rather than real signatures.

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::TreeIO;
#use Time::SoFar qw( runtime runinterval figuretimes );

my $DEBUG = 0;

&usage() if $#ARGV<1 or $#ARGV>2;

my $seqsFile = $ARGV[0];
my $treeFile = $ARGV[1];
my $outfile = ($ARGV[2]) ? $ARGV[2] : ">-"; # write to STDOUT by default

die "ERROR: cannot find $seqsFile ($!)" unless (-e $seqsFile);
die "ERROR: cannot find $treeFile ($!)" unless (-e $treeFile);

my $inSeqs = Bio::SeqIO->newFh(-file => $ARGV[0], -format => "Fasta");
my $inTree = Bio::TreeIO->new(-file => $ARGV[1], -format => "newick");
my $output = Bio::SeqIO->newFh(-format => "Fasta");

my $sequences = &readSequences($inSeqs); # returned value refers to an array
my $seqIDs = &getSequenceIDs($sequences);
my $tree;
my $root;
# we can accommodate multiple input trees and needn't assume uniqueness 
while ($tree = $inTree->next_tree) {

# to do: add option to call Will Fischer's ladderize script

    # find pair of taxa that maximize pairwise distance 
    my @leaves = $tree->get_leaf_nodes;
    my $taxa = &getLeafIDs($tree);

    ### sanity check: ensure that each leaf has a corresponding sequence
    my $nLeaves = &isSequencePresentForEachLeaf($taxa,$seqIDs);# if $DEBUG;
    my $treeSize = $tree->number_nodes if $DEBUG;
    print STDERR "working with $nLeaves taxa and $treeSize nodes\n" if $DEBUG;

#    my @nodes = $tree->get_nodes();
    my @nodes = $tree->get_leaf_nodes;
#    print STDERR "working with $#nodes = ".$#nodes."\n" if $DEBUG;

    my $i = 0;
    my %sequence = %$sequences;
    my %isInTree;
    my $seq_id;
    foreach $seq_id (@$seqIDs) {
	$isInTree{$seq_id} = 0;
    }

    foreach my $l (@nodes) {
	my $leaf_id = $l->id;

	foreach $seq_id (@$seqIDs) {

	    if ($leaf_id eq $seq_id) {
		my $s = Bio::Seq->new ( -seq => $sequence{$seq_id}, -id => $seq_id);
		print $output $s;
# remove $seq_id from @seqIDs
		$isInTree{$seq_id} = 1;
		last;
	    }
	}
    }
    # all sequences found in tree have been printed
# now output any remaining sequences
# this was done on 10/20/2008 to show low-quality excluded sequences

    foreach $seq_id (@$seqIDs) {
	if ($isInTree{$seq_id} == 0) {
	    my $t = Bio::Seq->new ( -seq => $sequence{$seq_id}, -id => $seq_id);
	    print $output $t;
	}
    }
}

sub myway { 
#    return $a->bootstrap() - $b->bootstrap();
# this should do the same thing:
    if ($a->bootstrap() > $b->bootstrap()) {
	return 1;
    } elsif ($a->bootstrap() > $b->bootstrap()) {
	return -1;
    }
    else {
	return 0;
    }
}

sub getExternalNodeCount() {
# return number of leaves below given node or 0 if a leaf

    my ($n) = @_;
    my $node = $$n;

    return 0 if ($node->is_Leaf);

    my $nChildren = 0;
    foreach my $child ($node->each_Descendent) {
	$nChildren++ if ($child->is_Leaf);
    }

    return $nChildren;
}

sub readSequences() {
    my ($input) = @_;
    my %sequences;

    while ( my $s = <$input> ) {
	my $seqID = $s->id;
	$sequences{$seqID} = $s->seq;
    }
    return \%sequences;
}

sub getSequenceIDs() {
    my ($sequences) = @_;

    my @seqIDs = keys %$sequences;
    return \@seqIDs;
}

sub getLeafIDs() {
    my ($tree) = @_;

    my @nodes = $tree->get_nodes();
    my @leafIDs;
    my $nLeaves=0;

    foreach my $myNode (@nodes) {

	my @children = $myNode->get_Descendents();
	$nLeaves++ if $#children == -1;

	push @leafIDs, $myNode->id() if $#children==-1 && defined($myNode->id);
    }
    my $nDefinedLeaves = $#leafIDs+1;
    # NB: 0-length branches will yield nodes that are apparently missing

    print STDERR "found $nLeaves leaves ($nDefinedLeaves with ids defined)\n"
	if $DEBUG;
    return \@leafIDs;
}

sub isSequencePresentForEachLeaf() {

    my ($leafIDs,$seqIDs) = @_;

    foreach my $Lid (@$leafIDs) {

	my $isSequencePresent = 0;

	foreach my $Sid (@$seqIDs) {

	    if ($Lid eq $Sid) {
		$isSequencePresent = 1;
		last;
	    }
	}

	die "$0 ERROR: could not find a sequence for $Lid" 
	    unless $isSequencePresent;
    }

    my @Lids = @$leafIDs;
    return $#Lids+1;
}

sub usage() {

    print "USAGE: $0 fasta_file tree_file [output_fasta_file]\n";
    print "WHERE:\n";
    print "  fasta_file contains a sequence for each leaf in tree_file\n";
    print "  tree_file is a newick-formatted tree file\n";
    print "  output_fasta_file (optional) will contain reordered alignment\n";
    print "  (if unspecified, this goes to STDOUT.)\n";
    exit;
}
