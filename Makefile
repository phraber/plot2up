## 2up - version 0.0.1alpha
# to do: 
#  automagically scale leaf label size
#  add optional support for amino-acid sequences
#  add option for colored leaves/branches
#  add axis labels to pixel plot

# Edit this string to identify sequence for use at top of tree & highlighter plots:
TREEROOT = Ref-B.FR.83.HXB2_LAI_IIIB_BRU.K03455

# Edit this to change PNG file size and resolution
S_FACTOR = 20

BINDIR = /Users/phraber/plot2up/v1/bin
SOURCES = $(wildcard *.fasta)
TARGETS = $(SOURCES:.fasta=.pdf)

default: $(TARGETS)

all: default

# ref.fasta is 399-nt from gp41 ectodomain of HIV-1 subtype reference alignment
test: ref.pdf

clean:
	rm *.trim *.trim.table *.rtree *.phyml *.rfasta *.png *.tre subnames.* *.phyml_phyml_tree.txt *.bak *~ *.pdf *.out

%.pdf: %.png
	cat plot2up.r | R --vanilla --args $* $(S_FACTOR) \
	`head -1 $*.phyml | perl -ne '/(\d+)\ (\d+)/ and print $$2."\n";'` > $*.out

%.png: %.rfasta
	$(BINDIR)/pixel --SCALE=$(S_FACTOR) $*.rfasta
	mv $*.rfasta.png $*.png

%.rfasta: %.rtree
	$(BINDIR)/reorder.pl $*.fasta $*.rtree | $(BINDIR)/dotifyFasta.pl > $@

%.rtree: %.tre
	nw_reroot $< $(TREEROOT) | $(BINDIR)/ladderize.pl > $*.rtree

%.tre: %.phyml
	phyml -i $< -d nt -q -m GTR -f m -v e -c 4 -a e -t e -o tlr -b 0
	chmod 755 subnames.$*.fasta
	./subnames.$*.fasta $<_phyml_tree.txt
	rm $<_phyml_tree.txt.bak
	rm $<_phyml_stats.txt
	rm subnames.$*.fasta
	rm $*.trim
	rm $*.trim.table
	mv $<_phyml_tree.txt $@

%.phyml: %.fasta
	$(BINDIR)/trim_fasta_names -f $< > $*.trim
	length=`($(BINDIR)/padfasta $*.trim | cf - Table | perl -pe 's/_*\t/ /' > $*.trim.table) 2>&1 | perl -ne '/(\d+) positions/ and print $$1."\n";'`; \
	ntax=`wc -l <  $*.trim.table`; \
	( echo $$ntax $$length; cat $*.trim.table) > $@

.PRECIOUS: %.fasta %.phyml %.pdf %.png %.r %.rfasta %.rtree %.tre
