plot2up, version 0.0.1a, 23 August 2013
PURPOSE: Generates a PDF paired tree & pixel plot for each fasta file in cwd.
  USAGE: type 'make'
AUTHORS: Peter Hraber, Will Fischer, Bette Korber

I. DEPENDENCIES

Currently assuming Mac OS X environment, with interactive console.

A. Executables

  perl
        The perl interpreter (or a link to it) should be in /usr/bin.

  make or gmake
        Apple/OS X Developer command-line tools. 
	See http://www.macports.org/install.php

#  phyml
#	Not strictly necessary but helpful for all but the largest alignments.

	If you choose not to use phyml to make trees, you will have to
	customize the Makefile and/or ensure the Newick tree file
	contains the same prefix as the Fasta alignment file.  Be
	aware that make might over-write this file.

  R	http://cran.r-project.org/bin/macosx
	To ensure integrity of the dowloaded file, run 'md5' command on
	it and compare with published value.

#  Newick utilities
#  	Currently used only for nw_reroot but has great general utility.
#  	http://cegg.unige.ch/newick_utils

B. C library

  libgd
	This is the trickiest component to build and install from source.
	Using macports to obtain the pre-built 'gs' will install this for you.
	Install the library before installing the perl module below.
	See http://www.macports.org/install.php

C. PERL module/s

  bioperl
	See http://www.bioperl.org/wiki/Installing_BioPerl_on_Unix
  GD
	perl -MCPAN -e 'install "GD"'  (Easy, right?)

D. R packages

  png
    http://cran.r-project.org/package=png

  ape
    http://cran.r-project.org/package=ape

  To install these from within R using install.packages(), firewall rules 
  require that the proxy server be set as follows:
      export HTTP_PROXY=http://proxyout.lanl.gov:8080
      export HTTPS_PROXY=https://proxyout.lanl.gov:8080
      or within R, Sys.setenv("HTTP_PROXY=http://proxyout.lanl.gov:8080"), etc.

  Otherwise, manually download and install them from the command line, e.g:
    R CMD INSTALL ape_3.0-9.tgz
  You will have to download and install the other required packages:
    http://cran.r-project.org/package=nlme
    http://cran.r-project.org/package=gee
    http://cran.r-project.org/package=expm

II. Known Issues

- Files are mistreated when inserted into MS Office ppt or doc files.
- Goal is to replace png files generated by pixel with pure-R raster images.
- This version assumes nucleotide sequences.

END
