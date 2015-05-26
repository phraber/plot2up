#!/usr/local/bin/perl -w

# "ladderize" Newick format trees

# $_ = '(1,(((6,5),(3,4)),2));';
# to    ((((3,4),(5,6)),2),1);

# Newick trees end with a semicolon at the end of the line
$/ = ";\n";

OPTIONS: while (($_ = $ARGV[0]) =~ /^-/) {
	shift;
	/h/ and $no_header++,   next;  # remove fastDNAml header
	/n/ and $no_newlines++, next;  # don't interpolate newlines in trees
	/b/ and $no_brlens++,   next;  # remove branch-lengths from trees
	/s|t/ and $no_brlens++,
	          $no_header++,
	          $no_newlines++, next;  # fully Strip trees to Topology only
	/q/ and $add_quotes++,  next;  # add quotes around taxon names
	/r/ and $reverse++,  next;  # reverse order

	warn "option $_ not (currently) supported (skipping)\n";
	next OPTIONS;
}


while (<>) {

	# strip branch-lengths (if so instructed)
	# (do it before sort)
	if ($no_brlens) {
		#s/: *-*\.{0,1}\d+(\.\d+(e-\d+){0,1}){0,1}//go;
		#s/:[-\.]*\d+(\.\d+){0,1}(e-\d+){0,1}//go;
		strip_brlens(\$_);
	}
	
        # trim out all newlines
	s/\n//go;
	s/\A\[.*\] //o if $no_header;
	
	# find locations of opening and closing parentheses
	@opens  = &offsets($_,'(');
	@closes = &offsets($_,')');
	
	# find the matching open paren for each closing paren
	# and sort contents (separated by commas)
	while ($close = shift(@closes)) {
		$i = $#opens;
		
		# count down in opening (s until finding one whose offset
		# is less than the current closing )
		while ($opens[$i] > $close) {
			$i--;
		}
		# match is stuff inside parens
		my($start,$end) = ($opens[$i]+1,$close-$opens[$i]-1);
		$match = substr($_,$start,$end);

		# remove matched open paren from list
		splice(@opens,$i,1);

		# sort paren contents
		# (and change commas to #s to hide them from
		#  subsequent sorts)
		if (! $reverse) {
			$newmatch = join('#',sort by_clade (split(/,/,$match)));
		}
		else {
			$newmatch = join('#',sort by_clade_reverse (split(/,/,$match)));
		}

		# replace unsorted with sorted
		substr($_,$start,$end) = $newmatch;

	}
	# turn #s back into commas
	tr/#/,/;
	
	# add quotes around taxon names (if so instructed)
	if ($add_quotes) {
		s/([^\'\(\),:;]+)/'$1'/g;
		s/\'\'/'/go;
		s/:' *([\d\.]+)'/:$1/g;
	}

	
	# print
	unless ($no_newlines) {
		# add newline after fastDNAml comments
		s/\A(\[.*\] )/$1\n/o;
		print &wrap_commas($_);
	}
	else {
		print "$_\n";
	}
}

sub strip_brlens { # do substitution on string passed by reference
	my $tree = $_[0];
	#$$tree =~ s/:[-\.]*\d+(\.\d+){0,1}(e-\d+){0,1}//go;
	$$tree =~ s/: *(-*\d+\.e*-*\d*)//go;
}

sub by_clade { # remember, we converted commas to #s for each
               # clade being compared here;
               # sort first by number of taxa,
	       # next by branch length,
	       # finally alphabetically
	
	# initialize branch-length values (digits at the end)
	my(%num);
	#($num{$a} = $a) =~ s/.*: *(\d+(\.\d+){0,1})/$1/;
	#($num{$b} = $b) =~ s/.*: *(\d+(\.\d+){0,1})/$1/;
	#$b =~ /: *(-{0,1}\d*\.{0,1}\d+)$/o and $num{$b} = $1;
	$a =~ /: *(-{0,1}\d*\.{0,1}\d+(e-{0,1}\d+){0,1})$/o and $num{$a} = $1;
	$b =~ /: *(-{0,1}\d*\.{0,1}\d+(e-{0,1}\d+){0,1})$/o and $num{$b} = $1;

	# compare # counts (number of #s = number of taxa), 
	# then branch lengths, 
	# then names
	$a =~ tr/#/#/ <=> $b =~ tr/#/#/
	or
	$num{$a} <=> $num{$b}
	or
	$a cmp $b;
}

sub by_clade_reverse { # remember, we converted commas to #s for each
               # clade being compared here;
               # sort first by number of taxa,
	       # next by branch length,
	       # finally alphabetically
	
	# initialize branch-length values (digits at the end)
	my(%num);
	#($num{$a} = $a) =~ s/.*: *(\d+(\.\d+){0,1})/$1/;
	#($num{$b} = $b) =~ s/.*: *(\d+(\.\d+){0,1})/$1/;
	#$b =~ /: *(-{0,1}\d*\.{0,1}\d+)$/o and $num{$b} = $1;
	$a =~ /: *(-{0,1}\d*\.{0,1}\d+(e-{0,1}\d+){0,1})$/o and $num{$a} = $1;
	$b =~ /: *(-{0,1}\d*\.{0,1}\d+(e-{0,1}\d+){0,1})$/o and $num{$b} = $1;

	# compare # counts (number of #s = number of taxa), 
	# then branch lengths, 
	# then names
	$b =~ tr/#/#/ <=> $a =~ tr/#/#/
	or
	$num{$b} <=> $num{$a}
	or
	$b cmp $a;
}

sub offsets { 
	# given a string and a character,
	# return a list of offset values
	# for each occurrence in the string
	# of the character

	my($string,$char) = (@_);
	my(@offsets,$offset) = ();

	my $pos = -1;
	while ( ($pos = index($string,$char,$pos)) > -1 ) {
		push(@offsets,$pos);
		$pos++;
	}
	return(@offsets);
}

sub de_exponiate { # convert exponential notation (PAUP-style -- is it general?) to decimal
		   # (e.g. 6.02e-05 becomes 0.0000602) 
	# expect a reference to a string (convert the original)
	my $sr = shift;
	#$$sr =~ s/(\d)\.(\d+)e-(\d+)/'0.'.('0'x($3 - 1)).$1.$2/ge;
	$$sr =~ s/(\d\.\d+)e(-{0,1}\d+)/$1 * 10**$2/ge;
	return 1;
}

sub wrap_commas {
	# insert newlines (after commas only)
	# to wrap lines near 60 characters
	# (or another number, if given a second
	# argument)
	
	my($string,$line_wrap) = (@_);
	my($line) = '';
	my($formatted) = '';
	my($length) = 0;
	
	$line_wrap = 60 unless $line_wrap;
	
	@words = split(/,/,$string);
	while (@words) {
		my($length) = 0;
		my(@outwords) = ();
		
		# add on (comma-delimited) words 
		# until length is too great or
		# there are no words left
		#
		while (($length < $line_wrap) and @words) {
			push(@outwords,shift(@words));
			$length = length(join('x',@outwords));
		}

		# take off the last word 
		# if length really is too great 
		# and there is more than one word
		#
		if ($length > $line_wrap and scalar(@outwords) > 1) {
			unshift(@words,pop(@outwords));
		}

		# store the line
		$formatted .= join(',',@outwords,"\n");
	}
	# remove final comma
	$formatted =~ s/,\n\Z/\n/o;
	return $formatted;
}
