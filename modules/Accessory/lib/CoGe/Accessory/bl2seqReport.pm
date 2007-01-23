package CoGe::Accessory::bl2seqReport;
use strict;
use warnings;
use Data::Dumper;


###############################################################################
# bl2seqReport
###############################################################################
sub new {
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $this  = {};
	bless ($this, $class);
	my $file = "";
	if ( @_ ) {
		($file) = shift;
	} else {
		die "bl2seqReport error: new needs a file name!\n";
	}
	
	$this->{FILE} = undef;
	$this->{FH} = undef;
	$this->{LASTLINE} = undef;
	$this->{REPORT_DONE} = undef;

	# init
	$this->{FILE} = $file;
	open( FH, "< $file" ) or
		die "bl2seqReport error: $file wouldn't open!\n";
	$this->{FH} = \*FH;
	if ( $this->_parseHeader() ) {
		$this->{REPORT_DONE} = 0;
	} else {
		$this->{REPORT_DONE} = 1;
	}
	return $this;
}

sub DESTROY {
	# this object, until it's rewritten, misbehaves without this DESTROY
	# subroutine
	my $this = shift;
#	print STDERR "bl2seqReport DESTROY\n";
	close( $this->{FH} );
	foreach my $key ( keys %$this ) {
		$this->{$key} = undef;
	}
#	print STDERR Dumper( $this );
	$this = undef;
}



sub query    {shift->{QUERY}}
sub qname    {shift->{QUERY}}
sub sbjct    {shift->{SBJCT}}
sub subject    {shift->{SBJCT}}
sub sname    {shift->{SBJCT}}
sub qlength    {shift->{QLENGTH}}
sub slength    {shift->{SLENGTH}}
sub file    {shift->{FILE}}

sub _parseHeader {
	my ($this) = shift;
	my $FH = $this->{FH};
	while(<$FH>) {
		if ($_ =~ /^Query=\s+(.+)/)    { # get the query info
			my $query = $1;
			while(<$FH>) {
				last if $_ !~ /\S/;
				$query .= $_;
			}
			$query =~ s/\s+/ /g;
			$query =~ /(.*)\s+\((.*)\s+letters\)/;
			$this->{QUERY} = $1;
			my $qlength = $2;
			$qlength =~ s/,//g; #drop ",", if present
			$this->{QLENGTH} = $qlength;
		} elsif ( $_ =~ /^Query=/ ) {
			# this happens when the original fasta sequence
			# used in the bl2seq program contained a single
			# word on the header line and only happens to
			# the Query, not the Subject!
			my $query = "UNK-ACCN ";
			while(<$FH>) {
				last if $_ !~ /\S/;
				$query .= " " . $_;
			}
			$query =~ s/\s+/ /g;
			$query =~ /(.*)\s+\((.*)\s+letters\)/;
			$this->{QUERY} = $1;
			my $qlength = $2;
			$qlength =~ s/,//g; #drop ",", if present
			$this->{QLENGTH} = $qlength;
		} elsif ($_ =~ /^>(.*)/) { # get the subject info
	
			my $sbjctdef = $1;
			while(<$FH>) {
				last if $_ !~ /\S/;
				$sbjctdef .= $_;
			}
			$sbjctdef =~ s/\s+/ /g;
			$sbjctdef =~ /(.*)\s+Length =\s+(.*)\s*/;
			$this->{SBJCT} = $1;
			my $slength = $2;
			$slength =~ s/(.*)\s+/$1/;
			$this->{SLENGTH} = $slength;
		} elsif ($_ =~ / Score/) {
			$this->{LASTLINE} = $_;
			return 1;
		} elsif ($_ =~ /^Lambda/) {
			$this->{LASTLINE} = $_;
			return 0; # there's nothing in the report
		}
	}
}

sub _fastForward {
	my ($this) = shift;
	return 0 if $this->{REPORT_DONE};
	return 1 if $this->{LASTLINE} =~ / Score/;
	if ($this->{LASTLINE} =~ /^Lambda/) {
		$this->{REPORT_DONE} = 1;
		return 1;
	}
	my $FH = $this->{FH};
	while(<$FH>) {
		if ($_ =~ /^Lambda/) {
			$this->{LASTLINE} = $_;
			$this->{REPORT_DONE} = 1;
			return 1;
		}
	}
	warn "Possible parse error in _fastForward in bl2seqReport\n";
}

sub nextHSP {
	my $this = shift;
	$this->_fastForward() or return 0;
	return 0 if $this->{REPORT_DONE};

	
	############################
	# get and parse scorelines #
	############################
	my $scoreline = $this->{LASTLINE};
	my $FH = $this->{FH};
	my ($bits, $score,$p) = $scoreline =~
		/Score =\s+(\S+) bits \((\d+)\), Expect.* = (\S+)/;
	$p =~ s/,//g;
	my $identityline = <$FH>;
	return undef if not defined $identityline;
	my ($match, $length) = $identityline =~ /Identities = (\d+)\/(\d+)/;
	my ($positive) = $scoreline =~ /Positives = (\d+)/;
	$positive = $match if not defined $positive;
	my $strandline = <$FH>;
	return undef if not defined $strandline;
	my $strand = "";
	if ($strandline =~ /Strand/)
	  {
	    my ($strandtop, $strandbot) = $strandline =~ /Strand = (\S+) \/ (\S+)/;
	

	    if ( defined $strandtop and $strandtop eq "Plus" ) { $strand = "+" }
	    else { $strand = "-" }
	    if ( defined $strandbot and $strandbot eq "Plus" ) { $strand .= "+" }
	    else { $strand .= "-" }
	  }
	elsif ($strandline =~ /Frame = (\S+)/)
	  {
	    $strand = $1;
	  }
	die "parse error\n" if not defined $score;

	#######################
	# get alignment lines #
	#######################
	my(@hspline) = ();;
	while(<$FH>) {
		if ($_ !~ /\S/) {next;} # blank line, skip
		elsif ($_ =~ /^\s*Score/)     {$this->{LASTLINE} = $_; last}
		elsif ($_ =~ /^>|^Lambda/)   {
			$this->{LASTLINE} = $_;
			$this->{REPORT_DONE} = 1;
			last;
		}
		else {
			push @hspline, $_;           # store the query line
			my $l1 = <$FH>;              # either alignment line or sbjct line
			if ($l1 =~ /^Sbjct/) {
				push @hspline, "";  # dummy line, this is a -noseq option
				push @hspline, $l1; # so store a fake alignment and real sbjct
				next;
			}
			push @hspline, $l1;                 # grab/store the alignment line
			my $l2 = <$FH>; push @hspline, $l2; # grab/store the sbjct line
		}
	}
	
	#########################
	# parse alignment lines #
	#########################
	my ($ql, $sl, $as) = ("", "", "");
	my ($qb, $qe, $sb, $se) = (0,0,0,0);
	my (@QL, @SL, @AS) = (); # for better memory management
			
	for(my $i=0;$i<@hspline;$i+=3) {
		#warn $hspline[$i], $hspline[$i+2];
		$hspline[$i]   =~ /^Query:\s+(\d+)\s*(\S+)\s+(\d+)/;
		$ql = $2; $qb = $1 unless $qb; $qe = $3;
		
		my $offset = index($hspline[$i], $ql);
		$as = substr($hspline[$i+1], $offset, CORE::length($ql))
			if $hspline[$i+1];
		
		$hspline[$i+2] =~ /^Sbjct:\s+(\d+)\s*(\S+)\s+(\d+)/;
		$sl = $2; $sb = $1 unless $sb; $se = $3;
		
		push @QL, $ql; push @SL, $sl; push @AS, $as;
	}

	##################
	# the HSP object #
	##################
	$ql = join("", @QL);
	$sl = join("", @SL);
	$as = join("", @AS);
	my $qgaps = $ql =~ tr/-/-/;
	my $sgaps = $sl =~ tr/-/-/;
	my $hsp = bl2seqReport::HSP::new( $score, $bits, $match, $positive,
	$length, $p, $qb, $qe, $sb, $se, $ql, $sl, $as, $qgaps, $sgaps,
	$strand);
	return $hsp;
}

###############################################################################
# bl2seqReport::HSP
###############################################################################
package bl2seqReport::HSP;
use Data::Dumper;
sub new {
	my $hsp = bless {};
	($hsp->{SCORE}, $hsp->{BITS},
		$hsp->{MATCH}, $hsp->{POSITIVE}, $hsp->{LENGTH},$hsp->{P},
		$hsp->{QB}, $hsp->{QE}, $hsp->{SB}, $hsp->{SE},
		$hsp->{QL}, $hsp->{SL}, $hsp->{AS}, $hsp->{QG}, $hsp->{SG}, 
		$hsp->{STRAND} ) = @_;
	$hsp->{PERCENT} = int(1000 * $hsp->{MATCH}/$hsp->{LENGTH})/10;
	if ( $hsp->{P} =~ /^e/ ) {
		$hsp->{P} = "1" . $hsp->{P};
	}
	return $hsp;
}
sub score           {shift->{SCORE}}
sub bits            {shift->{BITS}}
sub percent         {shift->{PERCENT}}
sub match           {shift->{MATCH}}
sub positive        {shift->{POSITIVE}}
sub length          {shift->{LENGTH}}
sub P               {shift->{P}}
sub p               {shift->{P}}
sub P_val               {shift->{P}}
sub p_val               {shift->{P}}
sub pval               {shift->{P}}
sub eval               {shift->{P}}
sub queryBegin      {shift->{QB}}
sub queryEnd        {shift->{QE}}
sub sbjctBegin      {shift->{SB}}
sub sbjctEnd        {shift->{SE}}
sub queryAlignment  {shift->{QL}}
sub sbjctAlignment  {shift->{SL}}
sub alignmentString {shift->{AS}}
sub queryGaps       {shift->{QG}}
sub sbjctGaps       {shift->{SG}}
sub qbegin      {shift->{QB}}
sub qend        {shift->{QE}}
sub sbegin      {shift->{SB}}
sub send        {shift->{SE}}
sub qalign  {shift->{QL}}
sub salign  {shift->{SL}}
sub alignment_string {shift->{AS}}
sub qgaps       {shift->{QG}}
sub sgaps       {shift->{SG}}
sub qb              {shift->{QB}}
sub qe              {shift->{QE}}
sub sb              {shift->{SB}}
sub se              {shift->{SE}}
sub qa              {shift->{QL}}
sub sa              {shift->{SL}}
sub as              {shift->{AS}}
sub qg              {shift->{QG}}
sub sg              {shift->{SG}}
sub strand          {shift->{STRAND}}

sub DESTROY {
	my $this = shift;
	foreach my $key ( keys %$this ) {
		$this->{$key} = undef;
	}
}


1;
__END__
