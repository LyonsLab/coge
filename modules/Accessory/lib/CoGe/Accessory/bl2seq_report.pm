package CoGe::Accessory::bl2seq_report;

use strict;
use warnings;
use CoGe::Accessory::bl2seq_report::HSP;
use base qw(Class::Accessor);

use Data::Dumper;


BEGIN
  {
    use vars qw($VERSION);
    $VERSION = "0.01";
    __PACKAGE__->mk_accessors qw(hsps hsp_count eval_cutoff);
  }


###############################################################################
# bl2seqReport
###############################################################################
sub new {
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self  = {};
	bless ($self, $class);
	my $file = "";
	if ( @_ ) {
		($file) = shift;
	} else {
		die "bl2seqReport error: new needs a file name!\n";
	}
	
	$self->{FILE} = undef;
	$self->{FH} = undef;
	$self->{LASTLINE} = undef;
	$self->{REPORT_DONE} = undef;
	$self->hsp_count(0);
	# init
	$self->{FILE} = $file;
	open( FH, "< $file" ) or
		die "bl2seqReport error: $file wouldn't open!\n";
	$self->{FH} = \*FH;
	if ( $self->_parseHeader() ) {
		$self->{REPORT_DONE} = 0;
	} else {
		$self->{REPORT_DONE} = 1;
	}
	my @hsps;
	while (my $hsp = $self->nextHSP)
	  {
	    push @hsps, $hsp;
	  }
	close( $self->{FH} );
	$self->hsps(\@hsps);
	return $self;
}

#sub DESTROY {
	# this object, until it's rewritten, misbehaves without this DESTROY
	# subroutine
#	my $self = shift;
#	print STDERR "bl2seqReport DESTROY\n";
	#close( $self->{FH} );
#	foreach my $key ( keys %$self ) {
#		$self->{$key} = undef;
#	}
#	print STDERR Dumper( $self );
#	$self = undef;
#}



sub query    {shift->{QUERY}}
sub qname    {shift->{QUERY}}
sub sbjct    {shift->{SBJCT}}
sub subject    {shift->{SBJCT}}
sub sname    {shift->{SBJCT}}
sub qlength    {shift->{QLENGTH}}
sub slength    {shift->{SLENGTH}}
sub file    {shift->{FILE}}

sub _parseHeader {
	my ($self) = shift;
	my $FH = $self->{FH};
	while(<$FH>) 
	  {
	    if ($_ =~ /^Query=\s+(.+)/)    { # get the query info
	      my $query = $1;
	      while(<$FH>) {
		last if $_ !~ /\S/;
		$query .= $_;
	      }
	      $query =~ s/\s+/ /g;
			$query =~ /(.*)\s+\((.*)\s+letters\)/;
			$self->{QUERY} = $1;
			my $qlength = $2;
			$qlength =~ s/,//g; #drop ",", if present
			$self->{QLENGTH} = $qlength;
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
			$self->{QUERY} = $1;
			my $qlength = $2;
			$qlength =~ s/,//g; #drop ",", if present
			$self->{QLENGTH} = $qlength;
		} elsif ($_ =~ /^>(.*)/) { # get the subject info
	
			my $sbjctdef = $1;
			while(<$FH>) {
				last if $_ !~ /\S/;
				$sbjctdef .= $_;
			}
			$sbjctdef =~ s/\s+/ /g;
			$sbjctdef =~ /(.*)\s+Length =\s+(.*)\s*/;
			$self->{SBJCT} = $1;
			my $slength = $2;
			$slength =~ s/(.*)\s+/$1/;
			$self->{SLENGTH} = $slength;
		} elsif ($_ =~ / Score/) {
			$self->{LASTLINE} = $_;
			return 1;
		} elsif ($_ =~ /^Lambda/) {
			$self->{LASTLINE} = $_;
			return 0; # there's nothing in the report
		}
	}
}

sub _fastForward {
	my ($self) = shift;
	return 0 if $self->{REPORT_DONE};
	return 1 if $self->{LASTLINE} =~ / Score/;
	if ($self->{LASTLINE} =~ /^Lambda/) {
		$self->{REPORT_DONE} = 1;
		return 1;
	}
	my $FH = $self->{FH};
	while(<$FH>) {
		if ($_ =~ /^Lambda/) {
			$self->{LASTLINE} = $_;
			$self->{REPORT_DONE} = 1;
			return 1;
		}
	}
	warn "Possible parse error in _fastForward in bl2seqReport\n";
}

sub nextHSP {
	my $self = shift;
	$self->_fastForward() or return 0;
	return 0 if $self->{REPORT_DONE};
	
	############################
	# get and parse scorelines #
	############################
	my $scoreline = $self->{LASTLINE};
	my $FH = $self->{FH};
	my ($bits, $score,$p) = $scoreline =~
		/Score =\s+(\S+) bits \((\d+)\), Expect.* = (\S+)/;
	$p =~ s/,//g;
	my $identityline = <$FH>;
	return undef if not defined $identityline;
	
	my ($match, $length) = $identityline =~ /Identities = (\d+)\/(\d+)/;
	my ($positive) = $identityline =~ /Positives = (\d+)/;
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
		elsif ($_ =~ /^\s*Score/)     {$self->{LASTLINE} = $_; last}
		elsif ($_ =~ /^>|^Lambda/)   {
			$self->{LASTLINE} = $_;
			$self->{REPORT_DONE} = 1;
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
	my $hsp_count = $self->hsp_count();
	$hsp_count++;
	$self->hsp_count($hsp_count);
	my $qgaps = $ql =~ tr/-/-/;
	my $sgaps = $sl =~ tr/-/-/;
	my $hsp = new CoGe::Accessory::bl2seq_report::HSP
	  ({
	    score=>$score,
	    bits=>$bits,
	    match=>$match,
	    positive=>$positive,
	    length=>$length,
	    pval=>$p,
	    query_start=>$qb,
	    query_stop=>$qe,
	    subject_start=>$sb,
	    subject_stop=>$se,
	    query_alignment=>$ql,
	    subject_alignment=>$sl,
	    alignment=>$as,
	    query_gaps=>$qgaps,
	    subject_gaps=>$sgaps,
	    strand=>$strand,
	    number=>$hsp_count,
	   });
	return 0 if defined $self->eval_cutoff && $hsp->eval > $self->eval_cutoff;
	return $hsp;
}

1;

__END__
