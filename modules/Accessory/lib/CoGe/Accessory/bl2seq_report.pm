package CoGe::Accessory::bl2seq_report;

use strict;
use warnings;
use CoGe::Accessory::parse_report::HSP;
use base qw(Class::Accessor);
use Data::Dumper;

BEGIN
  {
    use vars qw($VERSION);
    $VERSION = "0.01";
  }
__PACKAGE__->mk_accessors('file', 'lastline', 'report_done', 'hsps',
	'hsp_count', 'eval_cutoff', 'query', 'subject', 'qlength', 'slength');

sub qname    {shift->query(@_)}
sub sbjct    {shift->subject(@_)}
sub sname    {shift->subject(@_)}

###############################################################################
# bl2seqReport
###############################################################################
sub new {
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $opts = shift;
	$opts = {} unless $opts;
	my $self = bless ({%$opts}, $class);
	$self->hsp_count(0);
	# init
	$self->process_file();
	$self->hsps([]) unless $self->hsps;
	return $self;
      }

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    my @hsps;
    return $self unless $file;
    $self->file($file);
    open (IN, $file) || warn "can't open $file for reading: $!";
    my $data = join ("", <IN>);
    close IN;
    my $header;
    my $hsps = [];
    return $self if !$data || $data =~ /no hits found/i;
    ($header, $data) = $data =~ /^(.*?Length.*?\d+)\s+(Score.*$)/sx;
    $self->_parseHeader($header);
    return $self unless $data;  #no data?
    foreach my $hsp_data (split /Score/, $data)
      {
	next unless $hsp_data;
#	print $hsp_data;
#	next if $hsp_data =~ /Query=/;
       my $hsp = $self->_processHSP("Score".$hsp_data);
       push @hsps, $hsp if $hsp;
      }

    #	close( $self->{FH} );
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

sub _parseHeader
  {
    my ($self) = shift;
    my $header = shift;
    return unless $header;

# mdb removed 1/28/14, issue 288
#    $header =~ s/.*Query/Query/s;
#    my ($query, $subject) = split(/\n\n+/,$header);#=~ /^(.*letters\))\s+(.*)$/s;
#    warn "Problem parsing header from file: ".$self->file."\n header:\n$header\n" unless $query && $subject;
#    if ($query =~ /Query=\s+(.+)/)
#      {
#	my $qname = $1;
#	$self->query($qname);
#      }
#    else
#      {
#	warn "problem parsing bl2seq report query name from $query\n";
#	$self->query("UNKNOWN");
#      }
#    if ($query =~ /Length=\s*(\d+)/)
#      {
#	my $qlength = $1;
#	$self->qlength($qlength);
#      }
#    else
#      {
#	warn "problem parsing bl2seq report query length from $query\n";
#	$self->qlength("ERROR");
#      }
#    if ($subject =~ /Subject=\s+(.+)/)
#      {
#	my $sname = $1;
#	$self->subject($sname);
#      }
#    else
#      {
#	warn "problem parsing bl2seq report subject name from $subject\n";
#	$self->subject("UNKNOWN");
#      }
#    if ($subject =~ /Length=\s*(\d+)/i)
#      {
#	my $slength = $1;
#	$self->slength($slength);
#      }
#    else
#      {
#	warn "problem parsing bl2seq report subject length from $subject\n";
#	$self->slength("ERROR");
#      }

    # mdb added 1/28/14, issue 288
	my ($qname, $qlen, $sname, $slen) = $header =~ /.*Query=\s?(\S+).+Length=(\d+).+Subject=\s?(\S+).+Length=(\d+)/s;
	print STDERR "$qname $qlen $sname $slen\n";

	if (not defined $qname) {
		warn "problem parsing bl2seq report query name\n";
		$qname = 'UNKNOWN';
	}
	if (not defined $sname) {
		warn "problem parsing bl2seq report subject name\n";
		$sname = 'UNKNOWN';
	}
	if (not defined $qlen) {
		warn "problem parsing bl2seq report query length\n";
		$qlen = 'ERROR';
	}
	if (not defined $slen) {
		warn "problem parsing bl2seq report subject length\n";
		$slen = 'ERROR';
	}

	$self->query($qname);
	$self->qlength($qlen);
	$self->subject($sname);
	$self->slength($slen);

    return;
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

sub _processHSP {
	my $self = shift;
	my $data = shift;
	my $info;
	($info,$data) = split /Query/, $data,2;
	$data = "Query".$data;
	############################
	# get and parse scorelines #
	############################
	my ($score, $bits,$p) = $info =~
		/Score =\s+(\S+) bits \((\d+)\),\s+Expect.* += +(\S+)/;
	$p =~ s/,//g;
	my ($match, $length) = $info =~ /Identities = (\d+)\/(\d+)/;
	my ($positive) = $info =~ /Positives = (\d+)/;
	$positive = $match if not defined $positive;
	my $strand = "";
	if ($info =~ /Strand=(\S+)\/(\S+)/)
	  {
	    my ($strandtop, $strandbot) = ($1, $2);
	    if ( defined $strandtop and $strandtop eq "Plus" ) { $strand = "+" }
	    else { $strand = "-" }
	    if ( defined $strandbot and $strandbot eq "Plus" ) { $strand .= "+" }
	    else { $strand .= "-" }
	  }
	elsif ($info =~ /Frame = (\S+)/)
	  {
	    $strand = $1;
	  }

#	print join ("\t", $score, $bits, $p, $match, $length, $positive, $strand),"\n";
	#######################
	# get alignment lines #
	#######################
	my(@hspline) = ();;
	foreach (split /\n/, $data)
	  {
	  if ($_ !~ /\S/)
	    {next;} # blank line, skip
	  elsif ($_ =~ /(^>)|(^Lambda)|(^\s+Database:)/)
	    {
	      last;
	    }
	  else {
	    push @hspline, $_;           # store the query line
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
		$hspline[$i]   =~ /^Query\s+(\d+)\s*(\S+)\s+(\d+)/;
		$ql = $2; $qb = $1 unless $qb; $qe = $3;

		my $offset = index($hspline[$i], $ql);
		$as = substr($hspline[$i+1], $offset, CORE::length($ql))
			if $hspline[$i+1];

		$hspline[$i+2] =~ /^Sbjct\s+(\d+)\s*(\S+)\s+(\d+)/;
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
	($qb, $qe) = ($qe,$qb) if $qb > $qe;
	($sb, $se) = ($se,$sb) if $sb > $se;
	my $hsp = new CoGe::Accessory::parse_report::HSP
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
