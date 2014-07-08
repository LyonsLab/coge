package CoGe::Accessory::blast_report;

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
	'hsp_count', 'eval_cutoff', 'query', 'subject', 'qlength', 'slength',
	'limit');

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
	$self->query({});
	$self->subject({});
	# init
	$self->process_file();
	return $self;
      }

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    my $limit = shift || $self->limit;
    my @hsps;
    return $self unless $file;
    open (IN, $file) || die "can't open $file for reading: $!";
    my $header;
    my $hsps = [];
    my $oldris = $/;
    $/ = "Query=";
    loop: while (my $qsection = <IN>)
      {
#	next if $qsection =~ /Reference:/; #skip header
	$qsection = "Query=".$qsection;
	$qsection =~ s/Query=$//;
	my ($query, $qlength) = $self->_parseQuery($qsection);
	next unless $query;
	$self->query->{$query}=0 unless $self->query->{$query};
	foreach my $ssection (split />/,$qsection)
	  {
	    next if $ssection =~ /Sequences producing/;
	    $ssection = ">".$ssection;
	    my ($subject, $slength) = $self->_parseSubject($ssection);
	    next unless $subject;
	    $self->subject->{$subject}=0 unless $self->subject->{$subject};
	    foreach my $hsp_data (split /Score/, $ssection)
	      {
		next unless $hsp_data;
		next if $hsp_data =~ /^>/;
		$self->query->{$query}++;
		$self->query->{$subject}++;
		my $hsp = $self->_processHSP(data=>"Score".$hsp_data, query_name=>$query, query_length=>$qlength, subject_name=>$subject, subject_length=>$slength);
		$hsp->query_coverage($hsp->length/$qlength) if $qlength;
		$hsp->quality($hsp->query_coverage*$hsp->percent_id);
		push @hsps, $hsp if $hsp;
		last loop if $limit && $self->hsp_count > $limit;
	      }
	  }
      }
    @hsps = sort {$b->score <=> $a->score} @hsps;
    my $count = 1;
    foreach my $hsp (@hsps)
      {
	$hsp->number($count);
	$count++;
      }
    close IN;
    #	close( $self->{FH} );
    $/ = $oldris;
    $self->hsps(\@hsps);
    return $self;
  }

sub _parseQuery
  {
    my $self = shift;
    my $data = shift;
    my ($query, $length) = $data =~/Query=(.*?)\nLength=(\d+)/sx;
    return unless $query;
    $query =~ s/\n/ /gs;
    $query =~ s/ +/ /gs;
    $query =~ s/\s\s+/ /g;
    $query =~ s/^\s+//;
    $query =~ s/\s+$//;
    $length =~ s/,//g;
    return ($query, $length);
  }

sub _parseSubject
  {
    my $self = shift;
    my $data = shift;
    my ($subject, $length) = $data =~ />(.*?)\nLength=(\S+)/sx;
    return unless $subject;
    $subject =~ s/\n//g;
    $subject =~ s/ +/ /g;
    $length =~ s/,//g;
    return ($subject,$length);
  }

sub _processHSP {
	my $self = shift;
	my %opts = @_;
	my $data = $opts{data};
	$data =~ s/Query=$//;
	my $query_name = $opts{query_name};
	my $subject_name = $opts{subject_name};
	my $query_length = $opts{query_length};
	my $subject_length = $opts{subject_length};
	my $info;
	($info,$data) = split /Query/, $data,2;
	$data = "Query".$data;
#	print STDERR $info,"\n!!!\n",$data,"\nXXXXX\n";

	############################
	# get and parse scorelines #
	############################
	my ($score, $bits,$p) = $info =~
		/Score =\s*(\S+) bits\s\((\d+)\),\s+Expect.*\s*=\s*(\S+)/;
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
	  next unless $hspline[$i] =~ /^Query/;
		$hspline[$i]   =~ /^Query\s+(\d+)\s*(\S+)\s+(\d+)/;
		$ql = $2;
		$qb = $1 unless $qb;
		$qe = $3;

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
	    query_name=>$query_name,
	    subject_name=>$subject_name,
	    query_length=>$query_length,
	    subject_length=>$subject_length,
	   });
	return 0 if defined $self->eval_cutoff && $hsp->eval > $self->eval_cutoff;
	return $hsp;
      }

1;

__END__
