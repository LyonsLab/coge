package CoGe::Accessory::chaos_report;

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
__PACKAGE__->mk_accessors qw(file hsps hsp_count query subject qlength slength gotname);

###############################################################################
# chaos_report -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->hsp_count(0) unless $self->hsp_count;;
	$self->process_file();
	return $self;
      }

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    my $tmp;
    my @data;
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
    {
     $tmp .= $_;
     if (/perc/)
     {push @data,$tmp;$tmp="";}
    }
    close IN;
    foreach my $alignments (@data)
    {
      my $hsp = $self->get_hsp($alignments);
      push @hsps,$hsps if $hsps;
    }
    $self->hsps(\@hsps);
    return $self;
}


sub _parseReport {
	my $self = shift;
	my $data = shift;
	my $hsp_count = $self->hsp_count();
	$hsp_count++;
	$self->hsp_count($hsp_count);
	my ($locs,$scores) = ("","");
	my ($sq1,$align,$sq2) = ("","","");
	my ($score,$nmatches,$nga,$ngb,$nletters,$perc) = (0,0,0,0,0,0);
        my ($name1,$start1,$stop1,$name2,$start2,$stop2,$score_plus,$strand) = (0,0,0,0,0,0,0,0);
	foreach my $things (split /\n/,$data)
	{ 
         next unless $things;
         $locs = $things if $things =~/;/;
      	 $scores = $things if $things =~/,/;
         push @aligns,$things if $things =~/([ATGC]|:)/;
	}
	for (my $i=0;$i<scalar(@$ref);$i+=3)
	{
   	  $sq1 .= $ref->[$i];
   	  $align .= $ref->[$i+1];
   	  $sq2 .= $ref->[$i+2];
   	  }
   	$align =~s/:/|/g;
	while ($align !~ /^\|/)
        	{$sq1 =~s/^.//;$sq2 =~s/^.//;$align =~s/^.//;}
  	while ($align !~ /\|$/)
       		{$sq1 =~s/.$//;$sq2 =~s/.$//;$align =~s/.$//;}
  	if ($locs =~/^(\w+\.?\d*)\s(\d+)\s(\d+).\s(\w+\.?\d*)\s(\d+)\s(\d+).+(\d+\.\d+)\s.(.)/)
    	{
    		($name1,$start1,$stop1,$name2,$start2,$stop2,$score_plus,$strand) = ($1,$2,$3,$4,$5,$6,$7,$8);
   	}
	if ($scores =~/^\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+).+(0?\.\d+)/)
    	{
      		($score,$nmatches,$nga,$ngb,$nletters,$perc) = ($1,$2,$3,$4,$5,$6);
    	}
  	$self->_getNames($name1,$name2) unless $self->gotname();
  	$perc *=100;
   	my $hsp = new CoGe::Accessory::parse_report::HSP
        ({
   	  score=>$score,
    	  match=>$nmatches,
    	  length=>$nletters,
    	  query_start=>$start1,
    	  query_stop=>$stop1,
    	  subject_start=>$start2,
    	  subject_stop=>$stop2,
    	  query_alignment=>$sq1,
    	  subject_alignment=>$sq2,
    	  alignment=>$align,
    	  query_gaps=>$nga,
    	  subject_gaps=>$ngb,
    	  percent_id=>$perc,
    	  number=>$count,
        });
    return $hsp;
 }
 
sub _getNames {
	my $self = shift;
	my $query = shift;
	my $subject = shift;
	$self->query($query);
	$self->subject($subject);
	$self->gotname(1);
	return $self;
}
	

1;

__END__
