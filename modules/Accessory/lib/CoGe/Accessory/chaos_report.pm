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
__PACKAGE__->mk_accessors('file', 'hsps', 'hsp_count', 'query', 'subject',
	'qlength', 'slength', 'gotname', 'DEBUG');

###############################################################################
# chaos_report -- Josh Kane  UC Berkeley
###############################################################################
sub new {
	my $proto = shift;
	my $opts = shift;
	$opts = {} unless $opts;
	my $class = ref($proto) || $proto;
	my $self = bless ({%$opts}, $class);
	$self->hsp_count(0) unless $self->hsp_count;
	$self->process_file();
	return $self;
      }

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    my $tmp;
    my (@data,@hsps);
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
      my $hsp = $self->_parseReport($alignments);
      push @hsps,$hsp if $hsp;
    }
    $self->hsps(\@hsps);
    return $self;
}

sub _parseReport {
	my $self = shift;
	my $data = shift;
	my $hsp_count = $self->hsp_count();
	$hsp_count++;
	my ($locs,$scores) = ("","");
	my ($sq1,$align,$sq2) = ("","","");
	my ($score,$nmatches,$nga,$ngb,$nletters,$perc) = (0,0,0,0,0,0);
        my ($name1,$start1,$stop1,$name2,$start2,$stop2,$score_plus,$strand) = (0,0,0,0,0,0,0,0);
        my @aligns;
	foreach my $things (split /\n/,$data)
	{
         next unless $things;
	 if ($things =~/;/)
	   {
	     $locs = $things;
	   }
	 elsif ($things =~ /,/)
	   {
	     $scores = $things;
	   }
	 else
	   {
	     push @aligns,$things;
	   }
	}
	for (my $i=0;$i<scalar(@aligns);$i+=3)
	{
   	  $sq1 .= $aligns[$i];
   	  $align .= $aligns[$i+1];
   	  $sq2 .= $aligns[$i+2];
   	  }
   	$align =~s/:/|/g;
  	if ($locs =~/^(.*?)\s+(\d+)\s(\d+).\s(.*?)\s+(\d+)\s(\d+)\D+(\d+\.\d+)\s.(.)/)
    	{
    		($name1,$start1,$stop1,$name2,$start2,$stop2,$score_plus,$strand) = ($1,$2,$3,$4,$5,$6,$7,$8);
    		print STDERR "($name1,$start1,$stop1,$name2,$start2,$stop2,$score_plus,$strand)\n" if $self->DEBUG;
   	}
	if ($scores =~/^\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+).+(0?\.\d+)/)
    	{
      		($score,$nmatches,$nga,$ngb,$nletters,$perc) = ($1,$2,$3,$4,$5,$6);
    	}
    	while ($align && $align !~ /^\|/)
        	{$sq1 =~s/^.//;$sq2 =~s/^.//;$align =~s/^.//;$start1++;$start2++;}
  	while ($align && $align !~ /\|$/)
       		{$sq1 =~s/.$//;$sq2 =~s/.$//;$align =~s/.$//;$stop1--;$stop2--;}
  	$strand = $strand =~/-/ ? "-1" : "1";
  	($start1,$stop1) = ($stop1,$start1) if $start1 > $stop1;
  	($start2,$stop2) = ($stop2,$start2) if $start2 > $stop2;
  	$self->_getNames($name1,$name2) unless $self->gotname();
  	$perc = sprintf("%.4f", $perc) * 100;
	#explicitly count gaps -- output from chaos is screwy
	$nga = $sq1 =~ tr/-/-/;
	$ngb = $sq2 =~ tr/-/-/;

   	my $hsp = new CoGe::Accessory::parse_report::HSP
        ({
   	  score=>$score_plus,
    	  match=>$nmatches,
    	  length=>$nletters,
    	  query_start=>$start1,
    	  query_stop=>$stop1,
    	  subject_start=>$start2,
    	  subject_stop=>$stop2,
    	  query_alignment=>$sq1,
    	  subject_alignment=>$sq2,
    	  alignment=>$align,
    	  strand=>$strand,
    	  query_gaps=>$nga,
    	  subject_gaps=>$ngb,
    	  percent_id=>$perc,
    	  number=>$hsp_count,
	  pval=>"N/A",
        });
    $self->hsp_count($hsp_count) if $hsp;
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
