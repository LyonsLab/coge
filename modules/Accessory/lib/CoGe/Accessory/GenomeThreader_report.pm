package CoGe::Accessory::GenomeThreader_report;

use strict;
use warnings;
use CoGe::Accessory::parse_report::HSP;
use base qw(Class::Accessor);
use Data::Dumper;
#use CoGeX;

#this only works when using protein sequences:
#./bin/gth -genomic eric_example/at_genomic.faa -protein eric_example/at_protein.faa -scorematrix /opt/apache/CoGe/data/blast/matrix/aa/BLOSUM62
__PACKAGE__->mk_accessors('file', 'genes', 'hsps', 'query', 'subject');

sub new {
        my $proto = shift;
        my $class = ref($proto) || $proto;
        my $opts = shift;
        $opts = {} unless $opts;
        my $self = bless ({%$opts}, $class);
	$self->genes([]);
	$self->hsps([]);
        # init
        $self->process_file() if $self->file();
        return $self;
      }

sub process_file
  {
    my $self = shift;
    my $file = shift || $self->file;
    unless (-r $file)
      {
	warn ("GenomeThreader_report.pm:  Can't read '$file' or no file was specified.  $!\n");
	return;
      }
    $/ = "********************************************************************************";
    open (IN, $file);
    my $count = 0;
    while (<IN>)
      {
	$self->process_block($_, $count);
	$count++;
      }
    close IN;
  }

sub process_block
  {
    my $self = shift;
    my $block = shift;
    my $num = shift;
    return unless $block && $block =~ /Protein Sequence:/;
    my ($query, $subject, $subj_strand, $subj_start, $subj_stop, $alignment);
    my @exons;
    my $alignment_flag=0;
    foreach my $line (split /\n/, $block)
      {
	if ($line =~ /^Protein Sequence/)
	  {
	    ($query) = $line =~ /Protein Sequence:.*description=(.*)/;
	  }
	elsif ($line =~ /^Genomic Template/)
	  {
	    ($subj_strand, $subj_start, $subj_stop, $subject) = $line =~ /Genomic Template:.*?strand=(.*?),.*from=(\d+).*to=(\d+).*description=(.*)/;
	  }
	elsif ($line =~ /^\s+Exon/)
	  {
#	    $num++;
	    my ($exon_num, $gen_start, $gen_stop) = $line =~ /Exon\s+(\d+)\s+(\d+)\s+(\d+)/;
	    ($gen_start, $gen_stop) = ($gen_stop, $gen_start) if $gen_start > $gen_stop;
	    my ($prot_start, $prot_stop) = $line =~ /Protein\s+(\d+)\s+(\d+)/;
	    push @exons, {
			 num=>$num,
			 genomic_start=>$gen_start,
			 genomic_stop=>$gen_stop,
			 protein_start=>$prot_start,
			 protein_stop=>$prot_stop,
			};
	  }
	elsif ($line =~ /^Alignment/)
	  {
	    $alignment_flag=1;
	    next;
	  }
	elsif ($line =~ /Predicted gene locations/)
	  {
	    last;
#	    $alignment_flag=0;  #turn of alignment flag if a query has more than one possible gene structure
	  }

	if ($alignment_flag)
	  {
	    $alignment .= $line."\n";
	  }
      }
    $subj_strand = $subj_strand =~ /-/ ? "-1" : 1;
#    my $query_info = parse_query_for_feature_info($query);
    foreach my $exon (@exons)
      {
	push @{$self->hsps}, new CoGe::Accessory::parse_report::HSP({
								     strand=>$subj_strand,
								     number=>$exon->{num},
								     subject_start=>$exon->{genomic_start},
								     subject_stop=>$exon->{genomic_stop},
								     query_start=>$exon->{protein_start},
								     query_stop=>$exon->{protein_stop},
								     alignment=>$alignment,
								     query_name=>$query,
								     subject_name=>$subject,
								     link_as_gene=>1,
								    });
      }
#    push @{$self->hsps}, {
#			  query=>$query,
#			  subject=>$subject,
#			  strand=>$subj_strand,
#			  subject_start=>$subj_start,
#			  subject_stop=>$subj_stop,
#			  alignment=>$alignment,
#			  exons=>\@exons,
#			 };
    return $num;
  }

# sub parse_query_for_feature_info
#   {
#     my $self = shift;
#     my $header = shift;
#     my %data;
#     my ($fid) = $header =~ /fid:(\d+)/;
#     return \%data unless $fid;
#     my $coge = CoGeX->dbconnect();
#     my $feat = $coge->resultset('Feature')->find($fid);
#     $data{strand} = $feat->strand;
#     foreach my $loc ($feat->locations)
#       {0
# 	push $data{locs},
#       }

#   }

1;
