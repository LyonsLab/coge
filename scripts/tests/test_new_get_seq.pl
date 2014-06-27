#!/usr/bin/perl -w

use strict;
use Benchmark;
use CoGeX;

my $coge = CoGeX->dbconnect();

my $dsgid = shift || 3068;

my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
my $faa_file = $dsg->file_path;
my $path = $faa_file;
$path =~ s/\/[^\/]*$//;
$path .= "/chr";
print $path,"\n";
my $t0 = new Benchmark;
my $i =0;
foreach my $feat ($coge->resultset('Feature')->search(
						      {
						       dataset_group_id=>$dsgid,
						       feature_type_id=>3,
						      },
						      {
						       join =>[{dataset=>'dataset_connectors'}],
#						       prefetch=>['locations', 'feature_names'],
						      }
						     ))
  {
    my ($name) = $feat->names;
    my ($seq2) = get_seq(feat=>$feat);
#    print $name,"\n";
#    my ($seq1) = $feat->genomic_sequence;
#     if ($seq1 eq $seq2)
#       {next;}
#     print join ("\t", $name,$feat->strand, $feat->chromosome, $feat->start, $feat->stop),"\n";
#     print ">".$seq1,"\n";
#     print ">".$seq2,"\n";
    $i++;
    print $i,"\n" unless $i%100;
  }
my $t1 = new Benchmark;
my $time = timestr(timediff($t1,$t0));
print "Time to run:  $time\n";

sub get_seq
  {
    my %opts = @_;
    my $feat = $opts{feat};
    my $tmp;
    open (IN, "$path/".$feat->chromosome);
    seek(IN, $feat->start-1, 0);
    read(IN, $tmp, $feat->stop-$feat->start+1);
    close IN;
    my @seqs;
    foreach my $loc (sort {$a->start <=> $b->start} $feat->locations)
      {
	my $sub_seq = substr($tmp, $loc->start-$feat->start, $loc->stop-$loc->start+1);
	if ($feat->strand =~ /-/)
	  {
	    unshift @seqs, reverse_complement($sub_seq);
	  }else{
	    push @seqs, $sub_seq;
	  }
      }
    my $seq = join( "", @seqs );
    return $seq;
  }

sub reverse_complement
  {
    my $seq = shift;# || $self->genomic_sequence;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/;
    return $rcseq;
  }
