#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Graphics::Feature::GAGA;
use CoGe::Graphics::Feature::AminoAcid;
use CoGe::Genome;
use Data::Dumper;

my $form = new CGI;
my $db = new CoGe::Genome;
my $c = new CoGe::Graphics::Chromosome;
my $start = $form->param('start') || $form->param('x') || shift || 6940000;#6944100;#1;
my $stop = $form->param('stop') || shift;# || 6948000;#6949600;#190000;
$stop = $start unless $stop;
my $di = $form->param('di') || shift || 6;
my $chr = $form->param('chr') || shift || 1;
my $iw = $form->param('iw') || 20000;
my $mag = $form->param('m') || $form->param('mag');
my $file = $form->param('file');#|| "./tmp/pict.png";
unless ($start && $stop && $di && $chr)
  {
    print STDERR "missing needed parameters: Start: $start, Stop: $stop, Info_id: $di, Chr: $chr\n";
  }

my $chr_length = $db->get_genomic_sequence_obj->get_last_position($di);
$c->chr_length($chr_length);
$c->iw($iw);
$c->max_mag((80));
$c->DEBUG(0);
$c->feature_labels(1);
$c->fill_labels(1);
$c->draw_chromosome(1);
$c->draw_ruler(1);
$c->set_region(start=>$start, stop=>$stop);
if ($mag)
  {
    $c->mag($mag);
  }
else
  {
    $c->mag($c->mag-1);
  }
$start = $c->_region_start;
$stop= $c->_region_stop;
#process nucleotides
my $seq = uc($db->get_genomic_sequence_obj->get_sequence(start=>$start, end=>$stop, chr=>$chr, info_id=>$di)); 
my $seq_len = length $seq;
my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw);
$chrs = 1 if $chrs < 1;
my $pos = 0;
$start = 1 if $start < 1;

while ($pos < $seq_len)
  {
    my $subseq = substr ($seq, $pos, $chrs);
    my $rcseq = $subseq;
    $rcseq =~ tr/ATCG/TAGC/;
    next unless $subseq && $rcseq;
    my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$subseq, strand=>1, start =>$pos+$start});
    my $f2 = CoGe::Graphics::Feature::GAGA->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
    $c->add_feature($f1, $f2);
    $pos+=$chrs;
  }

#process features
foreach my $feat($db->get_feature_obj->get_features_in_region(start=>$start, end=>$stop, info_id=>$di, chr=>$chr))
  {
    my $f;
    if ($feat->type->name =~ /Gene/i)
      {
	$f = CoGe::Graphics::Feature::Gene->new();
	$f->color([255,0,0,50]);
	foreach my $loc ($feat->locs)
	  {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	  }
	$f->order(1);
      }
    elsif ($feat->type->name =~ /CDS/i)
      {
	$f = CoGe::Graphics::Feature::Gene->new();
	$f->color([0,255,0, 50]);
	foreach my $loc ($feat->locs)
	  {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	  }
	$f->order(3);
      }
    elsif ($feat->type->name =~ /mrna/i)
      {
	$f = CoGe::Graphics::Feature::Gene->new();
	$f->color([0,0,255, 50]);
	foreach my $loc ($feat->locs)
	  {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	  }
	$f->order(2);
      }
    elsif ($feat->type->name =~ /rna/i)
      {
	$f = CoGe::Graphics::Feature::Gene->new();
	$f->color([200,200,200, 50]);
	foreach my $loc ($feat->locs)
	  {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	    $f->add_type(1);
	  }
	$f->order(2);
      }
	#Do we have any protein sequence we can use?
	foreach my $seq ($feat->sequences)
	  {
	    next unless $seq->seq_type->name =~ /prot/i;
	    my ($pseq) = $seq->sequence_data;
	    $pseq = reverse $pseq if $f->strand =~ /-/;
	    my (@segs) = sort {$a->[0] <=> $b->[0]} @{$f->segments};  
            my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw);
	    $chrs = 1 if $chrs < 1;
	    my $len = length $pseq;
	    my $pos = 0; #protein sequence position
	    my $i = 0;   #array index of segments
	    my ($start, $stop) = @{$segs[$i]};
	    my $fstart = $start; #feature start
	    my $sseq; #sub sequence of protein sequence
	    foreach my $aa (split //, $pseq)
	      {
		if ($fstart >= $stop || ($sseq && length ($sseq) >= $chrs) )
		  {
#		    print  $start."-".$stop."  ".$pos."  ",$fstart,"  ".$chrs." ".$sseq, "\n";
		    my $ao = CoGe::Graphics::Feature::AminoAcid->new({aa=>$sseq, start=>$fstart, strand => $f->strand, order=>4}); #create aminoacid object
		    $c->add_feature($ao);
		    $fstart += (length($sseq))*3;
		    if ($fstart-1 >= $stop)
		      {
		        $i++; #imrement array index of segments
			#if the segment breaks a codon, we need to figure out how much to back up
			#the new start of the amino acid object
			my $sum = 0;
			for my $j (0 .. ($i-1))
			  {
			    my $tmp = $segs[$j];
			    $sum += ($tmp->[1]-$tmp->[0]+1);
			  }
			$fstart = (($sum)%3);
			$fstart = 3-$fstart if $fstart;
		        ($start, $stop) = @{$segs[$i]} if $segs[$i]; #assign new segment start and stop positions
			$fstart += $start;
		        $pos = 0; #reset the position index
		      }
		    $sseq = undef;
		  }
	        $pos++; #protein sequence position
		$sseq .= $aa;
	      }
	    if ($sseq) #make sure to create the final protein object if needed
	      {
		my $ao = CoGe::Graphics::Feature::AminoAcid->new({aa=>$sseq, start=>$fstart, strand => $f->strand, order=>4}); #create aminoacid object
		$c->add_feature($ao);
	      }
	  }
    my ($name) = map {$_->name} $feat->names;
    $f->label($name);
    $f->type($feat->type->name);
    $c->add_feature($f);

    
  }
if ($file)
  {
    $c->generate_png(file=>$file);
  }
else
  {
    print "Content-type: image/png\n\n";
    $c->generate_png();
  }


#foreach my $i (1..10)
# { 
#   $c->mag($i);
#   $c->generate_png(file=>"tmp/test$i.png");
# }
