#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Genome;
use Data::Dumper;

my $form = new CGI;
my $db = new CoGe::Genome;
my $c = new CoGe::Graphics::Chromosome;
my $start = $form->param('start') || $form->param('x') || shift || 19666;#1;
my $stop = $form->param('stop') || shift || 20451;#190000;
$stop = $start unless $stop;
my $di = $form->param('di') || shift || 8;
my $chr = $form->param('chr') || shift || 3;
my $iw = $form->param('iw') || 1600;
my $mag = $form->param('m') || $form->param('mag');
my $file = $form->param('file');# || "./tmp/pict.png";
unless ($start && $stop && $di && $chr)
  {
    print STDERR "missing needed parameters: Start: $start, Stop: $stop, Info_id: $di, Chr: $chr\n";
  }

my $chr_length = $db->get_genomic_sequence_obj->get_last_position($di);
$c->chr_length($chr_length);
$c->iw($iw);
$c->max_mag((80));
$c->DEBUG(0);
$c->feature_labels(0);
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
    elsif ($feat->type->name =~ /rna/i)
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
    my ($name) = map {$_->name} $feat->names;
    $f->label($name);
    $f->type($feat->type->name);

    $c->add_feature($f);
  }
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
#    print STDERR $subseq,"\n";
    my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$subseq, strand=>1, start =>$pos+$start});
    my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
    $c->add_feature($f1, $f2);
    $pos+=$chrs;
    #    print "working on position ",$i+$start,"\n";
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
