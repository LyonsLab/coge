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
my $start = $form->param('start') || 25550;#1;
my $stop = $form->param('stop') || 27400;#190000;
$stop = $start unless $stop;
my $di = $form->param('di') || 8;
my $chr = $form->param('chr') || 3;

unless ($start && $stop && $di && $chr)
  {
    print STDERR "missing needed parameters: Start: $start, Stop: $stop, Info_id: $di, Chr: $chr\n";
  }

my $chr_length = $db->get_genomic_sequence_obj->get_last_position($di);
$c->chr_length($chr_length);
$c->iw(1600);
$c->max_mag((80));
#$c->num_mag(20);
$c->DEBUG(1);
$c->labels(0);
$c->fill_labels(1);
#print Dumper $c->mag_scale;
#exit;
#$stop = $chr_length;
$c->set_region(start=>25550, stop=>27400);
#$c->set_point(($stop-$start)/2);
#$c->set_point(25550);
print Dumper $c;
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
my %trans = (A=>'T',
             T=>'A',
	     C=>'G',
	     G=>'C',
	    );
my @seq;
foreach my $chr (split //, $db->get_genomic_sequence_obj->get_sequence(start=>$start, end=>$stop, chr=>$chr, info_id=>$di))
  {
     $chr = uc($chr);
     my $rc = $trans{$chr} || $chr;
     push @seq, [$chr, $rc];
  }
my $ i = 0;
foreach my $chr (@seq)
  {
    my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$chr->[0], strand=>1, start =>$i+$start});
    my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$chr->[1], strand=>-1, start =>$i+$start});
    $c->add_feature($f1, $f2);
    $i++;
    print "working on position ",$i+$start,"\n";
  }
foreach my $i (1..10)
 { 
   $c->mag($i);
   $c->generate_png(file=>"tmp/test$i.png");
 }
