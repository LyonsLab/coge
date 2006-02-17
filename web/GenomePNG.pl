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
my $start = $form->param('start') || $form->param('x') ||6932550;
my $stop = $form->param('stop');# || 6948000;#6949600;#190000;
$stop = $start unless $stop;
my $di = $form->param('di') || 6;
my $chr = $form->param('chr') ||$form->param('chromosome') || 1;
my $iw = $form->param('iw') || $form->param('width') || $form->param('tile size')|| $form->param('tile_size') || 256;
my $mag = $form->param('m') || $form->param('z') || $form->param('mag') || $form->param('magnification') || 4;
my $file = $form->param('file');# || "./tmp/pict.png";
unless ($start && $stop && $di && $chr)
  {
    print STDERR "missing needed parameters: Start: $start, Stop: $stop, Info_id: $di, Chr: $chr\n";
  }

my $chr_length = $db->get_genomic_sequence_obj->get_last_position($di);
$c->chr_length($chr_length);
$c->mag_scale_type("constant_power");
$c->iw($iw);
$c->max_mag((10));
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
#let's add the max top and bottom tracks to the image to keep it constant
my $f1= CoGe::Graphics::Feature->new({start=>1, order => 4, strand => 1});
$f1->iw(5);
$f1->ih(5);
$f1->gd->fill(2,2,$f1->get_color(255,255,255));
$c->add_feature($f1);
my $f2= CoGe::Graphics::Feature->new({start=>1, order => 4, strand => 1});
$f2->iw(5);
$f2->ih(5);
$f2->gd->fill(2,2,$f1->get_color(255,255,255));
$c->add_feature($f2);
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
	draw_prots(genomic_feat=>$feat, c=>$c, chrom_feat=>$f);
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
    next unless $f;
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

sub draw_prots
  {
    my %opts = @_;
    my $feat = $opts{genomic_feat};
    my $c = $opts{c};
    my $f = $opts{chrom_feat};
    #Do we have any protein sequence we can use?
    foreach my $seq ($feat->sequences)
      {
	next unless $seq->seq_type->name =~ /prot/i;
	my ($pseq) = $seq->sequence_data;
	my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw)/3;
	$chrs = 1 if $chrs < 1;
	my $pos = 0;
	while ($pos <= length $pseq)
	  {
	    my $aseq = substr($pseq, $pos, $chrs);
	    foreach my $loc ($seq->get_genomic_locations(start=>$pos+1, stop=>$pos+$chrs))
	      {
		my $ao = CoGe::Graphics::Feature::AminoAcid->new({aa=>$aseq, start=>$loc->start, stop=>$loc->stop, strand => $f->strand, order=>4});
		$ao->skip_overlap_search(1);
		$c->add_feature($ao);
		delete $loc->{__Changed}; #silence the warning from Class::DBI
	      }
	    
	    $pos+=$chrs;
	  }
      }
  }

