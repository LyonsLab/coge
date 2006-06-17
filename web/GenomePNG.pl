#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Graphics::Feature::GAGA;
use CoGe::Graphics::Feature::Exon_motifs;
use CoGe::Graphics::Feature::AminoAcid;
use CoGe::Graphics::Feature::Domain;
use CoGe::Genome;
use Data::Dumper;
use Benchmark;


my $t0 = new Benchmark;
my $form = new CGI;
print STDERR $form->self_url(-full=>1),"\n";
my $db = new CoGe::Genome;
my $c = new CoGe::Graphics::Chromosome;
my $start = $form->param('start') || $form->param('x') ||0;#28520458;
my $stop = $form->param('stop');# || 6948000;#6949600;#190000;
$start = 1 unless defined $start;
#$stop = $start unless $stop;
my $di = $form->param('di');
my $version = $form->param('v') || $form->param('version');
my $org_id = $form->param('o') || $form->param('org') ||$form->param('organism') || $form->param('org_id') || $form->param('oid');
my $chr = $form->param('chr') ||$form->param('chromosome');
my $iw = $form->param('iw') || $form->param('width') || $form->param('tile size')|| $form->param('tile_size') || 256;
my $mag = $form->param('m') || $form->param('mag') || $form->param('magnification');
my $z = $form->param('z');
my $file = $form->param('file');# || "./tmp/pict.png";
my $start_pict = $form->param('start_pict');
my $simple = $form->param('simple');
my $chr_start_height = $form->param('csh') || 200;
my $chr_mag_height = $form->param('cmh') || 5;
my $feat_start_height = $form->param('fsh') || 10;
my $feat_mag_height = $form->param('fmh') || 2;
my $BENCHMARK = $form->param('bm') || 0;


##some limits to keep from blowing our stack
use vars qw($MAX_FEATURES $MAX_NT);
$MAX_FEATURES = $iw*10; #one feature per pixel
$MAX_NT       = $iw*100000; #50,000 nt per pixel


unless (defined $start && $di && $chr)
  {
    print STDERR "missing needed parameters: Start: $start, Stop: $stop, Info_id: $di, Chr: $chr\n";
  }

if ($di) #there can be additional information about a chromosome for a particular version of the organism that is not in the same data_information.  Let's go find the organism_id and version for the specified data information item.
  {
    my $dio = $db->get_data_info_obj->search({data_information_id=>$di})->next;
    $org_id = $dio->org->id unless $org_id;
    $version = $dio->version unless $version;
  }

my %dids; #we will store a list of data_information objects that 
my $ta = new Benchmark;
if ($org_id)
  {

    foreach my $did ( $db->get_data_info_obj->search({organism_id=>$org_id, version=>$version}))
      {
	$dids{$did} = 1;
      }
  }
my $tb = new Benchmark;
my $finddid_time = timestr(timediff($tb, $ta));

$dids{$di}=1 if $di;
unless ($chr) {($chr) = $db->get_genomic_seq_obj->search({data_information_id=>$di})->next->chr if $di};
my @dids = keys %dids;
my $tc = new Benchmark;
foreach my $did (@dids)
  {

    my ($tstart, $tstop) = initialize_c(di=>$did,
					chr=>$chr,
					iw=>$iw,
					z=>$z,
					mag=>$mag,
					start=>$start,
					stop=>$stop,
					db=>$db,
					c=>$c,
					start_pict=>$start_pict,
					csh=>$chr_start_height,
					cmh=>$chr_mag_height,
					fsh=>$feat_start_height,
					fmh=>$feat_mag_height,
		 ) unless $c->chr_length;
    
    if ($c->chr_length)
      {
        $start = $tstart;
	$stop = $tstop;
        process_nucleotides(start=>$start, stop=>$stop, chr=>$chr, di=>$did, db=>$db, c=>$c) unless $simple;
	last;
      }
  }
my $td = new Benchmark;
my $init_c_time = timestr(timediff($td, $tc));

unless ($c->chr_length)
  {
    warn "error initializing the chromosome object.  Failed for valid chr_length\n";
    exit(0);
  }
my $t1 = new Benchmark;
foreach my $did (keys %dids)
  {
    my $taa = new Benchmark;
    process_features(start=>$start, stop=>$stop, chr=>$chr, di=>$did, db=>$db, c=>$c) unless $simple;
    my $tab = new Benchmark;
    my $feat_time = timestr(timediff($tab, $taa));
    print STDERR " processing did $did:   $feat_time\n" if $BENCHMARK;
  }
my $t2 = new Benchmark;
generate_output(file=>$file, c=>$c);	
my $t3 = new Benchmark;
my $init_time = timestr(timediff($t1, $t0));
my $feature_time = timestr(timediff($t2, $t1));
my $draw_time = timestr(timediff($t3, $t2));

print STDERR qq{
BENCHMARKING GenomePNG.pl

  Finding DID:          $finddid_time
  Init chr obj:         $init_c_time
Total Initialization:   $init_time
Processing features:    $feature_time
Image generation:       $draw_time


} if $BENCHMARK;

1;
sub initialize_c
  {
    my %opts = @_;
    my $di = $opts{di};
    my $chr = $opts{chr};
    my $iw = $opts{iw};
    my $z = $opts{z};
    my $mag = $opts{mag};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $db = $opts{db};
    my $c = $opts{c};
    my $csh=$opts{csh};
    my $cmh=$opts{cmh};
    my $fsh=$opts{fsh};
    my $fmh=$opts{fmh};
    my $start_pict = $opts{'start_pict'} || 'left';
    my $chr_length = $db->get_genomic_sequence_obj->get_last_position(di=>$di, chr=>$chr);
    return unless $chr_length;
    $c->chr_length($chr_length);
    $c->mag_scale_type("constant_power");
    $c->iw($iw);
    $c->max_mag((10));
    $c->DEBUG(0);
    $c->feature_labels(1);
    $c->fill_labels(1);
    $c->draw_chromosome(1);
    $c->draw_ruler(1);
    $c->chr_start_height($csh);
    $c->chr_mag_height($cmh);
    $c->feature_start_height($fsh);
    $c->feature_mag_height($fmh);
    if (defined $z) #the $z val is used by the program for making tiles of genomic views.
                #by convention, a z value of 0 means maximum magnification which is
        	#opposite the convention used in chromosome.pm.  Thus, we need
        	#to reformat the z value appropriately
      {
         my ($max) = sort {$b <=> $a} keys %{$c->mag_scale};
         $mag = $max-$z;
         $mag = 1 if $mag < 1;
         $mag = $max if $mag > $max;
         $c->start_picture($start_pict);
      }
    
    if ($mag)
      {
        $c->mag($mag);
      }
    else
      {
        $c->mag($c->mag-1);
      }
    $c->set_region(start=>$start, stop=>$stop);
    $start = $c->_region_start;
    $stop= $c->_region_stop;
    #let's add the max top and bottom tracks to the image to keep it constant
    my $f1= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => 1});
    $f1->merge_percent(0);
    $c->add_feature($f1);
    my $f2= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => -1});
    $f2->merge_percent(0);
    $c->add_feature($f2);
    return ($start, $stop);
}
sub process_nucleotides
  {
    my %opts = @_;
    my $start = $opts{start};
    my $stop = $opts{stop};
    if (abs ($stop-$start) > $MAX_NT)
      {
	warn "exceeded nucleotide limit of $MAX_NT (requested: ".abs($stop-$start).")\n";
      }
    my $chr = $opts{chr};
    my $di = $opts{di};
    my $db = $opts{db};
    my $c = $opts{c};
#    print STDERR "IN $0: start: $start, stop: $stop\n";
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
        my $rcseq = substr ($seq, $pos, $chrs);
        $rcseq =~ tr/ATCG/TAGC/;
        next unless $subseq && $rcseq;
        my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$subseq, strand=>1, start =>$pos+$start});
	my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
	#my $f2 = CoGe::Graphics::Feature::Exon_motifs->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
        #my $f2 = CoGe::Graphics::Feature::GAGA->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
        
        $c->add_feature($f1, $f2);
        $pos+=$chrs;
      }
  }

sub process_features
  {
    #process features
    my %opts = @_;
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $chr = $opts{chr};
    my $di = $opts{di};
    my $db = $opts{db};
    my $c = $opts{c};
    my $accn = $opts{accn};
    my $print_names = $opts{print_names};
    my $feat_count = $db->get_feature_obj->count_features_in_region(start=>$start, end=>$stop, info_id=>$di, chr=>$chr);
    if ($feat_count > $MAX_FEATURES)
      {
	warn "exceeded maximum number of features $MAX_FEATURES. ($feat_count requested)\nskipping.\n";
	return;
      }
    foreach my $feat($db->get_feature_obj->get_features_in_region(start=>$start, end=>$stop, info_id=>$di, chr=>$chr))
      {
        my $f;
#	print STDERR Dumper $feat;
#	print STDERR "!",join ("\t", map {$_->name} $feat->names, map {$_->start."-".$_->stop} $feat->locs),"\n";
	
        if ($feat->type->name =~ /Gene/i)
          {
#	    next;
#	    print STDERR $feat->names->next->name,": ", $feat->id,"\n";
	    $f = CoGe::Graphics::Feature::Gene->new();
	    $f->color([255,0,0,50]);
	    foreach my $loc ($feat->locs)
	      {
		$f->add_segment(start=>$loc->start, stop=>$loc->stop);
		$f->strand($loc->strand);
	      }
	    $f->order(1);
	    $f->overlay(1);
	    $f->mag(0.5);
          }
        elsif ($feat->type->name =~ /CDS/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,255,0, 50]);
        	$f->order(1);
		$f->overlay(3);
		draw_prots(genomic_feat=>$feat, c=>$c, chrom_feat=>$f);
		if ($accn)
		  {
		    foreach my $name (@{$feat->{QUALIFIERS}{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
		      }
		  }
          }
        elsif ($feat->type->name =~ /mrna/i || $feat->type->name =~ /exon/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,0,255, 50]);
        	$f->order(1);
		$f->overlay(2);
		$f->mag(0.75);
          }
        elsif ($feat->type->name =~ /rna/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([200,200,200, 50]);
        	$f->order(1);
		$f->overlay(2);
		if ($accn)
		  {
		    foreach my $name (@{$feat->{QUALIFIERS}{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
		      }
		  }

          }
        elsif ($feat->type->name =~ /functional domains/i)
          {
	    $f = CoGe::Graphics::Feature::Domain->new();
	    foreach my $loc ($feat->locs)
	      {
	        $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	        $f->strand($loc->strand);
	      }
	    $f->order(2);
          }
        next unless $f;
        foreach my $loc ($feat->locs)
	  {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	  }
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} map {$_->name} $feat->names;
	$f->label($name) if $print_names;
        $f->type($feat->type->name);
        $c->add_feature($f);
    }
  }

sub generate_output
  {
    my %opts = @_;
    my $file = $opts{file};
    my $c = $opts{c};
    if ($file)
      {
        $c->generate_png(file=>$file);
      }
    else
      {
        print "Content-type: image/png\n\n";
        $c->generate_png();
      }
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
		my $ao = CoGe::Graphics::Feature::AminoAcid->new({aa=>$aseq, start=>$loc->start, stop=>$loc->stop, strand => $loc->strand, order=>2});
		$ao->skip_overlap_search(1);
		$c->add_feature($ao);
		delete $loc->{__Changed}; #silence the warning from Class::DBI
	      }
	    
	    $pos+=$chrs;
	  }
      }
  }

