#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Graphics::Feature::GAGA;
use CoGe::Graphics::Feature::Exon_motifs;
use CoGe::Graphics::Feature::AminoAcid;
use CoGe::Graphics::Feature::Domain;
use CoGe::Graphics::Feature::Block;
use CoGe::Genome;
use CoGe::Graphics;
use File::Spec::Functions;
#use CoGe::Accessory::Tile::Cache;
use Getopt::Long;
use POSIX;
use Benchmark;

use vars qw($BASE_DIR);
$BASE_DIR = "/opt/apache/CoGe/_cache_";

my (@orgs, $version,$iw, $min_zoomi, $max_chr_length, $min_chr_length, @skip_oids, $max_zoomi, $help, $overwrite, $unprocess_ds, $zoom_step, $unprocess_zoom, @datasets, $start_pos, $stop_pos);

GetOptions(
	   "org|oid|o=s"=>\@orgs,
	   "dataset|ds=s"=>\@datasets,
	   "v=s" => \$version,
	   "iw=s"=>\$iw,
	   "min_zoom|minz=s"=>\$min_zoomi,
	   "max_zoom|maxz=s"=>\$max_zoomi,
	   "max_chr_length=s"=>\$max_chr_length,
	   "min_chr_length=s"=>\$min_chr_length,
	   "skip_oid=s"=>\@skip_oids,
	   "unprocessed_orgs|uo|unprocessed_datasets|ud"=>\$unprocess_ds,
	   "unprocessed_zoom|uz"=>\$unprocess_zoom,
	   "h|help"=>\$help,
	   "overwrite|ow=s"=>\$overwrite,
	   "zoom_step|zs=s"=>\$zoom_step,
	   "start_pos=i"=>\$start_pos,
	   "stop_pos=i"=>\$stop_pos,
	   );

help() if $help;

$start_pos = 0 unless defined $start_pos;
$overwrite = 1 unless defined $overwrite;
$min_zoomi = 5 unless $min_zoomi;
$zoom_step = 4 unless defined $zoom_step;
$max_chr_length = 0 unless $max_chr_length;
my %skip_oids = map {$_,1} @skip_oids;

$iw = 256 unless $iw;
my $db = new CoGe::Genome;
my $graphic = new CoGe::Graphics;

##some limits to keep from blowing our stack
use vars qw($MAX_FEATURES $MAX_NT);
$MAX_FEATURES = $iw*2; #one feature per pixel
$MAX_NT       = $iw*100000; #50,000 nt per pixel

$graphic->MAX_FEATURES($iw*2);
$graphic->MAX_NT($iw*100000);
my $orgs = get_orgs(orgs=>\@orgs, datasets=>\@datasets);

foreach my $org (@$orgs)
  {
    next if $skip_oids{$org};
    print "working on: ".$org->name,"\n";
    foreach my $ds ($org->dataset)
      {
	process_dataset(ds=>$ds);
      }
  }

foreach my $ds (@datasets)
  {
    process_dataset(ds=>$ds);
  }

sub process_dataset
  {
    my %opts = @_;
    my $ds = $db->get_dataset_obj->resolve_dataset($opts{ds});
    
    next if $version && $ds->version ne $version;
    my $dsname = $ds->name;
    $dsname .= ": ".$ds->desc if $ds->desc;
    print "Processing dataset $dsname. . .\n";
    if ( $unprocess_ds && -d "$BASE_DIR/ds__".$ds->id )
      {
	print "Skipping dataset because image directory exists and we are only working on unprocessed datasets.\n";
	next;
      }
    my ($max_zoom_tmp, $chr_len, $chr) = find_max_z(ds=>$ds);
    $max_zoomi = $max_zoom_tmp unless $max_zoomi;
    next unless $chr_len;
    $stop_pos = $chr_len unless defined $stop_pos;

    next if ($max_chr_length && $chr_len > $max_chr_length);
    next if ($min_chr_length && $chr_len < $min_chr_length);
    print "Total number of bp for chr $chr: $chr_len\n";
    next unless defined $max_zoomi;
    for (my $max_zoom =$max_zoomi; $max_zoom >= $min_zoomi; $max_zoom = $max_zoom - $zoom_step)
      {
	my $min_zoom = $max_zoom-$zoom_step+1;
	next if $max_zoom < $min_zoomi;
	my $max_chars = 10 * 2**$max_zoom;
	$max_chars = $chr_len if $max_chars > $chr_len;
	#march through windows on the chromosome of size max_chars
	foreach (my $c_start=$start_pos; $c_start < $stop_pos; $c_start+=$max_chars)
	  {
	    print "Working on chromosome chunk $c_start - ", $c_start+$max_chars,"\n";
	    my $seq = uc($db->get_genomic_sequence_obj->get_sequence(start=>$c_start, end=>$c_start+$max_chars, chr=>$chr, dataset=>$ds)); 
	    unless ($seq)
	      {
		die "no nucleotide sequence for start=>$c_start, end=>$c_start+$max_chars, chr=>$chr, dataset=>$ds\n";
	      }
	    my $seq_len = length $seq;
	    my $c = new CoGe::Graphics::Chromosome;
#	    $c->DEBUG(1);
	    $c->overlap_adjustment(0);
	    $graphic->initialize_c(
				   ds => $ds->id,
				   chr => $chr,
				   iw => $iw,
				   start => $c_start,
				   stop => $c_start+$max_chars,
				   db=>$db,
				   z=>$max_zoom, 
				   c=>$c,
				   csh=>200,
				   cmh=>5,
				   fsh=>10,
				   fmh=>2,
				  );
# 	    my $c = initialize_c(
# 				 ds => $ds->id,
# 				 chr => $chr,
# 				 iw => $iw,
# 				 start => $c_start,
# 				 stop => $c_start+$max_chars,
# 				 db=>$db,
# 				 z=>$max_zoom,
# 				);
	    my @cds_feats;
	    foreach my $ds2 ($ds->get_associated_datasets)
	      {
		#	    next;
#		my $cds_feats = process_features(start=>$c_start, stop=>$c_start+$max_chars-1, chr=>$chr, ds=>$ds2->id, db=>$db,c=>$c);
		my $cds_feats = $graphic->process_features(start=>$c_start, stop=>$c_start+$max_chars-1, chr=>$chr, ds=>$ds2->id, db=>$db,c=>$c);
		push @cds_feats, @$cds_feats;
	      }
	    #go through all the zoom levels for this region
	    foreach (my $z=$max_zoom; $z >= $min_zoom; $z--)# (0..$max_zoom..0)
	      {
		print "BASE_DIR: $BASE_DIR/ds__".$ds->id."/iw__$iw/z__".$z."\n";
		
		if ($unprocess_zoom && -d "$BASE_DIR/ds__".$ds->id."/chr__$chr/iw__$iw/z__".$z)
		  {
		    print "Skipping zoom level $z because image directory exists and unprocessed_zoom flag is set to true.\n";
		    next;
		  }

		my $chars = 10 * 2**$z;
		my $tot = ceil($max_chars/$chars);
		print "$chars characters per tile at zoom level $z\n";
		print "Total number of images to be generated: ", $tot,"\n";
		$c->delete_features('aa');
		my $ta = new Benchmark;
		foreach my $feat (@cds_feats)
		  {
		    $graphic->draw_prots(genomic_feat=>$feat, c=>$c);
		  }
		my $tb = new Benchmark;
		my $count = 0;		
		my $pr_time = timestr(timediff($tb, $ta));
		print qq{
Time to process proteins           $pr_time
};
		my $seq_pos = 0;
		#go through each window at this zoom level for this chromosomal window
		my $stop = $c_start+$max_chars > $chr_len ? $chr_len : $c_start+$max_chars;
		foreach (my $i=$c_start; $i< $stop; $i+=$chars)
		  {
		    $count++;    
		    my $clen = $chars;
		    if ($i+$chars > $c_start+$max_chars)
		      {
			$clen = ($max_chars) % $chars;
		      }
		    elsif ($i == 0)
		      {
			$clen = $chars-1;
		      }
		    else
		      {
			$clen = $chars;
		      }
		    my $subseq = substr($seq, $seq_pos, $clen);
		    
		    $seq_pos+=$clen;
		    
#		    print "$i - $stop\n";
		    print "Image $count of $tot";
		    print "\tclen: $clen, chars: $chars, max_chars: $max_chars, seq_len: $seq_len\n";
		    print "\tstart: ",$i,"($seq_pos), end: ", $i+$clen-1,", length: ", length $subseq,"\n";
		    print "\tlength subseq: ", length $subseq,"\n";
		    print "\t(range: ".($i)."-".($i+$chars-1).")\n\n";
		    generate_image(db=>$db, c=>$c, chr=>$chr, ds=>$ds, iw=>$iw, z=>$z, x=>$i, stop=>$i+$chars-1, seq=>$subseq);
		  }
	      }
	  }
      }
  }

sub generate_image
  {
    my %opts = @_;
    my $c = $opts{c};
    my $chr = $opts{chr};
    my $ds = $opts{ds};
    my $iw = $opts{iw};
    my $z = $opts{z};
    my $x = $opts{x};# || $opts{start};
    my $db = $opts{db};
    my $stop = $opts{stop};

    my $seq = $opts{seq};
    my $cache_file = get_file_name(chr=>$chr,ds=>$ds, iw=>$iw, z=>$z, x=>$x);    
    print "Cache file: $cache_file\n";
    print "File exists! Overwrite flag set to $overwrite\n" if -r $cache_file;
    if (-r $cache_file && $overwrite == 0)
			  {
			    print " Skipping generating image: file exists and overwrite is turned off.\n";
			    return;
			  }
    my $t0 = new Benchmark;
    $c->set_region(start=>$x, stop=>$stop);
    $graphic->process_nucleotides(start=>$x, stop=>$stop, chr=>$chr, ds=>$ds->id, db=>$db, c=>$c, seq=>$seq);
    my $t1 = new Benchmark;
#    my ($s, $e) = ($c->_region_start, $c->_region_stop);
    $c->_gd(0);
    $c->generate_region();
    my $t2 = new Benchmark;
    my $img = $c->gd->png;
    my $t3 = new Benchmark;
    #local( *IMG,$/ );
    open( IMG, ">$cache_file") or die "cant open $cache_file: $!\n"; 
    binmode IMG;
    print IMG $img; 
    close(IMG);
    
    my $t4 = new Benchmark;
    my $nt_time = timestr(timediff($t1, $t0));
    my $gr_time = timestr(timediff($t2, $t1));
    my $gi_time = timestr(timediff($t3, $t2));
    my $si_time = timestr(timediff($t4, $t3));
    print qq{
Time to process nucleotides:       $nt_time
Time to generate region:           $gr_time
Time to generate image:            $gi_time
Time to store image:               $si_time

 ----================================------

};
    if (-r $cache_file)
      {
	#`touch $cache_file`;
	`chmod 666 $cache_file`;
      }
  }

sub find_max_z
  {
    my %opts = @_;
    my $ds = $opts{ds};
    my $go = $ds->genomic_sequences->next;
    return unless $go;
    my $chr_len = $go->get_last_position(ds=>$ds);
    return unless $chr_len;
    my $max_zoom = ceil (log10($chr_len/10)/(log10(2)));
#    print "Max zoom: $max_zoom\t";
    my $feat_count = $db->get_feature_obj->count_features_in_region(start=>1, stop=>$chr_len, dataset=>$ds->id, chr=>$go->chr);
    #assume features are evenly spread across chromosome (I know this is not correct, but. . .)
    while ($feat_count > $MAX_FEATURES)
      {
	$max_zoom--;
	$feat_count /= 2;
      }
#    print "Returned zoom: $max_zoom\n";
    return ($max_zoom, $chr_len, $go->chr);
  }

sub initialize_c
  {
    my %opts = @_;
    my $ds = $opts{ds};
    my $chr = $opts{chr};
    my $iw = $opts{iw};
    my $z = $opts{z};
    my $mag = $opts{mag};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $db = $opts{db};
    my $csh=$opts{csh} || 200;
    my $cmh=$opts{cmh} || 5;
    my $fsh=$opts{fsh} || 10;
    my $fmh=$opts{fmh} || 2;
    my $start_pict = $opts{'start_pict'} || 'left';
    my $c = new CoGe::Graphics::Chromosome;

    my ($gen_seq) = $db->get_genomic_seq_obj->search({dataset_id=>$ds});
    return unless $gen_seq && $gen_seq->chr eq $chr;
    my $chr_length =
    $db->get_genomic_sequence_obj->get_last_position(ds=>$ds);
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
    $stop = $c->_region_stop;
    #let's add the max top and bottom tracks to the image to keep it constant
    my $f1= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => 1, image_width=>5, image_height=>5});
    $f1->merge_percent(0);
    $c->add_feature($f1);
    my $f2= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => -1,image_width=>5, image_height=>5});
    $f2->merge_percent(0);
    $c->add_feature($f2);
    return ($c);
}
sub process_nucleotides
  {
    my %opts = @_;
    my $start = $opts{start};
    my $chr = $opts{chr};
    my $ds = $opts{ds};
    my $db = $opts{db};
    my $c = $opts{c};
    my $seq = $opts{seq};
    #process nucleotides
    $c->delete_features('fill');
    my $seq_len = length $seq;
    my $chrs = int (($c->find_size_for_magnification)/$c->iw);
#    print "\tst:", $c->_region_stop, "-", $c->_region_start,"\n";
#    print "c start - stop:  ".$c->_region_start." - ".$c->_region_stop."  iw: ".$c->iw."\n";
#    print "\tNT chars: $chrs\n";
    $chrs = 1 if $chrs < 1;
    my $pos = 0;
    $start = 1 if $start < 1;
    my $count = 0;
    while ($pos < $seq_len)
      {
	$count++;
        my $subseq = substr ($seq, $pos, $chrs);

#	print "$count: $subseq\tsl: $seq_len\t ssl: ". length $subseq,"\n";
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
    my $ds = $opts{ds};
    my $db = $opts{db};
    my $c = $opts{c};
    my $accn = $opts{accn};
    my $print_names = $opts{print_names};
    my $size = $stop-$start;
    my $feat_count = $db->get_feature_obj->count_features_in_region(start=>$start, end=>$stop, dataset=>$ds, chr=>$chr);
    my @cds_feats;  #place to hold CDS features with protein sequences for later generation of protein sequence images;
    my $count = 0;
    foreach my $feat ($db->get_feature_obj->get_features_in_region(start=>$start-$size, end=>$stop+$size, dataset=>$ds, chr=>$chr))
      {
#	print $feat->genbank_location_string,"\n";
	$count++;
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
	    $f->overlay(1);
	    $f->mag(0.5);
          }
        elsif ($feat->type->name =~ /CDS/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,255,0, 50]);
        	$f->order(1);
		$f->overlay(3);
		push @cds_feats, $feat;
		if ($accn)
		  {
		    foreach my $name (@{$feat->{QUALIFIERS}{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
		      }
		  }
          }
        elsif ($feat->type->name =~ /mrna/i)
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
        elsif ($feat->type->name =~ /CNS/i)
          {
            $f = CoGe::Graphics::Feature::Block->new();
            $f->order(1);
            my $color = [255, 100, 255];
            $f->color($color);
          }
	else
	  {
	    print "Skipping feature of type: ",$feat->type->name,"\n";
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
    return \@cds_feats;
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
    #Do we have any protein sequence we can use?
    foreach my $seq ($feat->sequences)
      {
	next unless $seq->seq_type->name =~ /prot/i;
	my ($pseq) = $seq->sequence_data;
#	print "\tAdding protein sequence of length: ",length($pseq),". . .";

	my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw)/3;
	$chrs = 1 if $chrs < 1;
#	my $chrs = 1;
	my $pos = 0;
	while ($pos <= length $pseq)
	  {
	    my $aseq = substr($pseq, $pos, $chrs);
	    foreach my $loc ($seq->get_genomic_locations(start=>$pos+1, stop=>$pos+$chrs))
	      {
		my $ao = CoGe::Graphics::Feature::AminoAcid->new({aa=>$aseq, start=>$loc->start, stop=>$loc->stop, strand => $loc->strand, order=>2});
		$ao->skip_overlap_search(1);
		$ao->type('aa');
		$ao->mag(0.75);
		$c->add_feature($ao);
		delete $loc->{__Changed}; #silence the warning from Class::DBI
	      }
	    
	    $pos+=$chrs;
	  }
#	print "Done!\n";
      }
  }

sub get_file_name
    {
      my %opts = @_;
      my ($basedir,$dir) = get_dir_array(%opts);
      my $reldir = catfile(@$dir) . '.png';
      my $basepath = catfile(@$basedir);

      my $fn = catfile($basepath,@$dir) . '.png' ;
      ($fn) = $fn =~ /(.*)/;
      my $data;
      if(! -e $fn){
       make_dir($basepath,$dir);
       }
    return $fn;
    }

sub make_dir {
    my ($dirstr,$dir) = @_;
    pop @$dir;
    $dirstr = catfile($dirstr);
    umask 0;
    foreach my $subdir (@$dir){
      ($dirstr) = catdir($dirstr,$subdir) =~ /(.*)/;
      my $status = mkdir($dirstr,0777) if !-e $dirstr;
    }
  }
                                       
sub get_max {
    my $z = shift;
    my $MAX_FILES_PER_DIR = 1000;
    my $upt = 10 * 2 ** $z;
    return int(log($MAX_FILES_PER_DIR * $upt)/log(10));
 }

sub grep_match {
    my ($needle,$haystack) = @_;
    return grep { $needle eq $_ } @$haystack;
  }


sub get_dir_array {
    my %query_pairs = @_;
    my $SPATIAL_PARS = ['x'];
    my $ZOOM_PARS = 'z';
    #my $base_dir = "/opt/apache/CoGe";#$ENV{DOCUMENT_ROOT};
    my @basedir = split(/[\/\\]/,$BASE_DIR);
    my @dir;# = '_cache_';
    my @keyvals = qw(ds chr iw z x);
    my $MAX = get_max($query_pairs{$ZOOM_PARS});
    foreach my $key(@keyvals ){
         next if  !$key;
         my $val = $query_pairs{$key};
         my @vals;
         # 123456 becomes /x__123/456/  if $MAX == 3
         # 123456 becomes /x__123456/   if $MAX  > 5
         if(grep_match($key,$SPATIAL_PARS)){
             if($val !~ /(\D|-)+/){
	       $val =~ s/-/n/g;
	       $val = scalar reverse($val);
	       while($val =~  /(\d{1,$MAX}n?)/g){
		 unshift(@vals,scalar reverse($1));
	       }
	     }
	   }
	 @vals = ($val) unless @vals;
	 my $first = shift @vals;
	 push(@dir,$key . '__' . $first);
	 map { push(@dir,$_) } @vals;
       }
    return (\@basedir,\@dir); 
 }

sub get_orgs
  {
    my %opts = @_;
    my $orgs = $opts{orgs};
    my $ds = $opts{datasets};
    my @out_orgs;
    foreach my $org (@orgs)
      {
	my $o = $db->get_org_obj->resolve_organism($org);
	push @out_orgs, $o if $o;
      }
    unless (@out_orgs || $ds)
      {
	@out_orgs = $db->get_org_obj->retrieve_all();
      }
    return \@out_orgs;

  }

sub help
  {
    print qq{
Welcome to $0

The program fills the CoGe system's image cache for chromosome images used by GeLo.pl when viewing chromosomal regions.  

Options:

 org       =>  organism name or database ids that will be used for generating images
               if none are specified, the program will go through all the organisms in the database
               More than one of these may be specified

 dataset   =>  dataset name or database id that will be used for generating images.  This can be
               used as an alternative to specifying an organism

 v         =>  version of organism/dataset to use.  If not is specified, then all will be used.

 iw        => width of the images to generate (default 256)

 min_zoom  => the smallest zoom for which images will be generated (default: 6)

 max_zoom  => the largest zoom for which images will be generated (default is determined by the
              number of features present at various zoom levels.

 start_pos => the starting position from which to begin generating images (default: 0)

 stop_pos  => the stopping position from which to end generating images (default: end position of dataset)

 max_chr_length => skip processing genomes whose size is larger than this number (default: no limit)

 min_chr_length => skip processing genomes whose size is smaller than this number (default: no limit)

 skip_oid  => organisms to skip based on their db id.

 overwrite | ow => overwrites existing image (default: 1);

 unprocessed_orgs | uo | unprocessed_datasets | ud   =>  skip any dataset that already has associated images.
 
 unprocessed_zoom | uz    => skip generating images if there are already images for the requested zoom level.

 zoom_step | zs   => the number of zoom step to process at once.  This will cause the virual chromosome to 
                     be flushed periodically and free up memory.  This is needed to keep from running out of
                     memory in some cases.  DEFAULT: 4


};
    exit;
  }
