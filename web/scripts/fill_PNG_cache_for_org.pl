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

my (@orgs, $version,$iw, $min_zoomi, $max_chr_length, $min_chr_length, @skip_oids, $max_zoomi, $help, $overwrite, $unprocess_ds, $zoom_step, $unprocess_zoom, @datasets, $start_pos, $stop_pos, $base_dir);

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
	   "base_dir" =>\$base_dir,
	   );

help() if $help;

$BASE_DIR = $base_dir if $base_dir;
$start_pos = 0 unless defined $start_pos;
$overwrite = 1 unless defined $overwrite;
$min_zoomi = 5 unless defined $min_zoomi;
$zoom_step = 4 unless defined $zoom_step;
$max_chr_length = 0 unless $max_chr_length;
my %skip_oids = map {$_,1} @skip_oids;

$iw = 256 unless $iw;
my $db = new CoGe::Genome;
my $graphic = new CoGe::Graphics;
$graphic->DEBUG(0);
##some limits to keep from blowing our stack
use vars qw($MAX_FEATURES $MAX_NT);
$MAX_FEATURES = $iw*2; #one feature per pixel
$MAX_NT       = $iw*100000; #50,000 nt per pixel

$graphic->MAX_FEATURES($iw*100);
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
    print "Generating images for region ($start_pos - $stop_pos)\n";
    next unless defined $max_zoomi;
    print "Processing zoom levels: $max_zoomi - $min_zoomi, zoom-step: $zoom_step\n";
    for (my $max_zoom =$max_zoomi; $max_zoom >= $min_zoomi; $max_zoom = $max_zoom - $zoom_step-1)
      {
	my $min_zoom = $max_zoom-$zoom_step;
	$min_zoom = 0 if $min_zoom < 0;
	$min_zoom = $min_zoomi if $min_zoom < $min_zoomi;
	next if $max_zoom < $min_zoomi;
	my $max_chars = 10 * 2**$max_zoom;
	$max_chars = $chr_len if $max_chars > $chr_len;
	#march through windows on the chromosome of size max_chars
	foreach (my $c_start=$start_pos; $c_start <= $stop_pos; $c_start+=$max_chars)
	  {
	    my $c_stop = $c_start+$max_chars-1;
	    print "________________________________________________________________\n";
	    print "Working on chromosome region $c_start - ", $c_stop,"\n";
	    print "----------------------------------------------------------------\n\n";
	    
	    my $seq = uc($db->get_genomic_sequence_obj->get_sequence(start=>$c_start, end=>$c_stop, chr=>$chr, dataset=>$ds)); 

	    unless ($seq)
	      {
		die "no nucleotide sequence for start=>$c_start, end=>$c_start+$max_chars, chr=>$chr, dataset=>$ds\n";
	      }
	    my $seq_len = length $seq;
	    my $gfx = new CoGe::Graphics::Chromosome;
#	    $gfx->DEBUG(1);
	    $gfx->overlap_adjustment(0);
	    $gfx->benchmark(1);
	    $gfx->skip_duplicate_features(1);
	    $graphic->initialize_c(
				   ds => $ds->id,
				   chr => $chr,
				   iw => $iw,
				   start => $c_start,
				   stop => $c_start+$max_chars,
				   db=>$db,
				   z=>$max_zoom, 
				   c=>$gfx,
				   csh=>200,
				   cmh=>5,
				   fsh=>10,
				   fmh=>2,
				  );
	    my @cds_feats;
	    my $taa = new Benchmark;
	    foreach my $ds2 ($ds->get_associated_datasets)
	      {
		my $cds_feats = $graphic->process_features(start=>$c_start, stop=>$c_start+$max_chars-1, chr=>$chr, ds=>$ds2->id, db=>$db,c=>$gfx);

		push @cds_feats, @$cds_feats;
	      }
	    my $tbb = new Benchmark;		
	    my $ft_time = timestr(timediff($tbb, $taa));
	    print qq{        **Time to process features           $ft_time
};
	    #go through all the zoom levels for this region
	    my $aa_length = 0;
	    foreach (my $z=$max_zoom; $z >= $min_zoom; $z--)# (0..$max_zoom..0)
	      {
#		print "BASE_DIR: $BASE_DIR/ds__".$ds->id."/iw__$iw/z__".$z."\n";
		
		if ($unprocess_zoom && -d "$BASE_DIR/ds__".$ds->id."/chr__$chr/iw__$iw/z__".$z)
		  {
		    print "Skipping zoom level $z because image directory exists and unprocessed_zoom flag is set to true.\n";
		    next;
		  }

		my $chars = 10 * 2**$z;
		my $tot = ceil($max_chars/$chars);
		print "     Zoom level $z: ($chars characters per tile).  ";
		print "Total number of images to be generated: ", $tot,"\n";
		my $ta = new Benchmark;
		unless ($aa_length == 1)
		  {
		    $gfx->delete_features('aa');
		    foreach my $feat (@cds_feats)
		      {
			next if $feat->stop < $c_start || $feat->start > $c_stop;
			$aa_length =$graphic->process_proteins(genomic_feat=>$feat, c=>$gfx);
		      }
		  }
		my $tb = new Benchmark;		
		my $pr_time = timestr(timediff($tb, $ta));
		print qq{        **Time to process proteins           $pr_time
};
		my $img_count = 0;
		my $seq_pos = 0;
		#go through each window at this zoom level for this chromosomal window
		my $stop = $c_start+$max_chars > $chr_len ? $chr_len : $c_start+$max_chars;
		foreach (my $i=$c_start; $i< $stop; $i+=$chars)
		  {
		    my $j = $i+$chars-1;
		    $img_count++;    
		    my $clen;
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
		    print "\tImage $img_count of $tot";
		    print " ($i-$j), image length: $clen, chars: $chars, seq_len: $seq_len, length subseq: ", length $subseq,"\n";
#		    print "\tstart: ",$i,", end: ", $i+$clen-1,", length: ", length $subseq,"\n";
		    #print "\tlength subseq: ", length $subseq,"\n";
		    #print "\t(range: ".($i)."-".($i+$chars-1).")\n\n";
		    generate_image(db=>$db, c=>$gfx, chr=>$chr, ds=>$ds, iw=>$iw, z=>$z, x=>$i, stop=>$j, seq=>$subseq);
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
    print "\tCache file: $cache_file.";
    print " File exists! Overwrite flag: $overwrite" if -r $cache_file;
    print "\n";
    if (-r $cache_file && $overwrite == 0)
			  {
			    print " Skipping generating image: file exists and overwrite is turned off.\n";
			    return;
			  }
    
    my ($ta, $tb);
    if ($x < 1)
      { #let's make the end of the chromosome
	$ta = new Benchmark;
	my $start = 0-(10*2**$z);
	$c->set_region(start=>$start, stop=>-1);
	$c->_gd(0);
	$c->generate_region;
	my $img = $c->gd->png;
	my $file = get_file_name(chr=>$chr,ds=>$ds, iw=>$iw, z=>$z, x=>$start);
	open( IMG, ">$file") or die "cant open $file: $!\n"; 
	binmode IMG;
	print IMG $img; 
	close(IMG);
	`chmod 666 $file`;
	$tb = new Benchmark;
      }
    $c->delete_features('nt'); #clear the old nucleotides before processing anew
    my $t0 = new Benchmark;
    $c->set_region(start=>$x, stop=>$stop);
    $graphic->process_nucleotides(start=>$x, stop=>$stop, chr=>$chr, ds=>$ds->id, db=>$db, c=>$c, seq=>$seq);
    
    my $t1 = new Benchmark;
#    my ($s, $e) = ($c->_region_start, $c->_region_stop);
    $c->_gd(0);
    $c->generate_region();
    my $feats = $c->_features;
    my $t2 = new Benchmark;
    my $img = $c->gd->png;
    my $t3 = new Benchmark;
    #local( *IMG,$/ );
    open( IMG, ">$cache_file") or die "cant open $cache_file: $!\n"; 
    binmode IMG;
    print IMG $img; 
    close(IMG);
    
    my $t4 = new Benchmark;
    my $end_time = timestr(timediff($tb, $ta)) if ($ta && $tb);
    my $nt_time = timestr(timediff($t1, $t0));
    my $gr_time = timestr(timediff($t2, $t1));
    my $gi_time = timestr(timediff($t3, $t2));
    my $si_time = timestr(timediff($t4, $t3));
    print qq{        **Time to generate chr_end image:    $end_time
} if $end_time;
    print qq{
        **Time to process nucleotides:       $nt_time
        **Time to generate region:           $gr_time
        **Time to generate image:            $gi_time
        **Time to store image:               $si_time
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

 base_dir   => the directory in which the images will be deposited.


};
    exit;
  }
