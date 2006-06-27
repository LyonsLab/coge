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
use File::Spec::Functions;
#use CoGe::Accessory::Tile::Cache;
use Getopt::Long;
use POSIX;
use Benchmark;


my (@org_ids, $version,$iw, $min_zoom);

GetOptions("oid=s"=>\@org_ids,
	   "v=s" => \$version,
	   "iw=s"=>\$iw,
	   "min_zoom|mz=s"=>\$min_zoom,
#	   "z|zoom=s" => \$zoom,
#	   "u|url=s" => \$url
	   );
$min_zoom = 0 unless $min_zoom;

#$url = "synteny.cnr.berkeley.edu" unless $url;
#$url = "http://" .$url unless $url =~ /http:\/\//;
#$url .= "/CoGe/tiler.pl";
$iw = 256 unless $iw;
my $db = new CoGe::Genome;

##some limits to keep from blowing our stack
use vars qw($MAX_FEATURES $MAX_NT);
$MAX_FEATURES = $iw*5; #one feature per pixel
$MAX_NT       = $iw*100000; #50,000 nt per pixel

unless (@org_ids)
  {
    @org_ids = $db->get_org_obj->retrieve_all();
  }
foreach my $org_id (@org_ids)
  {
    my ($org) = ref ($org_id) =~ /org/i ? $org_id :$db->get_org_obj->retrieve($org_id);
    print "working on: ".$org->name,"\n";
    foreach my $di ($org->data_information)
      {
	next if $version && $di->version ne $version;
	$version = $di->version unless $version;
	my ($max_zoom, $chr_len, $chr) = find_max_z(di=>$di);
	print "Total number of bp for chr $chr: $chr_len\n";
	next unless defined $max_zoom;
	#    $max_zoom=10;
	my $max_chars = 10 * 2**$max_zoom;
	#march through windows on the chromosome of size max_chars
	foreach (my $c_start=0; $c_start <= $chr_len; $c_start+=$max_chars)
	  {
	    print "Working on chromosome chunk $c_start - ", $c_start+$max_chars,"\n";
	    my $seq = uc($db->get_genomic_sequence_obj->get_sequence(start=>$c_start, end=>$c_start+$max_chars, chr=>$chr, info_id=>$di)); 
	    my $seq_len = length $seq;
	    my $c = initialize_c(
				 di => $di->id,
				 chr => $chr,
				 iw => $iw,
				 start => $c_start,
				 stop => $c_start+$max_chars,
				 db=>$db,
				 z=>$max_zoom,
				);
	    my @cds_feats;
	    foreach my $di2 ($di->get_associated_data_infos)
	      {
		#	    next;
		my $cds_feats = process_features(start=>$c_start, stop=>$c_start+$max_chars-1, chr=>$chr, di=>$di2->id, db=>$db,c=>$c);
		push @cds_feats, @$cds_feats;
	      }
	    #go through all the zoom levels for this region
	    foreach (my $z=$max_zoom; $z >= $min_zoom; $z--)# (0..$max_zoom..0)
#	    foreach (my $z=0; $z >= 0; $z--)# (0..$max_zoom..0)
	      {
		my $chars = 10 * 2**$z;
		my $tot = ceil ($max_chars/$chars);
		print "$chars characters per tile at zoom level $z\n";
		print "Total number of images to be generated: ", $tot,"\n";
		#	next;
		$c->delete_features('aa');
		my $ta = new Benchmark;
		foreach my $feat (@cds_feats)
		  {
		    draw_prots(genomic_feat=>$feat, c=>$c);
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
		    my $t0 = new Benchmark;
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
		    my $cache_file = get_file_name(chr=>$chr,di=>$di, iw=>$iw, z=>$z, x=>$i);
		    print "Cache file: $cache_file\n";
		    print "\tclen: $clen, chars: $chars, max_chars: $max_chars, seq_len: $seq_len\n";
		    print "\tstart: ",$i,"($seq_pos), end: ", $i+$clen-1,", length: ", length $subseq,"\n\n";
		    print "\tlength subseq: ", length $subseq,"\n";

 		    $seq_pos+=$clen;
 		    #next;

		    $c->set_region(start=>$i, stop=>$i+$chars-1);
		    process_nucleotides(start=>$c_start+$i, chr=>$chr, di=>$di->id, db=>$db, c=>$c, seq=>$subseq);

		    my $t1 = new Benchmark;

		    my ($s, $e) = ($c->_region_start, $c->_region_stop);
		    print "Image $count of $tot";
		    print "\t(range: ".($i)."-".($i+$chars-1).")\n";
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

};
		    if (-r $cache_file)
		      {
			#`touch $cache_file`;
			`chmod 666 $cache_file`;
		      }
		  }
	      }
	  }
      }
  }

sub initialize_db
  {
  }

sub find_max_z
  {
    my %opts = @_;
    my $di = $opts{di};
    my $go = $di->genomic_sequences->next;
    return unless $go;
    my $chr_len = $go->get_last_position(di=>$di);
    return unless $chr_len;
    my $max_zoom = ceil (log10($chr_len/10)/(log10(2)));
#    print "Max zoom: $max_zoom\t";
    my $feat_count = $db->get_feature_obj->count_features_in_region(start=>1, stop=>$chr_len, info_id=>$di->id, chr=>$go->chr);
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
    my $di = $opts{di};
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

    my ($gen_seq) = $db->get_genomic_seq_obj->search({data_information_id=>$di});
    return unless $gen_seq && $gen_seq->chr eq $chr;
    my $chr_length =
    $db->get_genomic_sequence_obj->get_last_position(di=>$di);
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
    my $di = $opts{di};
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
    my $di = $opts{di};
    my $db = $opts{db};
    my $c = $opts{c};
    my $accn = $opts{accn};
    my $print_names = $opts{print_names};
    my $feat_count = $db->get_feature_obj->count_features_in_region(start=>$start, end=>$stop, info_id=>$di, chr=>$chr);
    my @cds_feats;  #place to hold CDS features with protein sequences for later generation of protein sequence images;
#    if ($MAX_FEATURES && $feat_count > $MAX_FEATURES)
#      {
#	warn "exceeded maximum number of features $MAX_FEATURES. ($feat_count requested)\nskipping.\n";
#	return;
#      }
    my $count = 0;
    foreach my $feat ($db->get_feature_obj->get_features_in_region(start=>$start, end=>$stop, info_id=>$di, chr=>$chr))
      {
	$count++;
	print "processing feature $count / $feat_count.  feature type ",$feat->type->name,"\n";
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
#		draw_prots(genomic_feat=>$feat, c=>$c);
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
	    print STDERR "Skipping feature of type: ",$feat->type->name,"\n";
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
	print "\tAdding protein sequence of length: ",length($pseq),". . .";

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
		$c->add_feature($ao);
		delete $loc->{__Changed}; #silence the warning from Class::DBI
	      }
	    
	    $pos+=$chrs;
	  }
	print "Done!\n";
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
      $dirstr = catfile($dirstr,shift @$dir);
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
#    my $qstr = $ENV{QUERY_STRING};
    my %query_pairs = @_;
    my $SPATIAL_PARS = ['x'];
    my $ZOOM_PARS = 'z';
    my $base_dir = "/opt/apache/CoGe";#$ENV{DOCUMENT_ROOT};
    my @basedir = split(/[\/\\]/,$base_dir);
    my @dir = '_cache_';;# = split(/[\/\\]/,$ENV{SCRIPT_NAME});
#    shift @dir;
#    $dir[$#dir] = '_cache_';
#    my @keyvals = map { split('=', $_) } split(/&/,$qstr);
#    my %query_pairs = @keyvals;
#    @keyvals = map { $_ % 2 == 0 && !grep_match($keyvals[$_],$SPATIAL_PARS)? $keyvals[$_] : '' } 0.. $#keyvals;
#    @keyvals = sort @keyvals;
    my @keyvals = qw(di chr iw z x);
#    push(@keyvals,@$SPATIAL_PARS);
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
                                           
