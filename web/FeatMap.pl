#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
use URI::Escape;
use CoGeX;
use CoGeX::Feature;
use POSIX;
use File::Temp;
use File::Basename;
use CoGe::Accessory::Restricted_orgs;
use CoGe::Graphics::GenomeView;
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $coge $cogeweb );


$TEMPDIR = "/opt/apache/CoGe/tmp/FeatMap";
$TEMPURL = "/CoGe/tmp/FeatMap";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;

$coge = CoGeX->dbconnect();

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       generate_feat_info=>\&generate_feat_info,
		      );
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();


sub gen_html
  {
    my $html;
    unless ($USER)
      {
	$html = login();
      }
    else
      {
	my ($body) = gen_body();
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	$template->param(TITLE=>'Feature Map');
	$template->param(HELP=>'FeatMap');
	# print STDERR "user is: ",$USER,"\n";
	my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);
	$template->param(LOGO_PNG=>"FeatMap-logo.png");
	$template->param(LOGON=>1) unless $USER->user_name eq "public";
	$template->param(DATE=>$DATE);
	#$template->param(LOGO_PNG=>"CoGeBlast-logo.png");
	$template->param(BOX_NAME=>'CoGe: Feature Map');
	$template->param(BODY=>$body);
	$html .= $template->output;
      }
  }
  
  sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatMap.tmpl');
    my $form = $FORM;
    my $no_values;
    my $sort_by_type = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $feat_list = [];
    foreach my $item ($form->param('fid'))
      {
	push @$feat_list, $item if $item =~ /^\d+$/;
      }
    my $dsid = $form->param('dsid') if $form->param('dsid');
    my $anticodon = $form->param('anticodon') if $form->param('anticodon');
    my $chr = $form->param('chr') if $form->param('chr');
    my $clade = $form->param('clade') if $form->param('clade');
    my $ftid = $form->param('ftid') if $form->param('ftid');
    push @$feat_list, @{get_fids_from_dataset(dsid=>$dsid, ftid=>$ftid, chr=>$chr)} if $dsid;
    my $chromosome_data = generate_chromosome_images(feature_list=>$feat_list, ftid=>$ftid, anticodon=>$anticodon, clade=>$clade);
    my $table = generate_table_data(feature_list=>$feat_list, ftid=>$ftid);
    if ($chromosome_data)
      {
	$template->param(CHROMOSOME_LOOP=>$chromosome_data);
	$template->param(FEAT_TABLE=>$table);
	return $template->output;
      }
    else
      {
	return "No feature ids were specified.";
      }
  }
  
  sub generate_table_data
  {
  	my %opts = @_;
    my $feat_list = $opts{feature_list};
    #print STDERR Dumper \$feat_list;
    my $ftid = $opts{ftid};
    return unless @$feat_list;
  	my @table;
  	$feat_list = [map {$coge->resultset("Feature")->find($_)} @$feat_list];
    $feat_list = [sort {$a->organism->name cmp $b->organism->name || $a->type->name cmp $b->type->name || $a->chromosome cmp $b->chromosome|| $a->start <=> $b->start}@$feat_list];
    foreach my $feat(@$feat_list)
    {
      unless ($feat)
	{
#	  warn "feature id $featid failed to return a valid feature object\n";
	  next;
	}
	if ($ftid) 
	{
	  next unless $feat->type->id eq $ftid;
	}
      my ($name) = $feat->names;
      push @table,{
		   CHECKBOX=>$feat->id,
		   FEAT_NAME=>qq{<span class="link" title="Click for Feature Information" onclick="show_feature_info('}.$feat->id.qq{')">}.$name."</span>",
		   LOC=>$feat->start,
		   CHR=>$feat->chr,
		   ORG=>$feat->organism->name,
		  };
    }
    return \@table;
  }
  
  sub get_fids_from_dataset
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $ftid = $opts{ftid};
    my $chr = $opts{chr};
    my $search = {dataset_id=>$dsid};
    my $join={};
    if ($ftid)
      {
	$search->{feature_type_id}=$ftid;
      }
    if ($chr)
      {
	$search->{chromosome}=$chr;
      }
    my @ids = map{$_->id}$coge->resultset('Feature')->search($search);
    return \@ids;
  }

  sub generate_chromosome_images
  {
    my %opts = @_;
    my $feat_list = $opts{feature_list};
    #print STDERR Dumper \$feat_list;
    my $ftid = $opts{ftid};
    my $clade = $opts{clade};
    my $anticodon = $opts{anticodon};
    return unless @$feat_list;
    
    $cogeweb = initialize_basefile(prog=>"FeatMap");
    my $width = $opts{width} || 1200;
    my $image_filename = $cogeweb->basefile;
    my $height = ($width / 16) <= 64 ? ($width / 16) : 64;
    
    my $scale = $opts{scale} || 'linear';
    my %data;
    my $filename;
    my (@data,@no_data);
    $feat_list = [map {$coge->resultset("Feature")->find($_)} @$feat_list];
    $feat_list = [sort {$a->organism->name cmp $b->organism->name || $a->type->name cmp $b->type->name || $a->chromosome cmp $b->chromosome|| $a->start <=> $b->start}@$feat_list];
    foreach my $feat(@$feat_list)
      {
	unless ($feat)
	  {
	    #	  warn "feature id $featid failed to return a valid feature object\n";
	    next;
	  }
	if ($ftid) 
	  {
	    next unless $feat->type->id eq $ftid;
	  }
	
	my $org = $feat->organism->name;
	push @{$data{$org}{feats}},$feat;
	$filename = $image_filename."_*.png";
      }
	
	
    foreach my $org (sort keys %data)
      {
	#first, initialize graphic
	$org =~ s/\s+$//;
	$data{$org}{image} =new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>$width, chromosome_height=>$height}) unless $data{$org}{image};
	foreach my $feat (@{$data{$org}{feats}})
	  {
	    my ($chr) = $feat->chr;
	    
	    my ($dsid) = $feat->dataset->id;
	    #add chromosome to graphic
	    unless ($data{$org}{chr}{$chr})
	      {
		next unless $dsid;
		my $ds = $coge->resultset('Dataset')->find($dsid);
		my $last_pos = $ds->last_chromosome_position($chr);
		$data{$org}{image}->add_chromosome(name=>"Chr: $chr",
						   end=>$last_pos,
						  );
		$data{$org}{chr}{$chr}=1;
	      }
	    my $up = $feat->strand =~ /-/ ? 0 : 1;
	    my ($name) = $feat->names;
	    $data{$org}{image}->add_feature(name=>$feat->id,
					    start=>$feat->start,
					    stop=>$feat->stop,
					    chr=>"Chr: $chr",
					    imagemap=>qq/class="imagemaplink" title="/.$name.qq/ (/.$feat->type->name.qq/)" onclick="show_feature_info('/.$feat->id.qq/')"/,
					    up=>$up,
					    color=>color_feat_by_type(feat=>$feat,anticodon=>$anticodon,clade=>$clade),
					   );
	  }
      }
    my $count = 1;
    foreach my $org (sort keys %data)
    {
	if ($data{$org}{image})
	  {
	    my $image_file;
	    my $image_map;
	    unless ($data{$org}{skip})	    
	      {
		my $x;
		$image_file = $cogeweb->basefile."_$count.png";
		($x, $image_file) = check_taint($image_file);
		
		$data{$org}{image}->image_width($width);
		$data{$org}{image}->chromosome_height($height);
		$data{$org}{image}->show_count(1);
		$data{$org}{image}->generate_png(filename=>$image_file);
		$image_map = $data{$org}{image}->generate_imagemap(mapname=>$cogeweb->basefilename."_".$count);
		my $map_file = $cogeweb->basefile."_$count.map";
		($x, $map_file) = check_taint($map_file);
		open (MAP, ">$map_file");
		print MAP $image_map;
		close MAP;
	    }
	    else
	      {
		my $x;
		$image_file = $data{$org}{image}."_$count.png";
		$image_map = get_map($cogeweb->basefile."_$count.map");
		#print STDERR $image_map_large,"\n";
		($x, $image_file) = check_taint($image_file);
	      }
	    
	    $image_file =~ s/$TEMPDIR/$TEMPURL/;
	   
	    push @data,  {ORG_NAME=>$org,CHR_IMAGE=>"<img src=$image_file ismap usemap='".$cogeweb->basefilename."_"."$count' border=0>$image_map"};
	    $count++;
	  }
	else
	  {
	    push @no_data,  {ORG_NAME=>"No Hits: <a href=".$data{$org}{file}. " target=_new>$org</a>"};
	  }
      }


    return  return [@data,@no_data];
  }
  
  sub get_map
  {
    my $file = shift;
    my $map;
    open (IN, $file) || die "$!";
    while (<IN>)
      {
	$map .= $_;
      }
    close IN;
    return $map;
  }
  
sub generate_feat_info
{
	my $fid = shift;
	#print STDERR $fid, "\n"; 
	my $width = 1000;
	my $feat = $coge->resultset('Feature')->find($fid);
	my $hpp = $feat->annotation_pretty_print_html();
	my $cs = new CoGe::Graphics::Chromosome ();
	my $ds = $coge->resultset('Dataset')->find($feat->dataset->id);
	my $last_pos = $ds->last_chromosome_position($feat->chr);
    $cs->chr_length($last_pos);
    $cs->iw($width);
    $cs->ih(200);
    $cs->draw_chromosome(1);
    $cs->draw_ruler(1);
    $cs->draw_chr_end(0);
    $cs->mag(0);
    $cs->mag_off(1);
    $cs->minor_tick_labels(0);
    $cs->major_tick_labels(1);
    $cs->draw_hi_qual(0);
    $cs->padding(2);
    $cs->set_region(start=>$feat->start - 1000, stop=>$feat->stop + 1000, forcefit=>1);
    $cs->auto_zoom(0);
    $cs->feature_height(10);
    $cs->overlap_adjustment(0);
    $cs->feature_labels(1);
    my $graphics = new CoGe::Graphics;
    $graphics->process_features(c=>$cs, layers=>{features=>{gene=>1, cds=>1, mrna=>1, rna=>1, cns=>1}}, ds=>$ds, chr=>$feat->chr, coge=>$coge);
    $cs->overlap_adjustment(1);
    my $sub_file = $cogeweb->basefile.".s.".$feat->id.".png";
    $cs->generate_png(file=>$sub_file);
    my ($name) = $feat->names;
    my $subject_link = qq{
Feature: $name<br>
<a title='Click for Interactive Genome View' href = 'GeLo.pl?chr=}.$feat->chr.qq{&ds=}.$feat->dataset->id.qq{&x=}.$feat->start.qq{&z=5' target=_new border=0><img src=$sub_file border=0></a>
};
    $subject_link =~ s/$TEMPDIR/$TEMPURL/;
    #print STDERR "made it this far!\n";
    return $hpp, $subject_link;
}

sub color_feat_by_type
{
	my %opts = @_;
	my $feat = $opts{feat};
	my $anticodon = $opts{anticodon};
	my $clade = $opts{clade};
	if ($anticodon) {
		my $anno = join (" ", map {$_->annotation} $feat->annotations);
        my ($anticodon) = $anno =~ /\(anticodon:\s+(\w+)\)/;
        if ($anticodon =~ /AGA/i) {return [0,0,255];}
        elsif ($anticodon =~ /GGA/i) {return [0,255,0];}
		elsif ($anticodon =~ /GCT/i) {return [255,0,0];}
		elsif ($anticodon =~ /TGA/i) {return [139,137,137];}
		elsif ($anticodon =~ /CGA/i) {return [255,127,36];}
		else {return [69,69,69];}
        
	}
	elsif($clade)
	{
	my $fid = $feat->id;
	if($fid == 4157016 || $fid == 4181732 || $fid == 4100298) {return [0,0,255];}
	elsif ($fid == 4140192) {return [0,255,0];}
	elsif ($fid == 4181648) {return [255,0,0];}
	elsif ($fid == 4140158) {return [139,137,137];}
	elsif ($fid == 4157026 || $fid == 4100322) {return [202,255,112];}
	elsif ($fid == 4118458) {return [255,215,0];}
	elsif ($fid == 4140080) {return [255,193,193];}
	elsif ($fid == 4181724) {return [139,58,58];}
	elsif ($fid == 4181644 || $fid == 4118320) {return [255,255,224];}
	elsif ($fid == 4099894 || $fid == 4100252 || $fid == 4100336) {return [238,0,238];}
	elsif ($fid == 4100262 || $fid == 4100168) {return [255,181,173];}
	elsif ($fid == 4099980) {return [135,206,250];}
	elsif ($fid == 4118318 || $fid == 4118442) {return [0,100,0];}
	elsif ($fid == 4118336 || $fid == 4118358) {return [102,205,170];}
	elsif ($fid == 4181618 || $fid == 4181692 || $fid == 4100260 || $fid == 4099954 || $fid == 4099962 || $fid == 4100136 || $fid == 4099922 || $fid == 4100018 || $fid == 4100014 || $fid == 4100280 || $fid == 4099904 || $fid == 4100286 || $fid == 4100182 || $fid == 4100178 || $fid == 4099916 || $fid == 4100040 || $fid == 4100126 || $fid == 4100222 || $fid == 4100200 || $fid == 4100034) {return [152,251,152];}
	elsif ($fid == 4099924) {return [176,196,222];}
	elsif ($fid == 4099986 || $fid == 4099982) {return [25,25,112];}
	elsif ($fid == 4181706) {return [123,104,238];}
	elsif ($fid == 4157038 || $fid == 4181604 || $fid == 4118342) {return [124,252,0];}
	elsif ($fid == 4140082) {return [34,139,34];}
	elsif ($fid == 4100026) {return [139,69,19];}
	elsif ($fid == 4157006 || $fid == 4140086) {return [205,133,63];}
	elsif ($fid == 4100176 || $fid == 4100158) {return [178,34,34];}
	elsif ($fid == 4100238 || $fid == 4099908) {return [0,255,255];}
	elsif ($fid == 4140210) {return [255,0,255];}
	elsif ($fid == 4140124) {return [0,255,255];}
	elsif ($fid == 4100312 || $fid == 4100002) {return [199,21,133];}
	else {return [148,0,211];}	
	}
	else	{
		my $feat_type = $feat->type->name;
		if ($feat_type =~ /gene/i) {return [0,0,255];} #blue
		elsif ($feat_type =~ /CDS/i) {return [0,255,0];} #green
		elsif ($feat_type =~ /mRNA/i) {return [255,0,0];} #red
		elsif ($feat_type =~ /tRNA/i) {return [139,137,137];} #light grey
		elsif ($feat_type =~ /RNA/i) {return [69,69,69];} #dark grey
		else {return [255,127,36];} #orange
	}

}
