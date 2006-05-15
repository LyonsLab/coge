#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use CoGe::Genome;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use POSIX;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $ACCN);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$ACCN = $FORM->param('accn');
$USER = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;
my $pj = new CGI::Ajax(
		       db_lookup=>\&db_lookup,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       accn_search=>\&accn_search,
		       clear_div=>\&clear_div,
		       get_anno=>\&get_anno,
		       show_location=>\&show_location,
		       show_express=>\&show_express,
		       gen_data=>\&gen_data,
		       get_dna_seq_for_feat => \&get_dna_seq_for_feat,
		       get_prot_seq_for_feat => \&get_prot_seq_for_feat,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);


sub get_prot_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->search(feature_id=>$featid);
    my ($seq) = $DB->get_protein_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    return $seq;
  }

sub get_dna_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->search(feature_id=>$featid);
    my $seq = $DB->get_genomic_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    return $seq;
  }

sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

sub clear_div
  {
    return;
  }

sub get_types
  {
    my ($accn, $info_ver) = @_;
    my $html;
    my $blank = qq{<input type="hidden" id="Type_name">};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $DB->get_features_by_name_and_version(name=>$accn, version=>$info_ver);
    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_anno(['accn_select','infoid','Type_name'],['anno'])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1);
  }


sub accn_search
  {
    my $accn = shift;

    my $blank = qq{<input type="hidden" id="accn_select">};
    return $blank unless length($accn) > 2;
    my $html;

    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {uc($_)} $DB->get_feature_name_obj->search_name($accn."%");
    $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select" SIZE="5" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank."Accession not found.\n" unless $html =~ /OPTION/;
    return $html;
  }

sub get_anno
  {
    my $accn = shift;
    return unless $accn;
    my $info_ver = shift;
    my $type = shift;
    my @feats;
    foreach my $feat ($DB->get_features_by_name_and_version(name=>$accn, version=>$info_ver))
      {
	push @feats, $feat if ($feat->type->name eq $type);
      }
    my $anno;
    $anno .= "<font class=small>Annotation count: ".scalar @feats."</font>\n<BR>\n" if scalar @feats;
    my $i = 0;
    foreach my $feat (@feats)
      {
	$anno .= join "\n<BR><HR><BR>\n", $feat->annotation_pretty_print_html();
	$anno .= qq{<DIV id="loc$i"><input type="button" value = "Click for chromosomal view" onClick="gen_data(['args__Generating chromosomal view image'],['loc$i']);show_location(['args__}.$feat->begin_location.qq{', 'args__}.$feat->end_location.qq{', 'args__}.$feat->chr.qq{', 'args__}.$feat->info->id.qq{'],['loc$i']);"></DIV>};
	$anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
	$anno .= qq{<DIV id="dnaseq$i"><input type="button" value = "Click for DNA sequence" onClick="gen_data(['args__retrieving sequence'],['dnaseq$i']);get_dna_seq_for_feat(['args__}.$feat->id.qq{'],['dnaseq$i']);"></DIV>};
	$anno .= qq{<DIV id="protseq$i"><input type="button" value = "Click for protein sequence" onClick="gen_data(['args__retrieving sequence'],['protseq$i']);get_prot_seq_for_feat(['args__}.$feat->id.qq{'],['protseq$i']);"></DIV>};

#	    qq{<DIV id="exp$i"></DIV>};
	$anno = "<font class=\"annotation\">No annotations for this entry</font>" unless $anno;
	$i++;
      }
    return ($anno);
  }

sub show_location
  {
    my %opts = @_;
    my $start = shift;#$opts{start};
    my $stop = shift; #$opts{stop};
    my $chr = shift; # $opts{chr};
    my $info_id = shift; # = $opts{info_id};
#    my $num= int(rand(100000));
#    my $file = "./tmp/pict$num.png";
#    system("GenomePNG.pl $start $stop $chr $info_id $file");
#    gen_pict(start=>$start, stop=>$stop, chr=>$chr, info_id=>$info_id,file=>$file);
#    return "<img src=$file>";
    my $z = ceil (log((10000+$stop-$start)/10)/log(2));
    my $link = qq{<img src="tiler.pl?start=$start&chr=$chr&di=$info_id&iw=2000&z=$z">\n};
    print STDERR $link;
    return $link;
  }


sub show_express
  {
    my %opts = @_;
    my $accn = shift;
    my $log = shift;
    my $div = shift;
    $log = 1 unless defined $log;
    $accn =~ s/\..*//;
    my $link = qq{<img src="expressiontree.pl?locus=$accn&label=1&rw=80&rh=8&name=1&legend=1&mean=1&log_trans=$log">\n};
    $log = $log ? 0 : 1;
    my $type = $log ? "log transformed" : "normal";
    print STDERR $link;
    $link .= qq{<br><input type="button" value = "Click for $type expression tree" onClick="gen_image([],['exp$div']);show_express(['args__}.$accn.qq{','args__}.$log.qq{','args__}.$div.qq{'],['exp$div']);">};
    return $link;
  }


sub gen_html
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatView.tmpl');
    $template->param(TITLE=>'Genome Location Feature Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(TABLE_NAME1=>"Feature Selection");
    $template->param(ACCN=>$ACCN);	
#    $template->param(ADDIIONAL_JS=>'accn_search_chain(0);');
    if ($ACCN)
      {

      }
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub get_data_source_info_for_accn
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="infoid">};
    return $blank unless $accn;
    my @feats = $DB->get_feats_by_name($accn);
    my %sources;
    foreach my $feat (@feats)
      {
	my $val = $feat->data_info;
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->data_source->name;
	my $org = $val->org->name;
	my $title = "$org: $sname, version $ver";
	$sources{$title} = $ver
      }
    my $html;
    $html .= qq{
<SELECT name = "infoid" id="infoid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
    my $count = 0;
    foreach my $title (reverse sort keys %sources)
      {
	my $ver = $sources{$title};
	$html .= qq{  <option value="$ver" >$title\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Version count: ".$count ."</font>\n<BR>\n".$html, 1);
  }

sub get_data_source_info_for_accn_old
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="infoid">};
    return $blank unless $accn;
    my @feats = $DB->get_feats_by_name($accn);
    my %sources;
    foreach my $feat (@feats)
      {
	$sources{$feat->data_info->id} = $feat->data_info;
      }
    my $html;
    $html .= "<font class=small>Data count: ".scalar keys (%sources) ."</font>\n<BR>\n";
    $html .= qq{
<SELECT name = "infoid" id="infoid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
#<SELECT name = "infoid" id="infoid" MULTIPLE SIZE="5" onChange="get_types(['accn_select', 'infoid'],['FeatType']);">
    my $count = 0;
    foreach my $id (sort {$b <=> $a} keys %sources)
      {
	my $val = $sources{$id};
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->data_source->name;
	my $org = $val->org->name;
	$html .= qq{  <option value="$id" >$org: $sname, version $ver\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ($html);
  }


sub gen_pict
  {
    my %opts = @_;
    my $start = $opts{'start'};
    my $stop = $opts{'stop'};
    $stop = $start unless $stop;
    my $di = $opts{'info_id'};
    my $chr = $opts{'chr'};
    my $file = $opts{'file'};
    print STDERR Dumper @_;
    return unless ($start && $stop && $di);
    my $c = new CoGe::Graphics::Chromosome;
    my $chr_length = $DB->get_genomic_sequence_obj->get_last_position($di);
    $c->chr_length($chr_length);
    $c->iw(1600);
    $c->max_mag((80));
    $c->DEBUG(0);
    $c->feature_labels(0);
    $c->fill_labels(1);
    $c->draw_chromosome(1);
    $c->draw_ruler(1);
    $c->set_region(start=>$start, stop=>$stop);
    $c->mag($c->mag-1);
    $start = $c->_region_start;
    $stop = $c->_region_stop;

    foreach	 my $feat ($DB->get_feature_obj->get_features_in_region(start=>$start, end=>$stop, info_id=>$di, chr=>$chr)) 
      {
	my $f;
	if ($feat->type->name =~ /Gene/i) {
	  $f = CoGe::Graphics::Feature::Gene->new();
	  $f->color([255,0,0,50]);
	  foreach	 my $loc ($feat->locs) {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	  }
	  $f->order(1);
	} elsif ($feat->type->name =~ /CDS/i) {
	  $f = CoGe::Graphics::Feature::Gene->new();
	  $f->color([0,255,0, 50]);
	  foreach	 my $loc ($feat->locs) {
	    $f->add_segment(start=>$loc->start, stop=>$loc->stop);
	    $f->strand($loc->strand);
	  }
	$f->order(3);
      } elsif ($feat->type->name =~ /rna/i) {
	$f = CoGe::Graphics::Feature::Gene->new();
	$f->color([0,0,255, 50]);
	foreach	 my $loc ($feat->locs) {
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
    my $seq = uc($DB->get_genomic_sequence_obj->get_sequence(start=>$start, end=>$stop, chr=>$chr, info_id=>$di)); 
    my $seq_len = length $seq;
    my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw);
    my $pos = 0;
    $start = 1 if $start < 1;
    while ($pos < $seq_len)
      {
	my $subseq = substr ($seq, $pos, $chrs);
	my $rcseq = $subseq;
	$rcseq =~ tr/ATCG/TAGC/;
#	print STDERR $subseq,"\t", $rcseq,"\n";
	my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$subseq, strand=>1, start =>$pos+$start});
	my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
	$c->add_feature($f1, $f2);
	$pos+=$chrs;
	#    print "working on position ",$i+$start,"\n";
      }
    $c->generate_png(file=>"/opt/apache/CoGe/$file");
  }
