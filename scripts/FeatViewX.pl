#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use CoGeX;
use POSIX;

$ENV{PATH} = "/opt/apache2/bpederse/cogex";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $ACCN $FID %restricted_orgs);

# set this to 1 to print verbose messages to logs
$DEBUG = 1;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$ACCN = $FORM->param('accn');
($USER) = 1; #CoGe::Accessory::LogUser->get_user();
if (!$USER || $USER =~ /public/i)
  {
    $restricted_orgs{papaya} = 1;
  }

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $DB = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#$DB->storage->debug(1);

my $pj = new CGI::Ajax(
		       db_lookup=>\&db_lookup,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       cogesearch=>\&cogesearch,
		       get_anno=>\&get_anno,
		       show_express=>\&show_express,
		       gen_data=>\&gen_data,
		      );
$pj->JSDEBUG(2);
$pj->DEBUG(2);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);

sub gen_data
  {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
  }

sub get_types
  {
    my ($accn, $dsid) = @_;
    my $blank = qq{<input type="hidden" id="Type_name">};
    my %seen;
    $DB->storage->debug(1);
    my $rs = $DB->resultset('Feature')->esearch({
            'feature_names.name' => { 'like' => '%' . $accn . '%' }
            ,'dataset.dataset_id' => $dsid
            }, {
                'order_by'  => ['feature_type.name']
               }

            );
    my @opts;
    while (my $feat = $rs->next()){
        my $ftype = $feat->feature_type->name;
        next if $seen{$ftype};
        $seen{$ftype}++;
        push(@opts,"<option>$ftype<option");
    }
    my $html = "<font class=small>Type count: "
         . scalar @opts
         ."</font>\n<BR>\n"
         . qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_anno(['accn_select','Type_name', 'dsid'],['anno'])" >\n}
         .  join ("\n", @opts)
         . "\n</SELECT>\n";

    $html =~ s/option/option SELECTED/;
    return $blank unless $html =~ /option/;
    return ($html, 1);
  }

sub cogesearch
  {
    my $accn = shift;
    my $type = shift;
    my $org = shift;
    my $blank = qq{<input type="hidden" id="accn_select">};
    return $blank unless (length($accn) > 2 || $type || $org);
    my %seen;
    my $rs = $DB->resultset('Feature')->esearch({
            'feature_names.name' => { 'like' =>  $accn . '%' }
            ,'feature_type.feature_type_id' => $type
          #  ,'organism.name'     => {'like' => "%$org%" }
            }, {
                'order_by'  => ['feature_type.name']
                , 'limit'   => 1000
               }

            );
    my @opts;
    while( my $feat= $rs->next()){
        next if $org and $feat->dataset->organism->organism_id != $org;
        my $fnames = $feat->feature_names();
        while (my $fname = $fnames->next()){
            my $name = $fname->name;
            next if $seen{$name};
            $seen{$name}=1;
            push(@opts,"<option>$name</option>");
        }
    }
    my $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n"
         . qq{<SELECT id="accn_select" SIZE="5" MULTIPLE onChange="source_search_chain(); " >\n}
         .  join ("\n", @opts)
         .  "\n</SELECT>\n";

    $html =~ s/<option/<option SELECTED/;
    return $blank."No results found.\n" unless $html =~ /option/;
    return $html;
  }

sub get_anno
  {
    my $accn = shift;
    return unless $accn;
    my $type = shift;
    my $dataset_id = shift;
    my @feats;
    my $rs = $DB->resultset('Feature')->search({
            'feature_names.name' => { 'like' =>  $accn . '%' }
            ,'dataset.dataset_id' => $dataset_id
            ,'feature_type.name' => $type
            }, {
                'join'      =>
                    ['feature_names','feature_type','locations','dataset' ]
                ,'prefetch' => ['feature_type','locations','dataset']
                ,'order_by'  => ['feature_type.name']
                , 'limit'   => 1000
               }

            );
    #TODO: move this after
    my $i = 0;
    my $anno = '';
    while( my $feat= $rs->next()){

        push(@feats,$feat->feature_type->name);
        $i++;
        my $featid = $feat->feature_id;
        my $chr = $feat->chromosome;
        my $rc = 0;
        my $pro = 0;
        my $ds = $feat->dataset->dataset_id;
        my $x = $feat->start;
        my $z = 10;
        $anno .= qq{<DIV id="loc$i"><input type="button" value = "Click for chromosomal view" onClick="window.open('GeLo.pl?chr=$chr&ds=$ds&INITIAL_CENTER=$x,0&z=$z');"></DIV>};
        $anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
        $anno .= qq{<DIV id="dnaseq$i"><input type="button" value = "Click for DNA sequence" onClick="window.open('SeqView.pl?featid=$featid&dsid=$ds&chr=$chr&featname=$accn&rc=$rc&pro=$pro','window', 'menubar=1, resizable=1, scrollbars=1, location=1,left=20, top=20, width=700, height=610');"></DIV>};
        #"gen_data(['args__retrieving sequence'],['dnaseq$i']);get_dna_seq_for_feat(['args__}.$feat->id.qq{'],['dnaseq$i']);"></DIV>};
        #$anno .= qq{<DIV id="rcdnaseq$i"><input type="button" value = "Click for DNA sequence reverse complement" onClick="gen_data(['args__retrieving sequence'],['rcdnaseq$i']);get_rcdna_seq_for_feat(['args__}.$feat->id.qq{'],['rcdnaseq$i']);"></DIV>};

    $anno .= "<font class=small>Annotation count: ". scalar @feats .  "</font>\n<BR>\n" if scalar @feats;
	$anno = "<font class=\"annotation\">No annotations for this entry</font>" unless $anno;
      }
      print STDERR $anno . "\n";
    return ($anno);
  }

#'
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
    my $html;
    unless ($USER)
      {
	$html = login();
      }
    else
      {
	my $template = HTML::Template->new(filename=>'/opt/apache/bpederse/cogex/tmpl/generic_page.tmpl');
	$template->param(LOGO_PNG=>"FeatView-logo.png");
	$template->param(TITLE=>'Feature Viewer');
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(bOX_NAME=>"Feature Selection");
	my $body = gen_body();
	$template->param(BODY=>$body);

	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/bpederse/cogex/tmpl/FeatView.tmpl');

    $template->param(ACCN=>$ACCN);
    $template->param(TYPE_LOOP=> [{TYPE=>"<OPTION VALUE=0>All</OPTION>"},map {{TYPE=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)}
    $DB->resultset('FeatureType')->search() ]);

    my @orgs;

    my $rs = $DB->resultset('Organism')->search({},{order_by => 'name'});

    #TODO: move this after
    my $rscount = $rs->count;
    my $anno = "<font class=small>Annotation count: ". $rscount .  "</font>\n<BR>\n" if $rscount;

    while( my $org =$rs->next())
      {
	push @orgs, $org unless $restricted_orgs{$org->name};
      }
      @orgs = sort @orgs;
    $template->param(ORG_LOOP=> [{ORG=>"<OPTION VALUE=0>All</OPTION>"},
        map {{ORG=>"<OPTION value=\"".$_->organism_id."\">".$_->name."</OPTION>"}}  @orgs]);
    my $html = $template->output;
    $html =~ s/(>gene<\/OPTION)/ SELECTED$1/;
    return $html;
  }

sub get_data_source_info_for_accn
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="dsid">};
    return $blank unless $accn;

    my $rs = $DB->resultset('Feature')->search({
            'feature_names.name' => { 'like' =>  $accn . '%' }
            }, {
                'join'      =>
                    ['feature_names','feature_type','locations','dataset']
                ,'prefetch' => ['feature_type','locations','dataset']
                ,'order_by'  => ['feature_type.name']
                , 'limit'   => 1000
               });
    my %sources;
    while( my $feat = $rs->next()){
        my $val = $feat->dataset;
        my $name = $val->name;
        my $ver = $val->version;
        my $desc = $val->description;
        my $sname = $val->data_source->name;
        my $ds_name = $val->name;
        my $org = $val->organism->name;
        my $title = "$org: $ds_name ($sname, v$ver)";
        next if $restricted_orgs{$org};
        $sources{$title} = $val->dataset_id;
    }
    my $html;
    $html .= qq{
<SELECT name = "dsid" id="dsid" MULTIPLE SIZE="5" onChange="get_types_chain();">
};
    my $count = 0;
    foreach my $title (sort keys %sources)
      {
	my $id = $sources{$title};
	$html .= qq{  <option value="$id" >$title\n};
	$html =~ s/option/option selected/ unless $count;
	$count++;
      }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Dataset count: ".$count ."</font>\n<BR>\n".$html, 1);
  }
