#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
use CoGe::Genome;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $DB);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$USER = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;

my $pj = new CGI::Ajax(
		       expand_list=>\&expand_list,
		       close_list=>\&close_list,
		       edit_list => \&edit_list,
		       delete_list => \&delete_list,
		       create_list => \&create_list,
		       delete_feature => \&delete_feature,
		       add_feature => \&add_feature,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);
#print gen_html();

sub gen_html
  {
    my ($body, $seq_names, $seqs) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');

    $template->param(TITLE=>'CoGe: Feature List Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"FeatListView-logo.png");
    $template->param(BODY=>$body);
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatListView.tmpl');
    $template->param(FRONT_PAGE => 1);
    $template->param(FEATURE_LIST=>get_feat_lists());
    return $template->output;
  }

sub get_feat_lists
  {
    #need to make sure to do a user validation at some point
    my @lists;
    foreach my $fl ($DB->get_feature_list_obj->retrieve_all())
      {
	my $f_count = scalar $fl->flc;
	push @lists, {LIST_NAME=>$fl->name, LIST_DESC=>$fl->desc, FEATURE_COUNT=>$f_count, LIST_ID=>$fl->id};
      }
    return \@lists;
  }

sub expand_list
  {
    my $lid = shift;
    my $fl = $DB->get_feature_list_obj->retrieve($lid);
    my @feats;
    foreach my $f ($fl->features)
      {
	push @feats, {FEAT_TYPE=>$f->type->name, FEAT_NAME => join ", ", sort map {"<a href = FeatView.pl?accn=".$_->name.">".$_->name."</a>"} $f->names};
      }
    my $tmpl = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatListView.tmpl');
#    $tmpl->param(LIST_ID=>$fl->id);
    $tmpl->param(EXPAND_LIST=>1);
    $tmpl->param(FEATURE_INFO=>\@feats);
    my $close = qq{
<input type="IMAGE" src="picts/open.png" onClick="close_list(['args__$lid'],['FLID$lid', 'showfeat$lid'])">
};
    return $tmpl->output, $close;
  }

sub edit_list
  {
    my $lid = shift;
    my $fl = $DB->get_feature_list_obj->retrieve($lid);
    my @feats;
    foreach my $f ($fl->features)
      {
	push @feats, {FEAT_ID=>$f->id, FEAT_ORG=>$f->org->name, FEAT_DATA=>$f->dataset->name."(v".$f->dataset->version.")", FEAT_TYPE=>$f->type->name, LIST_ID=>$lid, FEAT_NAME => join ", ", sort map {"<a href = FeatView.pl?accn=".$_->name.">".$_->name."</a>"} $f->names};
      }
    my $tmpl = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FeatListView.tmpl');
    $tmpl->param(EDIT_LIST=>1);
    $tmpl->param(LIST_NAME=>$fl->name);
    $tmpl->param(LIST_DESC=>$fl->desc);
    $tmpl->param(LIST_ID=>$fl->id);
    $tmpl->param(FEATURE_INFO=>\@feats);
    return $tmpl->output;
  }

sub close_list
  {
    my $lid = shift;
    my $fl = $DB->get_feature_list_obj->retrieve($lid);
    my $open = qq{<input type="IMAGE" src="picts/close.png" onClick="expand_list(['args__$lid'],['FLID$lid', 'showfeat$lid'])">};
    return scalar $fl->flc, $open;
  }

sub create_list
  {
    
  }

sub delete_list
  {
    my $lid = shift;
    my $fl = $DB->get_feature_list_obj->retrieve($lid);
    $fl->delete;
  }

sub delete_feature
  {
    my $fid = shift;
    my $lid = shift;
    foreach my $flc ($DB->get_feature_list_connector_obj->retrieve(feature_list_id=>$lid, feature_id=>$fid))
      {
	$flc->delete;
      }
    return $lid;
  }

sub add_feature
  {
    my $id = shift;
    my $lid = shift;
    return $lid unless $id =~ /^\d+$/;
    return $lid unless ref ($DB->get_feature_obj->retrieve($id)) =~ /Feature/;
    return $lid unless $lid =~ /^\d+$/;
    $DB->get_feature_list_connector_obj->find_or_create(feature_list_id=>$lid, feature_id=>$id);
    return $lid;
  }

sub search_feature
  {
    
  }
