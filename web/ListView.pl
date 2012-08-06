#! /usr/bin/perl -w

use strict;
use CGI;
#use CGI::Ajax;
use JSON::XS;
use lib '/home/mbomhoff/CoGe/Accessory/lib'; #FIXME remove
use lib '/home/mbomhoff/CoGeX/lib'; #FIXME 8/2/12 remove
use CoGe_dev::Accessory::LogUser;
use CoGe_dev::Accessory::Web;
use CoGeX_dev;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use URI::Escape;
use Data::Dumper;
use File::Path;

no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL);
$P = CoGe_dev::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX_dev->connect($connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL = $P->{URL};
$COGEDIR = $P->{COGEDIR};
$TEMPDIR = $P->{TEMPDIR}."ListView/";
mkpath ($TEMPDIR, 0,0777) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL}."ListView/";


my ($cas_ticket) =$FORM->param('ticket');
$USER = undef;
($USER) = CoGe_dev::Accessory::Web->login_cas(cookie_name=>$COOKIE_NAME, ticket=>$cas_ticket, coge=>$coge, this_url=>$FORM->url()) if($cas_ticket);
($USER) = CoGe_dev::Accessory::LogUser->get_user(cookie_name=>$COOKIE_NAME,coge=>$coge) unless $USER;

%FUNCTION = (
	     gen_html=>\&gen_html,
	     get_list_info=>\&get_list_info,
    );

dispatch();

sub dispatch
{
    my %args = $FORM->Vars;
    my $fname = $args{'fname'};
    if($fname)
    {
	#my %args = $FORM->Vars;
	#print STDERR Dumper \%args;
	if($args{args}){
	    my @args_list = split( /,/, $args{args} );
	    print $FORM->header, $FUNCTION{$fname}->(@args_list);
       	}
	else{
	    print $FORM->header, $FUNCTION{$fname}->(%args);
	}
    }
    else{
	print $FORM->header, gen_html();
    }
}

sub gen_html
  {
    my $html;    
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
    $template->param(HELP=>'/wiki/index.php?title=ListView');
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $template->param(TITLE=>qq{Managing Data});
    $template->param(PAGE_TITLE=>qq{ListView});
    $template->param(LOGO_PNG=>"ListView-logo.png");
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    my ($body, $list_name) = gen_body();
    $template->param(BODY=>$body);
    $template->param(BOX_NAME=>$list_name);
    $template->param(ADJUST_BOX=>1);
    $html .= $template->output;
  }


sub gen_body
  {
      my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'ListView.tmpl');
      $template->param(PAGE_NAME=>$FORM->url);
      $template->param(MAIN=>1);
      my $lid;
      $lid = $FORM->param('lid');
      my ($list_info, $list_name) = get_list_info(lid=>$lid);
      $template->param(LIST_INFO=>$list_info);
      $template->param(ADMIN_AREA=>1) if $USER->is_admin;
      my $open;
      $open = $FORM->param('open') if defined $FORM->param('open');
      my $box_open = $open ? 'true' : 'false';
      $template->param(EDIT_BOX_OPEN=>$box_open);
      return $template->output, $list_name;
  }

sub get_list_info
  {
    my %opts = @_;
    my $lid=$opts{lid};
    return "Must have valid list id.\n" unless ($lid);
    my ($list) = $coge->resultset('List')->find($lid);
    return "Unable to create list object for $lid\n" unless ($list);
    my $name = $list->name;
    $name .= ": ". $list->description if $list->description;
    return $list->annotation_pretty_print_html(), $name;
  }
