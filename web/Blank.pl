#! /usr/bin/perl -w

use strict;
use CGI;
#use CGI::Ajax;
use JSON::XS;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use URI::Escape;
use Data::Dumper;
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $coge $cogeweb %FUNCTION $COOKIE_NAME $FORM $URL $COGEDIR $TEMPDIR $TEMPURL);
$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};
$URL = $P->{URL};
$COGEDIR = $P->{COGEDIR};
$TEMPDIR = $P->{TEMPDIR}."GenomeList/";
mkpath ($TEMPDIR, 0,0777) unless -d $TEMPDIR;
$TEMPURL = $P->{TEMPURL}."GenomeList/";


my ($cas_ticket) =$FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(ticket=>$cas_ticket, coge=>$coge, this_url=>$FORM->url()) if($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name=>$COOKIE_NAME,coge=>$coge) unless $USER;

%FUNCTION = (
	     gen_html=>\&gen_html,
	     test_ajax=>\&test_ajax,
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
    $template->param(HELP=>'/wiki/index.php?title=BLANK');
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $template->param(TITLE=>qq{Blank Page});
    $template->param(PAGE_TITLE=>qq{Blank Page});
    $template->param(LOGO_PNG=>"Blank-logo.png");
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(BODY=>gen_body());
    $template->param(ADJUST_BOX=>1);
    $html .= $template->output;
  }


sub gen_body
  {
      my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'Blank.tmpl');
      $template->param(PAGE_NAME=>$FORM->url);
      
      return $template->output;
  }

sub test_ajax
  {
    my %opts = @_;
    my $output = "Test worked!\n";
    $output .= "<pre>"."Args: ".Dumper (\%opts)."</pre>";
  }

