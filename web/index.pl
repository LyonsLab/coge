#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Cookie;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use Digest::MD5 qw(md5_base64);
use CoGeX;

use vars qw($USER $FORM $DATE $update $coge);

$ENV{PATH} = "/opt/apache/CoGe";
$FORM = new CGI;
($USER) = CoGe::Accessory::LogUser->get_user();
#print STDERR Dumper $USER->user_name;

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       login=>\&login,
		      );
$update =0;
$coge = new CoGeX->dbconnect;
print $pj->build_html($FORM, \&gen_html);
#print gen_html();



sub gen_html
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Comparative Genomics Homepage');
    $template->param(HELP=>'CoGe');

    if ($FORM->param('logout') || !$USER)
      {
	$template->param(USER=>"Not logged in");
      }
    else
      {
	my $name = $USER->user_name;
	$name = $USER->first_name if $USER->first_name;
	$name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
	$template->param(USER=>$name);
	$template->param(LOGON=>1) unless $USER->user_name eq "public";

      }
    $template->param(DATE=>$DATE);
    $template->param(BOX_NAME=>"Welcome!");
    $template->param(LOGO_PNG=>"CoGe-logo.png");
    $template->param(BODY=>gen_body());

    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $tmpl = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/index.tmpl');
    my $html;
    my %cookies = fetch CGI::Cookie;
    my $disable = 0;
    if ($cookies{'anim_speed'})
     {
      $disable = $cookies{'anim_speed'}->value() eq "0" ? 1 : 0;
     }
    if ($update)
      {
	$tmpl->param(update=>1);
      }
    elsif ($USER && !$FORM->param('logout') && !$FORM->param('login'))
      {
#	$disable = 1;
	if ($disable)
	  {
	    $tmpl->param(DISABLE=>1);
	    $tmpl->param(ACTIONS=>[map {{ACTION=>$_->{ACTION}, DESC=>$_->{DESC}}} @{actions()}  ]);
	  }
	else
	  {
	    $tmpl->param(ANIMATE=>1);
	    $tmpl->param(ACTIONS=>actions());
	  }
	$tmpl->param('INTRO'=>1);
      }
#     else
#       {
# 	$tmpl->param(LOGIN=>1);
# 	my $url = $FORM->param('url') if $FORM->param('url');
# 	$url =~ s/:::/;/g if $url;
# 	$tmpl->param(url=>$url);
#       }

    if ($FORM->param('logout'))
      {
	my $session = md5_base64($USER->user_name.$ENV{REMOTE_ADDR});
	($session) = $coge->resultset('UserSession')->find({session=>$session});
	$session->delete if $session;
	$tmpl->param(READY=>"delete_cookie();");
	$tmpl->param(LOGIN=>1);
      }
    $tmpl->param(LOGIN=>1) if $FORM->param('login');
    $html .= $tmpl->output;
    return $html;
  }

sub actions
  {
    my @actions = (
		   {
		    ID=>1,
		    LOGO=>qq{<a href="./GEvo.pl"><img src="/CoGe/picts/carousel/GEvo-logo.png" width="227" height="75" border="0"></a>},
		    ACTION=>qq{<a href="./GEvo.pl">GEvo</a>},
		    DESC => qq{Compare sequences using a variety of sequence comparison algorithms to find regions of similarity and displays the results graphically.<a href="/CoGe/docs/help/GEvo/Overview.html">Learn more...</a>},
		    SCREENSHOT=>qq{<a href="./GEvo.pl"><img src="/CoGe/picts/preview/GEvo.png"border="0"></a>},
		   },
		   {
		    ID=>2,
		    LOGO=>qq{<a href="./FeatView.pl"><img src="/CoGe/picts/carousel/FeatView-logo.png" width="227" height="75" border="0"></a>},
		    ACTION=>qq{<a href="./FeatView.pl">FeatView</a>},
		    DESC => qq{Find and display information about a genomic feature (e.g. gene) by searching names and annotations.},
		    SCREENSHOT=>qq{<a href="./FeatView.pl"><img src="/CoGe/picts/preview/FeatView.png" width="400" height="241" border="0"></a>},
		   },
# 		   {
# 		    ID=>3,
# 		    LOGO=>qq{<a href="./MSAView.pl"><img src="/CoGe/picts/carousel/MSAView-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./MSAView.pl">MSAView: Multiple Sequence Alignment Viewer</a>},
# 		    DESC   => qq{Allows users to submit a multiple sequence alignment in FASTA format (if people would like additional formats, please request via e-mail) in order to quickly check the alignment, find conserved regions, etc.  This program also generates a consensus sequence from the alignment and displays some basic statistics about the alignment.},
# 		    SCREENSHOT=>qq{<a href="./MSAView.pl"><img src="/CoGe/picts/preview/MSAView.png"border="0"></a>},
# 		   },
# 		   {
# 		    ID=>4,
# 		    LOGO=>qq{<a href="./TreeView.pl"><img src="/CoGe/picts/carousel/TreeView-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./TreeView.pl">TreeView: Phylogenetic Tree Viewer</a>},
# 		    DESC   => qq{Allows users to submit a tree file and get a graphical view of their tree.  There is support for drawing rooted and unrooted trees, zooming and unzooming functions, and coloring and shaping nodes based on user specifications.},
# 		    SCREENSHOT=>qq{<a href="./FeatView.pl"><img src="/CoGe/picts/preview/TreeView.png"border="0"></a>},
# 		   },
		   {
		    ID=>5,
		    LOGO => qq{<a href="./GenomeView.pl"><img src="/CoGe/picts/carousel/GenomeView-logo.png" width="227" height="75" border="0"></a>},
		    ACTION => qq{<a href="./GenomeView.pl">GenomeView</a>},
		    DESC   => qq{Get an overview of organisms in CoGe's genomes database and provides a dynamic, interactive genome browser.},
		    SCREENSHOT => qq{<img src="/CoGe/picts/preview/GenomeView.png" border="0"></a>},
		   },
		   {
		    ID => 6,
		    LOGO => qq{<a href="./CoGeBlast.pl"><img src="/CoGe/picts/carousel/CoGeBlast-logo.png" width="227" height="75" border="0"></a>},
		    ACTION => qq{<a href="./CoGeBlast.pl">CoGeBlast</a>},
		    DESC   => qq{Blast sequences from the CoGe system against multiple genomes and view the results in an interactive graphical system.},
		    SCREENSHOT => qq{<a href="./CoGeBlast.pl"><img src="/CoGe/picts/preview/Blast.png" width="400" height="241" border="0"></a>},
		   },
# 		   {
# 		    ID => 7,
# 		    LOGO => qq{<a href="./docs/help/CoGe"><img src="/CoGe/picts/carousel/FAQ-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./docs/help/CoGe/">CoGe Faq</a>},
# 		    DESC   => qq{What is CoGe?  This document covers some of the basics about what CoGe is, how it has been designed, and other information about the system.},
# 		    SCREENSHOT => qq{<a href="./docs/help/CoGe"><img src="/CoGe/picts/preview/app_schema.png" border="0"></a>},
# 		   },
		  );
    return \@actions;
  }

sub login
  {
    my ($name, $pwd, $url) = @_;
    my $coge = CoGeX->dbconnect();
    my ($u) = $coge->resultset('User')->search({user_name=>$name});
    my $pwdc = $u->check_passwd(pwd=>$pwd) if $u;
    $url = $FORM->param('url') unless $url;
    $url = $FORM->url() unless $url;
    if ($pwdc)
      {
	my $session = md5_base64($name.$ENV{REMOTE_ADDR});
	my $sid = $coge->log_user(user=>$u,session=>$session);
	my $c = CoGe::Accessory::LogUser->gen_cookie(session=>$session);
	return ('true', $c, $url );
      }
    elsif ($name =~ /^public$/i)
      {
	my $c = CoGe::Accessory::LogUser->gen_cookie(session=>"public");
	return ('true', $c,  $url);
      }
    else
      {
	return ('false',undef, $url);
      }
  }
