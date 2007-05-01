#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Genome;

use vars qw($USER $UID $LAST_LOGIN $FORM $DATE);

$FORM = new CGI;
($USER, $UID, $LAST_LOGIN) = CoGe::Accessory::LogUser->get_user();

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       login=>\&login,
		      );
print $pj->build_html($FORM, \&gen_html);




sub gen_html
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Comparative Genomics Homepage');
    $template->param(HELP=>'CoGe');
    $template->param(HEAD=>'<SCRIPT language="JavaScript" type="text/javascript" src="./js/kaj.stable.js"></SCRIPT>');
    if ($FORM->param('logout') || !$USER)
      {
	$template->param(USER=>"Not logged in");
      }
    else
      {
	$template->param(USER=>$USER);
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
    if ($USER && !$FORM->param('logout'))
      {
	$tmpl->param(ACTIONS=>actions());
      }
    else
      {
	$tmpl->param(LOGIN=>1);
	my $url = $FORM->param('url') if $FORM->param('url');
	$url =~ s/:::/;/g if $url;
	$tmpl->param(url=>$url);
      }
    $html .= $tmpl->output;
    if ($FORM->param('logout'))
      {
	my $c = CoGe::Accessory::LogUser->gen_cookie(user_name=>$USER, uid=>$UID, session=>$LAST_LOGIN, exp=>"Thu, 01-Jan-1970 00:00:01 GMT");
	$html .= "<script language = 'javascript'>window.onload=delete_cookie('$c');</script>";
      }
    return $html;
  }

sub actions
  {
    my @actions = (
		   {
		    ACTION=>qq{<a href="./GEvo.pl">GEvo: Genome Evolution Analysis</a>},
		    DESC => qq{Allows uses to compare sequences using a variety of sequence comparison algorithms to find regions of similarity and displays the results graphically.  This program has several unqiue options to aid in the finding of evolutionary conserved non-coding sequences and large syntenic regions.},
		   },
		   {
		    
		    ACTION=>qq{<a href="./FeatView.pl">FeatView: Feature Viewer</a>},
		    DESC => qq{Allows users to find and display information about a genomic feature.  A <font class=oblique>feature</font> in terms of sequence or genomic data denotes some genomic region that has something known about it.  For example, if some genomic region is known to be a gene, and that gene makes an mRNA transcript, and that mRNA is translated into a functional protein, and that protein has several functional domains, it will be likely that there will be several genomic features in that region:  One to denote the gene, one or more to denote the mRNA, and one or more to denote the protein, and one or more to denote the functional domains.  },
		   },
		   {
		    ACTION => qq{<a href="./MSAView.pl">MSAView: Multiple Sequence Alignment Viewer</a>},
		    DESC   => qq{Allows users to submit a multiple sequence alignment in FASTA format (if people would like additional formats, please request via e-mail) in order to quickly check the alignment, find conserved regions, etc.  This program also generates a consensus sequence from the alignment and displays some basic statistics about the alignment.},
		   },
		   {
		    ACTION => qq{<a href="./TreeView.pl">TreeView: Phylogenetic Tree Viewer</a>},
		    DESC   => qq{Allows users to submit a tree file and get a graphical view of their tree.  There is support for drawing rooted and unrooted trees, zooming and unzooming functions, and coloring and shaping nodes based on user specifications.},
		   },
		   {
		    ACTION => qq{<a href="./GeLo.pl">GeLo: Genome Location Viewer</a>},
		    DESC   => qq{Allows users to get an overview of organisms in genomes database and provides a dynamic, interactive browser of genomic information.},
		   },
		   {
		    ACTION => qq{<a href="./CoGeBlast.pl">CoGe Blast</a>},
		    DESC   => qq{Blast sequences from the CoGe system against various genomic sequence sets.},
		   },
		   {
		    ACTION => qq{<a href="./help/CoGe/">CoGe Faq</a>},
		    DESC   => qq{What is CoGe?  This document covers some of the basics about what CoGe is, how it has been designed, and other information about the system.},
		   },
		  );
    return \@actions;
  }

sub login
  {
    my ($name, $pwd, $url) = @_;
    my $db = new CoGe::Genome;
    my ($u) = $db->get_user_obj->search(user_name=>$name);
    my $pwdc = $u->check_passwd(pwd=>$pwd) if $u;
    $url = $FORM->param('url') unless $url;
    $url = $FORM->url() unless $url;
    if ($pwdc)
      {
	my $sid = $db->get_user_session_obj->log_user(user=>$u);
	my $c = CoGe::Accessory::LogUser->gen_cookie(user_name=>$name, uid=>$u->id, session=>$sid);
	return (1, $c, $url );
      }
    elsif ($name =~ /^public$/i)
      {
	my $c = CoGe::Accessory::LogUser->gen_cookie(user_name=>$name);
	return (1, $c,  $url);
      }
    else
      {
	return (0,undef, $url);
      }
  }

