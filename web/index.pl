#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Cookie;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use Digest::MD5 qw(md5_base64);
use CoGeX;
use vars qw($P $USER $FORM $DATE $URL $update $coge);

$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = $P->{COGEDIR};
$URL = $P->{URL};
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
#print $FORM->header, gen_html();



sub gen_html
  {
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
    $template->param(TITLE=>'The Place to Compare Genomes');
    $template->param(PAGE_TITLE=>'ANKoCG');
    
    $template->param(HELP=>'/wiki/index.php');

    if ($FORM->param('logout') || !$USER)
      {
	$template->param(USER=>"public");
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
    my $welcome = "Welcome to CoGe!&nbsp&nbsp&nbsp";
    unless ($update)
      {
      	$welcome .= "<span class='small'>Organisms: ";
	$welcome .= commify($coge->resultset('Organism')->count());
	$welcome .= "&nbsp&nbsp&nbsp&nbsp   Genomes: ";
	$welcome .= commify($coge->resultset('DatasetGroup')->count());
	$welcome .= "&nbsp&nbsp&nbsp&nbsp   Nucleotides: ";
	my $seq_length = $coge->resultset('GenomicSequence')->get_column('sequence_length');
	$welcome .= commify($seq_length->sum);
	$welcome .= "</span>";
      }
    $template->param(BOX_NAME=>$welcome);
     $template->param(ADJUST_BOX=>1);
    $template->param(LOGO_PNG=>"CoGe-logo.png");
    $template->param(BODY=>gen_body());
    #$template->param(DEV=>1);
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $tmpl = HTML::Template->new(filename=>$P->{TMPLDIR}.'index.tmpl');
    my $html;
    my $disable = 1;
    $disable = 0 if $FORM->param('wheel');
    if ($update)
      {
	$tmpl->param(update=>1);
      }
    elsif ($USER && !$FORM->param('logout') && !$FORM->param('login'))
      {
	$tmpl->param(ACTIONS=>[map {{ACTION=>$_->{NAME}, DESC=>$_->{DESC}, LINK=>$_->{LINK}}} sort {$a->{ID} <=> $b->{ID}}@{actions()}  ]);
	$tmpl->param('INTRO'=>1);
	$tmpl->param('CREDITS'=>1);
      }
    my $url = $FORM->param('url') if $FORM->param('url');
    $url =~ s/:::/;/g if $url;
    $url = $URL.$url unless $url =~ /$URL/;
    $tmpl->param(url=>$url);

    if ($FORM->param('logout'))
      {
	my $session = md5_base64($USER->user_name.$ENV{REMOTE_ADDR});
	$session =~ s/\+/1/g;
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
		    ID=>5,
		    LOGO=>qq{<a href="./GEvo.pl"><img src="/CoGe/picts/carousel/GEvo-logo.png" width="227" height="75" border="0"></a>},
		    ACTION=>qq{<a href="./GEvo.pl">GEvo</a>},
		    LINK=>qq{./GEvo.pl},
		    DESC => qq{Compare sequences and genomic regions to discover patterns of genome evolution.  <a href ="GEvo.pl?prog=blastz;accn1=at1g07300;fid1=4091274;dsid1=556;chr1=1;dr1up=20000;dr1down=20000;gbstart1=1;gblength1=772;accn2=at2g29640;fid2=4113333;dsid2=557;chr2=2;dr2up=20000;dr2down=20000;gbstart2=1;rev2=1;num_seqs=2;autogo=1" target=_new>Example.</a>},
		    SCREENSHOT=>qq{<a href="./GEvo.pl"><img src="/CoGe/picts/preview/GEvo.png"border="0"></a>},
		    NAME=>"GEvo: High-resolution sequence analysis of genomic regions",
		   },
		   {
		    ID=>3,
		    LOGO=>qq{<a href="./FeatView.pl"><img src="/CoGe/picts/carousel/FeatView-logo.png" width="227" height="75" border="0"></a>},
		    ACTION=>qq{<a href="./FeatView.pl">FeatView</a>},
		    LINK=>qq{./FeatView.pl},
		    DESC => qq{Find and display information about a genomic feature (e.g. gene). <a href = "FeatView.pl?accn=at1g07300" target=_new>Example.</a>},
		    SCREENSHOT=>qq{<a href="./FeatView.pl"><img src="/CoGe/picts/preview/FeatView.png" width="400" height="241" border="0"></a>},
		    NAME=>qq{FeatView: Searching for genomic features by name},
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
		    ID=>1,
		    LOGO => qq{<a href="./OrganismView.pl"><img src="/CoGe/picts/carousel/OrganismView-logo.png" width="227" height="75" border="0"></a>},
		    ACTION => qq{<a href="./OrganismView.pl">OrganismView</a>},
		    LINK=>qq{./OrganismView.pl},
		    DESC   => qq{Search for organisms, get an overview of their genomic make-up, and visualize them using a dynamic, interactive genome browser. <a href="OrganismView.pl?org_name=k12" target=_new>Example.</a>},
		    SCREENSHOT => qq{<img src="/CoGe/picts/preview/OrganismView.png" border="0"></a>},
		    NAME=>qq{OrganismView:  Search for organisms and perform analyses on their genomes },
		   },
		   {
		    ID => 2,
		    LOGO => qq{<a href="./CoGeBlast.pl"><img src="/CoGe/picts/carousel/CoGeBlast-logo.png" width="227" height="75" border="0"></a>},
		    ACTION => qq{<a href="./CoGeBlast.pl">CoGeBlast</a>},
		    LINK=>qq{./CoGeBlast.pl},
		    DESC   => qq{Blast sequences against any number of organisms in CoGe.},
		    SCREENSHOT => qq{<a href="./CoGeBlast.pl"><img src="/CoGe/picts/preview/Blast.png" width="400" height="241" border="0"></a>},
		    NAME=>qq{CoGeBlast:  Blast sequences against any number of genomes of your choosing},
		   },
# 		   {
# 		    ID => 7,
# 		    LOGO => qq{<a href="./docs/help/CoGe"><img src="/CoGe/picts/carousel/FAQ-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./docs/help/CoGe/">CoGe Faq</a>},
# 		    DESC   => qq{What is CoGe?  This document covers some of the basics about what CoGe is, how it has been designed, and other information about the system.},
# 		    SCREENSHOT => qq{<a href="./docs/help/CoGe"><img src="/CoGe/picts/preview/app_schema.png" border="0"></a>},
# 		   }
		   {
		    ID=>4,
		    LOGO=>qq{<a href="./SynMap.pl"><img src="/CoGe/picts/SynMap-logo.png"  border="0"></a>},
		    ACTION=>qq{<a href="./SynMap.pl">SynMap</a>},
		    LINK=>qq{./SynMap.pl},
		    DESC => qq{Compare any two genomes to identify regions of synteny.  <a href="SynMap.pl?dsgid1=3068;dsgid2=8;D=20;g=10;A=5;w=0;b=1;ft1=1;ft2=1;dt=geneorder;ks=1;autogo=1" target=_mew>Example.</a>  <span class=small>(Powered by <a href=http://dagchainer.sourceforge.net/ target=_new>DAGChainer</a></span>)},
		    SCREENSHOT=>qq{<a href="./SynMap.pl"><img src="/CoGe/picts/preview/SynMap.png" border="0" width="400" height="320"></a>},
		    NAME=>qq{SynMap:  Whole genome syntenic dotplot anlayses},
		   },

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
	$session =~ s/\+/1/g;
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

sub commify 
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
