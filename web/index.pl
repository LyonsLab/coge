#!/usr/bin/perl -w

use AuthCAS;

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
use CGI::Log;



no warnings 'redefine';
use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $USER $FORM $DATE $URL $update $coge $COOKIE_NAME);

$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');
$ENV{PATH} = $P->{COGEDIR};
$URL = $P->{URL};
$FORM = new CGI;
#print STDERR Dumper $USER->user_name;

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_latest_genomes=>\&get_latest_genomes,
		      );
$update =0;

$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) =$FORM->param('ticket');
CoGe::Accessory::Web->login_cas(ticket=>$cas_ticket, coge=>$coge, this_url=>$FORM->url()) if($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name=>$COOKIE_NAME,coge=>$coge);
#logout is only called through this program!  All logouts from other pages are redirected to this page
CoGe::Accessory::Web->logout_cas(cookie_name=>$COOKIE_NAME, coge=>$coge, user=>$USER, form=>$FORM) if $FORM->param('logout');

#print $FORM->header, gen_html();
print $pj->build_html($FORM, \&gen_html);


sub gen_html
  {
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
    $template->param(TITLE=>'The Place to <span style="color: #119911">Co</span>mpare <span style="color: #119911">Ge</span>nomes');
    $template->param(PAGE_TITLE=>'ANKoCG');
    $template->param(HELP=>'/wiki/index.php');
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
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
	$welcome .= "&nbsp&nbsp&nbsp&nbsp   Genomic Features: ";
	$welcome .= commify($coge->resultset('Feature')->count());
	$welcome .= "&nbsp&nbsp&nbsp&nbsp   Annotations: ";
	$welcome .= commify($coge->resultset('Annotation')->count());

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
    if ($update)
      {
	$tmpl->param(update=>1);
      }
    elsif ($USER && !$FORM->param('logout') && !$FORM->param('login'))
      {
	$tmpl->param(ACTIONS=>[map {{ACTION=>$_->{NAME}, DESC=>$_->{DESC}, LINK=>$_->{LINK}}} sort {$a->{ID} <=> $b->{ID}}@{actions()}  ]);
	$tmpl->param('INTRO'=>1);
      }
    my $url = $FORM->param('url') if $FORM->param('url');
    if ($url)
     {
        $url =~ s/:::/;/g if $url;
        $tmpl->param(url=>$url);
     }
    $html .= $tmpl->output;
    return $html;
  }





sub actions
  {
    my @actions = (
		   {
		    ID=>5,
		    LOGO=>qq{<a href="./GEvo.pl"><img src="picts/carousel/GEvo-logo.png" width="227" height="75" border="0"></a>},
		    ACTION=>qq{<a href="./GEvo.pl">GEvo</a>},
		    LINK=>qq{./GEvo.pl},
		    DESC => qq{Compare sequences and genomic regions to discover patterns of genome evolution.  <a href ="GEvo.pl?prog=blastz;accn1=at1g07300;fid1=4091274;dsid1=556;chr1=1;dr1up=20000;dr1down=20000;gbstart1=1;gblength1=772;accn2=at2g29640;fid2=4113333;dsid2=557;chr2=2;dr2up=20000;dr2down=20000;gbstart2=1;rev2=1;num_seqs=2;autogo=1" target=_new>Example.</a>},
		    SCREENSHOT=>qq{<a href="./GEvo.pl"><img src="picts/preview/GEvo.png"border="0"></a>},
		    NAME=>"GEvo: High-resolution sequence analysis of genomic regions",
		   },
		   {
		    ID=>3,
		    LOGO=>qq{<a href="./FeatView.pl"><img src="picts/carousel/FeatView-logo.png" width="227" height="75" border="0"></a>},
		    ACTION=>qq{<a href="./FeatView.pl">FeatView</a>},
		    LINK=>qq{./FeatView.pl},
		    DESC => qq{Find and display information about a genomic feature (e.g. gene). <a href = "FeatView.pl?accn=at1g07300" target=_new>Example.</a>},
		    SCREENSHOT=>qq{<a href="./FeatView.pl"><img src="picts/preview/FeatView.png" width="400" height="241" border="0"></a>},
		    NAME=>qq{FeatView: Searching for genomic features by name},
		   },
# 		   {
# 		    ID=>3,
# 		    LOGO=>qq{<a href="./MSAView.pl"><img src="picts/carousel/MSAView-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./MSAView.pl">MSAView: Multiple Sequence Alignment Viewer</a>},
# 		    DESC   => qq{Allows users to submit a multiple sequence alignment in FASTA format (if people would like additional formats, please request via e-mail) in order to quickly check the alignment, find conserved regions, etc.  This program also generates a consensus sequence from the alignment and displays some basic statistics about the alignment.},
# 		    SCREENSHOT=>qq{<a href="./MSAView.pl"><img src="picts/preview/MSAView.png"border="0"></a>},
# 		   },
# 		   {
# 		    ID=>4,
# 		    LOGO=>qq{<a href="./TreeView.pl"><img src="picts/carousel/TreeView-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./TreeView.pl">TreeView: Phylogenetic Tree Viewer</a>},
# 		    DESC   => qq{Allows users to submit a tree file and get a graphical view of their tree.  There is support for drawing rooted and unrooted trees, zooming and unzooming functions, and coloring and shaping nodes based on user specifications.},
# 		    SCREENSHOT=>qq{<a href="./FeatView.pl"><img src="picts/preview/TreeView.png"border="0"></a>},
# 		   },
		   {
		    ID=>1,
		    LOGO => qq{<a href="./OrganismView.pl"><img src="picts/carousel/OrganismView-logo.png" width="227" height="75" border="0"></a>},
		    ACTION => qq{<a href="./OrganismView.pl">OrganismView</a>},
		    LINK=>qq{./OrganismView.pl},
		    DESC   => qq{Search for organisms, get an overview of their genomic make-up, and visualize them using a dynamic, interactive genome browser. <a href="OrganismView.pl?org_name=k12" target=_new>Example.</a>},
		    SCREENSHOT => qq{<img src="picts/preview/OrganismView.png" border="0"></a>},
		    NAME=>qq{OrganismView:  Search for organisms and perform analyses on their genomes },
		   },
		   {
		    ID => 2,
		    LOGO => qq{<a href="./CoGeBlast.pl"><img src="picts/carousel/CoGeBlast-logo.png" width="227" height="75" border="0"></a>},
		    ACTION => qq{<a href="./CoGeBlast.pl">CoGeBlast</a>},
		    LINK=>qq{./CoGeBlast.pl},
		    DESC   => qq{Blast sequences against any number of organisms in CoGe.},
		    SCREENSHOT => qq{<a href="./CoGeBlast.pl"><img src="picts/preview/Blast.png" width="400" height="241" border="0"></a>},
		    NAME=>qq{CoGeBlast:  Blast sequences against any number of genomes of your choosing},
		   },
# 		   {
# 		    ID => 7,
# 		    LOGO => qq{<a href="./docs/help/CoGe"><img src="picts/carousel/FAQ-logo.png" width="227" height="75" border="0"></a>},
# 		    ACTION => qq{<a href="./docs/help/CoGe/">CoGe Faq</a>},
# 		    DESC   => qq{What is CoGe?  This document covers some of the basics about what CoGe is, how it has been designed, and other information about the system.},
# 		    SCREENSHOT => qq{<a href="./docs/help/CoGe"><img src="picts/preview/app_schema.png" border="0"></a>},
# 		   }
		   {
		    ID=>4,
		    LOGO=>qq{<a href="./SynMap.pl"><img src="picts/SynMap-logo.png"  border="0"></a>},
		    ACTION=>qq{<a href="./SynMap.pl">SynMap</a>},
		    LINK=>qq{./SynMap.pl},
		    DESC => qq{Compare any two genomes to identify regions of synteny.  <a href="SynMap.pl?dsgid1=3068;dsgid2=8;D=20;g=10;A=5;w=0;b=1;ft1=1;ft2=1;dt=geneorder;ks=1;autogo=1" target=_mew>Example.</a>  <span class=small>(Powered by <a href=http://dagchainer.sourceforge.net/ target=_new>DAGChainer</a></span>)},
		    SCREENSHOT=>qq{<a href="./SynMap.pl"><img src="picts/preview/SynMap.png" border="0" width="400" height="320"></a>},
		    NAME=>qq{SynMap:  Whole genome syntenic dotplot anlayses},
		   },

		  );
    return \@actions;
  }

sub get_latest_genomes
  {
    my %opts = @_;
    my $limit = $opts{limit} || 20;
    my @db = $coge->resultset("DatasetGroup")->search({},
						      {
						       distinct=>"organism.name",
						       join=>"organism",
						       prefetch=>"organism",
						       order_by=>"dataset_group_id desc",
						       rows=>$limit,
						      }
						);
  #  ($USER) = CoGe::Accessory::LogUser->get_user();
    my $html = "<table class=small>";
    $html .= "<tr><th>".join("<th>",qw(Organism  &nbsp Length&nbsp(nt) &nbsp Related Link ));
    my @opts;
    my %org_names;
    foreach my $dsg (@db)
      {
	next if $dsg->restricted && $USER->has_access_to_genome(dsg);
	next if $org_names{$dsg->organism->name};
	$org_names{$dsg->organism->name}=1;
	my $orgview_link = "OrganismView.pl?oid=".$dsg->organism->id;
	my $entry = qq{<tr>};
#	$entry .= qq{<td><span class='ui-button ui-corner-all' onClick="window.open('$orgview_link')"><span class="ui-icon ui-icon-link"></span>&nbsp&nbsp</span>};

	$entry .= qq{<td><span class="link" onclick=window.open('$orgview_link')>};
	my $name = $dsg->organism->name;
	$name = substr($name,0,40)."..." if length($name) > 40;
	$entry.=$name;
	$entry .= qq{</span>};

#	$entry .= ": ".$dsg->name if $dsg->name;
	$entry .= "<td>(v".$dsg->version.")&nbsp";
	$entry .= "<td align=right>".commify($dsg->length)."<td>";
	my @desc = split/;/,$dsg->organism->description;
	while ($desc[0] && !$desc[-1]) {pop @desc;}
	$desc[-1] =~ s/^\s+//;
	$desc[-1] =~ s/\s+$//;
	my $orgview_search = "OrganismView.pl?org_desc=".$desc[-1];
	$entry .= qq{<td><span class="link" onclick="window.open('$orgview_search')">Search</span>};
	$entry .= qq{<td>};
	$entry .= qq{<img onClick="window.open('$orgview_link')" src = "picts/other/CoGe-icon.png" title="CoGe" class=link>};

	my $search_term = $dsg->organism->name;
	$entry .= qq{<img onclick="window.open('http://www.ncbi.nlm.nih.gov/taxonomy?term=$search_term')" src = "picts/other/NCBI-icon.png" title="NCBI" class=link>};
	$entry .= qq{<img onclick="window.open('http://en.wikipedia.org/w/index.php?title=Special%3ASearch&search=$search_term')" src = "picts/other/wikipedia-icon.png" title="Wikipedia" class=link>};
	$search_term =~ s/\s+/\+/g;
	$entry .= qq{<img onclick="window.open('http://www.google.com/search?q=$search_term')" src="picts/other/google-icon.png" title="Google" class=link>};
	$entry .= qq{</tr>};
	push @opts, $entry;#, "<OPTION value=\"".$item->organism->id."\">".$date." ".$item->organism->name." (id".$item->organism->id.") "."</OPTION>";
      }
    $html .= join "\n", @opts;
    $html .= "</table>";
    return $html;
    }

sub commify 
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
