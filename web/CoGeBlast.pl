#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use CoGe::Genome;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use LWP::Simple;
use LWP::Simple::Post qw(post post_xml);
use URI::Escape;
use POSIX;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $DB %restricted_orgs);

$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
if (!$USER || $USER =~ /public/i)
  {
    $restricted_orgs{papaya} = 1;
  }
$DB = new CoGe::Genome;
$FORM = new CGI;

#my %ajax = CoGe::Accessory::Web::ajax_func();

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seq=>\&get_seq,
		       cogesearch=>\&cogesearch,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       get_sequence=>\&get_sequence,
		       get_url=>\&get_url,
		       check_seq=>\&check_seq,
		       set_seq=>\&set_seq,
		       find_radio=>\&find_radio,
		       #%ajax,
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();


sub gen_html
  {
    my ($body) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'CoGe: BLAST');
    $template->param(HELP=>'BLAST');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"CoGeBlast-logo.png");
    #$template->param(BOX_NAME=>$title);
    $template->param(BODY=>$body);
    my $html;
    $html .= $template->output;
  }

sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    my $form = $FORM;
    my $featid = $form->param('featid') || 0;
    my $upstream = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    if ($featid)
     {
       #print STDERR "featid: ", $featid, ", up: ", $upstream, ", down: ", $downstream, "\n";
       my $seq = get_sequence($featid, $upstream, $downstream, 1);
       $template->param(SEQUENCE=>$seq);
     }

    $template->param(BLAST_FRONT_PAGE=>1);
    $template->param(TYPE_LOOP=> [{TYPE=>"<OPTION VALUE=0>All</OPTION>"},map {{TYPE=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} $DB->all_feature_types]);
    my @orgs;
    foreach my $org ($DB->all_orgs)
      {
	push @orgs, $org unless $restricted_orgs{$org->name};
      }
    $template->param(ORG_LOOP=> [{ORG=>"<OPTION VALUE=0>All</OPTION>"},map {{ORG=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} @orgs]);
    $template->output;
  }
  
sub get_types
  {
    my ($accn, $dsid) = @_;
    my $html;
    my $blank = qq{<input type="hidden" id="Type_name">};
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {$_->type->name} $DB->get_features_by_name_and_dataset_id(name=>$accn, id=>$dsid);

    $html .= "<font class=small>Type count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_sequence(['accn_select','Type_name', 'dsid'],['seq_box'])" >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $blank unless $html =~ /OPTION/;
    return ($html, 1);
  }


sub cogesearch
  {
    my $accn = shift;
    my $type = shift;
    my $org = shift;
    my $blank = qq{<input type="hidden" id="accn_select">};
#    print STDERR Dumper @_;
    return $blank unless length($accn) > 2 || $type || $org;
    my $html;
    my %seen;
    my @opts = sort map {"<OPTION>$_</OPTION>"} grep {! $seen{$_}++} map {uc($_)} $DB->get_feature_name_obj->power_search(accn=>$accn."%", type=>$type, org=>$org);
    if (@opts > 5000)
      {
	return $blank."Search results over 1000, please refine your search.\n";
      }
    $html .= "<font class=small>Name count: ".scalar @opts."</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select" SIZE="5" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank."No results found.\n" unless $html =~ /OPTION/;
    return $html;
  }
  
sub get_data_source_info_for_accn
  {
    my $accn = shift;
    my $blank = qq{<input type="hidden" id="dsid">};
    return $blank unless $accn;
    my @feats = $DB->get_feats_by_name($accn);
    my %sources;
    foreach my $feat (@feats)
      {
	my $val = $feat->dataset;
	my $name = $val->name;
	my $ver = $val->version;
	my $desc = $val->description;
	my $sname = $val->data_source->name;
	my $ds_name = $val->name;
	my $org = $val->org->name;
	my $title = "$org: $ds_name ($sname, v$ver)";
	#next if $restricted_orgs{$org};
	$sources{$title} = $val->id;
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

sub get_sequence
  {
    my $accn = shift;
    my $type = shift;
    my $dsid = shift;
    my $seqview = shift;
    my $seq;
    my $featid;
    my $upstream = 0;
    my $downstream = 0;
    #print STDERR "type: ", $type, ", dsid: ", $dsid, ", accn: ", $accn, "\n";
    unless ($seqview)
    {
     my @feats;
     foreach my $feat ($DB->get_features_by_name_and_dataset_id(name=>$accn, id=>$dsid))
      {
	push @feats, $feat if ($feat->type->name eq $type);
      }
     foreach my $feat (@feats) 
     {
      $featid = $feat->id;
     }
    }
    my $fasta;
    #print STDERR "dsid: ", $dsid, "\n";
    unless ($seqview) 
    {
     $fasta = generate_fasta(featid=>$featid, dsid=>$dsid, feat_name=>$accn);
    }
    #print STDERR "\nfeatid: ", $featid, ", up: ", $upstream, ", down: ", $downstream, "\n";
    if ($seqview)
     {
      my ($feat) = $DB->get_feat_obj->retrieve($accn);
      unless (ref($feat) =~ /Feature/i)
	{
	 return "Unable to retrieve Feature object for id: $accn $type";
	}
      else {	  
      $seq = $feat->genomic_sequence(upstream=>$type, downstream=>$dsid);
      $seq = join ("\n", wrap('','',$seq));
      $seq = ($fasta. $seq);
     }
    return $seq;
    }
    else
     {
      my ($feat) = $DB->get_feat_obj->retrieve($featid);
         unless (ref($feat) =~ /Feature/i)
	  {
	    return "Unable to retrieve Feature object for id: $accn $type";
	  }
     else {	  
      $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);
      $seq = join ("\n", wrap('','',$seq));
      $seq = ($fasta. $seq);
     }
     return $seq;
    }
    #print STDERR "accn select: ",$featid, ", type name: ",$type, ", dsid: ", $dsid, "\n";
  }
  
sub check_seq
  {
    my $blast_type = shift;
    my $accn = shift;
    my $dsid = shift;
    my $type = shift;
    #print STDERR "type: ", $type, ", dsid: ", $dsid, ", accn: ", $accn, "\n";
    my $featid;
    my @feats;
    foreach my $feat ($DB->get_features_by_name_and_dataset_id(name=>$accn, id=>$dsid))
      {
	push @feats, $feat if ($feat->type->name eq $type);
      }
    foreach my $feat (@feats) 
    {
      $featid = $feat->id;
    }
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
     #print STDERR "featid: ", $featid, "\n";
    my $seq;
    my $fasta = generate_fasta(featid=>$featid, dsid=>$dsid, feat_name=>$accn);
    #print STDERR "blast_type: ", $blast_type, ", blast_p: ", $blast_p, "\n";
    if ($blast_type == 1) {
      ($seq) = $DB->get_protein_sequence_for_feature($feat);
      #print STDERR "AAs: ", $seq, "\n";
      $seq = "No sequence available" unless $seq;
    }
    else {
      $seq = $feat->genomic_sequence();
    }
    $seq = join ("\n", wrap('','',$seq));
    $seq = ($fasta. $seq);
    #print STDERR "Seq is: ", $seq, "\n";
    return $seq;
  }

sub get_url
  {
    my $url = shift;
    if ($url eq "megablast") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&AUTO_FORMAT=Semiauto&CLIENT=web&DATABASE=Test%2Fgpipe%2F9606%2Fallcontig_and_rna&DESCRIPTIONS=100&ENTREZ_QUERY=All+organisms&EXPECT=10&FILTER=L&FILTER=R&FORMAT_BLOCK_ON_RESPAGE=None&FORMAT_ENTREZ_QUERY=All+organisms&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&GET_SEQUENCE=on&HITLIST_SIZE=100&LAYOUT=TwoWindows&MASK_CHAR=2&MASK_COLOR=1&MEGABLAST=on&MYNCBI_USER=4389532342&NEW_VIEW=on&NUM_OVERVIEW=100&PAGE=MegaBlast&PERC_IDENT=None,+1,+-2&PROGRAM=blastn&REPEATS=repeat_9606&SEARCH_NAME=mb&SERVICE=plain&SET_DEFAULTS=Yes&SET_DEFAULTS.x=47&SET_DEFAULTS.y=12&SHOW_LINKOUT=on&SHOW_OVERVIEW=on&USER_TYPE=2&WORD_SIZE=28&WWW_BLAST_TYPE=mapview&dbtype=hc&END_OF_HTTPGET=Yes&QUERY=";}
    elsif ($url eq "blastn") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&CLIENT=web&DESCRIPTIONS=100&ENTREZ_QUERY=%28none%29&EXPECT=10&FILTER=L&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&NCBI_GI=off&PAGE=Nucleotides&PROGRAM=blastn&SERVICE=plain&SET_DEFAULTS.x=34&SET_DEFAULTS.y=8&SHOW_OVERVIEW=on&END_OF_HTTPGET=yes&GET_SEQUENCE=yes&NEW_VIEW=yes&SEARCH_NAME=bn&QUERY=";}
    elsif ($url eq "blastp") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=250&ALIGNMENT_VIEW=Pairwise&CDD_SEARCH=on&CLIENT=web&COMPOSITION_BASED_STATISTICS=on&DATABASE=nr&DESCRIPTIONS=500&ENTREZ_QUERY=%28none%29&EXPECT=10&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&I_THRESH=0.005&MATRIX_NAME=BLOSUM62&NCBI_GI=on&PAGE=Proteins&PROGRAM=blastp&SERVICE=plain&SET_DEFAULTS.x=41&SET_DEFAULTS.y=5&SHOW_OVERVIEW=on&END_OF_HTTPGET=Yes&SHOW_LINKOUT=yes&GET_SEQUENCE=yes&SEARCH_NAME=bp&QUERY=";}
      elsif ($url eq "blastx") {
        return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&CLIENT=web&DATABASE=nr&DESCRIPTIONS=100&ENTREZ_QUERY=%28none%29&EXPECT=10&FILTER=L&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&GENETIC_CODE=1&NCBI_GI=on&PAGE=Translations&PROGRAM=blastx&SERVICE=plain&SET_DEFAULTS.x=37&SET_DEFAULTS.y=5&SHOW_OVERVIEW=on&UNGAPPED_ALIGNMENT=no&END_OF_HTTPGET=Yes&SHOW_LINKOUT=yes&GET_SEQUENCE=yes&SEARCH_NAME=bx&QUERY=";}
    elsif ($url eq "tblastn") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&CLIENT=web&DATABASE=nr&DESCRIPTIONS=100&ENTREZ_QUERY=%28none%29&EXPECT=10&FILTER=L&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&GENETIC_CODE=0&NCBI_GI=on&PAGE=Translations&PROGRAM=tblastn&SERVICE=plain&SET_DEFAULTS.x=23&SET_DEFAULTS.y=10&SHOW_OVERVIEW=on&UNGAPPED_ALIGNMENT=no&END_OF_HTTPGET=Yes&SHOW_LINKOUT=yes&GET_SEQUENCE=yes&SEARCH_NAME=tbn&QUERY=";}
    elsif ($url eq "tblastx") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&CLIENT=web&DATABASE=nr&DESCRIPTIONS=100&ENTREZ_QUERY=%28none%29&EXPECT=10&FILTER=L&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&GENETIC_CODE=1&NCBI_GI=on&PAGE=Translations&PROGRAM=tblastx&SERVICE=plain&SET_DEFAULTS.x=21&SET_DEFAULTS.y=9&SHOW_OVERVIEW=on&UNGAPPED_ALIGNMENT=yes&END_OF_HTTPGET=Yes&SHOW_LINKOUT=yes&SEARCH_NAME=tbx&QUERY=";}
    elsif ($url eq "psiblast") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=250&ALIGNMENT_VIEW=Pairwise&CLIENT=web&COMPOSITION_BASED_STATISTICS=on&DATABASE=nr&CDD_SEARCH=on&DESCRIPTIONS=500&ENTREZ_QUERY=%28none%29&EXPECT=10&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&I_THRESH=0.005&MATRIX_NAME=BLOSUM62&NCBI_GI=on&PAGE=Proteins&PROGRAM=blastp&RUN_PSIBLAST=on&SERVICE=plain&SET_DEFAULTS.x=36&SET_DEFAULTS.y=5&SHOW_OVERVIEW=on&END_OF_HTTPGET=Yes&SHOW_LINKOUT=yes&GET_SEQUENCE=yes&SEARCH_NAME=psi_bp&QUERY=";}
    else {
      return 1;}
  }
  
# sub aa_or_nu_seq 
# {
#   my $blast = shift;
#   print STDERR "blast: ", $blast, "\n";
#   my %blasts = (megablast=>0,
#   		blastn=>0,
#   		blastp=>1,
#   		blastx=>0,
#   		tblastn=>1,
#   		tblastx=>0,
#   		psiblast=>1
#   			   );
#   my $type = $blasts{$blast};
#   print STDERR "type: ", $type, "\n";
#   return $type;
# }

sub generate_fasta
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my $dsid = $opts{dsid};
    my $feat_name = $opts{feat_name};
    my $ds = $DB->get_dataset_obj->retrieve($dsid);
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    my $fasta = ">".$ds->org->name."(v.".$feat->version.") ".", Name: ".$feat_name.", Type: ".$feat->type->name.", Location: ".$feat->genbank_location_string.", Chromosome: ".$feat->chr.", Strand: ".$feat->strand."\n";
    return $fasta;
  }
  
# sub ncbi_blast
#   {
#     my $blank = shift;
#     my $url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&LAYOUT=TwoWindows&AUTO_FORMAT=Semiauto&ALIGNMENTS=50&ALIGNMENT_VIEW=Pairwise&CLIENT=web&DESCRIPTIONS=100&ENTREZ_QUERY=%28none%29&EXPECT=10&FILTER=L&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&NCBI_GI=off&PAGE=Nucleotides&PROGRAM=blastn&SERVICE=plain&SET_DEFAULTS.x=34&SET_DEFAULTS.y=8&SHOW_OVERVIEW=on&END_OF_HTTPGET=yes&GET_SEQUENCE=yes&NEW_VIEW=yes&SEARCH_NAME=bn";
#     my $seq = shift;
#     #print STDERR "url: ", $url, "\nseq: ", $seq, "\n";
#     my $html = post($url, 'QUERY=ataboy');
#     #print STDERR "html: ", $html, "\n";
#   }

sub find_radio
  {
    my %opts = @_;
    my $radio = $opts{'bn'};
    return $radio;
  }