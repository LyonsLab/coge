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
		       get_pretty_print=>\&get_pretty_print,
		       blast_param=>\&blast_param,
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
    my $chr = $form->param('chr') || 0;
    my $upstream = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $dsid = $form->param('dsid') || 0;
    my $feat_name = $form->param('featname');
    $template->param(JAVASCRIPT=>1);
    $template->param(BLAST_FRONT_PAGE=>1);
    $template->param(TYPE_LOOP=> [{TYPE=>"<OPTION VALUE=0>All</OPTION>"},map {{TYPE=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} $DB->all_feature_types]);
    my @orgs;
    foreach my $org ($DB->all_orgs)
      {
	push @orgs, $org unless $restricted_orgs{$org->name};
      }
    $template->param(ORG_LOOP=> [{ORG=>"<OPTION VALUE=0>All</OPTION>"},map {{ORG=>"<OPTION value=\"".$_->id."\">".$_->name."</OPTION>"}} sort {uc($a->name) cmp uc($b->name)} @orgs]);
    #my $param = blast_param('blast_type_n');
    #print STDERR $param;
    if ($featid)
    {
    	$template->param(FEATID=>$featid);
    	$template->param(UPSTREAM=>$upstream);
    	$template->param(DOWNSTREAM=>$downstream);
    	$template->param(DSID=>$dsid);
    	$template->param(FEATNAME=>$feat_name);
    	$template->param(SEQVIEW=>1);
    	my $seq = get_sequence($feat_name,$featid, $dsid, 0, 1, $upstream, $downstream);
        $template->param(SEQUENCE=>$seq);
    }
    if ($chr)
    {
    	$template->param(CHR=>$chr);
    	$template->param(UPSTREAM=>$upstream);
    	$template->param(DOWNSTREAM=>$downstream);
    	$template->param(DSID=>$dsid);
    	$template->param(SEQVIEW=>2);
    	my $seq = get_sequence(0,$chr, $dsid, 0, 2, $upstream, $downstream);
    	$template->param(SEQUENCE=>$seq);
    }
    
    
    $template->param(USER_NAME=>$USER);
    #$template->param(DEFAULT_PARAM=>$param);
    $template->param(REST=>1);
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
    $html .= qq{<SELECT id="Type_name" SIZE="5" MULTIPLE onChange="get_sequence(['accn_select','Type_name', 'dsid'],['seq_box']);get_pretty_print(['accn_select','Type_name', 'dsid'],['pretty_print_anno']);check_seq_for_blast();" >\n};
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
    my $blast_type = shift || 0;
    my $seqview = shift || 0;
    my $upstream = shift || 0;
    my $downstream = shift || 0;
    my $seq;
    my $featid;
    #print STDERR "type: ", $type, ", dsid: ", $dsid, ", accn: ", $accn, ", blast_type: ", $blast_type, "\n";
    my $fasta;
    if ($seqview == 2)
    {
    	$seq = $DB->get_genomic_sequence(start=>$upstream,
					 stop=>$downstream,
					 chr=>$type,
					 dataset_id=>$dsid);
	$fasta = generate_fasta_without_featid(chr=>$type, dsid=>$dsid, start=>$upstream, stop=>$downstream);
	$seq = join ("\n", wrap('','',$seq));
	$seq = ($fasta. $seq);
	return $seq;
    }
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
    else
     {$featid = $type;}
    #print STDERR "dsid: ", $dsid, "\n";
    #print STDERR "\nfeatid: ", $featid, ", up: ", $upstream, ", down: ", $downstream, "\n";
     $fasta = generate_fasta_with_featid(featid=>$featid, dsid=>$dsid, feat_name=>$accn);
    
#     if ($seqview)
#      {
#       my ($feat) = $DB->get_feat_obj->retrieve($featid);
#       unless (ref($feat) =~ /Feature/i)
# 	{
# 	 return "Unable to retrieve Feature object for id: $accn";
# 	}
#       else {	  
#       $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);
#       $seq = join ("\n", wrap('','',$seq));
#       $seq = ($fasta. $seq);
#      }
#     return $seq;
#     }
#     else
#      {

    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $accn $type";
    }
    else 
    {	  
      if ($blast_type eq  "blast_type_p") {
       ($seq) = $DB->get_protein_sequence_for_feature($feat);
       $seq = "No sequence available" unless $seq;
      }
      else {
       $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);}
      $seq = join ("\n", wrap('','',$seq));
      $seq = ($fasta. $seq);
    }
    return $seq;
    #print STDERR "accn select: ",$featid, ", type name: ",$type, ", dsid: ", $dsid, "\n";
  }
  
# sub check_seq
#   {
#     my $blast_type = shift;
#     my $accn = shift;
#     my $dsid = shift;
#     my $type = shift;
#     my $seqview
#     #print STDERR "Blast_type: ", $blast_type, ", Accn : ", $accn, ", dsid: ", $dsid, ", type is: ", $type, "\n";
#     #print STDERR "Hello, World!\n";
#     #print STDERR "type is: ", $type, ", dsid: ", $dsid, ", accn: ", $accn, "\n";
#     my $featid;
#     my @feats;
#     foreach my $feat ($DB->get_features_by_name_and_dataset_id(name=>$accn, id=>$dsid))
#       {
# 	push @feats, $feat if ($feat->type->name eq $type);
#       }
#     foreach my $feat (@feats) 
#     {
#       $featid = $feat->id;
#     }
#     my ($feat) = $DB->get_feat_obj->retrieve($featid);
#      #print STDERR "featid: ", $featid, "\n";
#     my $seq;
#     my $fasta = generate_fasta(featid=>$featid, dsid=>$dsid, feat_name=>$accn);
#     #print STDERR "blast_type: ", $blast_type, ", blast_p: ", $blast_p, "\n";
#     if ($blast_type eq  "blast_type_p") {
#       ($seq) = $DB->get_protein_sequence_for_feature($feat);
#       #print STDERR "AAs: ", $seq, "\n";
#       $seq = "No sequence available" unless $seq;
#     }
#     else {
#       $seq = $feat->genomic_sequence();
#     }
#     $seq = join ("\n", wrap('','',$seq));
#     $seq = ($fasta. $seq);
#     #print STDERR "Seq is: ", $seq, "\n";
#     return $seq;
#   }

sub get_url
  {
    my $url = shift;
    my $expect = shift;
    my $db = shift;
    #print STDERR "expect: ", $expect, "\n";
    if ($url eq "blastn") {
      return "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=on&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&EXPECT=$expect&DATABASE=$db";}
    elsif ($url eq "blastp") {
      return "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&EXPECT=$expect&DATABASE=$db";}
      elsif ($url eq "blastx") {
        return "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=blastx&BLAST_PROGRAMS=blastx&PAGE_TYPE=BlastSearch&EXPECT=$expect&DATABASE=$db";}
    elsif ($url eq "tblastn") {
      return "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastn&BLAST_PROGRAMS=tblastn&PAGE_TYPE=BlastSearch&EXPECT=$expect&DATABASE=$db";}
    elsif ($url eq "tblastx") {
      return "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastx&BLAST_PROGRAMS=tblastx&PAGE_TYPE=BlastSearch&EXPECT=$expect&DATABASE=$db";}
    else {
      return 1;}
  }
#   
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

sub generate_fasta_with_featid
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
  
sub generate_fasta_without_featid
  {
    my %opts = @_;
    my $chr = $opts{chr};
    my $dsid = $opts{dsid};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $ds = $DB->get_dataset_obj->retrieve($dsid);
    my $fasta = ">".$ds->org->name.", Location: ".$start."-".$stop.", Chromosome: ".$chr."\n";
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
    my $expect = $opts{'expect'};
    my $db = $opts{'db'};
    #print STDERR "radio: ",$radio,", expect: ",$expect, ", db: ",$db,"\n";
    return $radio;
  }

sub get_pretty_print
{
    my $accn = shift;
    my $type = shift;
    my $dsid = shift;
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
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $accn $type";
    }
    my $html = $feat->annotation_pretty_print_html();
    return $html;
}

sub blast_param
{
    my $seq_type = shift || "blast_type_n";
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    if ($seq_type eq "blast_type_n") {
        $template->param(BLAST_NU=>1);}
    else {
	$template->param(BLAST_PRO=>1);}
    my $html = $template->output;
    return $html;
}
