#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use CoGe::Genome;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $DB );

$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;
$FORM = new CGI;

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seq=>\&get_seq,
		       get_extend_box=>\&get_extend_box,
		       gen_title=>\&gen_title,
			);
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();
sub gen_html
  {
    #print "DNA Sequence Button Doesn't Work. Doesn't change Box title when button clicked. Need to Hide buttons relating to upstream/downstream changes. Need to make extend button disappear text boxes when clicked again. Deal with textarea/copy text issue; CoGe title.";
    my ($body) = gen_body();
    my $rc;
    my $pro; 
    my ($title) = gen_title(protein=>$pro, rc=>$rc);
    my $form = $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    #$template->param(TITLE=>'CoGe: Sequence Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"SeqView-logo.png");
    $template->param(BOX_NAME=>qq{<DIV id="box_name">$title</DIV>});
    $template->param(BODY=>$body);
    $template->param(POSTBOX=>gen_foot());
    $template->param(CLOSE=>1);
    $template->param(HEAD=>qq{<script src="js/kaj.stable.js"></script>});
    #print STDERR gen_foot()."\n";
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = $FORM;
    my $feat_id = $form->param('featid');
    my $chr = $form->param('chr');
    my $dsid = $form->param('dsid');
    my $feat_name = $form->param('name');
    my $rc = $form->param('rc');
    my $pro = $form->param('pro');   
    my $upstream = $form->param('upstream');
    my $downstream = $form->param('downstream');
    my $start = $form->param('start');
    my $stop = $form->param('stop');
    my $seq;
    unless ($feat_id)
    {
      my $strand = $form->param('strand');
      $strand = check_strand(strand=>$strand, rc=>$rc);
      $seq = get_seq(pro=>$pro,
		      rc=>$rc,
		      chr=>$chr,
		      dsid=>$dsid,
		      strand=>$strand, 
		      upstream=>$upstream, 
		      downstream=>$downstream,
		      start=>$start,
		      stop=>$stop,
		      );
    }
    else
    {
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    my $strand = $feat->strand;
    $strand = check_strand(strand=>$strand, rc=>$rc);
    $seq = get_seq(featid=>$feat_id, 
		      pro=>$pro,
		      rc=>$rc,
		      chr=>$chr,
		      dsid=>$dsid,
		      featname=>$feat_name,
		      strand=>$strand, 
		      upstream=>$upstream, 
		      downstream=>$downstream,
		      );
    }
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    $template->param(TEXT=>1);
    $template->param(SEQ=>$seq);
    #$seq = qq{<TABLE style="width: 628px; height: 300px; overflow: auto;"><TR><TD align=left valign="top">$seq</TD></TR></TABLE>};
    return qq{<DIV id="seq" style="width: 623px; height: 300px; overflow: auto;">$seq</div>};
  }
 
sub check_strand
{
    my %opts = @_;
    my $strand = $opts{'strand'};
    my $rc = $opts{'rc'};
    if ($rc)
      {
        if ($strand =~ /-/)
          {
            $strand = "1";
          }
        else
          {
            $strand = "-1";
          }
      }
    return $strand;
}

sub get_seq
  {
    my %opts = @_;
    my $feat_id = $opts{'featid'};
    my $pro = $opts{'pro'};
    my $rc = $opts{'rc'};
    my $chr = $opts{'chr'};
    my $dsid = $opts{'dsid'};
    my $feat_name = $opts{'featname'};
    my $upstream = $opts{'upstream'};
    my $downstream = $opts{'downstream'};
    my $strand = $opts{'strand'};
    my $start = $opts{'start'};
    my $stop = $opts{'stop'};
    my $ds = $DB->get_dataset_obj->retrieve($dsid);
    my $seq;
    my $fasta;
    unless ($feat_id)
    {
      $fasta = ">".$ds->org->name.", Location: ".$start.": ".$stop.", Chromosome: ".$chr.", Strand: ".$strand."\n";
      $fasta = qq{<FONT class="main"><i>$fasta</i></FONT>};
      $columns = 80;
      $fasta = join ("\n", wrap('','',$fasta));
    }
    else
    {
    #print Dumper \%opts;
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    $fasta = ">".$ds->org->name."(v.".$feat->version."), Type: ".$feat->type->name.", Location: ".$feat->genbank_location_string.", Chromosome: ".$chr.", Strand: ".$strand.", Name: ".$feat_name."\n";
    $fasta = qq{<FONT class="main"><i>$fasta</i></FONT>};
    $columns = 80;
    $fasta = join ("\n", wrap('','',$fasta));
    }
    unless ($pro)
    {
      $seq .= get_dna_seq_for_feat (featid=>$feat_id,
      				    dsid=>$dsid, 
      				    rc=>$rc, 
      				    upstream=>$upstream,
      				    downstream=>$downstream, 
      				    start=>$start,
      				    stop=>$stop,
      				    chr=>$chr);
    }
    else
    {
      $seq .= get_prot_seq_for_feat($feat_id);
    }
    #print length($seq);
    my $up = upstream_color(seq=>$seq, upstream=>$upstream);
    my $down = downstream_color(seq=>$seq, downstream=>$downstream);
    my $main = main_color(seq=>$seq, upstream=>$upstream, downstream=>$downstream);
    $seq = join("", $up, $main, $down);
    $seq = join("", $fasta, $seq);
    return $seq;
  }
  
sub gen_foot
  {
    my $form = $FORM;
    my $feat_id = $form->param('featid');
    my $chr = $form->param('chr');
    my $dsid = $form->param('dsid');
    my $feat_name = $form->param('name');
    my $rc = $form->param('rc');
    my $pro = $form->param('pro');
    my $upstream = $form->param('upstream');
    my $downstream = $form->param('downstream');
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    my $strand;
    unless ($feat_id)
    {$strand = $form->param('strand');}
    else
    {$strand = $feat->strand;}
    $strand = check_strand(strand=>$strand, rc=>$rc);
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    
#   my $button = qq{<input type="submit" value="$name">};
    #my $url = $form->url(-full=>1, -query=>1);
    
    #$url =~ s/rc=.//g;
    #$url =~ s/pro=.//g;
    #my $url2 = $url;
    #my $url3 = $url;
    #$url .= "&rc=0";
    #$url2 .="&rc=1";
    #$url3 .= "&pro=1";

    $template->param(BOTTOM_BUTTONS=>1);
    $template->param(BUTTON_LOOP=>[
                            	   {BUTTON_NUM=>1, 
                                    BUTTON_NAME=>'DNA Sequence', 
                                    FEATID=>$feat_id,
    			            PRO=>0,
    				    RC=>0,
    			     	    CHR=>$chr,
    				    DSID=>$dsid,
    				    FEATNAME=>$feat_name,
   				    STRAND=>$strand,
   				    RC=>0,
    				    PRO=>0},
                           	   {BUTTON_NUM=>2, 
                            	    BUTTON_NAME=>'Reverse Complement', 
                            	    FEATID=>$feat_id,
   				    PRO=>0,
    				    RC=>1,
    			     	    CHR=>$chr,
    				    DSID=>$dsid,	
    				    FEATNAME=>$feat_name,
   				    STRAND=>$strand,
   				    RC=>1,
   				    PRO=>0},
                                  {BUTTON_NUM=>3, 
                                    BUTTON_NAME=>'Protein Sequence',
                                    FEATID=>$feat_id,
   				    PRO=>1,
    			    	    CHR=>$chr,
    				    DSID=>$dsid,	
    				    FEATNAME=>$feat_name,
   				    STRAND=>$strand,
   				    RC=>0,
   				    PRO=>1},
    				  ]);
    $template->param(FEATID=>$feat_id);
    $template->param(PRO=>$pro);
    $template->param(RC=>$rc);
    $template->param(CHR=>$chr);
    $template->param(DSID=>$dsid);
    $template->param(FEATNAME=>$feat_name);
    $template->param(STRAND=>$strand);
    
    #$button .= qq{<DIV id="Proseq"><input type="button" value="$name2" onClick="window.location='$url2'"> };
    #return qq{<FORM ACTION="SeqView.pl?featid=$feat_id&dsid=$dsid&chr=$chr&name=$feat_name&rc=$rc">$button</FORM>};
   # return qq{<FORM METHOD = "GET", ACTION="$url">$button</FORM>};
   return $template->output;

  }
    
sub get_dna_seq_for_feat
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my $dsid = $opts{dsid};
    my $rc = $opts{rc};
    my $upstream = $opts{upstream};
    my $downstream = $opts{downstream};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $chr = $opts{chr};
    my $seq;
    unless ($featid)
    {
    $seq = $DB->get_genomic_sequence(start=>$start,
    					stop=>$stop,
    					chr=>$chr,
    					dsid=>$dsid);
    }
    else
    {
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    unless (ref($feat) =~ /Feature/i)
      {
        return "Unable to retrieve Feature object for id: $featid";
      }
    $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);
    if ($rc)
    {
      $seq = reverse_complement($seq);
    }
    }
    $columns = 60;
    $seq = join ("\n", wrap('','',$seq));
    return $seq;
  }
  
sub get_prot_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    my ($seq) = $DB->get_protein_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    $columns = 60;
    $seq = join ("\n", wrap('','',$seq));
    return $seq;
  }
  
sub reverse_complement
  {
    my $seq = shift;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCG/TAGC/; 
    return $rcseq;
  }
  
sub upstream_color
    {
      my %opts = @_;
      my $seq = $opts{seq};
      my $upstream = $opts{upstream};
      my $up = substr($seq, 0, $upstream);
      return qq{<FONT class="up">$up</FONT>};
    }

sub downstream_color
    {
      my %opts = @_;
      my $seq = $opts{seq};
      my $downstream = $opts{downstream};
      my $down = substr($seq, ((length $seq)-($downstream)), length $seq);
      return qq{<FONT class="down">$down</FONT>};
     }
     
sub main_color
    {
      my %opts = @_;
      my $seq = $opts{seq};
      my $upstream = $opts{upstream};
      my $downstream = $opts{downstream};
      my $main = substr($seq, $upstream, ((length $seq) - ($downstream)));
      return qq{<FONT class="main">$main</FONT>};
    }
    
sub gen_title
    {
      my %opts = @_;
      my $rc = $opts{rc};
      my $pro = $opts{pro};
      my $title;
      unless ($pro)
      {
       $title = $rc ? "Reverse Complement" : "DNA Sequence";
      }
      else
      {
       $title = "Protein Sequence";
      }
      return $title;
    }
