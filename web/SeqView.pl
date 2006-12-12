#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use CoGe::Genome;
use Text::Wrap qw($columns &wrap);

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
			);
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();
sub gen_html
  {
    my ($body) = gen_body();
    my $form = $FORM;
    my $rc = $form->param('rc');
    my $protein = $form->param('pro');
    my $title;
    unless ($protein)
    {
      $title = $rc ? "Reverse Complement" : "DNA Sequence";
    }
    else
    {
      $title = "Protein Sequence";
    }
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'CoGe: Sequence Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(BOX_NAME=>"$title");
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
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    my $strand = $feat->strand;
    $strand = check_strand(strand=>$strand, rc=>$rc);
    my $seq = get_seq(featid=>$feat_id, 
		      pro=>$pro,
		      rc=>$rc,
		      chr=>$chr,
		      dsid=>$dsid,
		      featname=>$feat_name,
		      strand=>$strand, 
		      upstream=>$upstream, 
		      downstream=>$downstream,
		      ); 
    return qq{<DIV id="seq">$seq</div>};
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
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    my $ds = $DB->get_dataset_obj->retrieve($dsid);
    my $seq = ">".$ds->org->name."(v.".$feat->version."), Type: ".$feat->type->name.", Location: ".$feat->genbank_location_string.", Chromosome: ".$chr.", Strand: ".$strand.", Name: ".$feat_name."\n";
    unless ($pro)
    {
      $seq .= get_dna_seq_for_feat (featid=>$feat_id, rc=>$rc, upstream=>$upstream, downstream=>$downstream);
    }
    else
    {
      $seq .= get_prot_seq_for_feat($feat_id);
    }
    return qq{<textarea readonly class="seq">$seq</textarea>};
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
    my $strand = $feat->strand;
    $strand = check_strand(strand=>$strand, rc=>$rc);
   
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    
#   my $button = qq{<input type="submit" value="$name">};
    my $url = $form->url(-full=>1, -query=>1);
    
    $url =~ s/rc=.//g;
    $url =~ s/pro=.//g;
    my $url2 = $url;
    my $url3 = $url;
    $url .= "&rc=0";
    $url2 .="&rc=1";
    $url3 .= "&pro=1";
    
    $template->param(BOTTOM_BUTTONS=>1);
    $template->param(BUTTON_LOOP=>[
                                   {BUTTON_NUM=>1, BUTTON_NAME=>'DNA Sequence', BUTTON_URL=>"window.location='$url'"},
                                   {BUTTON_NUM=>2, BUTTON_NAME=>'Reverse Complement', BUTTON_URL=>"window.location='$url2'"},
                                   {BUTTON_NUM=>3, BUTTON_NAME=>'Protein Sequence', BUTTON_URL=>"window.location='$url3'"},
				   {BUTTON_NUM=>4, BUTTON_NAME=>'Extend Sequence', BUTTON_URL=>"get_extend_box([],['extend'])"},
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

sub get_extend_box
  {
    return qq
      {
       <table><tr>
       <TD>UPSTREAM:<td><INPUT type="text" id="upstream" value="0">
       <TD>DOWNSTREAM:<td><INPUT type="text" id="downstream" value="0">
       </table>
      };
  }
    
sub get_dna_seq_for_feat
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my $rc = $opts{rc};
    my $upstream = $opts{upstream} || 0;
    my $downstream = $opts{downstream} || 0;
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    unless (ref($feat) =~ /Feature/i)
      {
        return "Unable to retrieve Feature object for id: $featid";
      }
    my $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);
    if ($rc)
    {
      $seq = reverse_complement($seq);
    }
    $columns = 76;
    $seq = join ("\n", wrap('','',$seq));
    return $seq;
  }
  
sub get_prot_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    my ($seq) = $DB->get_protein_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    $columns = 76;
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
