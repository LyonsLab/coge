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
use POSIX;

$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $DB $coge);

$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;
$FORM = new CGI;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seq=>\&get_seq,
		       gen_title=>\&gen_title,
		       find_feats=>\&find_feats,
		       parse_url=>\&parse_url,
			);
$pj->js_encode_function('escape');
#print $pj->build_html($FORM, \&gen_html);
print $FORM->header;
print gen_html();

sub gen_html
  {
    my $html;
    unless ($USER)
      {
	$html = login();
      }
    else
     {
    my $form = $FORM;
    my $feat_id = $form->param('featid');
    my $rc = $form->param('rc');
    my $pro;
    my ($title) = gen_title(protein=>$pro, rc=>$rc);
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Fasta Viewer');
    $template->param(HELP=>'FastaView');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"FastaView-logo.png");
    $template->param(BOX_NAME=>qq{<DIV id="box_name">$title</DIV>});
    $template->param(BODY=>gen_body());
#    $template->param(POSTBOX=>gen_foot());
    $html .= $template->output;
    }
    return $html;
  }

sub gen_body
  {
    my $form = $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FastaView.tmpl');
    my $seqs;
    foreach my $featid ($form->param('featid'))
      {
	my ($feat) = $coge->resultset('Feature')->find($featid);
	next unless $feat;
	$seqs .= $feat->fasta(col=>100);
      }
    $template->param(SEQ=>$seqs);
    return $template->output;
  }
  
sub reverse_complementq
  {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATCG/TAGC/;
    return $seq;
  }

sub get_prot_seq_for_feat
  {
    my $featid = shift;
    #print STDERR "featid: ", $featid, "\n";
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    my ($seq) = $DB->get_protein_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    #print $seq;
    $columns = 60;
    $seq = join ("\n", wrap('','',$seq));
    return $seq;
  }
 
sub gen_title
    {
      my %opts = @_;
      my $rc = $opts{'rc'} || 0;
      my $pro = $opts{'pro'};
      my $title;
      unless ($pro)
      {
       if ($rc == 2)
        {$title = "Six Frame Translation";}
       elsif ($rc == 1)
        {$title = "Reverse Complement";}
       else
        {$title = "DNA Sequence";}
      }
      else
      {
       $title = "Protein Sequence";
      }
      return $title;
    }

sub sixframe
	{
	  my %opts = @_;
	  my $seq = $opts{seq};
	  my $fasta = $opts{fasta};
	  my $key;
	  my $sixframe;
     	  my $sequence = $DB->get_feat_obj->frame6_trans(seq=>$seq);
          #print STDERR Dumper ($sequence);
          foreach $key (sort {abs($a) <=> abs($b) || $b <=> $a} keys %$sequence)
           {
      	     $seq = join ("\n", wrap('','',$sequence->{$key}));
      	     $sixframe .= qq/$fasta Frame $key\n$seq\n/;
           }
          return $sixframe;
        }
	
