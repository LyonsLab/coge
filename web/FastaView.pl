#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use POSIX;

$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $coge);

$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;

$coge = CoGeX->dbconnect();

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seqs=>\&get_seqs,
		       export_to_file=>\&export_to_file,
			);
$pj->js_encode_function('escape');
if ($FORM->param('text'))
    {
      print $FORM->header('text');
      print gen_html();
    }
else
    {
      print $pj->build_html($FORM, \&gen_html);
    }


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
    my $rc = $form->param('rc');
    my $prot = $form->param('prot');
    my $text = $form->param('text');
    my $textbox = $text ? 0 : 1;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Fasta Viewer');
    $template->param(PAGE_TITLE=>'FastaView');
    $template->param(HELP=>'/wiki/index.php?title=FastaView');
    my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"FastaView-logo.png");
    $template->param(BOX_NAME=>qq{<DIV id="box_name">Sequences:</DIV>});
    my @fids;
    push @fids, $form->param('featid') if $form->param('featid');
    push @fids, $form->param('fid') if $form->param('fid');
    my $gstid = $form->param('gstid') if $form->param('gstid');
    my $seqs = get_seqs(prot=>$prot, fids=>\@fids, textbox=>$textbox, gstid=>$gstid);
    if ($text)
      {
	return  $seqs;
      }
    $template->param(BODY=>gen_body(fids=>\@fids, seqs=>$seqs, gstid=>$gstid));
    $html .= $template->output;
    }
    return $html;
  }

sub get_seqs
  {
    my %opts = @_;
    my $fids = $opts{fids};
    my $prot = $opts{prot};
    my $textbox = $opts{textbox};
    my $name_only = $opts{name_only};
    my $gstid = $opts{gstid};
    my @fids = ref($fids) =~ /array/i ? @$fids : split/,/, $fids;
    my $seqs;
    foreach my $featid (@fids)
      {
	my ($fid, $gstidt);
	if ($featid =~ /_/)
	  {
	    ($fid, $gstidt) = split /_/, $featid;
	  }
	else
	  {
	    ($fid, $gstidt) = ($featid, $gstid);
	  }
	my ($feat) = $coge->resultset('Feature')->find($fid);
	next unless $feat;
	$seqs .= $feat->fasta(col=>100, prot=>$prot, name_only=>$name_only, gstid=>$gstidt);
      }
    $seqs = qq{<textarea id=seq_text name=seq_text class=backbox readonly ondblclick="this.select();" style="height: 400px; width: 750px; overflow: auto;">$seqs</textarea>} if $textbox;
    return $seqs;
  }

sub gen_body
  {
    my %opts = @_;
    my $seqs = $opts{seqs};
    my $fids = $opts{fids};
    my $gstid = $opts{gstid} || 1;
    $fids = join (",", @$fids) if ref($fids) =~ /array/i;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/FastaView.tmpl');
    $template->param(BOTTOM_BUTTONS=>1);
    $template->param(SEQ=>$seqs) if $seqs;
    $template->param(FIDS=>qq{<input type=hidden id=fids value=$fids><input type=hidden id=gstid value=$gstid>});
    return $template->output;
  }
	
