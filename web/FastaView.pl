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
no warnings 'redefine';



use vars qw($P $TEMPDIR $TEMPURL $FORM $USER $DATE $coge);
$P = CoGe::Accessory::Web::get_defaults('coge.conf');
$ENV{PATH} = $P->{COGEDIR};

$TEMPDIR = $P->{TEMPDIR};
$TEMPURL = $P->{TEMPURL};
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

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
       my $prot = $form->param('prot');
       my $text = $form->param('text');
       my $name_only = $form->param('no');
       my $upstream = $form->param('up') || 0;
       my $downstream = $form->param('down') || 0;
       my $textbox = $text ? 0 : 1;
       my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'generic_page.tmpl');
#       $template->param(TITLE=>'Fasta Viewer');
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

       my $seqs = get_seqs(prot=>$prot, fids=>\@fids, textbox=>$textbox, gstid=>$gstid, name_only=>$name_only, upstream=>$upstream, downstream=>$downstream);
       if ($text)
	 {
	   return  $seqs;
	 }
       $template->param(BODY=>gen_body(fids=>\@fids, seqs=>$seqs, gstid=>$gstid, prot=>$prot, up=>$upstream, down=>$downstream));
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
    my $upstream = $opts{upstream};
    my $downstream = $opts{downstream};
    #print STDERR Dumper \%opts;
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
	$seqs .= $feat->fasta(col=>80, prot=>$prot, name_only=>$name_only, gstid=>$gstidt, upstream=>$upstream, downstream=>$downstream);
      }
    $seqs = qq{<textarea id=seq_text name=seq_text class="ui-widget-content ui-corner-all backbox" ondblclick="this.select();" style="height: 400px; width: 750px; overflow: auto;">$seqs</textarea>} if $textbox;
    return $seqs;
  }

sub gen_body
  {
    my %opts = @_;
    my $seqs = $opts{seqs};
    my $fids = $opts{fids};
    my $gstid = $opts{gstid} || 1;
    my $prot = $opts{prot} || 0;
    my $up = $opts{up} || 0;
    my $down = $opts{down} || 0;
    $fids = join (",", @$fids) if ref($fids) =~ /array/i;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'FastaView.tmpl');
    $template->param(BOTTOM_BUTTONS=>1);
    $template->param(SEQ=>$seqs) if $seqs;
    $template->param(FIDS=>qq{<input type=hidden id=fids value=$fids><input type=hidden id=gstid value=$gstid>});
    $template->param(PROT=>$prot);
    $template->param(UP=>$up);
    $template->param(DOWN=>$down);
    return $template->output;
  }
	
