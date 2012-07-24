#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use POSIX;
use File::Path;

no warnings 'redefine';



use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $TEMPDIR $TEMPURL $FORM $USER $DATE $coge $COOKIE_NAME);
$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');
$ENV{PATH} = $P->{COGEDIR};

$TEMPDIR = $P->{TEMPDIR}."FastaView/";
$TEMPURL = $P->{TEMPURL}."FastaView/";
mkpath ($TEMPDIR, 0,0777) unless -d $TEMPDIR;

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) =$FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(cookie_name=>$COOKIE_NAME, ticket=>$cas_ticket, coge=>$coge, this_url=>$FORM->url()) if($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(cookie_name=>$COOKIE_NAME,coge=>$coge) unless $USER;

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seqs=>\&get_seqs,
		       gen_file=>\&gen_file,
			);
$pj->js_encode_function('escape');
if ($FORM->param('text'))
    {
      my $header = "Content-disposition: attachement; filename=CoGe_";#test.gff\n\n";
      $header .= int(rand(100000));
      $header .=".faa\n\n";

      print $header;
      print gen_html();
    }
else
    {
      print $pj->build_html($FORM, \&gen_html);
    }
sub gen_html
  {
    my $html;
    my $form = $FORM;
    my $prot = $form->param('prot');
    my $text = $form->param('text');
    my $name_only = $form->param('no');
    my $id_only = $form->param('io');
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
    
    my ($seqs, $seq_count, $feat_count, $warning) = get_seqs(prot=>$prot, fids=>\@fids, textbox=>$textbox, gstid=>$gstid, name_only=>$name_only, id_only=>$id_only, upstream=>$upstream, downstream=>$downstream);
    if ($text)
      {
	return  $seqs;
      }
    $template->param(BODY=>gen_body(fids=>\@fids, seqs=>$seqs, seq_count=>$seq_count, feat_count=>$feat_count, gstid=>$gstid, prot=>$prot, up=>$upstream, down=>$downstream, message=>$warning));
    $html .= $template->output;
  return $html;
  }

sub get_seqs
  {
    my %opts = @_;
    my $fids = $opts{fids};
    my $prot = $opts{prot};
    my $textbox = $opts{textbox};
    my $name_only = $opts{name_only};
    my $id_only = $opts{id_only};
    my $gstid = $opts{gstid};
    my $upstream = $opts{upstream};
    my $downstream = $opts{downstream};
    my @fids = ref($fids) =~ /array/i ? @$fids : split/,/, $fids;
    my %seen = ();
    @fids = grep {!$seen{$_}++} @fids;
 
    my $seqs;
    my $seq_count = 0;
    my $fid_count = 0;
    foreach my $item (@fids)
      {
        foreach my $featid (split /,/, $item)
          {
	    $fid_count++;
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
	    unless ($feat)
	      {
		$seqs .= ">Not found: $featid\n";
		next;
	      }
	    my ($dsg) = $feat->dataset->dataset_groups;
	    if ($dsg->restricted && !$USER->has_access_to_genome($dsg))
	      {
		$seqs .= ">Restricted: $featid\n";
		next;
	      }
	    my $tmp = $feat->fasta(col=>100, prot=>$prot, name_only=>$name_only, fid_only=>$id_only, gstid=>$gstidt, upstream=>$upstream, downstream=>$downstream);
	    $seq_count += $tmp =~ tr/>/>/;
	    $seqs .= $tmp;
	  }
      }
    my $warning;
    %seen = ();
    while ($seqs =~ /(>.*\n)/g)
      {
	$warning = "Warning: Duplicate sequence names" if $seen{$1};
	$seen{$1}=1;
      }
    $seqs = qq{<textarea id=seq_text name=seq_text class="ui-widget-content ui-corner-all backbox" ondblclick="this.select();" style="height: 400px; width: 800px; overflow: auto;">$seqs</textarea>} if $textbox;
    return $seqs, $seq_count, $fid_count, $warning;
  }

sub gen_file
  {
    my %opts = @_;
    my $fids = $opts{fids};
    my $prot = $opts{prot};
    my $textbox = $opts{textbox};
    my $name_only = $opts{name_only};
    my $id_only = $opts{id_only};
    my $gstid = $opts{gstid};
    my $upstream = $opts{upstream};
    my $downstream = $opts{downstream};
    my @fids = ref($fids) =~ /array/i ? @$fids : split/,/, $fids;
    my ($seqs) = get_seqs(prot=>$prot, fids=>\@fids, gstid=>$gstid, name_only=>$name_only, id_only=>$id_only, upstream=>$upstream, downstream=>$downstream);
    my $file = $TEMPDIR."Seqs_".int(rand(100000000)).".faa";
    open (OUT, ">".$file);
    print OUT $seqs;
    close OUT;
    my $url = $file;
    $url =~ s/$TEMPDIR/$TEMPURL/;
    $url =~ s/^\/[^\/]*//;
    $url = $P->{SERVER}.$url;
    return $url;
  }

sub gen_body
  {
    my %opts = @_;
    my $seqs = $opts{seqs};
    my $seq_count = $opts{seq_count};
    my $feat_count = $opts{feat_count};
    my $message = $opts{message};
    my $fids = $opts{fids};
    my $gstid = $opts{gstid} || 1;
    my $prot = $opts{prot} || 0;
    my $up = $opts{up} || 0;
    my $down = $opts{down} || 0;
    $fids = join (",", @$fids) if ref($fids) =~ /array/i;
    my $template = HTML::Template->new(filename=>$P->{TMPLDIR}.'FastaView.tmpl');
    $template->param(BOTTOM_BUTTONS=>1);
    $template->param(SEQ=>$seqs) if $seqs;
    $template->param(SEQ_COUNT=>$seq_count) if defined $seq_count;
    $template->param(FEAT_COUNT=>$feat_count) if defined $feat_count;
    $template->param(WARNING=>$message) if defined $message;
    $template->param(FIDS=>qq{<input type=hidden id=fids value=$fids><input type=hidden id=gstid value=$gstid>});
    $template->param(PROT=>$prot);
    $template->param(UP=>$up);
    $template->param(DOWN=>$down);
    return $template->output;
  }
	
