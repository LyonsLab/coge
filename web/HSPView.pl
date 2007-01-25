#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use GS::LogUser;
use HTML::Template;
use Data::Dumper;
use CoGe::Genome;
use CoGe::Accessory::bl2seq_report;

# for security purposes
$ENV{PATH} = "/opt/apache2/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM);
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
# set this to 1 to print verbose messages to logs
$DEBUG = 0;

$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$CGI::POST_MAX= 60 * 1024 * 1024; # 24MB
$CGI::DISABLE_UPLOADS = 0; 
($USER) = GS::LogUser->get_user();
my $pj = new CGI::Ajax(
		       go=>\&gen_html,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');

#print $pj->build_html($FORM, \&gen_html);
print "Content-Type: text/html\n\n";
print gen_html();
sub gen_html
  {
    my $form = shift || $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(LOGO_PNG=>"HSPView-logo.png");
    $template->param(TITLE=>'Blast HSP Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    my $bfile = $form->param('blast_report');
    if ($bfile && -r $bfile)
      {
	my $hsp_num = $form->param('hsp_num');
	$bfile =~ /([^\/]*$)/;
        $template->param(BOX_NAME=>qq{<a href=/CoGe/tmp/$1>$1</a>}. " HSP: $hsp_num") if $1;
        $template->param(BODY=>gen_body($bfile, $hsp_num));
      }
    else
      {
	$template->param(BODY=>"Can't read or find file $bfile");
      }
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $blast_file = shift;
    my $hsp_num = shift;
    my $blast = new CoGe::Accessory::bl2seq_report($blast_file);
    my $hsp;
    foreach my $item (@{$blast->hsps})
      {
	$hsp = $item if $item->number eq $hsp_num;
      }
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/HSPView.tmpl');
       $template->param(Query=>$blast->query);
    $template->param(Subject=>$blast->subject);
    $template->param(QGap=>$hsp->qgap);
    $template->param(SGap=>$hsp->sgap);
    my $qseq = $hsp->qseq;
    $qseq =~ s/-//g;
    my $sseq = $hsp->sseq;
    $sseq =~ s/-//g;
    $template->param(QIdent=>sprintf("%.1f",$hsp->match/length($qseq)*100)."%");
    $template->param(SIdent=>sprintf("%.1f",$hsp->match/length($sseq)*100)."%");
    $template->param(QSim=>sprintf("%.1f",$hsp->positive/length($qseq)*100)."%");
    $template->param(SSim=>sprintf("%.1f",$hsp->positive/length($sseq)*100)."%");
    $template->param(QLen=>length($qseq));
    $template->param(SLen=>length($sseq));
    $template->param(QMis=>length($qseq)-$hsp->positive);
    $template->param(SMis=>length($sseq)-$hsp->positive);
    $template->param(strand=>$hsp->strand);
    $template->param(eval=>$hsp->eval);
    $template->param(Match=>$hsp->match);
    $template->param(Sim=>$hsp->positive);
    my $qstart = $FORM->param('qstart');
    my $qstop = $FORM->param('qstop');
    my $qchr = $FORM->param('qchr');
    my $qds = $FORM->param('qds');
    my $qstrand = $FORM->param('qstrand') || 1;
    my $qrc = $qstrand =~ /-/ ? "1" : "0";
    my $qseqview_link = "SeqView.pl?start=$qstart&stop=$qstop&chr=$qchr&dsid=$qds&strand=$qstrand&rc=$qrc" if $qstart && $qstop && $qchr && $qds;
    my $qseq_out = "<pre>".seqwrap($qseq)."</pre>";
    $qseq_out .= "<font class=small><a href = $qseqview_link target=_new>Open Sequence in SeqView</a></font>" if $qseqview_link;
      
    $template->param(qseq=>$qseq_out);
    $template->param(sseq=>"<pre>".seqwrap($sseq)."</pre>");
    my @qln = split /\n/,seqwrap($hsp->qseq,100);
    my @sln = split /\n/,seqwrap($hsp->sseq,100);
    my @aln = split /\n/,seqwrap($hsp->align,100);
    my $align;
    for (my $i=0; $i<@aln; $i++)
      {
	$align .= join ("\n", $qln[$i],$aln[$i],$sln[$i])."\n\n";
      }
    $template->param(align=>"<pre>".$align."</pre>");
    my $html;
    $html .= $template->output;
    return $html;
  }

sub seqwrap {
	my $seq = shift;
	my $wrap = shift || 50;
	my @chars = split( //, $seq );
	my $newseq = "";
	my $counter = 0;
	foreach my $char ( @chars ) {
	if ( $counter < $wrap ) {
		$newseq .= $char;
	} else {
		$newseq .= "\n" . $char;
		$counter = 0;
	}
		$counter++;
	}
	return($newseq);
}
