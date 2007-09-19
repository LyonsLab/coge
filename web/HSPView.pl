#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use GS::LogUser;
use HTML::Template;
use Data::Dumper;
use DBI;

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
    $template->param(TITLE=>'HSP Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    my $report_file = $form->param('blast_report') || $form->param('report');
    my $dbfile = $form->param('db');
    if ($report_file && -r $report_file)
      {
	my $hsp_num = $form->param('hsp_num') || $form->param('num');
	$report_file =~ /([^\/]*$)/;
        $template->param(BOX_NAME=>qq{<a href=/CoGe/tmp/$1>$1</a>}. " HSP: $hsp_num") if $1;
        $template->param(BODY=>gen_body(report_file=>$report_file, hsp_num=>$hsp_num, db_file=>$dbfile));
      }
    else
      {
	$template->param(BODY=>"Can't read or find file $report_file");
      }
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my %opts = @_;
    my $report_file = $opts{report_file};
    my $db_file = $opts{db_file};
    my $hsp_num = $opts{hsp_num};
#    my $report;
#    $report= new CoGe::Accessory::bl2seq_report({file=>$report_file}) if $report_file =~ /bl2seq/i;
#    $report= new CoGe::Accessory::blastz_report({file=>$report_file}) if $report_file =~ /blastz/i;
#    $report= new CoGe::Accessory::lagan_report({file=>$report_file}) if $report_file =~ /lagan/i;
#    $report= new CoGe::Accessory::chaos_report({file=>$report_file}) if $report_file =~ /chaos/i;
#    $report= new CoGe::Accessory::dialgn_report({file=>$report_file}) if $report_file =~ /dialign/i;
#    print STDERR Dumper $report;

    my $hsps = get_info_from_db(db_file=>$db_file, hsp_num=>$hsp_num, report_file=>$report_file);
    my ($qname, $sname) = $hsps->[0]{hsp}=~ />?\(?(.*?)-(.*)\)?<?/;

    my $hsp;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/HSPView.tmpl');
    $template->param(Query=>$qname);
    $template->param(Subject=>$sname);
    my $qgap = $hsps->[0]{alignment} =~ tr/-/-/;
    my $sgap = $hsps->[1]{alignment} =~ tr/-/-/;

    $template->param(QGap=>$qgap);
    $template->param(SGap=>$sgap);
    my $qseq = $hsps->[0]{alignment};
    $qseq =~ s/-//g;
    my $sseq = $hsps->[1]{alignment};
    $sseq =~ s/-//g;
    $template->param(QLoc=>$hsps->[0]{start}."-".$hsps->[0]{stop});
    $template->param(SLoc=>$hsps->[1]{start}."-".$hsps->[1]{stop});
    $template->param(QIdent=>sprintf("%.1f",$hsps->[0]{match}/length($qseq)*100)."%") if $qseq;
    $template->param(SIdent=>sprintf("%.1f",$hsps->[1]{match}/length($sseq)*100)."%") if $sseq;
    $template->param(QSim=>sprintf("%.1f",$hsps->[0]{match}/length($qseq)*100)."%") if $qseq;
    $template->param(SSim=>sprintf("%.1f",$hsps->[1]{match}/length($sseq)*100)."%") if $sseq;
#    $template->param(QSim=>sprintf("%.1f",$hsp->positive/length($qseq)*100)."%");
#    $template->param(SSim=>sprintf("%.1f",$hsp->positive/length($sseq)*100)."%");
    $template->param(QLen=>length($qseq));
    $template->param(SLen=>length($sseq));
    $template->param(QMis=>length($qseq)-$hsps->[0]{match}-$qgap);
    $template->param(SMis=>length($sseq)-$hsps->[1]{match}-$sgap);
    $template->param(strand=>$hsps->[0]{orientation});
    $template->param(eval=>$hsps->[0]{eval});
    $template->param(Match=>$hsps->[0]{match});
    $template->param(Sim=>$hsps->[0]{match});

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
    my @qln = split /\n/,seqwrap($hsps->[0]{alignment},100);
    my @sln = split /\n/,seqwrap($hsps->[1]{alignment},100);
    my @aln = split /\n/,seqwrap(get_alignment($hsps->[0]{alignment},$hsps->[1]{alignment}),100);
#    my @aln = split /\n/,seqwrap($hsp->align,100);
    my $align;
    for (my $i=0; $i<@aln; $i++)
      {
	$align .= join ("\n", $qln[$i],$aln[$i],$sln[$i])."\n\n";
      }
    $template->param(align=>"<pre>".$align."</pre>");
#    $template->param(newstuff=>get_info_from_db(db_file=>$db_file, hsp_num=>$hsp_num, report_file=>$report_file));
    my $html;
    $html .= $template->output;
    return $html;
  }

sub get_alignment
  {
    my @seq1 = split //, shift;
    my @seq2 = split //, shift;
    my $aln;
    for (my $i = 0; $i < @seq1; $i++)
      {
	$aln .= $seq1[$i] eq $seq2[$i] ? "|" : " ";
      }
    return $aln;
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


sub get_info_from_db
  {
    my %opts = @_;

    my $db = $opts{db_file};
    my $hsp_num = $opts{hsp_num};
    my $report_file = $opts{report_file};
    my ($base) = $db =~ /^(.*?)\./;
    my ($set1, $set2) = $report_file=~ /_(\d+)-(\d+)/;
    my $output;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$TEMPDIR/$db","","");
    my $statement = qq{SELECT annotation, pair_id, id from image_data where name like "$hsp_num-%" and (image = "$base}."_".qq{$set1.png" or image = "$base}."_".qq{$set2.png")};
#    print STDERR $statement,"\n";
    my $sth = $dbh->prepare($statement);
    $sth->execute();
    my %hsps;
    my %ids;
    while(my $item = $sth->fetchrow_arrayref)
      {
	my $split;
	$split = "<br>" if $item->[0] =~/<br>/;
	$split = "\n" if $item->[0] =~/\n/;
	$split = "\\n" if $item->[0] =~/\\n/;
	$split = "<br\/>" if $item->[0] =~/<br\/>/;
	$split = "&#10;" if $item->[0] =~/&#10;/;
	 my %data;
	foreach my $item (split /$split/i, $item->[0])
	  {
	    my ($tmp1, $tmp2) = split/:/,$item,2;
	    $tmp2 =~ s/^\s+//;
	    $tmp2 =~ s/\s+$//;
	    $data{$tmp1}=$tmp2;
	  }
	$data{HSP} =~ s/<.*?>//g;
	$data{HSP} =~ s/^\d+\s*//;
	my ($start, $stop, $orientation) = $data{Location} =~ /(\d+)-(\d+)\s+\(?(.*)\)?/;

	$hsps{$item->[2]}= {
		     hsp=>$data{HSP},
		     start=>$start,
		     stop=>$stop,
		     orientation=>$orientation,
		     alignment=>$data{Sequence},
		     match=>$data{Match},
		     length=>$data{Length},
		     identity=>$data{Identity},
		     eval=>$data{E_val},
		     location=>$data{Location},
		    };
	$ids{$item->[1]}++;
	$ids{$item->[2]}++;
      }
    my @hsps;
    #uberlame workaround ,but hey, it works
    foreach my $id (keys %ids)
      {
	next unless $ids{$id}==2;
	push @hsps, $hsps{$id};
      }
    return \@hsps;
  }
