#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CGI;
use Image::Size;
use HTML::Template;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CGI::Ajax;

$ENV{PATH} = "/opt/apache/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw ($FORM $USER $DATE $TMPDIR $TMPURL);

$FORM = new CGI;
$TMPDIR = "/opt/apache/CoGe/tmp/GEvo/";
$TMPURL = "/CoGe/tmp/GEvo";
($USER) = CoGe::Accessory::LogUser->get_user();
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));


my $pj = new CGI::Ajax (
			gen_html => \&gen_html,
		       );

$pj->js_encode_function('escape');

print $pj->build_html($FORM, \&gen_html);

#print $FORM->header, gen_html();


sub gen_html
  {
#    my $form = $FORM;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(BODY=>gen_body());
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(TITLE=>"GEvo direct:  reviewing past results.");
#    print $form->header;
#    print STDERR $template->output;
    return $template->output;
  }

sub gen_body
    {
      my $form = $FORM;
      my $name = $form->param('name');
      my $gobe_version = `svnversion /opt/apache/CoGe/gobe/flash`;
      $gobe_version =~ s/\n//g;;

      my %files;
      my $tiny;
      
      open (CMD, "/bin/ls $TMPDIR/$name"."* |");
      while (<CMD>)
	{
	  my $touch = "touch $_";
	  my $x;
	  ($x, $touch) = check_taint($touch);
	  `$touch`;
	  foreach (split /\n/)
	    {
	      if (/\.anno/)
		{
		  push @{$files{anno}},$_;
		}
	      elsif (/\.faa/)
		{
		  push @{$files{faa}},$_;
		}
	      elsif (/\.png/)
		{
		  push @{$files{png}},$_;
		}
	      elsif (/\.log/)
		{
		  open (IN, $_);
		  while (my $line = <IN>)
		    {
		      if ($line =~ /tiny url: (.*)/i)
			{
			  $tiny= $1;
			  last;
			}
		    }
		  close IN;
		  
		  push @{$files{log}},$_;
		}
	      elsif (/\.sqlite/)
		{
		  push @{$files{sqlite}},$_;
		}
	      else
		{
		  push @{$files{report}},$_;
		}
	    }
	}
      close CMD;
#print STDERR Dumper \%files;
      my $h=0;
      my $w=0;
      my $seq_num = 0;
      foreach my $img (@{$files{png}})
	{
	  my ($x, $y) = imgsize($img);
	  $h += $y+0.1*$y;
	  $w = $x;
	  $seq_num++;
	}
      $w+=400;
      my $html;
      $html .= qq{<DIV id=flash_viewer></DIV>};
      $html .= qq{<table>};
      $html .= qq{<tr valign=top><td class = small>Alignment reports};
      my $i = 1;
      foreach my $report (@{$files{report}})
	{
	  $report =~ s/$TMPDIR/$TMPURL/;
	  $html .= "<div><font class=small><A HREF=\"$report\" target=_new>View alignment output $i</A></font></DIV>\n";
	  $i++;
	}

      $html .= qq{<td class = small>Fasta files};
      $i=1;
      foreach my $item (@{$files{faa}})
	{
	  $item =~ s/$TMPDIR/$TMPURL/;
	  $html .= "<div><font class=small><A HREF=\"$item\" target=_new>Fasta file $i</A></font></DIV>\n";
	  $i++;
	}
      
      $html .= qq{<td class = small><a href = "http://baboon.math.berkeley.edu/mavid/gaf.html">GAF</a> annotation files};
      $i=1;
      foreach my $item (@{$files{anno}})
	{
	  $item =~ s/$TMPDIR/$TMPURL/;
	  $html .= "<div><font class=small><A HREF=\"$item\" target=_new>Annotation file $i</A></font></DIV>\n";
	  $i++;
	}
      $html .= qq{<td class = small>SQLite db};
      my $dbname = $files{sqlite}[0];
      $dbname =~ s/$TMPDIR/$TMPURL/;
      $html .= "<div class=small><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
      $html .= qq{<td class = small>Log File};
      my $logfile = $files{log}[0];
      $logfile =~ s/$TMPDIR/$TMPURL/;
      $html .= "<div class=small><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
      $html .= qq{<td class = small>GEvo Link<div class=small><a href=$tiny target=_new>$tiny<br>(See log file for full link)</a></div>};
      $html .= qq{</table>};
      my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo_direct.tmpl');
      $template->param('STUFF'=>$html);
      $template->param('WIDTH'=>$w);
      $template->param('HEIGHT'=>$h);
      $template->param('HEIGHT'=>$h);
      $template->param('GOBE_VERSION'=>$gobe_version);
      $template->param('SEQ_NUM'=>$seq_num);
      $template->param('BASE_NAME'=>$name);
      return $template->output;
    }

#$html = qq{
#<script src="/CoGe/gobe/static/swfobject.js" ></script>
#$html
#<DIV id="flash_viewer"></DIV>
#<SCRIPT language="JavaScript">
#    var so = new SWFObject("/CoGe/gobe/flash/gobe.swf?$gobe_version", "gobe", $w, $h, "9", "#FFFFFF");
#    so.useExpressInstall('/CoGe/gobe/flash/expressinstall.swf');
#    so.addVariable('n', $seq_num);
#    so.addVariable('base_url', '/CoGe/');
#    so.addVariable('img_url', '/CoGe/tmp/GEvo/');
#    so.addVariable('base_name', '$name');
#    so.addVariable('freezable', getQueryParamValue('freezable') || 'false'); 
#    so.addVariable('pad_gs', jQuery('#pad_gs') ? jQuery('#pad_gs').val()  : getQueryParamValue('pad_gs') || '0'); 
#    so.addVariable('gsid', getQueryParamValue('gsid') || getQueryParamValue('genespace_id') || '0'); 
#    so.write("flash_viewer");
#</script>
#};


