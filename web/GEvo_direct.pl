#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CGI;
use Image::Size;
use HTML::Template;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
$ENV{PATH} = "/opt/apache/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
my $form = new CGI;
my $name = $form->param('name');
my $tmpdir = "/opt/apache/CoGe/tmp/GEvo/";
my $gobe_version = `svnversion /opt/apache/CoGe/gobe/flash`;
$gobe_version =~ s/\n//g;;
my ($USER) = CoGe::Accessory::LogUser->get_user();
my $DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

my %files;
my $tiny;

open (CMD, "/bin/ls $tmpdir/$name"."* |");
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
	    while (<IN>)
	      {
		if (/tiny url: (.*)/i)
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

my $html;
$html .= qq{<DIV id=flash_viewer></DIV>};
$html .= qq{<table>};
$html .= qq{<tr valign=top><td class = small>Alignment reports};
my $i = 1;
foreach my $report (@{$files{report}})
  {
    
    $html .= "<div><font class=small><A HREF=\"$report\" target=_new>View alignment output $i</A></font></DIV>\n";
    $i++;
  }

$html .= qq{<td class = small>Fasta files};
$i=1;
foreach my $item (@{$files{faa}})
  {
    $html .= "<div><font class=small><A HREF=\"$item\" target=_new>Fasta file $i</A></font></DIV>\n";
    $i++;
  }

$html .= qq{<td class = small><a href = "http://baboon.math.berkeley.edu/mavid/gaf.html">GAF</a> annotation files};
$i=1;
foreach my $item (@{$files{anno}})
  {
    $html .= "<div><font class=small><A HREF=\"$item\" target=_new>Annotation file $i</A></font></DIV>\n";
    $i++;
  }
$html .= qq{<td class = small>SQLite db};
my $dbname = $files{sqlite}[0];
$html .= "<div class=small><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
$html .= qq{<td class = small>Log File};
my $logfile = $files{log}[0];
$html .= "<div class=small><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
$html .= qq{<td class = small>GEvo Link<div class=small><a href=$tiny target=_new>$tiny<br>(See log file for full link)</a></div>};
$html .= qq{</table>};

$html = qq{
<script src="/CoGe/gobe/static/swfobject.js" ></script>
$html
<DIV id="flash_viewer"></DIV>
<SCRIPT language="JavaScript">
    var so = new SWFObject("/CoGe/gobe/flash/gobe.swf?$gobe_version", "gobe", $w, $h, "9", "#FFFFFF");
    so.useExpressInstall('/CoGe/gobe/flash/expressinstall.swf');
    so.addVariable('n', $seq_num);
    so.addVariable('base_url', '/CoGe/');
    so.addVariable('img_url', '/CoGe/tmp/GEvo/');
    so.addVariable('base_name', '$name');
    so.addVariable('freezable', getQueryParamValue('freezable') || 'false'); 
    so.addVariable('pad_gs', jQuery('#pad_gs') ? jQuery('#pad_gs').val()  : getQueryParamValue('pad_gs') || '0'); 
    so.addVariable('gsid', getQueryParamValue('gsid') || getQueryParamValue('genespace_id') || '0'); 
    so.write("flash_viewer");
</script>
};
my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');

$template->param(BODY=>$html);
$template->param(USER=>$USER);
$template->param(DATE=>$DATE);
$template->param(TITLE=>"GEvo direct:  reviewing past results.");
print $form->header;
print $template->output;


# <html>
# <head>
# <script src="static/swfobject.js"></script>
# <script>
# /* e.g.  ('up', 1, 330) becomes: byId('drup1').value = 330; */
# /* so up is the 'left', down is the 'right' bar.   */
# function set_genespace(updown, idx, value){
#     var eid = 'dr' + updown + idx;
#     document.getElementById(eid).value = value;
# }

# function call_flash(){
#     document.getElementById('flashcontent').flash_hi();
# }
# </script>
# </head>
# <body>
# <input  name="drup1" id="drup1" size="10" value="10000" type="text">
# <input  name="drdown1" id="drdown1"  size="10" value="10000" type="text">
# <input  name="drup2" id="drup2" size="10" value="10000" type="text">
# <input  name="drdown2" id="drdown2" size="10" value="10000" type="text">
# <input  name="drup3" id="drup3" size="10" value="10000" type="text">
# <input  name="drdown3" id="drdown3" size="10" value="10000" type="text">

# <input  name="pad_gs" id="pad_gs" size="10" value="5000" type="text">

#     <div id="flashcontent" style="z-index=0;position:absolute">
#         <strong>You need to upgrade your Flash Player</strong>
#     </div>

#     <script type="text/javascript">
#         // <![CDATA[
#         var width = getQueryParamValue('w') || "1400"; 
#         var height = getQueryParamValue('h') || "1400"; 
#         var so = new SWFObject("flash/gobe.swf", "gobe", width, height, "9", "#FFFFFF");
#         so.useExpressInstall('flash/expressinstall.swf');
#         so.addVariable('n', getQueryParamValue('n') || 2); 
#         so.addVariable('freezable', getQueryParamValue('freezable') || 'false'); 
#         so.addVariable('pad_gs', getQueryParamValue('pad_gs') || '10000'); 
#         so.addVariable('gsid', getQueryParamValue('gsid') || getQueryParamValue('genespace_id') || '0'); 
#         so.addVariable('base_url', getQueryParamValue('base_url')|| '/CoGe/gobe/');
#         so.addVariable('img_url', getQueryParamValue('img_url')|| '/CoGe/gobe/tmp/');
#         so.addVariable('base_name', getQueryParamValue('base_name')|| 'GEvo_xMpPI9pS');
#         so.write("flashcontent");
        
#         // ]]>
#     </script>
# </body>
# </html>
