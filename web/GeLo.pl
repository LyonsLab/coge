#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use Bio::TreeIO;
use File::Temp;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $NEATO $DOT);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$NEATO = "/usr/bin/neato";
$DOT = "/usr/bin/dot";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$USER = CoGe::Accessory::LogUser->get_user();
print "Content-Type: text/html\n\n";
print gen_html($FORM);

sub gen_html
  {
    my $form = shift;
    my ($body, $seq_names, $seqs) = gen_body($form);
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');

    $template->param(TITLE=>'CoGe: Genome Location Viewer');
    $template->param(HEAD=>qq{
<link rel="stylesheet" type="text/css" href="css/tiler.css">
<script src=js/Dom.js> </script>
<script src=js/common.js> </script>
<script src=js/extendEvent.js> </script>
<script src=js/panner.js> </script>
<script src=js/tiler.js> </script>
<script src=js/tilerConfig.js> </script>
<script src=js/kaj.stable.js> </script>
});
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"GeLo-logo.png");
    $template->param(BODY=>$body);
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift;
    my $chr = $form->param('chr') || 1;
    my $di = $form->param('di') || 6;
    my $z = $form->param('z') || 7;
    my $x = $form->param('x') || 1;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GeLo.tmpl');
    $template->param(FRONT_PAGE => 1) unless ($chr && $di);
    $template->param(CHR=>$chr);
    $template->param(DI=>$di);
    $template->param(Z=>$z);
    $template->param(X=>$x);
    return $template->output;
  }
