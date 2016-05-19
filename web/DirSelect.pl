#! /usr/bin/perl -w

# NOTE: this file shares a lot of code with LoadExperiment.pl & LoadAnnotation.pl, replicate changes when applicable.

use strict;
use CGI;
use CoGeX;
use CoGe::Accessory::Web;
use HTML::Template;
use Data::Dumper;
no warnings 'redefine';

use vars qw(
  $CONF $USER $DB $FORM $LINK
);

$FORM = new CGI;
( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM
);

my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . 'widgets/FileSelect.tmpl' );
$template->param(
    DIR_SELECT => 1,
    DISABLE_IRODS_GET_ALL => 1
);
print $FORM->header, $template->output;
