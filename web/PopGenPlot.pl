#! /usr/bin/perl -w
use strict;
use CoGe::Algos::PopGen::Results qw( get_plot_data );
use CoGe::Accessory::Web;
use CGI;
use Data::Dumper;
use File::Spec::Functions qw( catfile );
use HTML::Template;
use JSON::XS;

use vars qw($CONF $USER $FORM $DB $PAGE_TITLE $LINK);

$PAGE_TITLE = 'PopGenPlot';

$| = 1;    # turn off buffering

$FORM = new CGI;

( $DB, $USER, $CONF, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my ($x, $y) = get_plot_data catfile($ENV{COGE_HOME}, 'test.txt'), $FORM->Vars->{'type'}, $FORM->Vars->{'column'};

my $template = HTML::Template->new( filename => $CONF->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
$template->param(
    x => encode_json($x),
    y => encode_json($y),
);
print $FORM->header, $template->output;
