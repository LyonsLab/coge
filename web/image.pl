#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;

my $FORM = new CGI;
my $id   = $FORM->param('id');
exit unless $id;

my ( $coge, $USER, $P ) = CoGe::Accessory::Web->init;
my $img = $coge->resultset('Image')->find($id);

$| = 1;    # turn off buffering
print "Content-type: image/png\n\n", $img->image;
