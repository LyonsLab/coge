#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;

my $FORM = new CGI;
my $id   = $FORM->param('id');
my ( $coge, $USER, $P ) = CoGe::Accessory::Web->init;
if (!$id) {
	$id	 = $USER->image_id;
}
exit unless $id;

my $img = $coge->resultset('Image')->find($id);

$| = 1;    # turn off buffering
print "Cache-Control: public\nContent-type: image/png\n\n", $img->image; # mdb added caching 11/24/15
