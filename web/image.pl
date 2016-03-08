#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::Web;

my $form = new CGI;
my $id = $form->param('id');
my ( $coge, $user ) = CoGe::Accessory::Web->init;
print STDERR "image.pl: ", $user->name, "\n";
if (!$id) {
	$id	 = $user->image_id;
}
exit unless $id;

my $img = $coge->resultset('Image')->find($id);
exit unless $img;

$| = 1;    # turn off buffering
print "Cache-Control: public\nContent-type: image/png\n\n", $img->image; # mdb added caching 11/24/15
