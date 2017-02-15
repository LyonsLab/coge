#! /usr/bin/perl -w

use strict;

use CGI;
use JSON::XS;

my $form = new CGI;
my $filename = '' . $form->param('image');
my $tmpfilename = $form->tmpFileName( $filename );
chmod 0666, $tmpfilename;
my $fh = $form->upload('image');
my $size = 0;
$size = -s $fh if $fh;

my $resp = encode_json({
    filename  => $filename,
    tmpfilename => $tmpfilename,
    size => $size
});
print $form->header, $resp;
