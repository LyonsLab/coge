#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::LogUser;
use CoGeX;
$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $coge);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
($USER) = CoGe::Accessory::LogUser->get_user();
$coge=CoGeX->dbconnect();

my $id = $FORM->param('id');
exit unless $id;
my $img = $coge->resultset('Image')->find($id);
print "Content-type: image/png\n\n";
print $img->image;
