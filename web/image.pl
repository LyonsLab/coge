#!/usr/bin/perl -w

use strict;
use CGI;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
no warnings 'redefine';

use vars qw($P $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS $connstr $DATE $DEBUG $USER $FORM $coge);
$P = CoGe::Accessory::Web::get_defaults($ENV{HOME}.'coge.conf');
$ENV{PATH} = $P->{COGEDIR};

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
($USER) = CoGe::Accessory::LogUser->get_user();
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS );

my $id = $FORM->param('id');
exit unless $id;
my $img = $coge->resultset('Image')->find($id);
print "Content-type: image/png\n\n";
print $img->image;
