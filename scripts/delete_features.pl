#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $GO $conf_file $coge $gid $restricted $P $ftid $user $pwd);

GetOptions ( "debug=s" => \$DEBUG,
             "go=s"    => \$GO,
             "conf_file|cf=s" => \$conf_file,
             "gid=i"=>\$gid, #genome id
	     "ftid=i"=>\$ftid, #feature type id
	     "u=s"=>\$user,
	     "p=s"=>\$pwd,
           );

$P = CoGe::Accessory::Web::get_defaults($conf_file);

unless ($P && $P->{DBNAME}) {usage();}

#my $TEMPDIR = $P->{TEMPDIR}."copy_genome";
#mkpath( $TEMPDIR,    1, 0777 );

my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};

$DBUSER = $user if $user;
$DBPASS = $pwd if $pwd;

my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->dbconnect(db_connection_string=>$connstr, db_name=>$DBUSER, db_passwd=>$DBPASS );

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

unless ($coge && $gid)
  {
    usage();
  }

my $genome = $coge->resultset('Genome')->find($gid);
foreach my $feat ($genome->features({feature_type_id=>$ftid}))
  {
    $feat->delete if $GO;
  }

sub usage
  {
    print qq{
Usage: $0 -cf <coge config file> -gid <genome id> -ftid <feature type id> -go 1
};
    exit;
  }
