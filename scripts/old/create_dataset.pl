#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my ($name, $desc, $link, $data_source_id, $version);

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

GetOptions(
	   "name=s"=>\$name,
	   "desc=s"=>\$desc,
	   "link=s"=>\$link,
	   "source_id"=>\$data_source_id,
	   "version=s"=>\$version
	  );

my $dataset = $coge->resultset('Dataset')->create(
					       {
						name                => $name,
						description         => $desc,
						link                => $link,
						data_source_id      => $data_source_id,,
						version             => $version,
					       }
					      );

if ($dataset->id)
  {
    print "Created dataset.  Database id:  ".$dataset->id,"\n";
  }
else
  {
    print "No database id associated with creating object.  Possible error.  Please check into this.\n";
  }
