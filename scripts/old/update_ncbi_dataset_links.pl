#!/usr/bin/perl -w

use strict;
use CoGeX;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
#my $coge = CoGeX->dbconnect();

foreach my $ds ($coge->resultset('Dataset')->search({data_source_id=>21}))
  {
    my $link = $ds->link;
    ($link) = $link =~ /([^\/]*$)/;
    $link =~ s/\..*//;
    $link = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&dopt=gbwithparts&list_uids=".$link;
    $ds->link($link);
    $ds->update;
  }
