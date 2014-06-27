#!/usr/bin/perl -w

use strict;
use Data::Dumper;
#use CoGeX;
use DBI;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $query = qq{select feature_name_id,name from feature_name where name like "% "};
my $sth = $coge->prepare($query);
$sth->execute();
my $updateq = qq{update feature_name set name = ? where feature_name_id = ?};
my $update = $coge->prepare($updateq);
while (my $item = $sth->fetchrow_arrayref)
  {
    print "Updating: '",$item->[1],"'\n";
    $item->[1] =~ s/\s+$//g;
    $update->execute($item->[1],$item->[0]);
  }
#    $name =~ s/\s+$//g;
