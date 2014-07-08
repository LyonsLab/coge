#!/usr/bin/perl -w

use strict;
use CoGeX;

my $dsgid = shift;

unless ( $dsgid )
  {
    print qq{Usage:  $0 <dsgid>

Will generate a dump of all features for a dataset_group in bed format
};
    exit;
  }

my $coge = CoGeX->dbconnect();
my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
foreach my $ds ($dsg->datasets)
  {
    my $rs = $coge->resultset('Feature');
    $rs->result_class('DBIx::Class::ResultClass::HashRefInflator');
    foreach my $ref ($rs->search({dataset_id => $ds->id}))
      {
	print join ("\t", $ref->{chromosome}, $ref->{start}, $ref->{stop}, $ref->{feature_id}),"\n";
      }
  }
