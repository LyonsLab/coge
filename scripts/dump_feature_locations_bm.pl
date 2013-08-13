#!/usr/bin/perl -w
use strict;
use MyDB;

# set up our database
my $dbo = new MyDB(
    dbhost => '',
    dbuser => '',
    dbpw   => '',
    dbname => ''
);

open( OUT, ">x.txt" ) or die;

# get all the features
my $features = $dbo->select( 'feature', ['*'] );

# go through each feature, get start, stop, chromosome, and strand and
# print to file
foreach my $f (@$features) {
    my $locations = $dbo->select(
        'location',
        [qw(start stop chromosome strand)],
        { feature_id => $f->{'feature_id'} },
        ['start asc']
    );
    print OUT $f->{'feature_id'}, "\t";
    print OUT $locations->[0]->{start},      "\t";
    print OUT $locations->[-1]->{stop},      "\t";
    print OUT $locations->[0]->{strand},     "\t";
    print OUT $locations->[0]->{chromosome}, "\n";

}
close(OUT);
