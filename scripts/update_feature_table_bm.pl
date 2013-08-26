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

open( IN, "</snb/bio3/work/temp/x.txt" ) or die;
while (<IN>) {
    chomp;
    my @temp = split( /\s+/, $_ );
    if ( scalar @temp == 5 ) {
        $dbo->update(
            'feature',
            {
                fstart      => $temp[1],
                fstop       => $temp[2],
                fstrand     => $temp[3],
                fchromosome => $temp[4]
            },
            { feature_id => $temp[0] }
        );
    }
}
close(IN);
