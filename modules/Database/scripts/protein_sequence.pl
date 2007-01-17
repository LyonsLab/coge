#! /usr/bin/perl -w

use strict;
use CoGeX;

#time
#mysql before update: real    0m54.935s
#mysql after  update: real    0m44.592s
#postgresql         : real    0m30.214s


#my $connstr = 'dbi:mysql:genomes:biocon:3306';
#my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

my $connstr = 'dbi:Pg:dbname=genomes;host=biocon;port=5432';
my $s = CoGeX->connect($connstr, 'bpederse', 'wsa47r' );

$s->storage->debug(1);

my $rs = $s->resultset('Feature')->search(
        { 
            'feature_type.name' =>   'CDS' ,
            'feature_names.name' => {like => 'At1g%','-not_like' => "%.%" }
        },
        {
            join => ['feature_names','feature_type'],
            prefetch => ['feature_names','feature_type']
        }

);


while (my $feat =$rs->next()){
    my $fn = $feat->feature_names;
    my $type = $feat->feature_type->name;
    print STDERR ">";
    map { print STDERR $_->name . ":". $type . "\t" } $fn->next();
    print STDERR "\n";

    # this prefetch avoids n calls where n is number of sequences #
    foreach my $seq ($feat->sequences({},{prefetch=>"sequence_type"})){
        print STDERR $seq->sequence_data if $seq->sequence_type->name =~ /protein/i;
    }
    print STDERR "\n\n";
}


