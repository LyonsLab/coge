#! /usr/bin/perl -w

use strict;
use CoGeX;

#time
#real    0m54.935s



my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$s->storage->debug(1);

my $rs = $s->resultset('Feature')->search(
        { 
            'feature_names.name' => {like => 'At1%','-not_like' => "%.%" },
            'feature_type.name' =>   'CDS' 
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


