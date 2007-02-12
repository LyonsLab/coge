#! /usr/bin/perl -w

use strict;
use CoGeX;


# 9min 9 sec
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );


$s->storage->debug(0);

my $org = $s->resultset('Dataset')->search(
    { 'organism.name' => {like => "%Poplar%" },
      'version'       => '1.1' 
    },
    { join => 'organism' }
);
my $did = $org->next()->dataset_id;


my $rs = $s->resultset('Feature')->search(
        { 
            'dataset_id' => $did,
#            'feature_names.name' =>  {like => 'P%%CDS'}, 
        },
        {
            join => ['feature_names','feature_type','sequences'],
            prefetch => ['feature_names','feature_type'],
            order_by => ['me.feature_id'],

        }
);

print "result count: " . $rs->count() . "\n";


my %FH;
foreach my $i (1..19){
    my $chr = $i < 10 ? "0". $i : $i;
    print $chr . "\n";
    open($FH{$chr},">/tmp/poplar/poplarchr$chr" . ".fasta");
}


while( my $feature = $rs->next()){
    my $name;
    #my $name = (grep { $_->name =~ /P\d+CDS/ } $feature->feature_names)[0]->name;
    map { print $_->name . "\t" } $feature->feature_names ;
    print "\n";
    next;
    print $name . "\n";
    my $chr = substr($name,1,2);
    my $fh = $FH{$chr};
    print $fh  ">" . $name . "\n";
    print $fh $feature->genome_sequence . "\n\n";
}
