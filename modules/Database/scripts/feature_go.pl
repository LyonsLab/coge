#!/usr/bin/perl -w

use CoGeX;
use strict;
use Data::Dumper;
$| = 1;

my $connstr = 'dbi:Pg:dbname=genomes;host=biocon;port=5432';
my $s = CoGeX->connect($connstr, 'bpederse', 'wsa47r' );


$connstr = 'dbi:mysql:genomes:biocon:3306';
$s = CoGeX->connect($connstr, 'cnssys', 'CnS' );

$s->storage->debug(1);


my $rs = $s->resultset('Feature')->search( { 
                    'feature_type.name' =>   'CDS' ,
                    'dataset.version'    => 6,
                    'feature_names.name' => ['-and' , 
                            { '-like' => 'At%'} ,
                            { '-not_like' => ['%.%'] }
                        ]
                    }, {
                        join => [ 'feature_names','feature_type','dataset'] ,
                        prefetch => [ 'feature_type','dataset'] , 
                    }, );

print "got resultset\n";

my %seen;
$s->storage->debug(0);
my %annotations;
while (my $feat = $rs->next()){
    my $fn = $feat->feature_names;
    my $type = $feat->feature_type->name;

    my $annos = $feat->annotations;
    my @names = grep { $_->name =~ /^At\dg\d{5}$/i } $fn->next();
    next unless @names;
    #print join("|", sort map { $_->name } @names) . "\n";
    my @gos = ();
    while (my $a = $annos->next()){
        next unless $a->annotation_type;
        next unless $a->annotation_type->annotation_type_group;
        next unless $a->annotation_type->annotation_type_group->name =~/^[Gg][Oo]/;
        #print $a->annotation_type->annotation_type_group->name  . "\n";
        push(@{$annotations{$names[0]->name}}, $a->annotation);
    }
    #map { print $_->annotation .  "\n" }  @gos;
}
foreach my $accn (keys %annotations){
    print $accn . "|" . join("|", @{$annotations{$accn}} ) . "\n";
}
   

