#! /usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;

my @chrs = qw/I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII XVIII XIX/;
my %cmap;
my $i = 1;
map { $cmap{$chrs[$_ -1 ]} = ($_ < 10) ? '0' . $_ : $_  } 1..19;

my %FH;

foreach my $chr (values %cmap){
    open($FH{$chr},">","/tmp/poplar_/poplarchr" .$chr . ".fasta");
}


# 9min 9 sec
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' );


$s->storage->debug(0);

#my $org = $s->resultset('Dataset')->search(
#    { 'organism.name' => {like => "%Poplar%" },
#      'version'       => '1.1' 
#    },
#    { prefetch => 'organism' }
#);
#
#my $did = $org->next()->dataset_id;

print "\n";
my $rs = $s->resultset('Feature')->esearch(
        { 
            'me.dataset_id' => 505,
            'feature_names.name' =>  {like => '%proteinId%'},
            'feature_type.name' => 'CDS'
        },
        {
            #-group_by => ['feature_names.name'], 
            join => 'feature_names'
        });


my $currentpid = '';
my $sequence = "";
my $currentgene = 0;
my $pid = 1;

while( my $feat = $rs->next()){
    my $chrnum = $cmap{$feat->chromosome};
    my $name = ($feat->names())[0];
    $sequence = $feat->genome_sequence();

    my $genenum = "0" x (6 - length($currentgene)) . $currentgene;

    my $fh = *STDOUT; # $FH{$chrnum};

    my $gname = "P" . $chrnum . "G" .  $genenum;

    #$feat->feature_na
    print $fh ">",$gname , "\n";
    print $feat->feature_type->name . "\n";
    print $fh length($sequence), "\n\n";
    ++$currentgene;

}
