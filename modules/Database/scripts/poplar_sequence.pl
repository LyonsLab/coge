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
    open($FH{$chr},">","/tmp/poplar/poplarchr" .$chr . ".fasta");
}


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
            'feature_names.name' =>  {like => '%proteinId%'}, 
        },
        {
            join => ['feature_names','locations'],
            prefetch => ['feature_names'],
            order_by => ['me.feature_id'],

        }
);

my $currentpid = '';
my $sequence = "";
my $currentgene = 1;


while( my $feat = $rs->next()){
    my $name = ($feat->feature_names())[0];
    $name = $name->name;
    my ($pid) = $name =~ /(proteinId\s\d+)/;
    $currentpid ||=$pid;
    if ($pid eq $currentpid){
        $sequence .= $feat->genome_sequence();
        next;
    }else{
        my $genenum = "0" x (6 - length($currentgene)) . $currentgene;
        my $chrnum = $cmap{(($feat->locations())[0])->chromosome};
        my $fh = $FH{$chrnum};
        my $gname = "P" . $chrnum . "G" .  $genenum;
        print $fh ">",$gname , "\n";
        print $fh $sequence, "\n\n";
        $sequence = $feat->genome_sequence();
        $currentpid = $pid;
        ++$currentgene;
    }

}

while( my $feature = $rs->next()){
    my $name;
    #my $name = (grep { $_->name =~ /P\d+CDS/ } $feature->feature_names)[0]->name;
    map { print $_->name . "\t" } $feature->feature_names ;
    print "\n";
    next;
#    print $name . "\n";
#    my $chr = substr($name,1,2);
#    my $fh = $FH{$chr};
#    print $fh  ">" . $name . "\n";
#    print $fh $feature->genome_sequence . "\n\n";
}
