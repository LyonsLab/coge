# -*- perl -*-
use strict;

use CoGeX;
use Getopt::Std;
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'bpederse', 'brent_cnr');
use Data::Dumper;
use Memoize;
use DB_File;

my %options;
getopt("otd", \%options);

my %org_hash = (
   rice        => 3
  ,arabidopsis => 1
  ,poplar      => 3387
  ,grape       => 351
  ,papaya      => 322
  ,maize       => 333
);

use Memoize;
use Memoize::Expire;
use Data::Dumper;
tie my %cache => 'Memoize::Expire', 'LIFETIME' => 20;
tie my %cachelist => 'Memoize::Expire', 'LIFETIME' => 10;

memoize(\&CoGeX::Dataset::genomic_sequences
          , NORMALIZE    => sub { print Dumper %cache; shift->genomic_sequence_id; }
          , SCALAR_CACHE => [ 'HASH', \%cache]
          , LIST_CACHE   => [ 'HASH', \%cachelist]);

tie my %persistent => 'DB_File', '/tmp/getgene2s.cache' , O_RDWR|O_CREAT, 0666;
memoize('get_feature_names_for_datasets', SCALAR_CACHE => [HASH => \%persistent], LIST_CACHE => 'MERGE', NORMALIZE => sub { join " ", @_ } );

my $organism = $options{o} or die "send in organism name i.e. -o rice.\n";
my $datasets = [sort map { $_->dataset_id  } @{$s->get_current_datasets_for_org($org_hash{$organism})}];
my $outdir   = ($options{d} or ".") . "/";


print STDERR "usings datasets: " . join(",", @$datasets) . " for $organism ...\n";

my $feature_names_ids = get_feature_names_for_datasets($datasets);
print STDERR "got " .  scalar(@$feature_names_ids) . " feature names\n";

if($options{t}){
    get_10kmers($organism, $datasets);
    exit();
}

get_accn_locs($organism, $datasets, $feature_names_ids);

sub get_accn_locs {
    my ($org, $datasets, $names_ids) = @_;
    my %seen;
    my %order;

    my %files;
    foreach my $name_id (@$names_ids){
        my ($name, $id, $start) = @$name_id;
        my $feat = $s->resultset('Feature')->search( {
                           'me.feature_id'     => $id 
                        ,  'me.chromosome' => {'NOT LIKE' => 'contig%' }
                        , 'feature_names.name' => $name
                        }, {
                             prefetch  =>  ['feature_type', 'feature_names']
                            ,order_by => 'feature_type.name'
                            ,limit    => 1
                        })->single(); 
        if(!$feat){ next; }
        my ($chr) = $feat->chromosome =~ /(\d+)/;
        if(length($chr) > 2){ next; } 
        #print STDERR $feat->feature_id . ", $name," .  $feat->chromosome . "," . $feat->feature_type->name . "\n";
        $chr = sprintf("%02i", $chr);
        if(!$files{$chr}){
            my $FH;
            my $filename = $outdir . $org . $chr . ".fasta";
            open($FH, ">", $filename);
            $files{$chr} = $FH;
            print STDERR "creating file $filename\n";
        }
        my $FH = $files{$chr};

        my $start = sprintf("%09i", $feat->start);
        my $stop  = sprintf("%09i", $feat->stop );

        my $header = $chr . "||" . $name . "||" . $start . "||" . $stop 
               . "||". $feat->strand . "||" . uc($feat->feature_type->name) 
               . "||". $feat->feature_id;

        $order{$header} = 1;
        print $FH  ">" . $header . "\n";

        print $FH $feat->genomic_sequence() . "\n";
    }
    print STDERR "creating file " . $org . ".order contained the list of genes in order ...\n";
    open(ORDER, ">", $outdir . $org . ".order");
    map { print ORDER $_ . "\n" } sort keys %order;
    map { close $_ } values %files;
    close(ORDER);
}



sub get_10kmers {
    my ($org, $datasets) = @_;
    my %files;
    my %order;
    foreach my $gs ( $s->resultset('GenomicSequence')->search(
                        { 'dataset_id'   => {'IN' => $datasets } 
                          , 'chromosome' => {'NOT LIKE' => 'contig%' }
                        }, { order_by => ['start'] } )) {
            my ($chr) = $gs->chromosome =~ /(\d+)/;
            if(length($chr) > 2){ next; } 
            my $file = $outdir . $org . "10kmers_chr" . $chr . ".fasta";
            my $FH;
            if (!$files{$file}) {
                print STDERR "creating file $file ...\n";
                open($FH, ">", $file);
                $files{$file} = $FH;
            }else{
                $FH = $files{$file};
            }
            my $header = sprintf("%02i",$chr) . "||NONE||" . sprintf("%09i", $gs->start());
            print $FH ">" . $header . "\n";
            print $FH $gs->_sequence_data() ."\n";

            $order{$header} = 1;
    }
    print STDERR "creating file " . $org . ".order contained the list of 10kmers ...\n";
    open(ORDER, ">", $outdir . $org . ".order");
    map { print ORDER $_ . "\n" } sort keys %order;
    close(ORDER);

}


sub get_feature_names_for_datasets {
    my $datasets = shift;
   
    my $rs = $s->resultset('FeatureName')->search( {
            'feature.dataset_id' => { 'IN' => $datasets }
            ,'me.name' => {'NOT REGEXP' => ',|\\-' }
            ,'feature.chromosome' => { 'NOT LIKE' => 'contig%' } # keep only chrs and super_contigs
            ,'feature_type.name' => { 'NOT LIKE' => '%contig%' }
            } , { 
               prefetch      =>  { 'feature' => 'feature_type' } 
               ,order_by => 'feature_type.name'  
               });
    my %seen;
    #my @names = grep { $_->[0] !~ /\W/ and !$seen{$_->[0]}++ } map { [uc($_->name), $_->feature_id] } $rs->all();
    my @names = grep { !$seen{$_->[0]}++ } map { [uc($_->name), $_->feature_id, $_->feature->start] } $rs->all();
    return [sort { $a->[2] cmp $b->[2] } @names];
}

