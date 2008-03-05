# -*- perl -*-
use strict;

use CoGeX;
my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'bpederse', 'brent_cnr');
#$s->storage->debug(1);

my $organism = $ARGV[0] or die "send in organism name i.e. $0 rice.\n";
print STDERR "making directory $organism";
mkdir($organism);

my ($genomic_sequence_type) = $s->resultset('GenomicSequenceType')->resolve('masked');

my ($org) = $s->resultset('Organism')->resolve($organism);

# MASKED !!!!!!!!!!!!!!!!!!
my $datasets = [sort map { $_->dataset_id } $org->current_datasets(genomic_sequence_type=>$genomic_sequence_type)];


print STDERR "usings datasets: " . join(", ", @$datasets) . " for $organism ...\n";


if (scalar($ARGV) > 1){
    print STDERR "fetching features...";
    my $feature_names_ids = get_feature_names_for_datasets($datasets, $organism);
    print STDERR "got " . scalar(@$feature_names_ids) . "\n";
    get_accn_locs($feature_names_ids);
}

foreach my $ds (@$datasets){
    my $ds = $s->resultset('Dataset')->resolve($ds);
    my %seen;
    open(FA, ">", $organism . ".fasta");
    foreach my $chr ($ds->get_chromosomes){
        # TODO: this will break with contigs.
        next if $chr =~ /^contig/;
        next if $chr =~ /random/;
        next if $chr =~ /scaffold/;

        print STDERR $chr . "\n" unless $seen{$chr};
        $seen{$chr} = 1;
        #  rice/chr01.fasta
        print FA "> $organism $chr\n";
        print FA $ds->get_genomic_sequence(chromosome => $chr);
    }
}
close FA;


sub get_accn_locs {
    my ($names_ids) = @_;
    my %seen;
    my %order;

    foreach my $name_id (@$names_ids){
        my ($name, $id, $start) = @$name_id;
        if($seen{$id}){ next; }
        my $feat = $s->resultset('Feature')->search( {
                           'me.feature_id'     => $id 
                        }, {
                             prefetch  =>  ['feature_type', 'locations']
                            ,order_by => 'feature_type.name'
                            ,limit    => 1
                        })->single(); 
        if(!$feat){ next; }
        my $chr = $feat->chromosome;
        my $chr_type = "chromosome";
        if ($chr =~ /super/i){
            $chr =~ s/contig//;
            $chr_type = "super";
        }
        elsif ($chr =~ /contig/i){
            $chr_type = "contig";
        }
        elsif ($chr =~ /random/i){
            #$chr_type = "random";
            next;

        }
        ($chr) = $feat->chromosome =~ /(\d+)/;
        $seen{$id}++;
        my $start = $feat->start;
        my $stop  = $feat->stop ;
        my $type  = uc($feat->feature_type->name);
        my $strand= $feat->strand;
        #print "$chr,$chr_type,$name,$type,$strand,$id";

        #print "$chr,$chr_type,NULL,$type,$strand,$id";
        print "$chr,$chr_type,$name,$type,$strand,$id";

        foreach my $l ($feat->locations()){
            print "," . $l->start, "," . $l->stop;
        }
        print "\n";
    }
}

sub get_feature_names_for_datasets {
    my $datasets = shift;
    my $notre = ',|\\-';
    my $org = shift;
    if( grep { $_ eq $org } ('rice', 'arabidopsis', 'grape', 'sorghum' )){
        $notre = ',|\\-|\\.';
    }
   
    my $rs = $s->resultset('FeatureName')->search( {
            'feature.dataset_id' => { 'IN' => $datasets }
            ,'me.name' => {'NOT REGEXP' => $notre }
            #,'me.primary_name' => 1
            #,'me.name' => { "LIKE" => 'Sb%g%'  }
            ,'feature_type.name' => { 'NOT LIKE' => '%contig%' }
            } , { 
               prefetch      =>  { 'feature' => 'feature_type' } 
               ,order_by => 'feature_type.name'  
               ,'distinct' => 'feature.feature_id'
               });

    my %seen;
    my @names;
    while(my $g = $rs->next()){
        #TODO: MAY BREAK WITH SOME.
        if($g->name =~/\s/){ next;}
        if($seen{$g->name}++){ next; }
        push(@names, [uc($g->name), $g->feature_id]);
    }
    return [sort { $a->[0] cmp $b->[0] } @names];
}

