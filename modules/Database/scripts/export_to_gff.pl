use strict;
use Data::Dumper;
use CoGeX;


my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'bpederse', 'brent_cnr');

my $organism = shift;
my $datasets = \@ARGV;


get_locs($datasets);
#get_sequence($datasets);


sub get_locs {
    my $SEP = "\t";
    my $datasets = shift;
    my $gene_rs = $s->resultset('Feature')->search( {
              'me.dataset_id' => { 'IN' => $datasets },
              'feature_type.name'  =>  'gene' 
            } , { 
               'prefetch'           => [ 'feature_type', 'feature_names'] 
              ,'order_by'           => [ 'me.chromosome', 'me.start']
            });

    #gff: chr  organization feature_type  start stop strand . name
    print "##gff-version\t3\n";
    while(my $g = $gene_rs->next()){
        my $gene_name = $g->feature_names()->next()->name;
        my $strand = $g->strand == 1 ? '+' : '-';
        my $attrs = "ID=$gene_name;Name=$gene_name";
        print join("\t", ($g->chr, 'ucb', $g->feature_type->name, $g->start, $g->stop, ".", $strand, ".", $attrs)) . "\n";
        my $sub_rs = $s->resultset('Feature')->search( {
              'me.dataset_id' => { 'IN' => $datasets },
              'feature_names.name'  =>  $gene_name
              , 'feature_type.name'  =>  { '!=' => 'gene' }
            } , { 
               'join'               => [ 'feature_names'],
               'prefetch'           => [ 'feature_type', 'locations'] 
              ,'order_by'           => [ 'me.chromosome', 'me.start']
            });
        while(my $f = $sub_rs->next()){
            my $locs = $f->locations({}, {'order' => 'start'});
            $attrs = "Parent=$gene_name;Name=$gene_name";
            while(my $loc = $locs->next()){
                print join("\t", ($f->chr, 'ucb', $f->feature_type->name, $loc->start, $loc->stop, ".", $strand, ".", $attrs)) . "\n";
        
            }

        }
    }
}


sub get_sequence {
    open(FA, ">", $organism . ".fasta");
    foreach my $ds (@$datasets){
        my $ds = $s->resultset('Dataset')->resolve($ds);
        my %seen;
        foreach my $chr ($ds->get_chromosomes){
            # TODO: this will break with contigs.
            #next if $chr =~ /^contig/;
            next if $chr =~ /random/;
            
            #$chr =~ s/scaffold/super/g; print STDERR "CHANGING scaffold => super\n";

            #next if $chr =~ /scaffold/;

            print STDERR $chr . "\n" unless $seen{$chr};
            $seen{$chr} = 1;
            #  rice/chr01.fasta
            print FA "> $chr\n";
            print FA $ds->get_genomic_sequence(chromosome => $chr) . "\n";
        }
    }
    close FA;
}
