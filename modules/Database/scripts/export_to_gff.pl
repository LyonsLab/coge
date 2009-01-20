use strict;
use Data::Dumper;
use CoGeX;


my $connstr = 'dbi:mysql:genomes:biocon:3306';
my $s = CoGeX->connect($connstr, 'bpederse', 'brent_cnr');

if(scalar(@ARGV) < 2){
    print "send in organism and datasets\n";
    print "$0 org [datasets]\n";
    exit();
}

my $organism = shift;
my $datasets = \@ARGV;

# gt sketch -seqid 1 -addintrons yes -start 1000 -style default.style -force -end 70000  out.png grape.gff3

my $chrs = get_locs($datasets);
my $schrs = get_sequence($datasets);

foreach my $k (keys %$chrs){
    if (! $schrs->{$k}){
        print STDERR "NOTFOUND in seqs:" . $k . "\n";
    }
}
foreach my $k (keys %$schrs){
    if (! $chrs->{$k}){
        print STDERR "NOTFOUND in locs:" . $k . "\n";
    }
}


sub get_locs {
    my $SEP = "\t";
    my $datasets = shift;

    my %chrs;
    foreach my $ds (@$datasets){
        my $dso = $s->resultset('Dataset')->resolve($ds);
        foreach my $chr ($dso->get_chromosomes){
            $chrs{$chr} = $dso->last_chromosome_position($chr);
        }
    }
    my @chrs = sort { $a cmp $b } keys %chrs;

    print "##gff-version\t3\n";
    foreach my $chr (@chrs){
        print "##sequence-region $chr 1 " . $chrs{$chr} . "\n";
    }
    my %seen = {};
    foreach my $chr (@chrs){
        
        my $gene_rs = $s->resultset('Feature')->search( {
                  'me.dataset_id' => { 'IN' => $datasets },
                  'me.chromosome' => $chr,
                  'feature_type.name'  =>  'gene' 
                } , { 
                   'prefetch'           => [ 'feature_type', 'feature_names'] 
                  ,'order_by'           => [ 'me.start']
                });

        #gff: chr  organization feature_type  start stop strand . name
        my %chrs;
        print STDERR $chr . "\n";
        print STDERR join(",", @$datasets) . "\n";
        while(my $g = $gene_rs->next()){
            my $gene_name = $g->feature_names()->next()->name;

            my $mrna_rs = $s->resultset('Feature')->search( {
                  'me.dataset_id' => { 'IN' => $datasets },
                  'me.chromosome' => $chr,
                  'feature_names.name' => $gene_name,
                  'feature_type.name'  =>  'mRNA' 
                } , { 
                    'join' => 'feature_names',
                   'prefetch'           => [ 'feature_type', 'feature_names'] 
                  ,'order_by'           => [ 'me.start']
                });

            #my $mrna = $mrna_rs->next();


            my $strand = $g->strand == 1 ? '+' : '-';
            my $clean_name = $gene_name;
            $clean_name =~ s/\s+/_/g;
            my $attrs = "ID=$clean_name;Name=$clean_name";
            my $gstr = join("\t", ($chr, 'ucb', $g->feature_type->name, $g->start, $g->stop, ".", $strand, ".", $attrs));
            if($seen{$gstr}){ next; }
            $seen{$gstr} = 1;
            print $gstr . "\n";
            my $parent = $clean_name;
            my $has_mrna = 0;
            while(my $f = $mrna_rs->next()){
                $attrs = "Parent=$parent;ID=$parent" . ".mRNA";
                my $has_mrna = 1;
                my $locs = $f->locations({}, {'order_by' => 'start'});
                while(my $loc = $locs->next()){
                    my $gstr = join("\t", ($f->chr, 'ucb', $f->feature_type->name, $loc->start, $loc->stop, ".", $strand, ".", $attrs));
                    if($seen{$gstr}){ next; }
                    $seen{$gstr} = 1;
                    print $gstr . "\n";
            
                }
            }
            if($has_mrna){
                $attrs = "Parent=$clean_name" . "mRNA";
            }
            else {
                $attrs = "Parent=$clean_name";
            }

            #print join("\t", ($chr, 'ucb', 'mRNA', $mrna->start, $mrna->stop, ".", $strand, ".", $attrs)) . "\n";
            $chrs{$g->chr} = 1;
            my $sub_rs = $s->resultset('Feature')->search( {
                  'me.dataset_id' => { 'IN' => $datasets },
                  'feature_names.name'  =>  $gene_name
                  , 'feature_type.name'  =>  { 'NOT IN' => ['gene', 'mRNA'] }
                } , { 
                   'join'               => [ 'feature_names'],
                   'prefetch'           => [ 'feature_type', 'locations'] 
                  ,'order_by'           => [ 'me.chromosome', 'me.start']
                });
            while(my $f = $sub_rs->next()){
                my $locs = $f->locations({}, {'order_by' => 'start'});
                while(my $loc = $locs->next()){
                    my $gstr = join("\t", ($f->chr, 'ucb', $f->feature_type->name, $loc->start, $loc->stop, ".", $strand, ".", $attrs));
                    if($seen{$gstr}){ next; }
                    print $gstr . "\n";
            
                }

            }
        }
    }
    return \%chrs;
}


sub get_sequence {
    open(FA, ">", $organism . ".fasta");
    my %chrs;
    foreach my $ds (@$datasets){
        my $ds = $s->resultset('Dataset')->resolve($ds);
        my %seen;
        foreach my $chr ($ds->get_chromosomes){
            # TODO: this will break with contigs.
            #next if $chr =~ /^contig/;
            next if $chr =~ /random/;
            
            #$chr =~ s/scaffold/super/g; print STDERR "CHANGING scaffold => super\n";

            #next if $chr =~ /scaffold/;

            #print STDERR $chr . "\n" unless $seen{$chr};
            $seen{$chr} = 1;
            $chrs{$chr} = 1;

            #  rice/chr01.fasta
            print FA "> $chr\n";
            print FA $ds->get_genomic_sequence(chromosome => $chr) . "\n";
        }
    }
    close FA;
    return \%chrs;
}
