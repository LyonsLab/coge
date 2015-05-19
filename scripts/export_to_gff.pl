#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my ( $coge, $fasta_name, $name_re, @datasets, $debug );

GetOptions(
    "fasta_name=s"                 => \$fasta_name,
    "name_search|name_re|name|n=s" => \$name_re,
    "dataset|ds=s"                 => \@datasets,
    "debug"                        => \$debug,
);
unless (@datasets) {
    print qq#
Welcome t0 $0

Usage:  $0 -dataset 24 -dataset NC_000001  -name_search regex_search -fasta_name output.faa
    or for multiple datasets
        $0 -ds 40504 -ds 40500 -ds 40495 -ds 40499 -ds 40503 -fasta_name arabidopsis_v9.fasta -name_re 'AT\\dG\\d{5}\$' > arabidopsis_v9.gff

Options:

 -dataset | -ds                        Datasets to retrieve, either by name or by database id

 -name_search | -name_re | -name | -n  OPTIONAL: Regular expression to search for specific feature names

 -fasta_name                           OPTIONAL:  create a fasta file of the sequences

 -debug                                OPTIONAL:  print debugging messages

#;
}

#$coge = CoGeX->dbconnect();
my $connstr = 'dbi:Pg:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

# gt sketch -seqid 1 -addintrons yes -start 1000 -style default.style -force -end 70000  out.png grape.gff3

my $chrs = get_locs( \@datasets );
my $schrs = get_sequence( ds => \@datasets, file_name => $fasta_name )
  if $fasta_name;

foreach my $k ( keys %$chrs ) {
    if ( !$schrs->{$k} ) {
        print STDERR "NOTFOUND in seqs:" . $k . "\n" if $debug;
    }
}
foreach my $k ( keys %$schrs ) {
    if ( !$chrs->{$k} ) {
        print STDERR "NOTFOUND in locs:" . $k . "\n" if $debug;
    }
}

sub get_locs {

    my $datasets = shift;
    my $SEP      = "\t";
    my %chrs;
    foreach my $ds (@$datasets) {
        my $dso = $coge->resultset('Dataset')->resolve($ds);
        foreach my $chr ( $dso->get_chromosomes ) {
            $chrs{$chr} = $dso->last_chromosome_position($chr);
        }
    }
    my @chrs = sort { $a cmp $b } keys %chrs;

    print "##gff-version\t3\n";
    foreach my $chr (@chrs) {
        print "##sequence-region $chr 1 " . $chrs{$chr} . "\n";
    }
    my %names = ();
    my %fids  = ();
    foreach my $chr (@chrs) {
        my %seen = ();

        #if ($i++ > 2){ print STDERR "*" x 1000 . "\nENDING EARLY";  last; }

        my $gene_rs = $coge->resultset('Feature')->search(
            {
                'me.dataset_id'     => { 'IN' => $datasets },
                'me.chromosome'     => $chr,
                'feature_type.name' => 'gene'
            },
            {
                'prefetch' => [ 'feature_type', 'feature_names' ],
                'order_by' => ['me.start']
            }
        );

        #gff: chr  organization feature_type  start stop strand . name
        my %chrs;
        print STDERR "dataset_ids: "
          . join( ",", @$datasets )
          . ";  chr: $chr\n"
          if $debug;
        while ( my $g = $gene_rs->next() ) {
            if ( $fids{ $g->feature_id } ) { next; }
            $fids{ $g->feature_id } = 1;
            my @gene_names;
            if ($name_re) {
                @gene_names =
                  grep { $_->name =~ /$name_re/i } $g->feature_names();
            }
            else {
                @gene_names = $g->feature_names();
            }
            my $gene_name;
            if ( scalar(@gene_names) == 0 ) {
                $gene_name = $g->feature_names()->next()->name;
                print STDERR "not from re:" . $gene_name . "\n" if $debug;
                print STDERR $name_re . "\n" if $debug;
            }
            else {
                foreach my $name (@gene_names) {
                    $gene_name = $name->name;
                    if ( !$names{$gene_name} ) { last; }
                    $gene_name = "";
                }

                $gene_name = $gene_names[0]->name;
            }
            if ( !$gene_name || $names{$gene_name} ) { next; }

            my $mrna_rs = $coge->resultset('Feature')->search(
                {
                    'me.dataset_id'      => { 'IN' => $datasets },
                    'me.chromosome'      => $chr,
                    'feature_names.name' => $gene_name,
                    'feature_type.name'  => 'mRNA'
                },
                {
                    'join'     => 'feature_names',
                    'prefetch' => [ 'feature_type', 'locations' ],
                    'order_by' => [ 'me.start', 'locations.start' ]
                }
            );

            #my $mrna = $mrna_rs->next();

            my $strand = $g->strand == 1 ? '+' : '-';
            my $clean_name = $gene_name;
            $names{$gene_name} = 1;
            $clean_name =~ s/\s+/_/g;
            my $attrs = "ID=$clean_name;Name=$clean_name;rname=$clean_name";

            my $gstr = join(
                "\t",
                (
                    $chr, 'ucb', $g->feature_type->name, $g->start, $g->stop,
                    ".", $strand, ".", $attrs
                )
            );
            if ( $seen{$gstr} ) { next; }
            $seen{$gstr} = 1;

            print $gstr . "\n";

            my $parent   = $clean_name;
            my $has_mrna = 0;
            while ( my $f = $mrna_rs->next() ) {
                if ( $fids{ $f->feature_id } ) { next; }
                $fids{ $f->feature_id } = 1;
                $attrs    = "Parent=$parent;ID=$parent" . ".mRNA;rname=$parent";
                $has_mrna = 1;

                foreach my $loc ( $f->locations() ) {
                    my $gstr = join(
                        "\t",
                        (
                            $f->chr, 'ucb', $f->feature_type->name, $loc->start,
                            $loc->stop, ".", $strand, ".", $attrs
                        )
                    );
                    if ( $seen{$gstr} ) { next; }
                    $seen{$gstr} = 1;
                    print $gstr . "\n";

                }
            }
            $attrs = "Parent=$clean_name";
            if ($has_mrna) {
                $attrs .= ".mRNA";
            }
            $attrs .= ";rname=$clean_name";    # keep it in another attr

#print join("\t", ($chr, 'ucb', 'mRNA', $mrna->start, $mrna->stop, ".", $strand, ".", $attrs)) . "\n";
            $chrs{ $g->chr } = 1;
            my $sub_rs = $coge->resultset('Feature')->search(
                {
                    'me.dataset_id'      => { 'IN' => $datasets },
                    'feature_names.name' => $gene_name,
                    'feature_type.name' => { 'NOT IN' => [ 'gene', 'mRNA' ] }
                },
                {
                    'join'     => ['feature_names'],
                    'prefetch' => [ 'feature_type', 'locations' ],
                    'order_by' => [ 'me.chromosome', 'me.start' ]
                }
            );
            while ( my $f = $sub_rs->next() ) {

                #my $locs = $f->locations({}, {'order_by' => 'start'});
                foreach my $loc ( $f->locations() ) {
                    my $gstr = join(
                        "\t",
                        (
                            $f->chr, 'ucb', $f->feature_type->name, $loc->start,
                            $loc->stop, ".", $strand, ".", $attrs
                        )
                    );
                    if ( $seen{$gstr} ) { next; }
                    $seen{$gstr} = 1;
                    print $gstr . "\n";

                }

            }
        }
    }
    return \%chrs;
}

sub get_sequence {
    my %opts      = @_;
    my $datasets  = $opts{ds};
    my $file_name = $opts{file_name};
    open( FA, ">", $file_name );
    my %chrs;
    foreach my $ds (@$datasets) {
        my $ds = $coge->resultset('Dataset')->resolve($ds);
        my %seen;
        foreach my $chr ( $ds->get_chromosomes ) {

            # TODO: this will break with contigs.
            next if $chr =~ /random/;
            $chrs{$chr} = 1;
            print FA ">$chr\n";
            print FA $ds->get_genomic_sequence( chromosome => $chr ) . "\n";
        }
    }
    close FA;
    return \%chrs;
}
