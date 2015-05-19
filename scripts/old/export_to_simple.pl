#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;

my ( $coge, $fasta_name, $name_re, @datasets, $dataset_group );

GetOptions(
    "fasta_name=s"        => \$fasta_name,
    "name_re=s"           => \$name_re,
    "dataset_group|dsg=s" => \$dataset_group,
);
unless ($dataset_group) {
    print qq#
Welcome t0 $0

Usage:  $0 -dsg 8120 -fasta_name brachy_v1.fasta > brachy_v1.bed

Options:

 -dataset_group | -dsg                 Datasets to retrieve, either by name or by database id
 -name_re                          -n  OPTIONAL: Regular expression to search for specific feature names
 -fasta_name                           OPTIONAL:  create a fasta file of the sequences

#;
}

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my $DSG = $coge->resultset('DatasetGroup')->resolve($dataset_group);
@datasets = map { $_->dataset_id } $DSG->datasets;

get_locs( \@datasets );
my $schrs = get_sequence( ds => \@datasets, file_name => $fasta_name )
  if $fasta_name;

sub get_name {
    my $g       = shift;
    my $name_re = shift;
    my @gene_names;
    if ($name_re) {
        @gene_names = grep { $_->name =~ /$name_re/i } $g->feature_names();
    }
    else {
        @gene_names = $g->feature_names();
    }
    my $gene_name;
    if ( scalar(@gene_names) == 0 ) {
        $gene_name = $g->feature_names()->next()->name;
        print STDERR "not from re:" . $gene_name . "\n";
        print STDERR $name_re . "\n";
    }
    else {
        $gene_name = $gene_names[0]->name;
    }
    return $gene_name;
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
    my @chrs  = sort { $a cmp $b } keys %chrs;
    my %names = ();
    my %fids  = ();

    foreach my $chr (@chrs) {
        my $rs = $coge->resultset('Feature')->search(
            {
                'me.dataset_id' => { 'IN' => $datasets },
                'me.chromosome' => $chr,
                'feature_type.name' =>
                  { 'NOT IN' => [ 'chromosome', 'scaffold' ] }
            },
            {
                'prefetch' => [ 'feature_type', 'feature_names' ],
                'order_by' => [ 'me.start',     '-feature_names.name' ]
            }
        );

        my %seen = ();
        while ( my $g = $rs->next() ) {
            if ( $fids{ $g->feature_id } ) { next; }
            $fids{ $g->feature_id } = 1;
            my $ftype     = $g->feature_type->name;
            my $gene_name = get_name( $g, $name_re );
            my $strand    = $g->strand == 1 ? '+' : '-';
            foreach my $loc ( $g->locations() ) {
                my $start = $loc->start;
                my $end   = $loc->stop;
                my $key   = $start . "|" . $end . "|" . $ftype;
                if ( $seen{$key} ) { next; }
                $seen{$key} = 1;
                print "$chr,$ftype,$start,$end,$strand,$gene_name\n";
            }
        }
    }
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
