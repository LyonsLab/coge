#!/usr/bin/perl -w
# -*- perl -*-
use strict;

use CoGeX;
use Getopt::Std;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $s = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
use Data::Dumper;
use DB_File;

# for orgs with lots of super contigs, only take this many
# (will take the first ones).
my $use_supers = 1;
my $MAX_CHR    = 99999;

my %options;
getopt( "otdcml", \%options );

my %ds_hash = ( 'rice' => [ 582 .. 593 ], );

my $organism = $options{o} or die "send in organism name i.e. -o rice.\n";
chomp $organism;
print STDERR $organism;
my $use_contigs = $options{c} or 0;
my ($org) = $s->resultset('Organism')->resolve($organism);
my $datasets;
my $chr_size_limit = $options{l};
if ( $options{m} ) {
    print STDERR "getting masked...\n";
    my ($genomic_sequence_type) =
      $s->resultset('GenomicSequenceType')->resolve('masked');
    $datasets = [
        sort map { $_->dataset_id } $org->current_datasets(
            genomic_sequence_type => $genomic_sequence_type
        )
    ];
}
else {
    print STDERR "getting NOT masked...\n";
    $datasets = [ sort map { $_->dataset_id } $org->current_datasets() ];
}

my $outdir = ( $options{d} or "." ) . "./";
chomp $outdir;

print STDERR "usings datasets: "
  . join( ",", @$datasets )
  . " for $organism ...\n";

if ( defined $options{t} ) {
    get_10kmers( $organism, $datasets, $use_contigs );
    exit();
}

my $feature_names_ids =
  get_feature_names_for_datasets( $datasets, $organism, $use_contigs );
print STDERR "got " . scalar(@$feature_names_ids) . " feature names\n";
get_accn_locs( $organism, $datasets, $feature_names_ids, $use_contigs );

sub get_accn_locs {
    my ( $org, $datasets, $names_ids, $use_contigs ) = @_;
    my %seen;
    my %order;

    my %files;
    foreach my $name_id (@$names_ids) {
        my ( $name, $id, $start ) = @$name_id;
        if ( $seen{$id} ) { next; }
        my $feat = $s->resultset('Feature')->search(
            {
                'me.feature_id'      => $id,
                'me.chromosome'      => { 'NOT LIKE' => 'contig%' },
                'feature_names.name' => $name
            },
            {
                prefetch => [ 'feature_type', 'feature_names' ],
                order_by => 'feature_type.name',
                limit    => 1
            }
        )->single();
        if ( !$feat ) { next; }
        if ( $feat->feature_type->name eq 'CNS' ) { next; }

        #print STDERR $feat->chromosome . "\n";
        if ( !$use_contigs && $feat->chromosome =~ /(^contig|random)/i ) {
            next;
        }
        if ( !$use_supers && $feat->chromosome =~ /super/i ) { next; }
        my ($chr) = $feat->chromosome =~ /(\d+)/;

        #print STDERR $chr . "\n";
        #if(length($chr) > 2){ next; }
        if ( $chr > $MAX_CHR ) { next; }
        $seen{$id}++;

#print STDERR $feat->feature_id . ", $name," .  $feat->chromosome . "," . $feat->feature_type->name . "\n";
        $chr = sprintf( "%04i", $chr );
        my $FH;
        if ( !$files{$chr} ) {
            my $filename = $outdir . $org . "chr" . $chr . ".fasta";
            open( $FH, ">", $filename );
            $files{$chr} = $FH;
            print STDERR "creating file $filename\n";
        }
        $FH = $files{$chr};

        my $start = sprintf( "%09i", $feat->start );
        my $stop  = sprintf( "%09i", $feat->stop );

        my $header =
            $chr . "||"
          . $start . "||"
          . $stop . "||"
          . $name . "||"
          . $feat->strand . "||"
          . uc( $feat->feature_type->name ) . "||"
          . $feat->feature_id;

        $order{$header} = 1;
        print $FH ">" . $header . "\n";

        print $FH $feat->genomic_sequence() . "\n";
    }
    print STDERR "creating file " . $org
      . ".order contained the list of genes in order ...\n";
    open( ORDER, ">", $outdir . $org . ".order" );
    map { print ORDER $_ . "\n" } sort keys %order;
    map { close $_ } values %files;
    close(ORDER);
}

sub get_10kmers {
    my ( $org, $datasets, $use_contigs ) = @_;
    my %files;
    my %order;
    foreach my $gs (
        $s->resultset('GenomicSequence')->search(
            { 'dataset_id' => { 'IN' => $datasets } },
            { order_by     => ['start'] }
        )
      )
    {
        if ( !$use_contigs && $gs->chromosome =~ /(^contig|random)/i ) { next; }
        my ($chr) = $gs->chromosome =~ /(\d+)/;

        #if(length($chr) > 2){ next; }
        if ( $chr > $MAX_CHR ) { next; }

        $chr = sprintf( "%04i", $chr );
        my $file = $outdir . $org . "tenkmers_chr" . $chr . ".fasta";
        my $header =
            $chr . "||"
          . sprintf( "%09i", $gs->start() ) . "||"
          . sprintf( "%09i", $gs->stop() )
          . "||10KMER";

        if (0) {
            my $FH;
            if ( !$files{$file} ) {
                if ( scalar( keys %files ) > 900 ) {
                    my $k = ( sort keys %files )[600];
                    delete $files{$k};
                }
                print STDERR "creating file $file ...\n";
                open( $FH, ">>", $file );
                $files{$file} = $FH;
            }
            else {
                $FH = $files{$file};
            }
            print $FH ">" . $header . "\n";
            print $FH $gs->sequence_data() . "\n";
        }
        else {
            print ">" . $header . "\n";
            print $gs->sequence_data() . "\n";
        }

        #$order{$header} = 1;
    }

#print STDERR "creating file " . $org . ".order contained the list of 10kmers ...\n";
#open(ORDER, ">", $outdir . $org . "tenkmers.order");
#map { print ORDER $_ . "\n" } sort keys %order;
#close(ORDER);

}

sub get_feature_names_for_datasets {
    my $datasets    = shift;
    my $notre       = ',|\\-';
    my $org         = shift;
    my $use_contigs = shift;
    if ( grep { $_ eq $org } ( 'rice', 'arabidopsis', 'grape' ) ) {
        $notre = ',|\\-|\\.';
    }

    my $rs = $s->resultset('FeatureName')->search(
        {
            'feature.dataset_id' => { 'IN'         => $datasets },
            'me.name'            => { 'NOT REGEXP' => $notre },
            'feature.chromosome' =>
              { 'NOT LIKE' => 'contig%' }    # keep only chrs and super_contigs
            ,
            'feature_type.name' => { 'NOT LIKE' => '%contig%' }
        },
        {
            prefetch   => { 'feature' => 'feature_type' },
            order_by   => 'feature_type.name',
            'distinct' => 'feature.feature_id'
        }
    );
    my %seen;
    my @names;
    while ( my $g = $rs->next() ) {
        if ( $seen{ $g->name }++ ) { next; }
        push( @names, [ uc( $g->name ), $g->feature_id, $g->feature->start ] );
    }
    return [ sort { $a->[2] <=> $b->[2] } @names ];
}
