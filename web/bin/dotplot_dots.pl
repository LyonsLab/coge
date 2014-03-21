#!/usr/bin/perl -w

use strict;
use GD;
use Getopt::Long;

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::SynMap_report;
use CoGe::Accessory::Utils qw( commify );
use DBI;
use Data::Dumper;
use DBI;
use POSIX;
use Sort::Versions;
use JSON::XS;
#encode_json(\%hash);

use vars
  qw($P $dagfile $alignfile $genomeid1 $genomeid2 $help $coge $CHR1 $CHR2 $basename $ks_db $assemble $GZIP $GUNZIP $URL $conffile %json_data);

GetOptions(
    "dagfile|d=s"         => \$dagfile,      #all dots
    "alignfile|a=s"       => \$alignfile,    #syntenic dots
    "genomeid1|genome1=i" => \$genomeid1,
    "genomeid2|genome2=i" => \$genomeid2,
    "help|h"              => \$help,
    "chr1|c1=s"           => \$CHR1,
    "chr2|c2=s"           => \$CHR2,
    "basename|b=s"        => \$basename,
    "ksdb|ks_db=s"        => \$ks_db,
#    "assemble=s"       => \$assemble,      #syntenic path assembly option
#    "axis_metrix|am=s" => \$axis_metric, #not needed, but send these data in JSON
#    "box_diags|bd=i"   => \$box_diags, #make sure to send blocks in JSON
    "config_file|cf=s" => \$conffile,
);

$P      = CoGe::Accessory::Web::get_defaults($conffile);
$GZIP   = $P->{GZIP};
$GUNZIP = $P->{GUNZIP};
$URL    = $P->{URL};

usage() if $help;
unless ( ( defined $dagfile && -r $dagfile )
    || -r $alignfile
    || -r "$alignfile.gz" )
{
    print qq{
Need to define an input file for dots.  Either the syntenic pairs file or the alignment file.
};
    usage();
}

if ( defined $dagfile and !( -r $dagfile || -r $dagfile . ".gz" ) ) {
    warn "dagfile specified but not present or readable: $!";
}

$dagfile = CoGe::Accessory::Web::gunzip($dagfile)
  if $dagfile;    # && $dagfile =~ /\.gz$/;
$alignfile = CoGe::Accessory::Web::gunzip($alignfile)
  if $alignfile;    # && $alignfile =~ /\.gz$/;

if ( $alignfile && -r $alignfile && $alignfile =~ /\.gz$/ ) {
    print "Problem decompressing $alignfile.\n";
    exit;
}

if ( $dagfile && -r $dagfile && $dagfile =~ /\.gz$/ ) {
    print "Problem decompressing $dagfile.\n";
    exit;
}

$basename = "test" unless $basename;
$coge = CoGeX->dbconnect($P);

my $synmap_report = new CoGe::Accessory::SynMap_report;

my ($genome1) = $coge->resultset('Genome')->find($genomeid1);
my ($genome2) = $coge->resultset('Genome')->find($genomeid2);
unless ($genome1) {
    warn "No genome found with dbid $genomeid1\n";
    return;
}
unless ($genome2) {
    warn "No genome found with dbid $genomeid2\n";
    return;
}

#get display order of chromosomes, get genome information
my ($org1info) = get_genome_info(
    genome => $genome1,
    chr    => $CHR1,
);

my ($org2info) = get_genome_info(
    genome => $genome2,
    chr    => $CHR2,
);

#not sure that these are needed if gene order is encoded in dag files
#get_gene_info( genomeid => $genomeid1, info => $org1info );
#get_gene_info( genomeid => $genomeid2, info => $org2info );

my $org1length = 0;
map { $org1length += $_->{length} } values %$org1info;
my $org2length = 0;
map { $org2length += $_->{length} } values %$org2info;
unless ( $org1length && $org2length ) {
    print STDERR qq{
Error:  one or both of the genomes has no effective length:
 Org1:  $org1length
 Org2:  $org2length
  };
    exit;
}

#add org_info to json_data now that order has been determined
add_genome_to_json(
    json_data => \%json_data,
    org_data  => $org1info,
    genomeid  => $genomeid1
);
add_genome_to_json(
    json_data => \%json_data,
    org_data  => $org2info,
    genomeid  => $genomeid2
);

my $pairs = get_pairs( file => $alignfile, chr1 => $CHR1, chr2 => $CHR2 )
  if $alignfile && $ks_db && -r $alignfile && -r $ks_db;

#get syntenic gene pairs for ks_data (if needed)
my $ksdata = get_ksdata(
    ks_db => $ks_db,
    pairs => $pairs
) if $ks_db && -r $ks_db;

#get dots for all matches
get_dots(
    file      => $dagfile,
    org1      => $org1info,
    org2      => $org2info,
    genomeid1 => $genomeid1,
    genomeid2 => $genomeid2,
    json_data => \%json_data
) if defined $dagfile && -r $dagfile;

#get_syntenic_dots
my $box_coords = get_dots(
    file           => $alignfile,
    org1           => $org1info,
    org2           => $org2info,
    genomeid1      => $genomeid1,
    genomeid2      => $genomeid2,
    ksdata         => $ksdata,
    json_data      => \%json_data,
    syntenic_pairs => 1,
);

$json_data{layers}{syntenic_blocks}{data}{rects}{$genomeid1}{$genomeid2} = $box_coords;
$json_data{layers}{syntenic_blocks}{style} = {
    strokeStyle => "rgba(0, 155, 0, 0.6)"
};

#write out JSON file of dots"
#print Dumper \%json_data;
open( OUT, ">" . $basename . ".json" ) || die "$!";
print OUT encode_json( \%json_data );
close OUT;

#CoGe::Accessory::Web::gzip($dagfile) if $dagfile && -r $dagfile;
#CoGe::Accessory::Web::gzip($alignfile) if $alignfile && -r $alignfile;
#generate_historgram of ks values if necessary

#This function appears to parse dagchainer output, generated in SynMap.pl, and draw the results to the GD graphics context.
sub get_dots {
    my %opts      = @_;
    my $file      = $opts{file};
    my $org1      = $opts{org1};
    my $org2      = $opts{org2};
    my $genomeid1 = $opts{genomeid1};
    my $genomeid2 = $opts{genomeid2};
    my $ksdata    = $opts{ksdata};
    my $json_data =
      $opts{json_data}; #EL: 9/12/2013: holder for data to be converted to a JSON string to test better dotplot viewer
    my $syntenic_pairs =
      $opts{syntenic_pairs}; #flag for whether these are syntenic gene pairs or general matches
    $syntenic_pairs = 0 unless defined $syntenic_pairs;

    my $data_label = $syntenic_pairs ? "syntenic_pairs" : "pairs";

    my $has_ksdata = keys %$ksdata ? 1 : 0;

    open( IN, $file )
      || die "Can't open $file: $!";    #this is where the problem lies!

    my ($boxes, $block_chr1, $block_chr2);

    #for storing bounds of syntenic blocks
    my ( $block_min_nt_1, $block_min_nt_2, $block_max_nt_1, $block_max_nt_2 );

    #for storing bounds of syntenic blocks
    my ($block_min_gene_1, $block_min_gene_2,  $block_max_gene_1, $block_max_gene_2);

    my $count = 0;
    while (<IN>) {
        chomp;
        if (/^#/) {
            push @{$boxes->{$block_chr1}{$block_chr2}},

       #coordinates is an array of values with: start1, stop1, start2, stop2
       [int($block_min_nt_1), int($block_max_nt_1), int($block_min_nt_2), int($block_max_nt_2)]
       #nucleotides=>[int($block_min_nt_1), int($block_max_nt_1), int($block_min_nt_2), int($block_max_nt_2)],
       #genes=>[int($block_min_gene_1), int($block_max_gene_1), int($block_min_gene_2), int($block_max_gene_2)],
            if defined $block_min_nt_1
                  && defined $block_min_nt_2
                  && defined $block_max_nt_1
                  && defined $block_max_nt_2
                  && defined $block_min_gene_1
                  && defined $block_min_gene_2
                  && defined $block_max_gene_1
                  && defined $block_max_gene_2;
            $block_min_nt_1   = undef;
            $block_min_nt_2   = undef;
            $block_max_nt_1   = undef;
            $block_max_nt_2   = undef;
            $block_min_gene_1 = undef;
            $block_min_gene_2 = undef;
            $block_max_gene_1 = undef;
            $block_max_gene_2 = undef;
            $block_chr1 = undef;
            $block_chr2 = undef;
        }

        next if /^#/;
        next unless $_;

        my @line = split /\t/;
        my $ks_vals;
        my @item1 = split /\|\|/, $line[1];
        my @item2 = split /\|\|/, $line[5];
        my $fid1  = $item1[6];
        my $fid2  = $item2[6];
        if ($has_ksdata) {
            $ks_vals = $ksdata->{$fid1}{$fid2} if $ksdata->{$fid1}{$fid2};
        }
        my ( $a, $chr1 ) = split /_/, $line[0], 2;
        my ( $b, $chr2 ) = split /_/, $line[4], 2;
        my $special = 0
          ; #stupid variable name so that if we are viewing a single chr to single chr comparison within the same organism, this will make collinear matches appear on the inverse section
        if ( $CHR1 && $CHR2 ) {
            if ( $CHR1 eq $chr2 && $CHR2 eq $chr1 ) {
                $special = 1;
                ( $chr1, $chr2 ) = ( $chr2, $chr1 );
            }
            else {
                next if $chr1 ne $CHR1;
                next if $chr2 ne $CHR2;
            }
        }
        #sometimes there will be data that is skipped, e.g. where chromosome="random";
        next unless $org1->{$chr1} && $org2->{$chr2};

        #absolute positions
        my($nt_min_1, $nt_max_1) = sort ($item1[1], $item1[2]);
        my($nt_min_2, $nt_max_2) = sort ($item2[2], $item2[2]);

        #relative position
        my ( $gene_order_1, $gene_order_2 ) = ( $item1[7], $item2[7] );
#        my $data_item = {
#			 chr1   => $chr1,
#			 chr2   => $chr2,
#			 feat1  => int($fid1),
#			 feat2  => int($fid2),
#			 #coordinates is an array of values with: start1, stop1, start2, stop2
#			 nucleotides=>[int($nt_min_1), int($nt_max_1), int($nt_min_2), int($nt_max_2)],
#			 genes=>[int($gene_order_1), int($gene_order_1), int($gene_order_2), int($gene_order_2)],
#			};
#        $data_item->{ks_data} = $ks_vals
#          if $ks_vals;    #check with Evan if missing data should still have key
#        push @{ $json_data->{$genomeid1}{$genomeid2}{$data_label} }, $data_item;

        my $data_item = [int($item1[1]), int($item1[2]), int($item2[1]), int($item2[2])];
#        $data_item->{ks_data} = $ks_vals
#          if $ks_vals;    #check with Evan if missing data should still have key
#          #push @{ $json_data->{data}->{$genomeid1}{$genomeid2}{$data_label}{$chr1}{$chr2}}, $data_item;
        push @{ $json_data->{layers}->{$data_label}{data}{lines}{$genomeid1}{$genomeid2}{$chr1}{$chr2}}, $data_item;

        if ($syntenic_pairs) {
            $json_data->{layers}->{$data_label}{style} = {
                strokeStyle => "rgb(0, 150, 0)"
            };
        }
        #syntenic blocks
        $block_min_nt_1 = $nt_min_1 unless $block_min_nt_1;
        $block_min_nt_1 = $nt_min_1 if $nt_min_1 < $block_min_nt_1;
        $block_min_nt_2 = $nt_min_2 unless $block_min_nt_2;
        $block_min_nt_2 = $nt_min_2 if $nt_min_2 < $block_min_nt_2;
        $block_max_nt_1 = $nt_max_1 unless $block_max_nt_1;
        $block_max_nt_1 = $nt_max_1 if $nt_max_1 > $block_max_nt_1;
        $block_max_nt_2 = $nt_max_2 unless $block_max_nt_2;
        $block_max_nt_2 = $nt_max_2 if $nt_max_2 > $block_max_nt_2;

        $block_min_gene_1 = $gene_order_1 unless $block_min_gene_1;
        $block_min_gene_1 = $gene_order_1 if $gene_order_1 < $block_min_gene_1;
        $block_min_gene_2 = $gene_order_2 unless $block_min_gene_2;
        $block_min_gene_2 = $gene_order_2 if $gene_order_2 < $block_min_gene_2;
        $block_max_gene_1 = $gene_order_1 unless $block_max_gene_1;
        $block_max_gene_1 = $gene_order_1 if $gene_order_1 > $block_max_gene_1;
        $block_max_gene_2 = $gene_order_2 unless $block_max_gene_2;
        $block_max_gene_2 = $gene_order_2 if $gene_order_2 > $block_max_gene_2;

        $block_chr1 = $chr1;
        $block_chr2 = $chr2;
        $count++;
    }

    close IN;
    push @{$boxes->{$block_chr1}{$block_chr2}},
       #coordinates is an array of values with: start1, stop1, start2, stop2
       [int($block_min_nt_1), int($block_max_nt_1), int($block_min_nt_2), int($block_max_nt_2)]
       #nucleotides=>[int($block_min_nt_1), int($block_max_nt_1), int($block_min_nt_2), int($block_max_nt_2)],
       #genes=>[int($block_min_gene_1), int($block_max_gene_1), int($block_min_gene_2), int($block_max_gene_2)],
      if defined $block_min_nt_1
          && defined $block_min_nt_2
          && defined $block_max_nt_1
          && defined $block_max_nt_2
          && defined $block_min_gene_1
          && defined $block_min_gene_2
          && defined $block_max_gene_1
          && defined $block_max_gene_2;

    return $boxes;
}

sub get_ksdata {
    my %opts  = @_;
    my $ks_db = $opts{ks_db};
    my $pairs = $opts{pairs};
    my %data;
    return \%data unless -r $ks_db;
    my $select = "select * from ks_data";
    my $dbh    = DBI->connect( "dbi:SQLite:dbname=$ks_db", "", "" );
    my $sth    = $dbh->prepare($select);
    $sth->execute();

    while ( my $data = $sth->fetchrow_arrayref ) {
        if ($pairs) {
            unless ( $pairs->{ $data->[1] }{ $data->[2] } ) {
                next;
            }
        }
        next unless $data->[3] && $data->[3] =~ /\d/;
        my %item = (
            KS => int($data->[3]),
            KN => int($data->[4]),
#            'KN_KS' => $data->[5], #don't think this is needed to be sent as it can be calculated
        );
        $data{ $data->[1] }{ $data->[2] } = \%item;

    }
    $sth->finish();
    undef $sth;
    $dbh->disconnect();
    return \%data;
}

sub get_pairs {
    my %opts = @_;
    my $file = $opts{file};
    my $chr1 = $opts{chr1};
    my $chr2 = $opts{chr2};
    my %data;
    open( IN, $file ) || die $!;
    while (<IN>) {
        chomp;
        next if /^#/;
        next unless $_;
        my @line  = split /\t/;
        my @item1 = split /\|\|/, $line[1];
        my @item2 = split /\|\|/, $line[5];
        next unless $item1[6] && $item2[6];
        if ($chr1) {
            next unless $item1[0] eq $chr1 || $item2[0] eq $chr1;
        }
        if ($chr2) {
            next unless $item1[0] eq $chr2 || $item2[0] eq $chr2;
        }
        $data{ $item1[6] }{ $item2[6] } = 1;
        $data{ $item2[6] }{ $item1[6] } = 1;
    }
    close IN;
    return \%data;
}

sub get_genome_info {
    my %opts   = @_;
    my $genome = $opts{genome};
    my $chr    = $opts{chr};
    my %data;
    foreach my $gs ( $genome->genomic_sequences ) {
        next if defined $chr && $chr ne $gs->chromosome;
        my $len = $gs->sequence_length;
        if ( $data{ $gs->chromosome } ) {
            warn "Duplicate chromosome:" . $gs->chromosome . "\n";
        }
        $data{ $gs->chromosome }{chr_length} = $len;
        $data{ $gs->chromosome }{length}     = $len;
    }
    return \%data;
}

sub add_genome_to_json {
    my %opts      = @_;
    my $json_data = $opts{json_data};
    my $org_data  = $opts{org_data};
    my $genomeid  = $opts{genomeid};
    my ($genome)  = $coge->resultset('Genome')->find($genomeid);
    my $dbh   = $coge->storage->dbh; #database handle for direct queries

    my %data;
    $data{orgId}   = int( $genome->organism->id );
    $data{orgName} = $genome->organism->name . " (v" . $genome->version . ")";
    $data{name}    = $genome->organism->name . " (v" . $genome->version . ")";
    my $order = 1;

    my $chrs = [];

    foreach my $chr ( keys %$org_data ) {

      #get number of genes in chromosome
      my $query = qq{
SELECT count(distinct(feature_id))
  FROM feature
  JOIN dataset_connector dc using (dataset_id)
 WHERE genome_id = $genomeid
   AND feature_type_id IN (3, 5, 8)
   AND feature.chromosome = '$chr'

};
        my ($gene_count) = $dbh->selectrow_array($query);



        push $chrs, {
				    name   => $chr,
				    nucleotides => int( $org_data->{$chr}{length} ),
#				    genes => int( $gene_count ),
#				    id     => int( $order++ )
				    , #eric changed from 'order' to 'id'.  Make sure that it isn't assumed to be an order.  See if this value can be dropped
		};
      }
    $data{chromosomes} = $chrs;
    $json_data{genomes}{$genomeid} = \%data;
}

#Print out info on script usage
sub usage {
    print qq{
Welcome to $0

This generates JSON formatted data for Syntenic Dotplots

General JSON data structure

data
   |
   -genomes: (array) genomes being compared
          | (hash)
          -genome_id1
                    | (hash)
                    -genome_id2
                              | (hash)
                              - id=genome_id
                              - orgId=organism_id
                              - orgName=organism_name
                              - name=organism_name
                              - chromosomes
                                          | (array)
                                          -name=chromosome_name
                                          -nucleotides=chromosome_length
                                          -genes=gene_length
                                          -id=arbitrary_id
                              - pairs:  Pairs of non-syntenic matches between genomes.  This is optional
                                    | (array)
                                    - chr1   = chromosome of match-mate 1
                                    - feat1  = feature_id of match-mate 1 (may be '0' if a genomic hit)
                                    - chr2   = chromosome of match-mate 2
                                    - feat2  = feature_id of match-mate 2 (may be '0' if a genomic hit)

                                    #location coordinates are an array of values: start1, stop1, start2, stop2
                                    - nucleotides = [start1, stop1, start2, stop2]
                                    - genes       = [start1, stop1, start2, stop2]

                              - syntenic_pairs:  Pairs of syntenic matches between genomes.
                                    | (array)
                                    - chr1   = chromosome of match-mate 1
                                    - feat1  = feature_id of match-mate 1 (may be '0' if a genomic hit)
                                    - chr2   = chromosome of match-mate 2
                                    - feat2  = feature_id of match-mate 2 (may be '0' if a genomic hit)

                                    #location coordinates are an array of values: start1, stop1, start2, stop2
                                    - nucleotides = [start1, stop1, start2, stop2]
                                    - genes       = [start1, stop1, start2, stop2]

                                    - ks_data: optional.
                                            | (hash)
                                            - KS  =  Synonymous mutation rate
                                            - KN  =  Nonsynonymous mutation rate
                                    - syntenic_blocks:
                                                    | (array)
                                                    #location coordinates are an array of values: start1, stop1, start2, stop2
                                                    - nucleotides = [start1, stop1, start2, stop2]
                                                    - genes       = [start1, stop1, start2, stop2]
           -<additional genome pairs>




dagfile      | d       path to dag file containing all the hits

alignfile    | a       path to .aligncoords file generated by
                       dagchainer containing just the diags
genomeid1       | genome1 | gid1    database id of genome on x-axis

genomeid2       | genome2 | gid1    database id of genome on y-axis

confile      | cf      CoGe configuration file for getting various parameters for connecting to database, etc.

chr1         | c1      only show data for this chromosome on the x axis

chr2         | c2      only show data for this chromosome on the y axis

basename     | b       base path and name for output files (test)

ks_db        | ksdb    specify a sqlite database with synonymous/nonsynonymous data
                       to color syntenic points

help         | h       print this message

};
    exit;
}

