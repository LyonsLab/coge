#!/usr/bin/perl -w
use v5.14;
use strict;

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Chromosomes;

use DBI;
use Data::Dumper;
use Getopt::Long;
use JSON::XS qw(encode_json);
use List::Util qw(reduce max min);
use POSIX;
use Sort::Versions;

our ($P, $dagfile, $alignfile, $genomeid1, $genomeid2, $help, $coge, $CHR1,
     $CHR2, $basename, $ks_db, $assemble, $GZIP, $GUNZIP, $URL, $config,
     %all_pairs, %base_data, %ks_data, $DEBUG);

GetOptions(
    "dagfile|d=s"      => \$dagfile,      #all dots
    "alignfile|a=s"    => \$alignfile,    #syntenic dots
    "genome1|gid1=i"   => \$genomeid1,
    "genome2|gid2=i"   => \$genomeid2,
    "help|h"           => \$help,
    "chr1|c1=s"        => \$CHR1,
    "chr2|c2=s"        => \$CHR2,
    "basename|b=s"     => \$basename,
    "ksdb|ks_db=s"     => \$ks_db,
    "debug|dbg=i"      => \$DEBUG,
    "config_file|cf=s" => \$config,
);

$GZIP   = $P->{GZIP};
$GUNZIP = $P->{GUNZIP};
$URL    = $P->{URL};

if($config) {
    $P          = CoGe::Accessory::Web::get_defaults($config);
    my $db      = $P->{DBNAME};
    my $host    = $P->{DBHOST};
    my $port    = $P->{DBPORT};
    my $user    = $P->{DBUSER};
    my $pass    = $P->{DBPASS};
    my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";

    $coge = CoGeX->connect( $connstr, $user, $pass );
}

unless ($coge) {
    die "error: unable to connect to the database";
}

usage() if $help;

my $has_alignfile = (defined $alignfile && -r $alignfile);
my $has_dagfile = (defined $dagfile && -r $dagfile );
my $has_ksfile = (defined $ks_db && -r $ks_db);

unless ($has_alignfile || $has_dagfile) {
    say STDERR "Need to define an input file for dots. Either the syntenic pairs file or the alignment file.";
    usage();
}

$basename = "test" unless $basename;

my ($genome1) = $coge->resultset('Genome')->find($genomeid1);
my ($genome2) = $coge->resultset('Genome')->find($genomeid2);

unless ($genome1) {
    die "error: No genome found with dbid $genomeid1\n";
}

unless ($genome2) {
    die "error: No genome found with dbid $genomeid2\n";
}

# get display order of chromosomes, get genome information
my ($org1info) = get_genome_info(genome => $genome1, chr => $CHR1);
my ($org2info) = get_genome_info(genome => $genome2, chr => $CHR2);

#not sure that these are needed if gene order is encoded in dag files
#get_gene_info( genomeid => $genomeid1, info => $org1info );
#get_gene_info( genomeid => $genomeid2, info => $org2info );

my $org1length = reduce { $a + $b } 0, (values $org1info);
my $org2length = reduce { $a + $b } 0, (values $org2info);

die "Organism 1 has an effective length of zero" unless $org1length;
die "Organism 2 has an effective length of zero" unless $org2length;

#add org_info to json_data now that order has been determined
add_genome_to_json(
    json_data => \%base_data,
    org_data  => $org1info,
    genomeid  => $genomeid1
);

add_genome_to_json(
    json_data => \%base_data,
    org_data  => $org2info,
    genomeid  => $genomeid2
);

my $pairs;

if ($has_ksfile && $has_alignfile) {
    $pairs = get_pairs( file => $alignfile, chr1 => $CHR1, chr2 => $CHR2 );
}

# get syntenic gene pairs for ks_data (if needed)
my $ksdata = get_ksdata(
    ks_db => $ks_db,
    pairs => $pairs
) if $has_ksfile;

my %ks_json;
$ks_json{datasets}{histogram}{kn}{layer} = "syntenic_pairs";
$ks_json{datasets}{histogram}{kn}{title} = "Ks Values";

$ks_json{datasets}{histogram}{ks}{layer} = "syntenic_pairs";
$ks_json{datasets}{histogram}{ks}{title} = "Kn Values";

$ks_json{datasets}{features}{layer} = "syntenic_pairs";
$ks_json{datasets}{features}{title} = "Features";

# get dots for all matches
get_dots(
    file      => $dagfile,
    org1      => $org1info,
    org2      => $org2info,
    genomeid1 => $genomeid1,
    genomeid2 => $genomeid2,
    ksdata    => $ksdata,
    ks_json   => \%ks_json,
    json_data => \%all_pairs
) if $has_dagfile;

# get_syntenic_dots
my $box_coords = get_dots(
    file           => $alignfile,
    org1           => $org1info,
    org2           => $org2info,
    genomeid1      => $genomeid1,
    genomeid2      => $genomeid2,
    json_data      => \%base_data,
    syntenic_pairs => 1,
);

$base_data{layers}{syntenic_blocks}{data}{rects}{$genomeid1}{$genomeid2} = $box_coords;
$base_data{layers}{syntenic_blocks}{style} = {
    strokeStyle => "rgba(25, 25, 25, 0.6)"
};

#write out JSON file of dots"
if (%base_data) {
    print Dumper \%base_data if $DEBUG;
    open( OUT, ">" . $basename . ".json" ) || die "$!";
    print OUT encode_json( \%base_data );
    close OUT;
}

if (%all_pairs) {
    print Dumper \%all_pairs if $DEBUG;
    open( OUT, ">" . $basename . ".all.json" ) || die "$!";
    print OUT encode_json( \%all_pairs);
    close OUT;
}

#generate_historgram of ks values if necessary
if (%ks_json) {
    print Dumper \%all_pairs if $DEBUG;
    open( OUT, ">" . $basename . ".datasets.json" ) || die "$!";
    print OUT encode_json( \%ks_json );
    close OUT;
}

sub parse_line {
    my @field = split(/\t/, shift);
    my @item1 = split(/\|\|/, $field[1]);
    my @item2 = split(/\|\|/, $field[5]);

    my ( undef, $chr1 ) = split(/_/, $field[0], 2);
    my ( undef, $chr2 ) = split(/_/, $field[4], 2);

    # absolute positions
    my ($s1, $e1, $s2, $e2) = (int($field[2]), int($field[3]), int($field[6]), int($field[7]));

    # relative position
    my ( $gene_order_1, $gene_order_2 ) = ( $item1[7], $item2[7] );

    return {
        chr1 => $chr1,
        chr2 => $chr2,
        fid1 => $item1[6],
        fid2 => $item2[6],
        s1   => $s1,
        e1   => $e1,
        s2   => $s2,
        e2   => $e2,
        g1   => $gene_order_1,
        g2   => $gene_order_2,
    };
}

# This function appears to parse dagchainer output, generated in SynMap.pl
sub get_dots {
    my %opts      = @_;
    my $file      = $opts{file};
    my $org1      = $opts{org1};
    my $org2      = $opts{org2};
    my $genomeid1 = $opts{genomeid1};
    my $genomeid2 = $opts{genomeid2};
    my $ksdata    = $opts{ksdata};
    my $ks_json   = $opts{ks_json};

    # EL: 9/12/2013: holder for data to be converted to a JSON string to test
    # better dotplot viewer
    my $json_data = $opts{json_data};

    #flag for whether these are syntenic gene pairs or general matches
    my $syntenic_pairs = $opts{syntenic_pairs} // 0; #/

    my $data_label = $syntenic_pairs ? "syntenic_pairs" : "pairs";
    my $has_ksdata = keys %$ksdata ? 1 : 0;

    #this is where the problem lies!
    open( IN, $file )
      || die "Can't open $file: $!";

    #for storing bounds of syntenic blocks
    my ($block_min_nt_1,   $block_min_nt_2,     $block_max_nt_1,
        $block_min_gene_1, $block_min_gene_2,   $block_max_gene_1,
        $block_max_gene_2, $boxes, $block_chr1, $block_max_nt_2,
        $block_chr2);

    my ($count, $block_count) = (0, 0);

    while (<IN>) {
        chomp;
        next if /^#/;
        next unless $_;

        my $fields = parse_line($_);

        my $fid1 = $fields->{fid1};
        my $fid2 = $fields->{fid2};
        my $ks_vals;

        if ($has_ksdata) {
            $ks_vals = $ksdata->{$fid1}{$fid2} if $ksdata->{$fid1}{$fid2};
        }

        my $chr1 = $fields->{chr1};
        my $chr2 = $fields->{chr2};

        # sometimes there will be data that is skipped, e.g. where
        # chromosome="random";
        next unless (defined $org1->{$chr1} && defined $org2->{$chr2});

        # Find max and min distance for entry 1
        my ($s1, $e1, $g1) = ($fields->{s1}, $fields->{e1}, $fields->{g1});
        my ($nt_min_1, $nt_max_1) = sort {$a <=> $b} ($s1, $e1);

        my ($s2, $e2, $g2) = ($fields->{s2}, $fields->{e2}, $fields->{g2});
        my ($nt_min_2, $nt_max_2) = sort {$a <=> $b} ($s2, $e2);

        # Complement coordinates if chromosome is complemented
        if ($org1->{$chr1}{rev}) {
            $s1 = int($org1->{$chr1}{length}) - $s1;
            $e1 = int($org1->{$chr1}{length}) - $e1;
        }
        if ($org2->{$chr2}{rev}) {
            $s2 = int($org2->{$chr2}{length}) - $s2;
            $e2 = int($org2->{$chr2}{length}) - $e2;
        }

#        my $data_item = {
#            chr1   => $chr1,
#            chr2   => $chr2,
#            feat1  => int($fid1),
#            feat2  => int($fid2),
#            #coordinates is an array of values with: start1, stop1, start2, stop2
#            nucleotides=>[int($nt_min_1), int($nt_max_1), int($nt_min_2), int($nt_max_2)],
#            genes=>[int($gene_order_1), int($gene_order_1), int($gene_order_2), int($gene_order_2)],
#           };
#        $data_item->{ks_data} = $ks_vals
#          if $ks_vals;    #check with Evan if missing data should still have key
#        push @{ $json_data->{$genomeid1}{$genomeid2}{$data_label} }, $data_item;

        if ($ks_vals) {
            $ks_json{datasets}{features}{data}{$count} = $ks_vals->{features};
            $ks_json{datasets}{histogram}{kn}{data}{$count} = $ks_vals->{kn};
            $ks_json{datasets}{histogram}{ks}{data}{$count} = $ks_vals->{ks};
        }

        # Save coordinates as data point
        my $data_item = [$s1, $e1, $s2, $e2];

#        $data_item->{ks_data} = $ks_vals
#          if $ks_vals;    #check with Evan if missing data should still have key
#          #push @{ $json_data->{data}->{$genomeid1}{$genomeid2}{$data_label}{$chr1}{$chr2}}, $data_item;
        $json_data->{layers}->{$data_label}{data}{lines}{$genomeid1}{$genomeid2}{$chr1}{$chr2}{$count} = $data_item;

        if ($syntenic_pairs) {
            $json_data->{layers}->{$data_label}{style} = {
                strokeStyle => "rgb(0, 150, 0)"
            };
        }

        #syntenic blocks in nucleotides
        $block_min_nt_1 = min grep { defined $_ } ($nt_min_1, $block_min_nt_1);
        $block_min_nt_2 = min grep { defined $_ } ($nt_min_2, $block_min_nt_2);
        $block_max_nt_1 = max grep { defined $_ } ($nt_max_1, $block_max_nt_1);
        $block_max_nt_2 = max grep { defined $_ } ($nt_max_2, $block_max_nt_2);

        #syntenic blocks in genes
        $block_min_gene_1 = min grep { defined $_ } ($g1, $block_min_gene_1);
        $block_min_gene_2 = min grep { defined $_ } ($g2, $block_min_gene_2);
        $block_max_gene_1 = max grep { defined $_ } ($g1, $block_max_gene_1);
        $block_max_gene_2 = max grep { defined $_ } ($g2, $block_max_gene_2);

        $block_chr1 = $chr1;
        $block_chr2 = $chr2;

        $count++;

    } continue {
        if (/^#/) {
            if (defined $block_min_nt_1
                && defined $block_min_nt_2
                && defined $block_max_nt_1
                && defined $block_max_nt_2
                && defined $block_min_gene_1
                && defined $block_min_gene_2
                && defined $block_max_gene_1
                && defined $block_max_gene_2) {

                #coordinates is an array of values with: start1, stop1, start2, stop2
                $boxes->{$block_chr1}{$block_chr2}{$block_count++} = [
                    int($block_min_nt_1),
                    int($block_max_nt_1),
                    int($block_min_nt_2),
                    int($block_max_nt_2),
                    #int($block_min_gene_1),
                    #int($block_max_gene_1),
                    #int($block_min_gene_2),
                    #int($block_max_gene_2),
                ];
            };

            # Reset chromosome-chromosome block
            $block_min_nt_1   = undef;
            $block_min_nt_2   = undef;
            $block_max_nt_1   = undef;
            $block_max_nt_2   = undef;
            $block_min_gene_1 = undef;
            $block_min_gene_2 = undef;
            $block_max_gene_1 = undef;
            $block_max_gene_2 = undef;
            $block_chr1       = undef;
            $block_chr2       = undef;
        }
    }

    close IN;

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
        my $fid1 = int($data->[1]);
        my $fid2 = int($data->[2]);

        $data{$fid1}{$fid2} = {
            ks       =>  $data->[3] + 0.0,
            kn       =>  $data->[4] + 0.0,
            features => { fid1 =>  $fid1, fid2 =>  $fid2, },
        };
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
        my @line  = split(/\t/);
        my @item1 = split(/\|\|/, $line[1]);
        my @item2 = split(/\|\|/, $line[5]);
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
#    foreach my $gs ( $genome->genomic_sequences ) {
#        next if defined $chr && $chr ne $gs->chromosome;
#        my $len = $gs->sequence_length;
#        if ( $data{ $gs->chromosome } ) {
#            warn "Duplicate chromosome:" . $gs->chromosome . "\n";
#        }
#        $data{ $gs->chromosome }{chr_length} = $len;
#        $data{ $gs->chromosome }{length}     = $len;
#    }
	my $c = CoGe::Core::Chromosomes->new($genome->id);
	while ($c->next) {
        next if defined $chr && $chr ne $c->name;
        my $len = $c->length;
        if ( $data{ $c->name } ) {
            warn "Duplicate chromosome:" . $c->name . "\n";
        }
        $data{ $c->name }{chr_length} = $len;
        $data{ $c->name }{length}     = $len;
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
#                   genes => int( $gene_count ),
#                   id     => int( $order++ )
                    , #eric changed from 'order' to 'id'.  Make sure that it isn't assumed to be an order.  See if this value can be dropped
        };
      }
    $data{chromosomes} = $chrs;
    $json_data->{genomes}{$genomeid} = \%data;
}

#Print out info on script usage
sub usage {
    print qq{
Welcome to $0

This generates JSON formatted data for Syntenic Dotplots

dagfile      | d       path to dag file containing all the hits

alignfile    | a       path to .aligncoords file generated by
                       dagchainer containing just the diags
genomeid1    | genome1 | gid1    database id of genome on x-axis

genomeid2    | genome2 | gid1    database id of genome on y-axis

confile      | cf      CoGe configuration file for getting various parameters for connecting to database, etc.

chr1         | c1      only show data for this chromosome on the x axis

chr2         | c2      only show data for this chromosome on the y axis

basename     | b       base path and name for output files (test)

ks_db        | ksdb    specify a sqlite database with synonymous/nonsynonymous data
                       to color syntenic points

help         | h       print this message

};
    exit 1;
}
