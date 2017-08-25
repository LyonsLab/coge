package CoGe::Core::Genome;

use strict;
use warnings;

use File::Spec::Functions;
use Sort::Versions;
use Data::Dumper;

use CoGe::Accessory::TDS qw(write read);
use CoGe::Accessory::Utils;
use CoGe::Core::Storage qw(get_genome_path);
use CoGe::Core::Favorites;
use CoGe::Accessory::IRODS qw($IRODS_METADATA_PREFIX);

BEGIN {
    our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw( has_statistic get_stats calc_gc calc_noncoding_gc
        get_wobble_histogram get_wobble_gc_diff_histogram get_feature_type_gc_histogram
        fix_chromosome_id read_fasta_index get_irods_metadata);
    @EXPORT_OK = qw(genomecmp genomecmp2 search_genomes);
}

my @LOCATIONS_PREFETCH = (
    { "feature_type_id" => 3 },
    {
        join => [
            'locations',
            { 'dataset' => { 'dataset_connectors' => 'genome' } }
        ],
        prefetch => [
            'locations',
            { 'dataset' => { 'dataset_connectors' => 'genome' } }
        ]
    }
);

sub genomecmp($$) {
    my ($a, $b) = @_;
    return genomecmp2($a, $b);
}

sub genomecmp2 {
    my ($a, $b, $favorites) = @_;

    my $namea = $a->name ? $a->name :  '';
    my $nameb = $b->name ? $b->name :  '';
    my $typea = $a->type ? $a->type->id : 0;
    my $typeb = $b->type ? $b->type->id : 0;
    my $favea = $favorites ? $favorites->is_favorite($a) : 0;
    my $faveb = $favorites ? $favorites->is_favorite($b) : 0;

    $b->certified <=> $a->certified
      || $faveb <=> $favea
      || $a->organism->name cmp $b->organism->name
      || versioncmp( $b->version, $a->version )
      || $typea <=> $typeb
      || $namea cmp $nameb
      || $a->id cmp $b->id;
}

sub search_genomes {
    my %opts = @_;
    my $db = $opts{db};
    my $search_term = $opts{search_term}; # id value, or keyword in name/description
    my $user = $opts{user}; # optional database user object
    my $sort = $opts{sort}; # optional boolean flag to sort results (not sorted by default for performance)
    return unless $db and $search_term;
    #my $include_deleted = $opts{include_deleted}; # optional boolean flag
    
    # Search genomes
    my $search_term2 = '%' . $search_term . '%';
    my @genomes = $db->resultset("Genome")->search(
        \[
            'genome_id = ? OR name LIKE ? OR description LIKE ?',
            [ 'genome_id', $search_term  ],
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );
    
    # Search organisms
    my @organisms = $db->resultset("Organism")->search(
        \[
            'name LIKE ? OR description LIKE ?',
            [ 'name',        $search_term2 ],
            [ 'description', $search_term2 ]
        ]
    );
    
    # Combine matching genomes and organisms, preventing duplicates
    my %unique;
    map { $unique{ $_->id } = $_ } @genomes;
    foreach my $organism (@organisms) {
        map { $unique{ $_->id } = $_ } $organism->genomes;
    }

    # Filter response
    my @filtered = sort genomecmp grep {
        !$_->restricted || (defined $user && $user->has_access_to_genome($_))
    } values %unique;

    # Sort
    if ($sort) {
        my $favorites;
        $favorites = CoGe::Core::Favorites->new(user => $user) if $user;
        @filtered = sort { genomecmp2($a, $b, $favorites) } @filtered;
    }
    
    return \@filtered;
}

sub get_wobble_histogram {
    my $genome = _get_genome_or_exit(@_);
    my $storage_path = _get_histogram_file($genome->id);

    my $data = read($storage_path);
    return $data->{wobble_histogram} if defined $data->{wobble_histogram};

    $data->{wobble_histogram} = _generate_wobble_content($genome);

    # Exit if generate failed
    unless(defined $data->{wobble_histogram}) {
        say STDERR "Genome::get_genome_wobble_content: generate wobble content failed!";
        exit;
    }

    say STDERR "Genome::get_genome_wobble_content: write failed!"
        unless write($storage_path, $data);

    # Return data
    return $data->{wobble_histogram};
}

sub get_feature_type_gc_histogram {
    my $genome = _get_genome_or_exit(shift);
    my $typeid = shift;

    unless (defined $typeid) {
        say STDERR "Genome::get_feature_type_gc_histogram: typeid is not defined!";
        exit;
    }

    my $key = 'feature_type_' . $typeid . '_gc_histogram';

    my $storage_path = _get_histogram_file($genome->id);
    my $data = read($storage_path);
    return $data->{$key} if defined $data->{$key};

    $data->{$key} = _generate_feature_type_gc($genome, $typeid);

    # Exit if generate failed
    unless(defined $data->{$key}) {
        say STDERR "Genome::get_genome_wobble_content: generate wobble content failed!";
        exit;
    }

    say STDERR "Genome::get_genome_wobble_content: write failed!"
        unless write($storage_path, $data);

    # Return data
    return $data->{$key};
}

sub get_wobble_gc_diff_histogram {
    my $genome = _get_genome_or_exit(@_);
    my $storage_path = _get_histogram_file($genome->id);

    my $data = read($storage_path);
    return $data->{wobble_gc_diff_histogram} if defined $data->{wobble_gc_diff_histogram};

    $data->{wobble_gc_diff_histogram} = _generate_wobble_gc_diff($genome);

    # Exit if generate failed
    unless(defined $data->{wobble_gc_diff_histogram}) {
        say STDERR "Genome::get_genome_wobble_content: generate wobble content failed!";
        exit;
    }

    say STDERR "Genome::get_genome_wobble_content: write failed!"
        unless write($storage_path, $data);

    # Return data
    return $data->{wobble_gc_diff_histogram};
}

sub has_statistic {
    my $genome = _get_genome_or_exit(shift);
    my $stat = shift;

    my $storage_path = _get_stats_file($genome->id);
    my $data = read($storage_path);

    return defined $data->{$stat};
}

sub _generate_wobble_content {
    my $genome = shift;
    my $gstid = $genome->type->id;
    my $wobble_content = {};

    my ($at, $gc, $n) = (0) x 3;

    foreach my $ds ($genome->datasets()) {
        foreach my $feat ($ds->features(@LOCATIONS_PREFETCH)) {
            my @gc = $feat->wobble_content( counts => 1 );
            $gc = $gc[0] if $gc[0] && $gc[0] =~ /^\d+$/;
            $at = $gc[1] if $gc[1] && $gc[1] =~ /^\d+$/;
            $n  = $gc[2] if $gc[2] && $gc[2] =~ /^\d+$/;

            my $total = 0;
            $total += $gc[0] if $gc[0];
            $total += $gc[1] if $gc[1];
            $total += $gc[2] if $gc[2];
            my $perc_gc = 100 * $gc[0] / $total if $total;

            $wobble_content->{$feat->id . '_' . $gstid} = {
                at => $at,
                gc => $gc,
                n => $n,
            };

            #skip if no values
            next unless $perc_gc;

            my $node = $wobble_content->{$feat->id . '_' . $gstid};
            $node->{percent_gc} = $perc_gc;
        }
    }

    return $wobble_content;
}

sub get_stats {
    my ($genome, $key, $calc_func) = @_;
    my $storage_path = _get_stats_file($genome->id);
    my $data = read($storage_path);
    return $data->{$key} if defined $data->{$key};

    $data->{$key} = $calc_func->($genome);
    write($storage_path, $data);
    return $data->{$key};
}

sub calc_gc {
    my $genome = shift;
    my $gstid = $genome->type->id;

    my ( $gc, $at, $n, $x ) = (0) x 4;

    foreach my $ds ($genome->datasets) {
        foreach my $chr ( $ds->chromosomes ) {
            my @gc = $ds->percent_gc( chr => $chr, seq_type => $gstid, count => 1 );
            $gc += $gc[0];
            $at += $gc[1];
            $n  += $gc[2];
            $x  += $gc[3];
        }
    }
    my $total = $gc + $at + $n + $x;
    return unless $total;

    return {
        total => $total,
        gc    => $gc / $total,
        at    => $at / $total,
        n     => $n  / $total,
        x     => $x  / $total,
    };
}

sub calc_noncoding_gc {
    my $genome = shift;
    my $gstid = $genome->type->id;

    my %seqs; # prefetch the sequences (slow for many seqs)
    map {
        $seqs{$_} = $genome->get_genomic_sequence( chr => $_, seq_type => $gstid )
    } $genome->chromosomes;

    my @datasets;
    push @datasets, $genome->datasets;
    foreach my $ds (@datasets) {
        foreach my $feat ($ds->features(@LOCATIONS_PREFETCH)) {
            foreach my $loc ( $feat->locs ) {
                if ( $loc->stop > length( $seqs{ $feat->chromosome } ) ) {
                    print STDERR "feature "
                      . $feat->id
                      . " stop exceeds sequence length: "
                      . $loc->stop . " :: "
                      . length( $seqs{ $feat->chromosome } ), "\n";
                }
                substr(
                    $seqs{ $feat->chromosome },
                    $loc->start - 1,
                    ( $loc->stop - $loc->start + 1 )
                ) = "-" x ( $loc->stop - $loc->start + 1 );
            }
        }
    }

    my ( $gc, $at, $n, $x ) = ( 0, 0, 0, 0 );

    foreach my $seq ( values %seqs ) {
        $gc += $seq =~ tr/GCgc//;
        $at += $seq =~ tr/ATat//;
        $n  += $seq =~ tr/nN//;
        $x  += $seq =~ tr/xX//;
    }

    my $total = $gc + $at + $n + $x;
    return unless $total;

    return {
        total => $total,
        gc    => $gc / $total,
        at    => $at / $total,
        n     => $n  / $total,
        x     => $x  / $total,
    };
}

#
# Private functions
#
sub _get_genome_or_exit {
    my $genome = shift;

    unless ($genome) {
        say STDERR "Genome::get_genome_wobble_content: genome not specified!";
        exit;
    }

    return $genome;
}

sub _get_histogram_file {
    catfile((get_genome_path(shift), "metadata/histograms.json"));
}

sub _get_stats_file {
    catfile((get_genome_path(shift), "metadata/stats.json"));
}

sub _generate_wobble_gc_diff {
    my $genome = shift;
    my $gstid = $genome->type->id;
    my $data = [];

    foreach my $ds ($genome->datasets) {
        foreach my $feat ($ds->features(@LOCATIONS_PREFETCH)) {
            my @wgc  = $feat->wobble_content();
            my @gc   = $feat->gc_content();
            my $diff = $gc[0] - $wgc[0] if defined $gc[0] && defined $wgc[0];
            push @$data, sprintf( "%.2f", 100 * $diff ) if $diff;
        }
    }

    return $data;
}

sub _generate_feature_type_gc {
    my ($genome, $typeid) = @_;
    my $gstid = $genome->type->id;
    my $gc_content = {};

    my (@items, @datasets);

    push @items, $genome;

    my %seqs; # prefetch the sequences with one call to genomic_sequence (slow for many seqs)
    foreach my $item (@items) {
        map {
            $seqs{$_} = $item->get_genomic_sequence( chr => $_, seq_type => $gstid )
        } $item->chromosomes;
    }

    my ($at, $gc, $n) = (0) x 3;

    my @params = (
        { "feature_type_id" => $typeid },
        {
            join => [
                'locations',
                { 'dataset' => { 'dataset_connectors' => 'genome' } }
            ],
            prefetch => [
                'locations',
                { 'dataset' => { 'dataset_connectors' => 'genome' } }
            ],
        }
    );

    foreach my $ds ($genome->datasets) {
        my @feats = $ds->features(@params);

        foreach my $feat (@feats) {
            my $seq = substr(
                $seqs{ $feat->chromosome },
                $feat->start - 1,
                $feat->stop - $feat->start + 1
            );

            $feat->genomic_sequence( seq => $seq );
            my @gc = $feat->gc_content( counts => 1 );

            $gc = $gc[0] if $gc[0] =~ /^\d+$/;
            $at = $gc[1] if $gc[1] =~ /^\d+$/;
            $n  = $gc[2] if $gc[2] =~ /^\d+$/;

            my $total = 0;
            $total += $gc[0] if $gc[0];
            $total += $gc[1] if $gc[1];
            $total += $gc[2] if $gc[2];

            my $perc_gc = 100 * $gc[0] / $total if $total;

            $gc_content->{$feat->id . '_' . $gstid} = {
                at => $at,
                gc => $gc,
                n => $n
            };

            #skip if no values
            next unless $perc_gc;
            my $node = $gc_content->{$feat->id . '_' . $gstid};
            $node->{percent_gc} = sprintf( "%.2f", $perc_gc );
        }
    }

    return $gc_content;
}

# Used by the load scripts:  load_genome.pl, load_experiment.pl, and load_annotation.pl
# mdb consolidated from load scripts into this module, 2/11/15 COGE-587
sub fix_chromosome_id { 
    my $chr = shift;        # chr id to fix
    my $genome_chr = shift; # optional hash ref of existing chromosome ids for genome
    return unless defined $chr;

    # Fix chromosome identifier
    $chr =~ s/^lcl\|//;           # remove leading 'lcl|'
    $chr =~ s/^gi\|//;            # remove leading 'gi|'
    $chr =~ s/chromosome//i;      # remove 'chromosome'
    $chr =~ s/^chr//i;            # remove leading 'chr'
    $chr = "0" if $chr =~ /^0+$/; # handle chromosome name '00' (or something like that) (EL added 2/13/14)
    $chr =~ s/^0+// unless $chr eq '0'; # remove leading 0s
    $chr =~ s/^_+//;              # remove underscores
    $chr =~ s/\s+/ /g;            # collapse whitespace to single space
    $chr =~ s/^\s//;              # remove leading whitespace
    $chr =~ s/\s$//;              # remove trailing whitespace
    $chr =~ s/\//_/g;             # replace '/' with '_' (mdb added 12/17/13 issue 266)
    $chr =~ s/\|$//;              # remove trailing pipes (mdb added 3/14/14 issue 332)
    $chr =~ s/\|/_/g;             # convert pipes to underscore (mdb added 8/13/15)
    $chr =~ s/\(/_/g;             # replace '(' with '_' (mdb added 2/11/15 COGE-587)
    $chr =~ s/\)/_/g;             # replace ')' with '_' (mdb added 2/11/15 COGE-587)
    $chr =~ s/_+/_/g;             # convert multiple underscores to single underscore (mdb added 8/13/15)
    return if ($chr eq '');

    # Convert 'chloroplast' and 'mitochondia' to 'C' and 'M' if needed
    if (defined $genome_chr) {
        if (   $chr =~ /^chloroplast$/i
            && !$genome_chr->{$chr}
            && $genome_chr->{"C"} )
        {
            $chr = "C";
        }
        if (   $chr =~ /^mitochondria$/i
            && !$genome_chr->{$chr}
            && $genome_chr->{"M"} )
        {
            $chr = "M";
        }
    }

    return $chr;
}

sub read_fasta_index {
    my $index_file = shift; # samtools .fai file
    
    my %contigs;
    open(my $fh, $index_file) or return;
    while (<$fh>) {
        my ($name, $size) = split("\t");
        unless ($name && $size) {
            print STDERR "read_fasta_index: error reading index file\n";
            return;    
        }
        
        $name = fix_chromosome_id($name); # mdb added 8/18/16 COGE-735
        $contigs{$name} = $size;
    }
    close($fh);
    
    return \%contigs;
}

sub get_irods_metadata {
    my $genome = shift;
    
    my %md = (
        $IRODS_METADATA_PREFIX.'link'              => "http://genomevolution.org",
        $IRODS_METADATA_PREFIX.'OrganismView-link' => "http://genomevolution.org/coge/OrganismView.pl?gid=".$genome->id,
        $IRODS_METADATA_PREFIX.'GenomeInfo-link'   => "http://genomevolution.org/coge/GenomeInfo.pl?gid=".$genome->id,
        $IRODS_METADATA_PREFIX.'genome-id'                => $genome->id,
        $IRODS_METADATA_PREFIX.'genome-organism-name'     => $genome->organism->name,
        $IRODS_METADATA_PREFIX.'genome-organism-taxonomy' => $genome->organism->description,
        $IRODS_METADATA_PREFIX.'genome-version'           => $genome->version,
        $IRODS_METADATA_PREFIX.'genome-sequence-type'     => $genome->type->info,
        $IRODS_METADATA_PREFIX.'genome-summary'           => $genome->info(hideRestrictedSymbol => 1),
        $IRODS_METADATA_PREFIX.'genome-certified'         => $genome->certified ? 'true' : 'false'
    );

    $md{$IRODS_METADATA_PREFIX.'genome-link'}            = $genome->link if ($genome->link);
    $md{$IRODS_METADATA_PREFIX.'genome-additional-info'} = $genome->message if ($genome->message);
    $md{$IRODS_METADATA_PREFIX.'genome-name'}            = $genome->name if ($genome->name);
    $md{$IRODS_METADATA_PREFIX.'genome-description'}     = $genome->description if ($genome->description);

    # Add sources
    my $i = 1;
    my @sources = $genome->source;
    foreach my $item (@sources) {
        my $source = $item->name;
        $source.= ": ".$item->description if $item->description;
        my $key = "genome-source";
        $key .= $i if scalar @sources > 1;
        $md{$IRODS_METADATA_PREFIX.$key} = $source;
        $md{$IRODS_METADATA_PREFIX.$key."-link"} = $item->link if $item->link;
        $i++;
    }

    # Add custom metadata
    foreach my $a ( $genome->annotations ) {
        my $group = (
                defined $a->type->group
            ? $a->type->group->name.','.$a->type->name
            : $a->type->name
        );

        $md{$group} = $a->info(hideRestrictedSymbol => 1);
    }

    return \%md;
}

1;
