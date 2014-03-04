package CoGe::Core::Genome;

use strict;
use warnings;

use File::Spec;

use CoGe::Accessory::Storage qw(get_genome_path);
use CoGe::Accessory::TDS qw(write read);

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw (
        get_wobble_histogram get_wobble_gc_diff_histogram get_feature_type_gc_histogram );
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
    File::Spec->catdir((get_genome_path(shift), "metadata/histograms.json"));
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

1;

