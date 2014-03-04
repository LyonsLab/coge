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
    @EXPORT = qw ( get_wobble_histogram );
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


1;

