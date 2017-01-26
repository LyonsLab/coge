package CoGe::Builder::Methylation::Metaplot;

use Moose;
extends 'CoGe::Builder::Methylation::Analyzer';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);
use String::ShellQuote;

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my %opts = @_;

    my ($bam_file, $genome, $experiment);
    if ($opts{data_files}) {
        $genome = $self->request->genome;
        ($bam_file) = @{$opts{data_files}};
    }
    else { # for when called from ExperimentView
        $genome     = $self->request->genome;
        $experiment = $self->request->experiment;
        $bam_file   = get_experiment_files($experiment->id, $experiment->data_type)->[0];
    }

    # Make sure genome is annotated as is required by metaplot script
    my $isAnnotated = $genome->has_gene_features;
    unless ($isAnnotated) {
        CoGe::Exception::Generic->throw(message => 'Genome must be annotated to generate metaplot');
    }

    #
    # Build the workflow
    #

    # Generate cached gff
    $self->add(
        $self->create_gff( #FIXME simplify this
            gid => $genome->id,
            output_file => get_gff_cache_path(
                gid => $genome->id,
                genome_name => sanitize_name($genome->organism->name),
                output_type => 'gff',
                params => {}
            )
        )
    );
    my $gff_file = $self->previous_output;
    
    # Generate metaplot
    $self->add(
        $self->metaplot(
            bam_file  => $bam_file,
            gff_file  => $gff_file,
            feat_type => 'gene'
        )
    );

    # Add metadata to results
    my $annotations = $self->generate_additional_metadata($self->previous_outputs->[1]);
    $self->add_to_previous(
        $self->add_metadata_to_results(
            annotations => $annotations,
            item_id     => ($experiment ? $experiment->id : undef),
            item_type   => 'experiment',
            locked      => 0
        )
    );
}

sub generate_additional_metadata {
    my $self = shift;
    my $metaplot_image_path = shift;
    my $metaplot_params = $self->params->{metaplot_params};

    my @annotations;
    push @annotations, "$metaplot_image_path|||metaplot|parameters: " . join(' ', map { $_.' '.$metaplot_params->{$_} } ('outer', 'inner', 'window'));
    
    return \@annotations;
}

sub metaplot {
    my $self = shift;
    my %opts = @_;
    my $bam_file    = $opts{bam_file};
    my $gff_file    = $opts{gff_file};
    my $feat_type   = $opts{feat_type};

    my $params = $self->params->{metaplot_params};
    my $outer       = $params->{outer}  // 2000;
    my $inner       = $params->{inner}  // 5000;
    my $window_size = $params->{window} // 100;

    my $cmd = catfile($self->conf->{SCRIPTDIR}, 'methylation', 'makeMetaplot.pl');
    $cmd = 'nice ' . $cmd;
    my $output_name = 'metaplot';

    return {
        cmd => $cmd,
        args => [
            ['-o',        $output_name, 0],
            ['-gff',      $gff_file,    0],
            ['-outside',  $outer,       0],
            ['-inside',   $inner,       0],
            ['-w',        $window_size, 0],
            ['-outRange', '',           0],
            ['-featType', $feat_type,   0],
            ['-cpu',      8,            0],
            ['-quiet',    '',           0], # disable frequent printing of feature ID
            ['-bam',      $bam_file,    0]
        ],
        inputs => [
            $bam_file,
            qq[$bam_file.bai],
            $gff_file
        ],
        outputs => [
            catfile($self->staging_dir, "$output_name.tab"),
            catfile($self->staging_dir, "$output_name.png")
        ],
        description => "Generating metaplot"
    };
}

sub add_metadata_to_results {
    my $self = shift;
    my %opts = @_;
    my $item_id     = $opts{item_id};
    my $item_type   = $opts{item_type};
    my $annotations = $opts{annotations}; # array ref
    my $locked      = $opts{locked} // 1;

    my $cmd = catfile($self->conf->{SCRIPTDIR}, "add_metadata_to_results.pl");

    my $log_file = catfile($self->staging_dir, "add_metadata_to_results", "log.txt");

    my $annotations_str = '';
    $annotations_str = join(';', @$annotations) if (defined $annotations && @$annotations);

    my $args = [
        ['-uid',         $self->user->id,               0],
        ['-wid',         $self->workflow->id,           0],
        ['-locked',      $locked,                       0],
        ['-annotations', shell_quote($annotations_str), 0],
        ['-config',      $self->conf->{_CONFIG_PATH},   0],
        ['-log',         $log_file,                     0]
    ];

    if ($item_id && $item_type) {
        push @$args, ['-item_id',   $item_id,   0];
        push @$args, ['-item_type', $item_type, 0];
    }

    return {
        cmd => $cmd,
        args => $args,
        inputs => [],
        outputs => [ $log_file ],
        description => "Adding metadata to results"
    };
}

1;