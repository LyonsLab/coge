package CoGe::Builder::Methylation::Metaplot;

use Moose;
extends 'CoGe::Builder::Methylation::Analyzer';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my $bam_file = shift;

    unless ($bam_file) { # for when called from ExperimentView
        my $experiment = self->request->experiment;
        $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];
    }

    my $genome = $self->request->genome;

    # Make sure genome is annotated as is required by metaplot script
    my $isAnnotated = $genome->has_gene_features;
    unless ($isAnnotated) {
        CoGe::Exception::Generic->throw(message => 'Genome must be annotated to generate metaplot');
    }

    #
    # Build the workflow
    #

    # Generate cached gff
    my $gff_task = create_gff_generation_job(
        gid => $genome->id,
        organism_name => $genome->organism->name
    );
    my $gff_file = $gff_task->{outputs}->[0];
    push @tasks, $gff_task;
    
    # Generate metaplot
    my $metaplot_task = create_metaplot_job(
        bam_file => $bam_file,
        gff_file => $gff_file,
        outer => $metaplot_params->{outer},
        inner => $metaplot_params->{inner},
        window_size => $metaplot_params->{window},
        feat_type => 'gene',
        staging_dir => $staging_dir
    );
    push @done_files, @{$metaplot_task->{outputs}};
    push @tasks, $metaplot_task;
    
    # Add metadata
    my $annotations = generate_additional_metadata($metaplot_params, $metaplot_task->{outputs}[1]);
    my $metadata_task = add_metadata_to_results_job(
        user => $user,
        wid => $wid,
        annotations => $annotations,
        staging_dir => $staging_dir,
        done_files => \@done_files,
        item_id => $experiment_id,
        item_type => 'experiment',
        locked => 0
    );
    push @tasks, $metadata_task;

    return {
        tasks => \@tasks,
        done_files => \@done_files
    };
}

sub generate_additional_metadata {
    my $metaplot_params = shift;
    my $metaplot_image_path = shift;
    
    my @annotations;
    push @annotations, "$metaplot_image_path|||metaplot|parameters: " . join(' ', map { $_.' '.$metaplot_params->{$_} } ('outer', 'inner', 'window'));
    
    return \@annotations;
}

sub create_metaplot_job {
    my %opts = @_;
    my $bam_file    = $opts{bam_file};
    my $gff_file    = $opts{gff_file};
    my $outer       = $opts{outer} // 2000;
    my $inner       = $opts{inner} // 5000;
    my $window_size = $opts{window} // 100;
    my $feat_type   = $opts{feat_type};
    my $staging_dir = $opts{staging_dir};

    my $cmd = catfile($CONF->{SCRIPTDIR}, 'methylation', 'makeMetaplot.pl');
    $cmd = 'nice ' . $cmd;
    my $output_name = 'metaplot';

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-o', $output_name, 0],
            ['-gff', $gff_file, 0],
            ['-outside', $outer, 0],
            ['-inside', $inner, 0],
            ['-w', $window_size, 0],
            ['-outRange', '', 0],
            ['-featType', $feat_type, 0],
            ['-cpu', 8, 0],
            ['-quiet', '', 0], # disable frequent printing of feature ID
            ['-bam', $bam_file, 0]
        ],
        inputs => [
            $bam_file,
            $bam_file . '.bai',
            $gff_file
        ],
        outputs => [
            catfile($staging_dir, "$output_name.tab"),
            catfile($staging_dir, "$output_name.png")
        ],
        description => "Generating metaplot"
    };
}

1;