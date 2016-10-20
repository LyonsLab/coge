package CoGe::Builder::Methylation::Metaplot;

use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Spec::Functions qw(catfile);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_workflow_paths);
use CoGe::Builder::CommonTasks qw(create_gff_generation_job add_metadata_to_results_job);

our $CONF = CoGe::Accessory::Web::get_defaults();

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw(Exporter);
    @EXPORT  = qw(build);
}

sub build {
    my $opts = shift;
    my $genome = $opts->{genome};
    my $user = $opts->{user};
    my $bam_file = $opts->{bam_file}; # path to sorted & indexed bam file
    my $experiment_id = $opts->{experiment_id}; # for when called from ExperimentView
#    my $metadata = $opts->{metadata};
#    my $additional_metadata = $opts->{additional_metadata};
    my $wid = $opts->{wid};
    my $methylation_params = $opts->{methylation_params};
    my $metaplot_params = $methylation_params->{metaplot_params};
    unless ($genome && $user && $bam_file && $wid && $methylation_params && $metaplot_params) {
        print STDERR " CoGe::Builder::Methylation::Metaplot ERROR, missing inputs: ",
            ($genome ? '' : ' genome '),
            ($user ? '' : ' user '),
            ($bam_file ? '' : ' bam_file '),
            ($wid ? '' : ' wid '),
            ($methylation_params ? '' : ' methylation_params '),
            ($metaplot_params ? '' : ' metaplot_params '),
            "\n";
        return;
    }

    # Setup paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);

    # Make sure genome is annotated as is required by metaplot script
    my $isAnnotated = $genome->has_gene_features;
    unless ($isAnnotated) {
        print STDERR "CoGe::Builder::Methylation::Metaplot ERROR, genome must be annotated to generate metaplot\n";
        return;
    }

    #
    # Build the workflow
    #
    my (@tasks, @done_files);
    
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