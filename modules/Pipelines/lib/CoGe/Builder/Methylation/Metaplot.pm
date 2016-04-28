package CoGe::Builder::Methylation::Metaplot;

use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Basename qw(basename);
use File::Spec::Functions qw(catdir catfile);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_workflow_paths);
use CoGe::Builder::CommonTasks qw(create_gff_generation_job);

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
    my $bam_file = $opts->{raw_bam_file}; # path to bam file -- important: this should be the unsorted version
                                          # see COGE-706 and http://seqanswers.com/forums/showthread.php?t=45192
#    my $metadata = $opts->{metadata};
#    my $additional_metadata = $opts->{additional_metadata};
    my $wid = $opts->{wid};
    my $methylation_params = $opts->{methylation_params};

    # Setup paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $wid);

    # Set metadata for the pipeline being used
#    my $annotations = generate_additional_metadata($read_params, $methylation_params);
#    my @annotations2 = CoGe::Core::Metadata::to_annotations($additional_metadata);
#    push @$annotations, @annotations2;

    # Make sure genome is annotated as is required by metaplot script
    my $isAnnotated = $genome->has_gene_features;
    return unless $isAnnotated;

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
        outside => 2000,
        inside => 5000,
        window_size => 100
    );
    push @done_files, @{$metaplot_task->{outputs}};
    push @tasks, $metaplot_task;

    return {
        tasks => \@tasks,
        done_files => \@done_files
    };
}

sub create_metaplot_job {
    my %opts = @_;
    my $bam_file    = $opts{bam_file};
    my $gff_file    = $opts{gff_file};
    my $outside     = $opts{outside};
    my $inside      = $opts{inside};
    my $window_size = $opts{window_size};

    my $cmd = catfile($CONF->{SCRIPTDIR}, 'methylation', 'makeMetaplot.pl');

    return {
        cmd => $cmd,
        script => undef,
        args => [
            ['-o', 'metaplot', 0],
            ['-gff', $gff_file, 0],
            ['-outside', $outside, 0],
            ['-inside', $inside, 0],
            ['-w', $window_size, 0],
            ['-outRange', '', 0],
            ['', '', $bam_file]
        ],
        inputs => [
            $bam_file,
            $bam_file . '.bai',
            $gff_file
        ],
        outputs => [
            'metaplot.tab',
            'metaplot.pdf'
        ],
        description => "Generating metaplot..."
    };
}

1;