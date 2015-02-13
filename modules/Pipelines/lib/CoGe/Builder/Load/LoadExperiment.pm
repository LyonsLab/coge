package CoGe::Builder::Tools::LoadExperiment;

use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Basename qw(fileparse basename dirname);
use File::Path qw(mkpath);
use File::Spec::Functions qw(catdir catfile);
use JSON qw(decode_json);
use URI::Escape::JavaScript qw(unescape);
use Switch;

use CoGe::Accessory::Jex;
use CoGe::Accessory::TDS qw(read);
use CoGe::Accessory::Utils qw(to_filename);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Workflow;
use CoGe::Core::Storage qw(get_workflow_paths);
use CoGe::Builder::CommonTasks;
#use CoGe::Builder::SNP::CoGeSNPs qw(build);
#use CoGe::Builder::SNP::Samtools qw(build);
#use CoGe::Builder::SNP::Platypus qw(build);
#use CoGe::Builder::SNP::GATK qw(build);

our $CONF = CoGe::Accessory::Web::get_defaults();

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK);
    require Exporter;

    $VERSION   = 0.1;
    @ISA       = qw(Exporter);
    @EXPORT    = qw(run build);
}

sub run {
    my $opts = shift;
    my $user = $opts->{user};

    # Connect to workflow engine and get an id
    my $jex = CoGe::Accessory::Jex->new( host => $CONF->{JOBSERVER}, port => $CONF->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Loading/analyzing experiment data', init => 1 );
    my $wid = $workflow->id;

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths( $user->name, $workflow->id );
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    my $tasks = build({
        result_dir => $result_dir,
        staging_dir => $staging_dir,
        wid => $wid,
        %{$opts},
    });
    $workflow->add_jobs($tasks);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

sub build {
    my $opts = shift;
    my $wid         = $opts->{wid};
    my $user        = $opts->{user};
    my $genome      = $opts->{genome};
    my $items       = $opts->{data};
    my $load_id     = $opts->{load_id};
    my $description = $opts->{description};
    my $result_dir  = $opts->{result_dir};
    my $staging_dir = $opts->{staging_dir};
    my $metadata    = $opts->{metadata};
    my $data        = $opts->{data};
    my $input_file  = $data->[0]->{path};
    my $file_type   = $data->[0]->{file_type};
    my $options     = $opts->{options};
    my $aligner     = ( $options && $options->{aligner} ? $options->{aligner}->{tool} : 'gsnap' );
    
    my @tasks;
    
    if ( $file_type eq 'fastq' || $file_type eq 'bam' ) {
        my $bam_file;
        
        if ( $file_type eq 'fastq' ) {
            # Add alignment workflow -- TODO make into a pipeline object
            (my $alignment_tasks, $bam_file, my $reheader_fasta, my $gff_file) = create_alignment_workflow(
                user => $user,
                wid => $wid,
                input_file => $input_file,
                genome => $genome,
                staging_dir => $staging_dir,
                result_dir => $result_dir,
                metadata => $metadata,
                options => $options
            );
            push @tasks, @$alignment_tasks;
        }
        elsif ( $file_type eq 'bam' ) {
            $bam_file = $opts->{files}[0]->{path};
        }
        else { # error
            die;
        }
        
        # Add expression workflow (if specified)
        if ( $options->{expression_params} ) {
            my @expression_tasks = CoGe::Builder::Expression::qTeller::build(
                user => $user,
                wid => $wid,
                genome => $genome,
                staging_dir => $staging_dir,
                result_dir => $result_dir,
                input_file => $bam_file,
                metadata => $metadata,
                options => $options
            );
            push @tasks, @expression_tasks;
        }
        
        # Add SNP workflow (if specified)
        if ( $options->{snp_params} ) {
            my $method = $options->{snp_params}->{method};
            my $conf = {
                user => $user,
                wid => $wid,
                genome => $genome,
                staging_dir => $staging_dir,
                result_dir => $result_dir,
                input_file => $bam_file,
                metadata => $metadata,
                options => $options
            };
            
            my @snp_tasks;
            switch ($method) { # TODO move this into subroutine
                case 'coge'     { @snp_tasks = CoGe::Builder::SNP::CoGeSNPs::build($conf); }
                case 'samtools' { @snp_tasks = CoGe::Builder::SNP::Samtools::build($conf); }
                case 'platypus' { @snp_tasks = CoGe::Builder::SNP::Platypus::build($conf); }
                case 'gatk'     { @snp_tasks = CoGe::Builder::SNP::GATK::build($conf); }
            }
            push @tasks, @snp_tasks;
        }
    }
    # Else, all other file types
    else {
        # Submit workflow to generate experiment
        push @tasks, create_load_experiment_job(
            user => $user,
            staging_dir => $staging_dir,
            result_dir => $result_dir,
            wid => $wid,
            gid => $genome->id,
            input_file => $input_file,
            metadata => $metadata
        );
    }

    print STDERR Dumper \@tasks, "\n";
    return \@tasks;
}

1;
