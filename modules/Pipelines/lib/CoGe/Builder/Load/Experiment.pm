package CoGe::Builder::Load::Experiment;

use Moose;

use Data::Dumper qw(Dumper);
use Switch;
use File::Spec::Functions qw(catfile);

use CoGe::Core::Storage qw(get_workflow_paths);
use CoGe::Builder::CommonTasks;
use CoGe::Builder::Common::Alignment qw(build);
use CoGe::Builder::Expression::qTeller qw(build);
use CoGe::Builder::SNP::CoGeSNPs qw(build);
use CoGe::Builder::SNP::Samtools qw(build);
use CoGe::Builder::SNP::Platypus qw(build);
use CoGe::Builder::SNP::GATK qw(build);

sub build {
    my $self = shift;
    
    # Validate inputs
    my $gid = $self->params->{gid};
    return unless $gid;
    my $data = $self->options->{source_data};
    return unless (defined $data && @$data);
    my $metadata = $self->params->{metadata};
    return unless $metadata;
    
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    # Get genome
    my $genome = $self->db->resultset('Genome')->find($gid);
    return unless $genome;
    # TODO add permissions check here -- or will it happen in Request::Genome?
    
    # Initialize workflow
    $self->workflow($self->jex->create_workflow(name => "Load Experiment", init => 1));
    return unless $self->workflow->id;
    
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $self->workflow->id);
    $self->workflow->logfile(catfile($result_dir, "debug.log"));
    
    # Build workflow steps
    my @tasks;
    my $file_type = $data->[0]->{file_type}; # type of first data file
    
    if ( $file_type eq 'fastq' || $file_type eq 'bam' ) {
        my $bam_file;
        
        if ( $file_type eq 'fastq' ) {
            # Add alignment workflow
            my ($alignment_tasks, $alignment_outputs) = CoGe::Builder::Common::Alignment::build(
                user => $self->user,
                wid => $self->workflow->id,
                input_files => $data,
                genome => $genome,
                metadata => $metadata,
                options => $self->options,
                cutadapt_params => $self->params->{cutadapt_params},
                alignment_params => $self->params->{alignment_params}
            );
            push @tasks, @$alignment_tasks;
            $bam_file = $alignment_outputs->{bam_file};
        }
        elsif ( $file_type eq 'bam' ) {
            $bam_file = $data->[0]->{path};
        }
        else { # error -- should never happen
            die "invalid file type";
        }
        
        # Add expression workflow (if specified)
        if ( $self->params->{expression_params} ) {
            my @expression_tasks = CoGe::Builder::Expression::qTeller::build(
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_file => $bam_file,
                metadata => $metadata,
                options => $self->options,
                params => $self->params->{expression_params}
            );
            push @tasks, @expression_tasks;
        }
        
        # Add SNP workflow (if specified)
        if ( $self->params->{snp_params} ) {
            my $method = $self->params->{snp_params}->{method};
            my $params = {
                user => $self->user,
                wid => $self->workflow->id,
                genome => $genome,
                input_file => $bam_file,
                metadata => $metadata,
                options => $self->options
            };
            
            my @snp_tasks;
            switch ($method) { # TODO move this into subroutine
                case 'coge'     { @snp_tasks = CoGe::Builder::SNP::CoGeSNPs::build($params); }
                case 'samtools' { @snp_tasks = CoGe::Builder::SNP::Samtools::build($params); }
                case 'platypus' { @snp_tasks = CoGe::Builder::SNP::Platypus::build($params); }
                case 'gatk'     { @snp_tasks = CoGe::Builder::SNP::GATK::build($params); }
                else            { die "unknown SNP method"; }
            }
            push @tasks, @snp_tasks;
        }
    }
    # Else, all other file types
    else {
        # Submit workflow to generate experiment
        push @tasks, create_load_experiment_job(
            user => $self->user,
            staging_dir => $staging_dir,
            result_dir => $result_dir,
            wid => $self->workflow->id,
            gid => $genome->id,
            input_file => $data->[0]->{path},
            metadata => $metadata
        );
    }

#    print STDERR Dumper \@tasks, "\n";
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

with qw(CoGe::Builder::Buildable);

1;