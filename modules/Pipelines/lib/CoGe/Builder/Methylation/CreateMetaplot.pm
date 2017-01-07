package CoGe::Builder::Methylation::CreateMetaplot;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);
use Switch;
use CoGe::Core::Storage qw(get_experiment_files);
use CoGe::Builder::Methylation::Metaplot qw(build);

sub get_name {
    #my $self = shift;
    return 'Create Metaplot';
}

sub build {
    my $self = shift;
    
    # Validate inputs #TODO add error messages to STDERR
    my $eid = $self->params->{eid} || $self->params->{experiment_id};
    my $methylation_params = $self->params->{methylation_params};
    unless ($eid && $methylation_params) {
        print STDERR "CoGe::Builder::Methylation::CreateMetaplot ERROR, missing 'eid' or 'methylation_params'\n";
        return;
    }
    my $metaplot_params = $self->params->{methylation_params}->{metaplot_params};
    unless ($metaplot_params) {
        print STDERR "CoGe::Builder::Methylation::CreateMetaplot ERROR, missing 'metaplot_params'\n";
        return;
    }    
    
    # Get experiment
    my $experiment = $self->db->resultset('Experiment')->find($eid);
    unless ($experiment) {
        print STDERR "CoGe::Builder::Methylation::CreateMetaplot ERROR, experiment '$eid' not found\n";
        return;
    }    
    my $genome = $experiment->genome;
    
    # Get input file
    my $bam_file = get_experiment_files($experiment->id, $experiment->data_type)->[0];

    #
    # Build workflow steps
    #
    my @tasks;
    
    # Add methylation analysis workflow
    my $metaplot_workflow = CoGe::Builder::Methylation::Metaplot::build({
        user => $self->user,
        wid => $self->workflow->id,
        genome => $genome,
        bam_file => $bam_file,
        experiment_id => $eid,
#        metadata => $metadata,
        methylation_params => $self->params->{methylation_params}
    });
    push @tasks, @{$metaplot_workflow->{tasks}};
        
    $self->workflow->add_jobs(\@tasks);
    
    return 1;
}

__PACKAGE__->meta->make_immutable;

1;