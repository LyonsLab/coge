package CoGe::Builder::Export::Genome;

use Moose;
with qw(CoGe::Builder::Buildable);

use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Gff;

sub get_name {
    return "Export Genome";
}

sub build {
    my $self = shift;
    
    my $request = { #FIXME better way to do this?
        params      => $self->params,
        requester   => $self->requester,
        db          => $self->db,
        user        => $self->user,
        conf        => $self->conf,
        workflow    => $self->workflow,
        staging_dir => $self->staging_dir,
        result_dr   => $self->result_dir,
        outputs     => $self->outputs
    };

    my $fasta = CoGe::Builder::Export::Fasta->new($request);
    $fasta->build();
    
    my $gff = CoGe::Builder::Export::Gff->new($request);
    $gff->build();
    
    return 1;
}

1;
