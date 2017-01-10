package CoGe::Builder::Export::Genome;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Gff;

sub get_name {
    return "Export Genome"; #TODO add genome id and outputs
}

sub build {
    my $self = shift;
    
    my $request = { #FIXME better way to do this? See Moose constructors
        params      => $self->params,
        db          => $self->db,
        user        => $self->user,
        conf        => $self->conf,
        workflow    => $self->workflow,
        staging_dir => $self->staging_dir,
        result_dr   => $self->result_dir,
        outputs     => $self->outputs
    };

    my $outputs = $self->params->{output_types};
    my %outputs;
    %outputs = map { lc($_) => 1 } @$outputs if $outputs;

    if (!$outputs || $outputs{'fasta'}) {
        my $fasta = CoGe::Builder::Export::Fasta->new( $request );
        $fasta->build();
    }

    if (!$outputs || $outputs{'gff'} || $outputs{'bed'} || $outputs{'tbl'}) {
        my $gff = CoGe::Builder::Export::Gff->new( $request );
        $gff->build();
    }

    return 1;
}

__PACKAGE__->meta->make_immutable;

1;
