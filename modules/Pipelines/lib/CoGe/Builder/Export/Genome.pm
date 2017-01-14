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
    
    my $outputs = $self->params->{output_types};
    my %outputs = map { lc($_) => 1 } @$outputs if $outputs;

    if (!$outputs || $outputs{'fasta'}) {
        my $fasta = CoGe::Builder::Export::Fasta->new($self);
        $fasta->build();
    }

    if (!$outputs || $outputs{'gff'} || $outputs{'bed'} || $outputs{'tbl'}) {
        my $gff = CoGe::Builder::Export::Gff->new($self);
        $gff->build();
    }
}

__PACKAGE__->meta->make_immutable;

1;
