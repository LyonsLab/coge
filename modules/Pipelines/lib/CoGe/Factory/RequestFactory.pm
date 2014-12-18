package CoGe::Factory::RequestFactory;

use Moose;
use CoGe::Request::Experiment;
use CoGe::Request::Genome;

has 'user'    => (
    is        => 'ro',
    required  => 1
);

has 'db'      => (
    is        => 'ro',
    required  => 1
);

has 'jex'     => (
    is        => 'ro',
    required  => 1
);

sub get {
    my ($self, $message) = @_;

    my $options = {
        db         => $self->db,
        jex        => $self->jex,
        user       => $self->user,
        options    => $message->{options},
        parameters => $message->{parameters}
    };

    if ($message->{type} eq "export_gff" ||
        $message->{type} eq "export_fasta")
    {
        return CoGe::Request::Genome->new($options);
    }
    elsif ($message->{type} eq "export_experiment") 
    {
        return CoGe::Request::Experiment->new($options);
    }
}

1;
