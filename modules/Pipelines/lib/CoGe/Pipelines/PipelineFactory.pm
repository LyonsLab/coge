package CoGe::Pipelines::PipelineFactory;

use Moose;

use CoGe::Builder::GffBuilder;

has 'conf' => (
    is => 'ro',
    required => 1
);

has 'user' => (
    is  => 'ro',
    required => 1
);

has 'jex' => (
    is => 'ro',
    required => 1
);

sub get {
    my ($self, $message) = @_;

    my $options = {
        params   => $message->{parameters},
        options  => $message->{options},
        jex      => $self->jex,
        user     => $self->user,
        conf     => $self->conf
    };

    my $builder;

    if ($message->{type} eq "gff_export") {
        $builder = CoGe::Builder::GffBuilder->new($options);
    } else {
        return;
    }

    # Construct the workflow
    $builder->build;

    # Fetch the workflow constructed
    return $builder->get;
}

1;
