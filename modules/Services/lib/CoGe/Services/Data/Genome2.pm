package CoGe::Services::Data::Genome2;
use Mojo::Base 'Mojolicious::Controller';

sub search {
    my $self = shift;
    $self->render(json => { success => Mojo::JSON->true });
}

sub fetch {
    my $self = shift;
    $self->render(json => { success => Mojo::JSON->true });
}

1;
