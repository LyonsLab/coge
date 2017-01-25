package CoGe::Services::API::Annotations;

use Mojo::Base 'Mojolicious::Controller';
use CoGe::Core::Metadata qw( get_type_groups );
use CoGe::Services::Auth;

sub fetch_type_groups {
    my $self = shift;
    my ($db) = CoGe::Services::Auth::init($self);
    my $type_groups = get_type_groups($db);
    $self->render(json => $type_groups) if $type_groups;
}

1;
