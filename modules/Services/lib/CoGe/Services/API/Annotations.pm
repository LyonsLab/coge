package CoGe::Services::API::Annotations;

use Data::Dumper;
use Mojo::Base 'Mojolicious::Controller';

use CoGe::Core::Metadata qw( get_type_groups search_annotation_types );
use CoGe::Services::Auth;

sub fetch_type_groups {
    my $self = shift;
    my ($db) = CoGe::Services::Auth::init($self);
    my $type_groups = get_type_groups($db);
    $self->render(json => $type_groups) if $type_groups;
}

sub search_types {
    my $self = shift;
    my $search_term = $self->param('search_term');
    my $type_group = $self->param('type_group');
    my ($db) = CoGe::Services::Auth::init($self);
    my $types = search_annotation_types($search_term, $type_group, $db);
    warn Dumper $types;
    $self->render(json => $types) if $types;
}

1;
