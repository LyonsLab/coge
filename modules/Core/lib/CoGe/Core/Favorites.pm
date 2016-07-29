package CoGe::Core::Favorites;
use strict;
use warnings;

use Moose;
use Data::Dumper;
use CoGe::Core::Notebook;

# Attributes -------------------------------------------------------------------
has 'user' => (
    is       => 'ro',
    required => 1
    #isa      =>
);

has 'notebook' => (
    is       => 'rw',
    required => 0
);

# Private methods --------------------------------------------------------------
sub BUILD { # called after constructor
    my $self = shift;
    
    if ($self->user->is_public) {
        # error
        return;
    }
    
    my $notebook;
    unless ($notebook = $self->_find_favorites_notebook()) {
        # error
    }
    $self->notebook($notebook);
}

sub _find_favorites_notebook {
    my $self = shift;
    
    # Search for favorites notebook
    foreach my $notebook ($self->user->lists) {
        if ($notebook->name eq 'Favorites') {
            return $notebook;
        }
    }
    
    # Create favorites notebook
    my $notebook = $self->_create_favorites_notebook();
    unless ($notebook) {
        # error
        return;
    }
    
    return $notebook;
}

sub _create_favorites_notebook {
    my $self = shift;
    
    my $db = $self->user->result_source->schema;
    return create_notebook(
        db => $db, user => $self->user, type => 'mixed', name => 'Favorites',
        descriptoin => 'This notebook contains genomes you have marked as your favorite.  It is created automatically by CoGe.');
}

sub _get_favorites {
    my $self = shift;
    
    my @favorites;
    push @favorites, @{$self->notebook->genomes} if $self->notebook->genomes;
    push @favorites, @{$self->notebook->experiments} if $self->notebook->experiments;
    #push @favorites, $notebook->features if $notebook->features;
    
    return \@favorites;
}

# Public methods ---------------------------------------------------------------
sub is_favorite {
    my $self = shift;
    my $item = shift; # genome/experiment DBIX object
    
    foreach ($self->notebook->child_connectors({ child_id => $item->id, child_type => $item->item_type })) {
        return 1;
    }
    
    return 0;
}

sub add {
    my $self = shift;
    my $item = shift; # genome/experiment DBIX object
    return if $self->is_favorite($item);
    
    my $db = $self->user->result_source->schema;
    add_items_to_notebook(db => $db, user => $self->user, notebook => $self->notebook, items => [ $item ]);
}

sub remove {
    my $self = shift;
    my $item = shift; # genome/experiment DBIX object
    
    my $db = $self->user->result_source->schema;
    remove_items_from_notebook(db => $db, user => $self->user, notebook => $self->notebook, items => [ $item ]);
}

sub toggle {
    my $self = shift;
    my $item = shift; # genome/experiment DBIX object
    
    if ($self->is_favorite($item)) {
        $self->remove($item);
        return 0;
    }
    else {
        $self->add($item);
        return 1;
    }
}

#__PACKAGE__->meta->make_immutable;
#no Moose;
#1;
#__END__
1;
