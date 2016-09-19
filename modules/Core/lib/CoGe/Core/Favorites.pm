package CoGe::Core::Favorites;
use strict;
use warnings;

use Moose;

# Attributes -------------------------------------------------------------------
has 'user' => (
    is       => 'ro',
    required => 1
    #isa      =>
);

# Private methods --------------------------------------------------------------
sub BUILD { # called after constructor
    my $self = shift;
    
    if (!$self->user || $self->user->is_public) {
        # error
        return;
    }
}

# Public methods ---------------------------------------------------------------
sub is_favorite {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    return unless $item;
    return unless ($self->user && !$self->user->is_public);
    
    my ($favorite) = $self->user->favorites({ child_id => $item->id, child_type => $item->item_type });
    return defined($favorite);
}

sub get {
    my $self = shift;
    my %opts = @_;
    my $onlyMine = $opts{onlyMine}; # optional flag to exclude items the user doesn't own
    my $notMine  = $opts{notMine};  # optional flag to exclude items the user owns
    return unless ($self->user && !$self->user->is_public);
    
    my @favorites;
    foreach ($self->user->favorites()) {
        my $item = $_->child;
        next if ($notMine  && $self->user->is_owner(item => $item));
        next if ($onlyMine && !$self->user->is_owner(item => $item));
        push @favorites, $item; 
    }
    
    return wantarray ? @favorites : \@favorites;
}

sub add {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    return if $self->is_favorite($item);
    
    my $favorite = $self->user->add_to_favorites({ child_id => $item->id, child_type => $item->item_type });
    return $favorite;
}

sub remove {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    return unless $self->is_favorite($item);
    
    my ($favorite) = $self->user->favorites({ child_id => $item->id, child_type => $item->item_type });
    $favorite->delete;
    return 1;
}

sub toggle {
    my $self = shift;
    my $item = shift; # genome/experiment/notebook DBIX object
    
    if ($self->is_favorite($item)) {
        $self->remove($item);
        return 0;
    }
    else {
        $self->add($item);
        return 1;
    }
}

1;
