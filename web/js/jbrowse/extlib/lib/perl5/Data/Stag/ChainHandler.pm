package Data::Stag::ChainHandler;

=head1 NAME

  Data::Stag::ChainHandler

=head1 SYNOPSIS


=cut

=head1 DESCRIPTION


=head1 PUBLIC METHODS -

=cut

use strict;
use base qw(Data::Stag::Base Data::Stag::Writer);

use vars qw($VERSION);
$VERSION="0.11";

sub init {
    my $self = shift;
    return;
}

sub subhandlers {
    my $self = shift;
    $self->{_subhandlers} = shift if @_;
    return $self->{_subhandlers};
}

sub blocked_event {
    my $self = shift;
    if (@_) {
        my $e = shift;
        $self->{_blocked_event} = $e;
        unless (ref($e)) {
            $e = [$e];
        }
        my %h = map {$_=>1} @$e;
        $self->blocked_event_h(\%h);
    }
    return $self->{_blocked_event};
}

sub blocked_event_h {
    my $self = shift;
    $self->{_blocked_event_h} = shift if @_;
    return $self->{_blocked_event_h};
}

sub is_blocked {
    my $self = shift;
    my $e = shift;
    if (!$e) {
        $self->throw("must pass arg to is_blocked");
    }
    my $is = $self->blocked_event_h->{$e};
    return $is;
}

sub start_event {
    my $self = shift;
    my $ev = shift;
    my $stack = $self->elt_stack;
    push(@$stack, $ev);

    my $sh = $self->subhandlers;
    my $is_blocked = grep {$self->is_blocked($_)} @$stack;
    if ($is_blocked) {
        $sh->[0]->start_event($ev);
    }
    else {
        foreach (@$sh) {
            $_->start_event($ev);
        }
    }
}

sub evbody {
    my $self = shift;
    my $ev = shift;
    my @args = @_;

    my $stack = $self->elt_stack;

    my $sh = $self->subhandlers;
    if (grep {$self->is_blocked($_)} @$stack) {
        $sh->[0]->evbody($ev, @args);
    }
    else {
        foreach (@$sh) {
            $_->evbody($ev, @args);
        }
    }
    
    return;
}

sub end_event {
    my $self = shift;
    my $ev = shift;

    my $stack = $self->elt_stack;
    $ev = pop @$stack;
    if (!$ev) {
        $self->throw("no event name on stack");
    }

    my $sh = $self->subhandlers;


    my $inside_blocked = grep {$self->is_blocked($_)} @$stack;
    if ($self->is_blocked($ev) &&
	!$inside_blocked) {

	# condition:
	# end of a blocked event, and we are 
	# not inside another blocked event
        my ($h, @rest) = @$sh;

        my @R = $h->end_event($ev);
        foreach my $handler (@rest) {
	    if (@R) {
		$handler->event(@$_) foreach @R;
		@$_ = () foreach @R;
	    }
        }
        my $tree = $h->tree;
        $tree->free;
    }
    else {

        if (grep {$self->is_blocked($_)} @$stack) {
            $sh->[0]->end_event($ev);
        }
        else {
            foreach (@$sh) {
                $_->end_event($ev);
            }
        }
    }
}

1;
