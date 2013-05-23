package Cache::Ref::Role::WithDoublyLinkedList;
BEGIN {
  $Cache::Ref::Role::WithDoublyLinkedList::AUTHORITY = 'cpan:NUFFIN';
}
BEGIN {
  $Cache::Ref::Role::WithDoublyLinkedList::VERSION = '0.04';
}
use MooseX::Role::Parameterized;

parameter name => (
    isa => "Str",
    required => 1,
);

parameter [qw(value_offset next_offset prev_offset)] => (
    isa => "Int",
    required => 1,
);

foreach my $method (qw(head tail shift pop unshift push)) {
    parameter "${method}_method" => (
        isa => "Str",
        is  => "ro",
        lazy => 1,
        default => sub {
            my $p = shift;
            $p->name . "_" . $method;
        }
    );
}

role {
    my $p = shift;

    my $name = $p->name;

    my ( $head_attr, $tail_attr ) = ("${name}_head", "${name}_tail");

    # technically a doubly linked list is just the inverse of itself
    # so really this is like including two inverted parameterized roles (the
    # notion of next/prev is reversed) for the endian methods
    for (
        {
            attr => $head_attr,
            rev_attr => $tail_attr,
            vo => $p->value_offset,
            no => $p->next_offset,
            po => $p->prev_offset,
            head_method => $p->head_method,
            shift_method => $p->shift_method,
            unshift_method => $p->unshift_method
        },
        {
            attr => $tail_attr,
            rev_attr => $head_attr,
            vo => $p->value_offset,
            no => $p->prev_offset,
            po => $p->next_offset,
            head_method => $p->tail_method,
            shift_method => $p->pop_method,
            unshift_method => $p->push_method,
        },
    ) {
        my ( $attr, $rev_attr, $value, $next, $prev ) = @{$_}{qw(attr rev_attr vo no po)};

        has $attr => ( is => "rw" );

        method $_->{head_method} => sub { (shift->$attr || return)->[$value] };

        method $_->{shift_method} => sub {
            my $self = shift;

            my $node = $self->$attr or return;

            if ( my $neighbor = $node->[$next] ) {
                $self->$attr($neighbor);
                $neighbor->[$prev] = undef;
            } else {
                # list is empty, clear both attrs
                $self->$attr(undef);
                $self->$rev_attr(undef);
            }

            return $node->[$value];
        };

        method $_->{unshift_method} => sub {
            my ( $self, @values ) = @_;

            my $head = $self->$attr;

            my @ret;

            foreach my $v ( reverse @values ) {
                # cons up a new list
                my $new = [];
                $new->[$value] = $v;
                $new->[$next] = $head;
                $head->[$prev] = $new;
                push @ret, $new;
                $head = $new;
            }

            $self->$attr($head);
            $self->$rev_attr($ret[0]) unless $self->$rev_attr;

            return ( @ret == 1 ? $ret[0] : reverse @ret );
        }
    }

    # these methods are per linked list

    my ( $next_offset, $prev_offset ) = ( $p->next_offset, $p->prev_offset );

    method "${name}_clear" => sub {
        my $self = shift;

        my $cur = $self->$head_attr;

        while ( $cur ) {
            my $next = $cur->[$next_offset];
            $cur = (); # FIXME not so general purpose
            #@{$cur}[$next_offset, $prev_offset] = (); # more general purpose? don't care...
            $cur = $next;
        }

        $self->$head_attr(undef);
        $self->$tail_attr(undef);
    };

    method "${name}_set_next" => sub {
        my ( $self, $node, $next ) = @_;
        $node->[$next_offset] = $next;
    };

    method "${name}_set_prev" => sub {
        my ( $self, $node, $prev ) = @_;
        $node->[$prev_offset] = $prev;
    };

    method "${name}_link_sequence" => sub {
        my ( $self, @nodes ) = @_;

        return unless @nodes;

        my $prev = shift @nodes;

        delete $prev->[$prev_offset];

        foreach my $node ( @nodes ) {
            $prev->[$next_offset] = $node; # $prev->next($l)
            $node->[$prev_offset] = $prev;
            $prev = $node;
        }

        delete $prev->[$next_offset];

        return;
    };

    method "${name}_splice" => sub {
        my ( $self, @nodes ) = @_;

        return unless @nodes;

        foreach my $node ( @nodes ) {
            # detach node from its current place in the list
            if ( $node->[$prev_offset] ) {
                $node->[$prev_offset][$next_offset] = $node->[$next_offset];
            } else {
                # $node is currently head, so unmark it as such
                $self->$head_attr($node->[$next_offset]);
            }

            if ( $node->[$next_offset] ) {
                $node->[$next_offset][$prev_offset] = $node->[$prev_offset];
            } else {
                # $node is currently tail, so unmark it as such
                $self->$tail_attr($node->[$prev_offset]);
            }

            delete @{ $node }[$prev_offset, $next_offset];
        }

        return
    };
};

# ex: set sw=4 et:

__PACKAGE__;


__END__
=pod

=encoding utf-8

=head1 NAME

Cache::Ref::Role::WithDoublyLinkedList

=head1 AUTHOR

Yuval Kogman

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2010 by Yuval Kogman.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

