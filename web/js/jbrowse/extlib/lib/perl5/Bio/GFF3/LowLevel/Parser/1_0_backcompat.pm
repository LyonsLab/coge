package Bio::GFF3::LowLevel::Parser::1_0_backcompat;
BEGIN {
  $Bio::GFF3::LowLevel::Parser::1_0_backcompat::AUTHORITY = 'cpan:RBUELS';
}
{
  $Bio::GFF3::LowLevel::Parser::1_0_backcompat::VERSION = '1.5';
}
# ABSTRACT: compatibility layer to support Bio::GFF3::LowLevel::Parser 1.0 API

use strict;
use warnings;


sub new {
    my $class = shift;
    my $self = {
        parser => Bio::GFF3::LowLevel::Parser->open( @_ ),
        item_buffer => [],
    };
    return bless $self, $class;
}


sub next_item {
    my ( $self ) = @_;
    my $item_buffer = $self->{item_buffer};

    # try to get more items if the buffer is empty
    $self->_buffer_items unless @$item_buffer;

    # return the next item if we have some
    return shift @$item_buffer if @$item_buffer;

    # if we were not able to get any more items, return nothing
    return;
}

sub _buffer_items {
    my ( $self ) = @_;
    my $item = $self->{parser}->next_item;
    unless( ref $item eq 'ARRAY' ) {
        push @{$self->{item_buffer}}, $item;
        return;
    }

    # convert all the features and child features back to non-arrayrefs
    push @{$self->{item_buffer}}, $self->_xform_1x( $item );
}

# take a 2.x feature arrayref, return a list of 1.x-compliant features
sub _xform_1x {
    my ( $self, $f ) = @_;
    return $f unless ref $f eq 'ARRAY';
    for my $line (@$f) {
        for my $attr ( 'child_features', 'derived_features' ) {
            $line->{$attr} = [
                map $self->_xform_1x( $_ ),
                grep $_ != $f,
                @{ $line->{$attr} }
             ];
        }
    }
    return @$f;
}

1;

__END__
=pod

=encoding utf-8

=head1 NAME

Bio::GFF3::LowLevel::Parser::1_0_backcompat - compatibility layer to support Bio::GFF3::LowLevel::Parser 1.0 API

=head1 FUNCTIONS

=head2 new( $file_or_filehandle, ... )

Constructor for 1.0 backcompat layer.  Do not use directly.

=head2 next_item

1.0 backcompat layer, wraps next_item() from the parser to transform features back to 1.0 format.

=head1 AUTHOR

Robert Buels <rmb32@cornell.edu>

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2012 by Robert Buels.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

