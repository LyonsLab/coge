package CoGeX::Result::Log;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::Log

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

__PACKAGE__->table("log");
__PACKAGE__->add_columns(
    "log_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "time",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
    "user_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "type",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 1 },
    "page",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 0,
        size          => 255
    },
    "description",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 255
    },
    "link",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 255
    },
    "status",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 1 },
    "comment",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 255
    },
    "parent_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "parent_type",
    { data_type => "TINYINT", default_value => undef, is_nullable => 1, size => 1 },
);
__PACKAGE__->set_primary_key("log_id");

__PACKAGE__->belongs_to( 'user' => "CoGeX::Result::User", 'user_id' );

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : generate a string of information about the log entry
 Returns   : a string
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub info {
    my $self = shift;
    my $user_name = ( $self->user ? $self->user->user_name : '' );
    return
        $self->time . ' '
      . $user_name . ' '
      . $self->page . ' '
      . $self->description . ' '
      . $self->link . ' '
      . $self->comment;
}

sub short_info {
    my $self = shift;
    return
        $self->time . ' | '
      . $self->page . ' | '
      . $self->description
      . ( $self->comment ? ' | ' . $self->comment : '' );
}

sub is_important {
    return shift->status == 1;
}

1;
