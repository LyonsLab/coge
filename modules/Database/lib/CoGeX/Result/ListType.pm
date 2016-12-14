package CoGeX::Result::ListType;

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# MDB 12/14/16:  PER COGE-800 THIS TABLE HAS BEEN DEPRECATED AND NEEDS TO BE REMOVED EVENTUALLY.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

CoGeX::Result::ListType

=cut

__PACKAGE__->table("list_type");

=head1 ACCESSORS

=head2 list_type_id

  data_type: 'bigint'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 description

  data_type: 'varchar'
  is_nullable: 1
  size: 1024

=cut

__PACKAGE__->add_columns(
    "list_type_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "name",
    { data_type => "varchar", is_nullable => 0, size => 255 },
    "description",
    { data_type => "varchar", is_nullable => 1, size => 1024 },
);
__PACKAGE__->set_primary_key("list_type_id");

__PACKAGE__->has_many( 'lists' => "CoGeX::Result::List", 'list_type_id' );

################################################ subroutine header begin ##

=head2 is_*

 Usage     :
 Purpose   : Test for particular list type
 Returns   : boolean
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub is_genome {
    return shift->list_type_id == 1;
}

sub is_experiment {
    return shift->list_type_id == 2;
}

sub is_owner {
    return shift->list_type_id == 3;
}

sub is_feature {
    return shift->list_type_id == 4;
}

sub is_mixed {
    return shift->list_type_id == 5;
}

sub is_other {
    return shift->list_type_id == 6;
}

1;

=head1 AUTHORS

 Eric Lyons

=head1 COPYRIGHT 2012

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
