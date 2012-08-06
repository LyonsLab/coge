package CoGeX_dev::Result::ListType;

use strict;
use warnings;

use base 'DBIx::Class::Core';


=head1 NAME

CoGeX_dev::Result::ListType

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
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "description",
  { data_type => "varchar", is_nullable => 1, size => 1024 },
);
__PACKAGE__->set_primary_key("list_type_id");

__PACKAGE__->has_many('lists'=>"CoGeX_dev::Result::List",'list_type_id');



################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : None
 Throws    : None
 Comments  : 



=cut

################################################## subroutine header end ##


1;
