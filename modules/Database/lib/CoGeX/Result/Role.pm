package CoGeX::Result::Role;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::Role

=cut

__PACKAGE__->table("role");

=head1 ACCESSORS

=head2 role_id

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
  size: 255

=cut

__PACKAGE__->add_columns(
  "role_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "description",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);
__PACKAGE__->set_primary_key("role_id");
# Role has many permissions
__PACKAGE__->has_many('role_permission_connectors'=>"CoGeX::Result::RolePermissionConnector",'role_id');
#__PACKAGE__->has_many('user_groups'=>'CoGeX::Result::UserGroup', 'role_id');
__PACKAGE__->has_many('user_connectors'=>'CoGeX::Result::UserConnector', 'role_id');


################################################ subroutine header begin ##

=head2 is_<ROLE>

 Usage     : 
 Purpose   : 
 Returns   :
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##
sub is_owner
 {
   return shift->name =~ /owner/i;
 }
 
sub is_editor
 {
   return shift->name =~ /editor/i;
 } 
 
sub is_reader
 {
   return shift->name =~ /reader/i;
 }
 
################################################ subroutine header begin ##

=head2 groups

 Usage     : $self->groups
 Purpose   : alias for $self->user_groups
 Returns   : DBIX::Class for user_group objects
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##
sub groups
 {
   return shift->user_groups(@_);
 }

################################################ subroutine header begin ##

=head2 permissions

 Usage     : 
 Purpose   : Returns permission objects
 Returns   : wantArray of permission objects
 Argument  : None
 Throws    : None
 Comments  : 

=cut

################################################## subroutine header end ##
sub permissions {

    my $self = shift;
    my @permissions=();

    foreach my $pc ($self->role_permission_connectors())
      {
        push@permissions, $pc->permission;
    }

    return wantarray ? @permissions : \@permissions;
}

1;


=head1 AUTHORS

 Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
