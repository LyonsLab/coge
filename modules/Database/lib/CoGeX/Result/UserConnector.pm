package CoGeX::Result::UserConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';
#use CoGeX;

=head1 NAME

CoGeX::Result::UserConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user_connector> table in the CoGe database.
The C<user_connector> table is used to associate C<experiment,genome,feature> records with C<user> records.

=head1 DESCRIPTION

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

# my $child_types = CoGeX::list_child_types();


__PACKAGE__->table("user_connector");
__PACKAGE__->add_columns(
	"user_connector_id",
	{ data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
	"parent_id", # user_id or user_group_id
	{ data_type => "INT", default_value => "", is_nullable => 1, size => 11 }, 
	"parent_type",
	{ data_type => "TINYINT", default_value => "", is_nullable => 0, size => 1 },	
	"child_id",
	{ data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
	"child_type",
	{ data_type => "TINYINT", default_value => "", is_nullable => 0, size => 1 },
	"role_id",
	{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },	
);
__PACKAGE__->set_primary_key("user_connector_id");

__PACKAGE__->belongs_to("experiment" => "CoGeX::Result::Experiment", "child_id");
__PACKAGE__->belongs_to("genome"     => "CoGeX::Result::Genome",     "child_id");
__PACKAGE__->belongs_to("feature"    => "CoGeX::Result::Feature",    "child_id");

__PACKAGE__->belongs_to("role" => "CoGeX::Result::Role", "role_id" );
__PACKAGE__->belongs_to("user" => "CoGeX::Result::User", { "foreign.user_id" => "self.parent_id" } );
__PACKAGE__->belongs_to("user_group" => "CoGeX::Result::UserGroup", { "foreign.user_group_id" => "self.parent_id" } );



################################################ subroutine header begin ##

=head2 type

 Usage     : 
 Purpose   : Alias to the child_type() method.
 Returns   : See child_type()
 Argument  : None
 Throws    : 
 Comments  : 
 
See Also   : 

=cut

################################################## subroutine header end ##

sub type
{
	return shift->child_type();
}

################################################ subroutine header begin ##

=head2 is_list

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : None
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

# sub is_list
# {
# 	return shift->child_type() == $child_types->{list};
# }

################################################ subroutine header begin ##

=head2 is_genome

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : None
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

# sub is_genome
# {
# 	return shift->child_type() == $child_types->{genome};
# }

################################################ subroutine header begin ##

=head2 is_feature

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : None
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

# sub is_feature
# {
# 	return shift->child_type() == $child_types->{feature};
# }


################################################ subroutine header begin ##

=head2 is_experiment

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : None
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

# sub is_experiment
# {
# 	return shift->child_type() == $child_types->{experiment};
# }


################################################ subroutine header begin ##

=head2 child

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : None
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

# sub child
# {
# 	my $self = shift;
	
# 	if ($self->is_experiment) {
# 		return $self->experiment;
# 	}
# 	elsif ($self->is_genome) {
# 		return $self->genome;
# 	}
# 	elsif ($self->is_feature) {
# 		return $self->feature;
# 	}
# 	elsif ($self->is_list) {
# 		return $self->child_list;	
# 	}
# 	else {
# 		die;
# 	}

# 	return;
# }


1;

=head1 BUGS

=head1 SUPPORT

=head1 AUTHORS

 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
