package CoGeX::Result::UserConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

use Data::Dumper;

=head1 NAME

CoGeX::Result::UserConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<user_connector> table in the CoGe database.
The C<user_connector> table is used to associate C<list,experiment,genome,feature> records with C<user>  and C<user_group>records.

=head1 DESCRIPTION

=head1 USAGE

 use CoGeX;

=head1 METHODS

=head1 AUTHORS

 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

my $node_types = CoGeX::node_types();

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

#FIXME hardcoded parent/child types below
__PACKAGE__->belongs_to("user"			=> "CoGeX::Result::User", 		{ "foreign.user_id" => "self.parent_id" } );
__PACKAGE__->belongs_to("parent_group" 	=> "CoGeX::Result::UserGroup",  { "foreign.user_group_id" => "self.parent_id" } );
__PACKAGE__->belongs_to("child_group" 	=> "CoGeX::Result::UserGroup",  { "foreign.user_group_id" => "self.child_id" } );
__PACKAGE__->belongs_to("experiment" 	=> "CoGeX::Result::Experiment",	{ "foreign.experiment_id" => "self.child_id" } );
__PACKAGE__->belongs_to("genome"     	=> "CoGeX::Result::Genome",     { "foreign.genome_id" => "self.child_id" } );
# mdb removed 8/21/15 COGE-648
#,{
#					 join=>['genomic_sequence_type', 'organism'],
#					 prefetch=>['genomic_sequence_type', 'organism'],
#					} );
__PACKAGE__->belongs_to("feature"    	=> "CoGeX::Result::Feature",    { "foreign.feature_id" => "self.child_id" } );
__PACKAGE__->belongs_to("list" 		 	=> "CoGeX::Result::List", 		{ "foreign.list_id" => "self.child_id" } );
# mdb removed 3/17/15 due to error after upgrade to ubuntu 14.04 and perl 5.18.2:
#   DBIx::Class::ResultSet::single(): single() can not be used on resultsets collapsing a has_many. Use find( \\%cond ) or next() instead at /usr/local/lib/perl/5.18.2/CoGeX/Result/UserConnector.pm line 226
#			{
#			 join =>['child_connectors', 'list_type'],
#			 prefetch => ['child_connectors','list_type'],
#			});
__PACKAGE__->belongs_to("role" 		 	=> "CoGeX::Result::Role", "role_id" );

################################################ subroutine header begin ##

=head2 debug

 Usage     : print debug info to stderr
 Purpose   :
 Returns   :
 Argument  : None
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub debug
{
	my $self = shift;
	print STDERR 'UserConnector: id=' . $self->id . ' parent=' . $self->parent_id . ':' . $self->parent_type . ' child=' . $self->child_id . ':' . $self->child_type . ' role_id=' . $self->role_id . "\n";
}

################################################ subroutine header begin ##

=head2 is_parent_XXXXXXXX

 Usage     :
 Purpose   :
 Returns   :
 Argument  : None
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub is_parent_user
{
	return shift->parent_type() == $node_types->{user};
}

sub is_parent_group
{
	return shift->parent_type() == $node_types->{group};
}

################################################ subroutine header begin ##

=head2 is_child_XXXXXXXX

 Usage     :
 Purpose   :
 Returns   :
 Argument  : None
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub is_child_list
{
	return shift->child_type() == $node_types->{list};
}

sub is_child_genome
{
	return shift->child_type() == $node_types->{genome};
}

sub is_child_feature
{
	return shift->child_type() == $node_types->{feature};
}

sub is_child_experiment
{
	return shift->child_type() == $node_types->{experiment};
}

sub is_child_group
{
	return shift->child_type() == $node_types->{group};
}

################################################ subroutine header begin ##

=head2 parent

 Usage     :
 Purpose   :
 Returns   :
 Argument  : None
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub parent
{
	my $self = shift;

	if ($self->is_parent_user) {
		return $self->user;
	}
	elsif ($self->is_parent_group) {
		return $self->parent_group;
	}
	else {
		die;
	}

	return;
}

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

sub child
{
	my $self = shift;
	
	if ($self->is_child_experiment) {
		return $self->experiment;
	}
	elsif ($self->is_child_genome) {
		return $self->genome;
	}
	elsif ($self->is_child_feature) {
		return $self->feature;
	}
	elsif ($self->is_child_list) {
		return $self->list;
	}
	elsif ($self->is_child_group) {
		return $self->child_group;
	}
	else {
		warn "Unknown type ".$self->child_type."\n";
	}

	return;
}

1;
