package CoGeX::Result::ListConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';
#use CoGeX;

=head1 NAME

CoGeX::Result::ListConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<list_connector> table in the CoGe database.
The C<list_connector> table is used to associate C<experiment,genome,feature> records with C<list> records.

=head1 DESCRIPTION

Has columns:
C<list_connector_id> (Primary Key)
Type: INT, Default: yes, Nullable: no, Size: 11
Primary identification key for table.

C<parent_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<list> table.

C<child_id>
Type: INT, Default: "", Nullable: no, Size: 11
Key for identifying the record in the C<experiment,genome,feature> table.

C<child_type>
Type: TINYINT, Default: "", Nullable: no, Size: 1
Child type indicator.

Belongs to CCoGeX::Result::Genome> via C<child_id>
Belongs to CCoGeX::Result::Experiment> via C<child_id>
Belongs to CCoGeX::Result::Feature> via C<child_id>
Belongs to CCoGeX::Result::List> via C<child_id>   -- for a list of lists!
Belongs to CCoGeX::Result::List> via C<list_id>    -- parent list

=head1 USAGE

 use CoGeX;

=head1 METHODS

=cut

my $node_types = CoGeX::node_types();

__PACKAGE__->table("list_connector");
__PACKAGE__->add_columns(
	"list_connector_id",
	{ data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
	"parent_id",
	{ data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
	"child_id",
	{ data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
	"child_type",
	{ data_type => "TINYINT", default_value => "", is_nullable => 0, size => 1 },
);
__PACKAGE__->set_primary_key("list_connector_id");

__PACKAGE__->belongs_to("experiment" => "CoGeX::Result::Experiment", "child_id");
__PACKAGE__->belongs_to("genome"     => "CoGeX::Result::Genome",     "child_id", {
					 join=>['genomic_sequence_type', 'organism'],
					 prefetch=>['genomic_sequence_type', 'organism'],
					});
__PACKAGE__->belongs_to("feature"    => "CoGeX::Result::Feature",    "child_id");

__PACKAGE__->belongs_to("child_list"   => "CoGeX::Result::List",  {'foreign.list_id' => 'self.child_id'}); # a list of lists
__PACKAGE__->belongs_to("parent_list"  => "CoGeX::Result::List",  {'foreign.list_id' => 'self.parent_id' } ); # parent list of a genome/experiment/feature/list

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

sub is_list
{
	return shift->child_type() == $node_types->{list};
}

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

sub is_genome
{
	return shift->child_type() == $node_types->{genome};
}

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

sub is_feature
{
	return shift->child_type() == $node_types->{feature};
}

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

sub is_experiment
{
	return shift->child_type() == $node_types->{experiment};
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

	if ($self->is_experiment) {
		return $self->experiment;
	}
	elsif ($self->is_genome) {
		return $self->genome;
	}
	elsif ($self->is_feature) {
		return $self->feature;
	}
	elsif ($self->is_list) {
		return $self->child_list;
	}
	else {
		warn "unknown child type " . $self->child_type;
	}

	return;
}

1;

=head1 AUTHORS

 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
