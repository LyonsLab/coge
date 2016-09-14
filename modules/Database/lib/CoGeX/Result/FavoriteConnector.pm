package CoGeX::Result::FavoriteConnector;

use strict;
use warnings;

use base 'DBIx::Class::Core';

use Data::Dumper;

=head1 NAME

CoGeX::Result::FavoriteConnector

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<favorite_connector> table in the CoGe database.
The C<favorite_connector> table is used to associate C<list,experiment,genome,feature> records with C<user> records.

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

__PACKAGE__->table("favorite_connector");
__PACKAGE__->add_columns(
	"favorite_connector_id",
	{ data_type => "INT", default_value => 1, is_nullable => 0, size => 11 },
	"user_id",
	{ data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
	"child_id",
	{ data_type => "INT", default_value => "", is_nullable => 0, size => 11 },
	"child_type",
	{ data_type => "TINYINT", default_value => "", is_nullable => 0, size => 1 },
);
__PACKAGE__->set_primary_key("favorite_connector_id");
__PACKAGE__->belongs_to("user"			=> "CoGeX::Result::User", 		"user_id" );
__PACKAGE__->belongs_to("experiment" 	=> "CoGeX::Result::Experiment",	{ "foreign.experiment_id" => "self.child_id" }, { where => { child_type  => $node_types->{experiment} } } );
__PACKAGE__->belongs_to("genome"     	=> "CoGeX::Result::Genome",     { "foreign.genome_id"     => "self.child_id" }, { where => { child_type  => $node_types->{genome} } } );
__PACKAGE__->belongs_to("feature"    	=> "CoGeX::Result::Feature",    { "foreign.feature_id"    => "self.child_id" }, { where => { child_type  => $node_types->{feature} } } );
__PACKAGE__->belongs_to("list" 		 	=> "CoGeX::Result::List", 		{ "foreign.list_id"       => "self.child_id" }, { where => { child_type  => $node_types->{notebook} } } );

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

sub debug {
	my $self = shift;
	# TODO
}

sub toggle {
    my $self = shift;
        
}

1;
