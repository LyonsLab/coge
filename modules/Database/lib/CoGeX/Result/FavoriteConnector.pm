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
__PACKAGE__->belongs_to("experiment" 	=> "CoGeX::Result::Experiment",	{ "foreign.experiment_id" => "self.child_id" } );
__PACKAGE__->belongs_to("genome"     	=> "CoGeX::Result::Genome",     { "foreign.genome_id"     => "self.child_id" } );
__PACKAGE__->belongs_to("feature"    	=> "CoGeX::Result::Feature",    { "foreign.feature_id"    => "self.child_id" } );
__PACKAGE__->belongs_to("list" 		 	=> "CoGeX::Result::List", 		{ "foreign.list_id"       => "self.child_id" } );

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

################################################ subroutine header begin ##

=head2 is_*

 Usage     :
 Purpose   :
 Returns   :
 Argument  : None
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub is_list {
    return shift->child_type() == $node_types->{list};
}

sub is_genome {
    return shift->child_type() == $node_types->{genome};
}

sub is_feature {
    return shift->child_type() == $node_types->{feature};
}

sub is_experiment {
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

sub child {
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
        return $self->list;
    }
    else {
        warn "unknown child type " . $self->child_type;
    }

    return;
}

1;
