package CoGeX::Result::ExperimentType;

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::ExperimentType

=cut

__PACKAGE__->table("experiment_type");

=head1 ACCESSORS

=head2 experiment_type_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0
  size: 11

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 description

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
			 "experiment_type_id",
			 { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
			 "name",
			 { data_type => "VARCHAR", is_nullable => 0, size => 255 }
			);
__PACKAGE__->set_primary_key("experiment_type_id");

__PACKAGE__->has_many('experiment_type_connectors'=>"CoGeX::Result::ExperimentTypeConnector","experiment_type_id");

################################################ subroutine header begin ##

=head2 experiments

 Usage     : $self->experiments
 Purpose   : pass through experiment_type_connector to fake a many-to-many connection with experiment
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub experiments
 {
   map {$_->experiment(@_)} shift->experiment_type_connectors();
 }

1;

=head1 AUTHORS

 Eric Lyons
 Matt Bomhoff

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
