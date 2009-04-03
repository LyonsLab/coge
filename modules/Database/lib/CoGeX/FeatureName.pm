package CoGeX::FeatureName;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::FeatureName

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<feature_name> table in the CoGe database.

=head1 DESCRIPTION

Has columns:
C<feature_name_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<name>
Type: VARCHAR, Default: "", Nullable: no, Size: 100


C<description>
Type: VARCHAR, Default: undef, Nullable: yes, Size: 255


C<feature_id>
Type: INT, Default: 0, Nullable: no, Size: 11


C<primary_name>
Type: TINYINT, Default: 0, Nullable: no, Size: 1


Belongs to C<CoGeX::Feature> via C<feature_id>

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "ResultSetManager", "Core");
__PACKAGE__->table("feature_name");
__PACKAGE__->add_columns(
  "feature_name_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "name",
  { data_type => "VARCHAR", default_value => "", is_nullable => 0, size => 100 },
  "description",
  {
    data_type => "VARCHAR",
    default_value => undef,
    is_nullable => 1,
    size => 255,
  },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "primary_name",
  { data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 1 },
);
__PACKAGE__->set_primary_key("feature_name_id");
__PACKAGE__->belongs_to("feature" => "CoGeX::Feature", "feature_id");



################################################ subroutine header begin ##

=head2 esearch

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub esearch : ResultSet {
    my $self = shift;
    my $join = $_[1]{'join'};
    my $prefetch = $_[1]{'prefetch'};

    map { push(@$prefetch, $_ ) } 
        ({ 'feature' => ['locations','feature_type'] });


    $_[1]{'join'} = $join;
    $_[1]{'prefetch'} = $prefetch;
    my $rs = $self->search(
         @_
    );
    return $rs;

}

1;


=head1 BUGS

=head1 SUPPORT

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

=cut
