  package CoGeX::ResultSet::FeatureAnnotation;

  use strict;
  use warnings;
  use base 'DBIx::Class::ResultSet';

################################################ subroutine header begin ##

=head2 esearch

 Usage     :
 Purpose   : Returns not only annotaion data, but related annotaion type and annotation type group.
 Returns   :
 Argument  :
 Throws    :
 Comments  : Extended SEARCH

See Also   :

=cut

################################################## subroutine header end ##

sub esearch
{
	my $self = shift;
	my $join = $_[1]{'join'};

	map { push(@$join, $_ ) } ('annotation_type');

	my $prefetch = $_[1]{'prefetch'};
	map { push(@$prefetch, $_ ) }
	     ('annotation_type',
	          { 'annotation_type' => 'annotation_type_group' }
	     );

	$_[1]{'join'} = $join;
	$_[1]{'prefetch'} = $prefetch;
	return $self->search( @_ );
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
