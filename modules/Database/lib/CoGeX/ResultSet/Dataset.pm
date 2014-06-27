  package CoGeX::ResultSet::Dataset;

  use strict;
  use warnings;
  use base 'DBIx::Class::ResultSet';

 ################################################ subroutine header begin ##

=head2 resolve

 Usage     :
 Purpose   :
 Returns   :
 Argument  :
 Throws    :
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub resolve {
    my $self = shift;
    my $info = shift;
    return $info if ref($info) =~ /Dataset/i;
    return $self->find($info) if $info =~ /^\d+$/;
    return $self->search({ 'name' => { '-like' => '%' . $info . '%'}},
			 ,{});
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
