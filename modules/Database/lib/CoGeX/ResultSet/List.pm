package CoGeX::ResultSet::List;

use strict;
use warnings;
use base 'DBIx::Class::ResultSet';

################################################ subroutine header begin ##

=head2 public_lists

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

sub public_lists
  {
    my $self = shift;
    my %opts = @_;
    return $self->search({restricted=>0});
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
