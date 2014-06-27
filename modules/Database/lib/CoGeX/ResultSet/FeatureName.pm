  package CoGeX::ResultSet::FeatureName;

  use strict;
  use warnings;
  use base 'DBIx::Class::ResultSet';

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

sub esearch {
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

=head1 AUTHORS

 Eric Lyons

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut
