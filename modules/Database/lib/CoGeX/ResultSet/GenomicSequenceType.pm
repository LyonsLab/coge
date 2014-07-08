  package CoGeX::ResultSet::GenomicSequenceType;

  use strict;
  use warnings;
  use base 'DBIx::Class::ResultSet';

 ################################################ subroutine header begin ##

=head2 resolve

 Usage     : resolve('<object or string here>')
 Purpose   :
 Returns   : A result set(?), or the GenomicSequenceType object.
 Argument  : One (no hash)
 Throws    :
 Comments  : Searches for results like (or not like?) the provided string.

See Also   :

=cut

################################################## subroutine header end ##

sub resolve {
    my $self = shift;
    my $info = shift;
    return $info if ref($info) =~ /GenomicSequenceType/;	#If $info is a reference to a GenomicSequenceType object, return it.
    return $self->find($info) if $info =~ /^\d+$/;		#If $info is all digits, query for it and return the result (some kind of reference number?)

	#Search for name with beginning (not?) like contents of $info
    my @res = $self->search({
			     'name' => { '-like' => $info . '%'},
			    }
			    ,{});

	#If previous search did not return a scalar, repeat, but search for $info anywhere in string
    @res = $self->search({
			     'name' => { '-like' => '%' . $info . '%'},
			    }
			    ,{}) unless scalar @res;

    return wantarray ? @res : \@res;	#Check context, return approprate data.
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
