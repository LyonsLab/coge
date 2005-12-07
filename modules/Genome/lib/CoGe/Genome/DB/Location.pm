package CoGe::Genome::DB::Location;
use strict;
use base 'CoGe::Genome::DB';

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->table('location');
    __PACKAGE__->columns(All=>qw{location_id start stop strand chromosome feature_id});
    __PACKAGE__->has_a(feature_id=>'CoGe::Genome::DB::Feature');
    __PACKAGE__->set_sql(delete_data_information=>qq{
DELETE location 
  FROM location
  JOIN feature using (feature_id)
 WHERE feature.data_information_id = ?
});

}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Location - Genome::DB::Location

=head1 SYNOPSIS

  use Genome::DB::Location
  blah blah blah


=head1 DESCRIPTION

Stub documentation for this module was created by ExtUtils::ModuleMaker.
It looks like the author of the extension was negligent enough
to leave the stub unedited.

Blah blah blah.


=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons
	CPAN ID: AUTHOR
	XYZ Corp.
	elyons@nature.berkeley.edu
	http://a.galaxy.far.far.away/modules

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

############################################# main pod documentation end ##


################################################ subroutine header begin ##

=head2 sample_function

 Usage     : How to use this function/method
 Purpose   : What it does
 Returns   : What it returns
 Argument  : What it wants to know
 Throws    : Exceptions and other anomolies
 Comments  : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.

See Also   : 

=cut

################################################## subroutine header end ##


sub new
{
    my ($class, %parameters) = @_;

    my $self = bless ({}, ref ($class) || $class);

    return ($self);
}

sub begin
  {
    my $self = shift;
    return $self->start();
  }

sub end
  {
    my $self = shift;
    return $self->stop();
  }

sub chr
  {
    my $self = shift;
    return $self->chromosome();
  }

sub feature
  {
    my $self = shift;
    return $self->feature_id();
  }

sub feat
  {
    my $self = shift;
    return $self->feature_id();
  }
sub id
  {
    my $self = shift;
    return $self->location_id();
  }

sub delete_data_information
  {
    my $self = shift;
    my $id = shift;
    my $sth = $self->sql_delete_data_information;
    print STDERR $id,"\n";
    return $sth->execute($id);
  }

################################################ subroutine header begin ##

=head2 combine_overlaps

 Usage     : @locs = $loc_obj->combine_overlaps(\@locs);
 Purpose   : Combines overlapping locations.  This can happen if there are
             multiple splice varients for a given gene.  For example, there 
             are three locations (start-end) 20-200, 175-210, and 10-150.  
             10-210 will be returned.  .
 Returns   : wantarray (an array or a ref to an array) or undef
 Argument  : an array or ref to an array of location objects
 Throws    : undef
 Comments  : This does not check anything other than the start and stop
           : of the location objects.  This routine also modifies Class::DBI
             objects and thus can generate warnings.  Be aware of updating
             those objects in the database! (not recommended)

See Also   : 

=cut

################################################## subroutine header end ##


sub combine_overlaps
  {
    my $self = shift;
    my @locs;
    foreach my $item (@_)
      {
	if (ref ($item) =~ /array/i)
	  {push @locs, @$item;}
	else
	  {push @locs, $item;}
      }
    return unless @locs;
    my $orig_count = scalar @locs;
    my %skip;
    my @new_locs;
    while (my $nloc = shift @locs)
      {
	next if $skip{$nloc->id};
	foreach my $loc (@locs)
	  {
	    if ($nloc->start >= $loc->start && $nloc->start <= $loc->stop)
	      {
		$skip{$loc->id} = 1;
		$nloc->start($loc->start);
		$nloc->discard_changes();
	      }
	    if ($nloc->stop <= $loc->end && $nloc->stop >= $loc->start)
	      {
		$skip{$loc->id} = 1;
		$nloc->stop($loc->stop);
		$nloc->discard_changes();
	      }
	  }
	push @new_locs, $nloc;
      }

    if ($orig_count != scalar @new_locs)
      {

	@new_locs = $self->combine_overlaps(@new_locs);
      }
    return wantarray ? @new_locs : \@new_locs;
  }

1; #this line is important and will help the module return a true value

