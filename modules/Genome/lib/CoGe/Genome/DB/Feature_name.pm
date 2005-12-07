package CoGe::Genome::DB::Feature_name;
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
    __PACKAGE__->table('feature_name');
    __PACKAGE__->columns(All=>qw{feature_name_id name description feature_id});
    __PACKAGE__->has_a(feature_id=>'CoGe::Genome::DB::Feature');
     __PACKAGE__->set_sql(delete_data_information=>qq{
DELETE feature_name 
  FROM feature_name
  JOIN feature using (feature_id)
 WHERE feature.data_information_id = ?
});

}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Feature_name - Genome::DB::Feature_name

=head1 SYNOPSIS

  use Genome::DB::Feature_name
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

sub desc
  {
    my $self = shift;
    return $self->description();
  }

sub id
  {
    my $self = shift;
    return $self->feature_name_id();
  }

sub delete_data_information
  {
    my $self = shift;
    my $id = shift;
    my $sth = $self->sql_delete_data_information;
    print STDERR $id,"\n";
    return $sth->execute($id);
  }


1; #this line is important and will help the module return a true value

