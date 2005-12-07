package CoGe::Genome::DB::Annotation_type;
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
    __PACKAGE__->table('annotation_type');
    __PACKAGE__->columns(All=>qw{annotation_type_id name description annotation_type_group_id});
    __PACKAGE__->has_a(annotation_type_group_id=>'CoGe::Genome::DB::Annotation_type_group');
    __PACKAGE__->has_many(annotations=>'CoGe::Genome::DB::Annotation');
 
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Annotation_type - Genome::DB::Annotation_type

=head1 SYNOPSIS

  use Genome::DB::Annotation_type
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

sub annotation_type_group
  {
    my $self = shift;
    return $self->annotation_type_group_id();
  }

sub anno_type_group
  {
    my $self = shift;
    return $self->annotation_type_group_id();
  }
  
sub type_group
  {
    my $self = shift;
    return $self->annotation_type_group_id();
  }
  
sub group
  {
    my $self = shift;
    return $self->annotation_type_group_id();
  }
  
sub annos
  {
    my $self = shift;
    return $self->annotations();
  }

sub desc
  {
    my $self = shift;
    return $self->description();
  }

sub id
  {
    my $self = shift;
    return $self->annotation_type_id();
  }

1; #this line is important and will help the module return a true value

