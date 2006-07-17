package CoGe::Genome::DB::Organism;
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
    __PACKAGE__->table('organism', 'organism_id');
    __PACKAGE__->columns(All=>qw{organism_id name description});
    __PACKAGE__->has_many(data_information=>'CoGe::Genome::DB::Data_information');

}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Organism - Genome::DB::Organism

=head1 SYNOPSIS

  use Genome::DB::Organism
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

sub data_info
  {
    my $self = shift;
    return $self->data_information();
  }

sub data_infos
  {
    my $self = shift;
    return $self->data_information();
  }

sub info
  {
    my $self = shift;
    return $self->data_information();
  }

sub infos
  {
    my $self = shift;
    return $self->data_information();
  }

sub datasets
  {
    my $self = shift;
    return $self->data_information();
  }

sub dataset
  {
    my $self = shift;
    return $self->data_information();
  }


sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->organism_id();
  }

################################################ subroutine header begin ##

=head2 get_features

 Usage     : my @features = $org->get_features(version=>$version, type=>"gene");
 Purpose   : Returns an array or array ref of CoGe::Genome::DB::Feature objects
             for the organism
 Returns   : an array or array ref of CoGe::Genome::DB::Feature objects
 Argument  : a hash of key->value pairs
             version => which version of data to use (if not specified this
                        routine will find the most current version of the data
                        and use that
             type    => the type of feature to return (e.g. "gene").  If this is
                        not specified it will return all features.
 Throws    : none
 Comments  : 
           : 

See Also   : CoGe::Genome::DB::Data_information->get_features_for_organism

=cut

################################################## subroutine header end ##


sub get_features
  {
    my $self = shift;
    my %opts = @_;
    my $db = new CoGe::Genome;
    return $db->get_data_info_obj->get_features_for_organism(org=>$self, %opts);
  }

1; #this line is important and will help the module return a true value

