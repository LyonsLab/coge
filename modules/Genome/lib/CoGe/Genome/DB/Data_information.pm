package CoGe::Genome::DB::Data_information;
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
    __PACKAGE__->table('data_information');
    __PACKAGE__->columns(All=>qw{data_information_id data_source_id organism_id name description link version});
    __PACKAGE__->has_a(data_source_id=>'CoGe::Genome::DB::Data_source');
    __PACKAGE__->has_a('organism_id'=>'CoGe::Genome::DB::Organism');
    __PACKAGE__->has_many(features=>'CoGe::Genome::DB::Feature');
    __PACKAGE__->has_many(genomic_sequences=>'CoGe::Genome::DB::Genomic_sequence');
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Data_source - Genome::DB::Data_source

=head1 SYNOPSIS

  use Genome::DB::Data_source
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


sub data_source
  {
    my $self = shift;
    return $self->data_source_id();
  }

sub source
  {
    my $self = shift;
    return $self->data_source_id();
  }

sub feats
  {
    my $self = shift;
    return $self->features();
    }

sub genomics_seqs
  {
    my $self = shift;
    return $self->genomics_sequences();
  }

sub seqs
  {
    my $self = shift;
    return $self->genomics_sequences();
  }

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->data_information_id();
  }

sub organism
  {
    my $self = shift;
    return $self->organism_id();
  }

sub org
  {
    my $self = shift;
    return $self->organism_id();
  }

sub species
  {
    my $self = shift;
    return $self->organism_id();
  }


1; #this line is important and will help the module return a true value

