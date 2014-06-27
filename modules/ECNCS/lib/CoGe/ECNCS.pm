package CoGe::ECNCS;
use strict;
use CoGe::ECNCS::DB::Algorithm;
use CoGe::ECNCS::DB::Algorithm_data;
use CoGe::ECNCS::DB::Algorithm_run;
use CoGe::ECNCS::DB::Author;
use CoGe::ECNCS::DB::Classification;
use CoGe::ECNCS::DB::Data_mask;
use CoGe::ECNCS::DB::Ecncs;
use CoGe::ECNCS::DB::Location;
use CoGe::ECNCS::DB::Spike;
use CoGe::ECNCS::DB::Status;
use Carp qw(cluck);
BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = '0.1';
    @ISA         = qw(Exporter);
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
}

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module.
## You better edit it!

=head1 NAME

CoGe::ECNCS - CoGe::ECNCS

=head1 SYNOPSIS

  use CoGe::ECNCS;
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

	HASH(0x813d9e0)
	CPAN ID: MODAUTHOR
	XYZ Corp.
	a.u.thor@a.galaxy.far.far.away
	http://a.galaxy.far.far.away/modules

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

perl(1).

=cut

#################### main pod documentation end ###################

#################### subroutine header begin ####################

=head2 sample_function

 Usage     : How to use this function/method
 Purpose   : What it does
 Returns   : What it returns
 Argument  : What it wants to know
 Throws    : Exceptions and other anomolies
 Comment   : This is a sample subroutine header.
           : It is polite to include more pod and fewer comments.

See Also   :

=cut

#################### subroutine header end ####################

sub new
{
    my ($class, %parameters) = @_;

    my $self = bless ({}, ref ($class) || $class);

    return $self;
}

sub get_algorithm_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Algorithm->new();
  }

sub get_alg_obj
  {
    my $self = shift;
    return $self->get_algorithm_obj();
  }

sub get_algorithm_data_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Algorithm_data->new();
  }

sub get_alg_data_obj
  {
    my $self = shift;
    return $self->get_algorithm_data_obj();
  }

sub get_algorithm_run_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Algorithm_run->new();
  }

sub get_alg_run_obj
  {
    my $self = shift;
    return $self->get_algorithm_run_obj();
  }

sub get_author_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Author->new();
  }

sub get_auth_obj
  {
    my $self = shift;
    return $self->get_author_obj();
  }

sub get_classification_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Classification->new();
  }

sub get_class_obj
  {
    my $self = shift;
    return $self->get_classification_obj();
  }

sub get_data_mask_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Data_mask->new();
  }

sub get_mask_obj
  {
    my $self = shift;
    return $self->get_data_mask_obj();
  }

sub get_ecncs_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Ecncs->new();
  }

sub get_cns_obj
  {
    my $self = shift;
    return $self->get_ecncs_obj();
  }

sub get_location_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Location->new();
  }

sub get_loc_obj
  {
    my $self = shift;
    return $self->get_location_obj();
  }

sub get_spike_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Spike->new();
  }

sub get_status_obj
  {
    my $self = shift;
    return CoGe::ECNCS::DB::Status->new();
  }

sub get_stat_obj
  {
    my $self = shift;
    return $self->get_status_obj();
  }

1;
# The preceding line will help the module return a true value
