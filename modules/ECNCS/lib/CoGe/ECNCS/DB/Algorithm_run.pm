package CoGe::ECNCS::DB::Algorithm_run;
use strict;
use base 'CoGe::ECNCS::DB';

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = '0.1';
    @ISA         = (@ISA, qw(Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
    __PACKAGE__->table('algorithm_run');
    __PACKAGE__->columns(All=>qw(algorithm_run_id algorithm_id algorithm_options status_id classification_id author_id note date));
    __PACKAGE__->has_a('algorithm_id'=>'CoGe::ECNCS::DB::Algorithm');
    __PACKAGE__->has_a('status_id'=>'CoGe::ECNCS::DB::Status');
    __PACKAGE__->has_a('classification_id'=>'CoGe::ECNCS::DB::Classification');
    __PACKAGE__->has_a('author_id'=>'CoGe::ECNCS::DB::Author');
    __PACKAGE__->has_many('data_masks'=>'CoGe::ECNCS::DB::Data_mask');
    __PACKAGE__->has_many('algorithm_data'=>'CoGe::ECNCS::DB::Algorithm_data');
    __PACKAGE__->has_many('ECNCS'=>'CoGe::ECNCS::DB::Ecncs');
  }

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

#################### main pod documentation begin ###################
## Below is the stub of documentation for your module.
## You better edit it!

=head1 NAME

CoGe::ECNCS::DB::algorithm_run - CoGe::ECNCS::DB::algorithm_run

=head1 SYNOPSIS

  use CoGe::ECNCS::DB::algorithm_run;
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
sub id
  {
    my $self = shift;
    return $self->algorithm_run_id();
  }

sub algorithm
  {
    my $self = shift;
    return $self->algorithm_id();
  }

sub algo
  {
    my $self = shift;
    return $self->algorithm_id();
  }

sub status
  {
    my $self = shift;
    return $self->status_id();
  }

sub classification
  {
    my $self = shift;
    return $self->classification_id();
  }

sub class
  {
    my $self = shift;
    return $self->classification_id();
  }

sub author
  {
    my $self = shift;
    return $self->author_id();
  }

sub masks
  {
    my $self = shift;
    return $self->masks();
  }

sub algo_data
  {
    my $self = shift;
    return $self->algorithm_data();
  }

sub data
  {
    my $self = shift;
    return $self->algorithm_data();
  }

sub ecncs
  {
    my $self = shift;
    return $self->ECNCS();
  }

sub cns
  {
    my $self = shift;
    return $self->ECNCS();
  }

sub results
  {
    my $self = shift;
    return $self->ECNCS();
  }

1;
# The preceding line will help the module return a true value
