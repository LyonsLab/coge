package CoGe::ECNCS::DB::Algorithm;
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
    __PACKAGE__->table('algorithm');
    __PACKAGE__->columns(All=>qw(algorithm_id name description link));
    __PACKAGE__->has_many(algorithm_runs=>'CoGe::ECNCS::DB::Algorithm_run');
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

CoGe::ECNCS::DB::algorithm - CoGe::ECNCS::DB::algorithm

=head1 SYNOPSIS

  use CoGe::ECNCS::DB::algorithm;
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

sub algo_runs
  {
    my $self = shift;
    return $self->algorithm_runs();
  }

sub runs
  {
    my $self = shift;
    return $self->algorithm_runs();
  }

sub id
  {
    my $self = shift;
    return $self->algorithm_id();
  }

1;
# The preceding line will help the module return a true value
