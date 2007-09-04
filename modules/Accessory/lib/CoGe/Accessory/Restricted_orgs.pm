package CoGe::Accessory::Restricted_orgs;

use strict;
use CoGeX;
use Data::Dumper;
use base 'Class::Accessor';
use CGI::Carp('fatalsToBrowser');
use CGI;
use DBIxProfiler;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $coge);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw (restricted_orgs);
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
    $coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#    $cogex->storage->debugobj(new DBIxProfiler());
#    $cogex->storage->debug(1);
    __PACKAGE__->mk_accessors qw(orgs);
 }

sub restricted_orgs
  {
    my ($self, %opts) = self_or_default(@_);
    my $user = $opts{user};
    my %orgs;
    if (!$user || $user =~ /public/i)
      {
	$orgs{papaya} = 1;
	$orgs{sorghum} = 1;
      }
    else 
      {
	delete $orgs{papaya};
	delete $orgs{sorghum};
      }
    return \%orgs;
  }

sub self_or_default { #from CGI.pm
    return @_ if defined($_[0]) && (!ref($_[0])) &&($_[0] eq 'CoGe::Accessory::Restricted_orgs');
    my $Q;

    unless (defined($_[0]) && 
            (ref($_[0]) eq 'CoGe::Accessory::Restricted_orgs' || UNIVERSAL::isa($_[0],'CoGe::Accessory::Restricted_orgs')) # slightly optimized for common case
            ) {
        $Q = CoGe::Accessory::Restricted_orgs->new unless defined($Q);
        unshift(@_,$Q);
    }
    return wantarray ? @_ : $Q;
}
1;
