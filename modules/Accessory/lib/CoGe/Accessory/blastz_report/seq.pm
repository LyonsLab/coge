package CoGe::Accessory::blastz_report::seq;
###############################################################################
# bl2seqReport::HSP
###############################################################################

use strict;
use base qw(Class::Accessor);

use Data::Dumper;

BEGIN
  {
    use vars qw($VERSION);
    $VERSION = "0.01";
    __PACKAGE__->mk_accessors(qw(file start stop length strand contig header));
  }

#ripped from class::Accessor
#################### subroutine header begin ####################

=head2 new

 Usage     : my $seq = new CoGe::Accessory::blastz_report::seq({
								       file=>$name,
								       start=>$start,
								       stop=>$stop,
								       length=>($stop-$start+1),
								       strand=>$strand,
								       contig=>$contig,
								      });
 Purpose   : create an object filled with accessor functions for sequence information from blastz
 Returns   : this object
 Argument  :
 Throws    :
 Comment   :
           :

See Also   :

=cut

#################### subroutine header end ####################

sub new {
    my($proto, $fields) = @_;
    my($class) = ref $proto || $proto;

    $fields = {} unless defined $fields;

    # make a copy of $fields.
    my $seq = bless {%$fields}, $class;
    return $seq;

  }

1;
