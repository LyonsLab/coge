package CoGe::Accessory::blastz_report::segment;
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
    __PACKAGE__->mk_accessors(qw(query_start query_stop subject_start subject_stop length strand query_alignment subject_alignment alignment number identity));
  }

#ripped from class::Accessor
#################### subroutine header begin ####################

=head2 new

 Usage     : my $seq = new CoGe::Accessory::blastz_report::segment({
								    query_start=>
								    query_stop=>
								    subject_start=>
								    subject_stop=>
								    length=>
								    strand=>
								    query_alignment=>
								    subject_aligment=>
								    alignment=>
								    identity=>
								    number=>
								   });
 Purpose   : create an object filled with accessor functions for an alignment segment from blastz
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
    my $seg = bless {%$fields}, $class;
    return $seg;

  }

sub query_begin      {shift->query_start(@_)}
sub query_end        {shift->query_stop(@_)}
sub subject_begin    {shift->subject_start(@_)}
sub subject_end      {shift->subject_stop(@_)}
sub qbegin      {shift->query_start(@_)}
sub qend        {shift->query_stop(@_)}
sub sbegin    {shift->subject_start(@_)}
sub send      {shift->subject_stop(@_)}
sub qstart      {shift->query_start(@_)}
sub qstop        {shift->query_stop(@_)}
sub sstart    {shift->subject_start(@_)}
sub sstop      {shift->subject_stop(@_)}
sub qb      {shift->query_start(@_)}
sub qe        {shift->query_stop(@_)}
sub sb    {shift->subject_start(@_)}
sub se      {shift->subject_stop(@_)}

sub qalign  {shift->query_alignment(@_)}
sub salign  {shift->subject_alignment(@_)}
sub qseq  {shift->query_alignment(@_)}
sub sseq  {shift->subject_alignment(@_)}
sub queryAlignment  {shift->query_alignment(@_)}
sub sbjctAlignment  {shift->subject_alignment(@_)}
sub qa  {shift->query_alignment(@_)}
sub sa  {shift->subject_alignment(@_)}

sub alignmentString {shift->alignment(@_)}
sub alignment_string {shift->alignment(@_)}
sub align {shift->alignment(@_)}
1;
