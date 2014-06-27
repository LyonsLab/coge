package CoGe::Accessory::parse_report::HSP;
###############################################################################
# parse_report::HSP
###############################################################################

use strict;
use base qw(Class::Accessor);

use Data::Dumper;

BEGIN
  {
    use vars qw($VERSION);
    $VERSION = "0.01";
    __PACKAGE__->mk_accessors('score', 'bits', 'percent_id', 'percent_sim',
    	'match', 'positive', 'length', 'pval', 'query_start', 'query_stop',
    	'subject_start', 'subject_stop', 'query_alignment', 'subject_alignment',
    	'alignment', 'query_gaps', 'subject_gaps', 'strand', 'number',
    	'qpercent_id', 'spercent_id', 'segments', 'contains_spike',
    	'query_name', 'subject_name', 'query_length', 'subject_length',
    	'query_coverage', 'quality', 'link_as_gene');
  }

#ripped from class::Accessor
sub new {
    my($proto, $fields) = @_;
    my($class) = ref $proto || $proto;

    $fields = {} unless defined $fields;

    # make a copy of $fields.
    my $hsp = bless {%$fields}, $class;
    unless ($hsp->percent_id())
     {$hsp->percent_id(int(1000 * $hsp->match/$hsp->length)/10) if $hsp->match && $hsp->length;}
    unless ($hsp->percent_sim())
     {$hsp->percent_sim(int(1000 * $hsp->positive/$hsp->length)/10) if $hsp->positive && $hsp->length;}
    $hsp->pval("1".$hsp->pval) if $hsp->pval && $hsp->pval =~ /^e/;
    $hsp->segments([]) unless $hsp->segments;
    return $hsp;

  }

sub P               {shift->pval(@_)}
sub p               {shift->pval(@_)}
sub P_val           {shift->pval(@_)}
sub p_val           {shift->pval(@_)}
sub eval            {shift->pval(@_)}

sub query      {shift->query_name(@_)}
sub subject     {shift->subject_name(@_)}
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

sub qpercent   {shift->qpercent_id(@_)}
sub spercent   {shift->spercent_id(@_)}

sub alignmentString {shift->alignment(@_)}
sub alignment_string {shift->alignment(@_)}
sub align {shift->alignment(@_)}

sub qgaps       {shift->query_gaps(@_)}
sub sgaps       {shift->subject_gaps(@_)}
sub qgap       {shift->query_gaps(@_)}
sub sgap       {shift->subject_gaps(@_)}
sub qg       {shift->query_gaps(@_)}
sub sg       {shift->subject_gaps(@_)}

1;
