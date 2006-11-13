package CoGe::Genome::DB::Genomic_sequence;
use strict;
use base 'CoGe::Genome::DB';
use Carp;
BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
     #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->table('genomic_sequence');
    __PACKAGE__->columns(All=>qw{genomic_sequence_id start stop chromosome sequence_data dataset_id});
    __PACKAGE__->has_a(dataset_id=>'CoGe::Genome::DB::Dataset');
    __PACKAGE__->set_sql(delete_dataset=>qq{
DELETE genomic_sequence 
  FROM genomic_sequence
 WHERE dataset_id = ?
});
    __PACKAGE__->set_sql(get_sequence=>qq{
SELECT *
  FROM genomic_sequence
 WHERE dataset_id = ?
   AND chromosome = ?
   AND (
       (start <= ? && stop >= ?)
    || (start <= ? && stop >= ?)
    || (start >  ? && stop <  ?)
   )
 ORDER BY start
});
    __PACKAGE__->set_sql(last_position=>qq{
SELECT stop
  FROM genomic_sequence
 WHERE dataset_id = ?
 ORDER BY stop DESC
 LIMIT 1
});
    __PACKAGE__->set_sql(last_chromosome_position=>qq{
SELECT stop
  FROM genomic_sequence
 WHERE dataset_id = ?
   AND chromosome = ?
 ORDER BY stop DESC
 LIMIT 1
});


  __PACKAGE__->set_sql(get_chromosomes=> qq{
SELECT DISTINCT chromosome
  FROM genomic_sequence
  WHERE dataset_id = ?
  ORDER BY  dataset_id asc
});

}
########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Genomic_sequence - Genome::DB::Genomic_sequence

=head1 SYNOPSIS

  use Genome::DB::Genomic_sequence
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

sub dataset
  {
    my $self = shift;
    return $self->dataset_id();
  }

sub data_information
  {
    my $self = shift;
    print STDERR "data_information subroutine is obselete. Please use dataset";
    return $self->dataset_id();
  }

sub information
  {
    my $self = shift;
    print STDERR "information subroutine is obselete. Please use dataset";
    return $self->dataset_id();
  }

sub data_info
  {
    my $self = shift;
    print STDERR "data_info subroutine is obselete. Please use dataset";
    return $self->dataset_id();
  }

sub info
  {
    my $self = shift;
    print STDERR "info subroutine is obselete. Please use dataset";
    return $self->dataset_id();
  }

sub data_information_id
  {
    my $self = shift;
    carp "!!!data_information_id is obselete. Please use dataset_id";
    return $self->dataset_id();
  }

sub begin
  {
    my $self = shift;
    return $self->start(@_);
  }

sub end
  {
    my $self = shift;
    return $self->stop(@_);
  }

sub chr
  {
    my $self = shift;
    return $self->chromosome(@_);
  }

sub seq
  {
    my $self = shift;
    return $self->sequence_data(@_);
  }

sub seq_data
  {
    my $self = shift;
    return $self->sequence_data();
  }

sub id
  {
    my $self = shift;
    return $self->genomic_sequence_id();
  }

sub delete_dataset
  {
    my $self = shift;
    my $id = shift;
    my $sth = $self->sql_delete_dataset;
#    print STDERR $id,"\n";
    return $sth->execute($id);

#    return $self->search_delete_data_source($id);
  }

################################################ subroutine header begin ##

=head2 get_sequence

 Usage     : $object->get_sequence(start   => $start, 
                                   stop    => $stop, 
                                   chr     => $chr,
                                   dataset => $dataset);

 Purpose   : gets the genomic sequence for the specified conditions
 Returns   : a string (containing the genomic sequence)
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             dataset => dataset object, or dataset database id, or dataset name
                        uses CoGe::Genome:DB::Dataset->resolve_dataset
             strand  => 1 or -1.  Default 1.
                        if negative strand is requested, the complement
                        of the dna seq will be returned
 Throws    : undef if no sequence is obtained
 Comments  : You must provide an info_id

See Also   : 

=cut

################################################## subroutine header end ##

sub get_sequence
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{'start'} || $opts{'START'} || $opts{begin} || $opts{BEGIN};
    $start = 1 unless $start;
    $start = 1 if $start < 1;
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    $stop = $start unless $stop;
    my $chr = $opts{chr} || $opts{CHR} || $opts{chromosome} || $opts{CHROMOSOME};
    my $ds_id = $opts{dataset} || $opts{dataset_id} || $opts{info_id} || $opts{INFO_ID} || $opts{data_info_id} || $opts{DATA_INFO_ID};
    my $coge = new CoGe::Genome;
    my $ds = $coge->get_dataset_obj->resolve_dataset($ds_id);
    my $strand = $opts{strand} || $opts{STRAND} || 1;
    my $sth = $self->sql_get_sequence();
    $sth->execute($ds->id, $chr, $start, $start, $stop, $stop, $start, $stop);
    my ($seq, $seq_start, $seq_stop); #(4, 1, 2)
    while (my $q = $sth->fetch())
      {
	$seq .= $q->[4];
	$seq_start = $q->[1] unless $seq_start;
	$seq_start = $q->[1] if $q->[1] < $seq_start;
	$seq_stop = $q->[2] unless $seq_stop;
	$seq_stop = $q->[2] if $q->[2] > $seq_stop;
      }
    $sth->finish;
    return undef unless $seq;
    my $trim_start = $start - $seq_start;
    my $trim_stop = $seq_stop - $stop;
    #trim end first
    $seq = substr($seq, 0, length($seq) - $trim_stop);
    #trip beginning
    $seq = substr($seq, $trim_start);
    if ($strand =~ "-")
      {
	$seq = reverse $seq;
	$seq =~ tr/atcgATCG/tagcTAGC/;
      }
    return $seq;
  }

################################################ subroutine header begin ##

=head2 get_last_position

 Usage     : my $last = $genome_seq_obj->get_last_position(dataset=>$dataset, chr=>$chr);
 Purpose   : gets the last genomic sequence position for a dataset given a chromosome
 Returns   : an integer that refers to the last position in the genomic sequence refered
             to by a dataset given a chromosome
 Argument  : chr=> chromosome
             dataset=>dataset object, name, or database id
                      if no dataset is specified, then it uses itself to find the related
                      dataset by $self->dataset
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_last_position
  {
    my $self = shift;
    my %opts = @_;
    my $ds = $opts{dataset} || $opts{ds} || $self->dataset_id;
    my $coge = new CoGe::Genome;
    $ds = $coge->get_dataset_obj->resolve_dataset($ds);
    my $chr = $opts{chr};
    my $id = ref ($ds) =~ /dataset/i ? $ds->id : $ds;
    my $sth = $chr ? $self->sql_last_chromosome_position() : $self->sql_last_position();
    $chr ? $sth->execute($id, $chr) : $sth->execute($id);
    my $q = $sth->fetch();
    $sth->finish;
    my $stop = $q->[0];
    return $stop;
  }
################################################ subroutine header begin ##

=head2 get_dataset_chromosome_sequence

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_dataset_chromosome_sequence
  {
    my $self = shift;
    my $ds = shift;
    my $chr = shift;
    my $id = ref ($ds) =~ /dataset/i ? $ds->id : $ds;
    my $stop = $self->get_last_position($id);
    return $self->get_sequence(start=>1, stop=>$stop, chr=>$chr, dataset_id=>$id);
  }
################################################ subroutine header begin ##

=head2 get_chromosome_for_dataset

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_chromosome_for_dataset
  {
    my $self = shift;
    my $ds = shift;
    my $id = ref ($ds) =~ /dataset/i ? $ds->id : $ds;
    my $sth = $self->sql_get_chromosomes();
    $sth->execute($id);
    my @chrs;
    while (my $row = $sth->fetchrow_arrayref())
      {
	push @chrs, $row->[0];
      }
    $sth->finish;
    return unless scalar @chrs;
    return wantarray ? @chrs : \@chrs;
  }
sub get_chromosome_for_data_information
  {
    my $self = shift;
    carp "get_chromosome_for_data_information is obselete.  Use get_chromosome_for_dataset\n";
    return $self->get_chromosome_for_dataset(@_);
  }

################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


1; #this line is important and will help the module return a true value

