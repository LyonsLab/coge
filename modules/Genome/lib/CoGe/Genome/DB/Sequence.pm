package CoGe::Genome::DB::Sequence;
use strict;
use base 'CoGe::Genome::DB';
use CoGe::Genome::DB::Location;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->table('sequence');
    __PACKAGE__->columns(Primary=>qw{sequence_id});
    __PACKAGE__->columns(Others=>qw{sequence_type_id sequence_data feature_id});
    __PACKAGE__->has_a(sequence_type_id=>'CoGe::Genome::DB::Sequence_type');
    __PACKAGE__->has_a(feature_id=>'CoGe::Genome::DB::Feature');

}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Sequence

=head1 SYNOPSIS

  use CoGe::Genome::DB::Sequence
  my $seq = CoGe::Genome::DB::Sequence->search(feature_id=>$feat_id);
  my $seq_type = $seq->seq_type();
  my $feat = $seq->feat();
  my ($feat_name) = $feat->name();
  print ">" . $feat_name->name . " " . $feat_name->desc . " " . $seq_type->name . "\n";
  print $seq->seq . "\n";


=head1 DESCRIPTION

The sequence table in the genomes database stores a sequence associated with a feature.
Since this database is designed to store genomic sequences and associated information, 
the sequences stored for a feature in the sequence table are usually translated protein 
sequences. (The actual DNA genomic sequence is stored in the genomic_sequence table of 
the database.)  This table is related to the feature and sequence_type tables.

This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.

The columns for this table are:
 sequence_id
 sequence_type_id
 sequence_data
 feature_id

Related objects that can be accessed through this object are:
 CoGe::Genome::DB::Sequence_type
 CoGe::Genome::DB::Feature



=head1 USAGE



=head1 BUGS



=head1 SUPPORT



=head1 AUTHOR

	Eric Lyons
	elyons@nature.berkeley.edu

=head1 COPYRIGHT

This program is free software licensed under the...

	The Artistic License

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

CoGe::Genome
CoGe::Genome::DB
CoGe::Genome::DB::Sequence_type
CoGe::Genome::DB::Feature
Class::DBI

perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

 sequence_id      =>  database entry id
 id               =>  alias for sequence_id

 sequence_type_id =>  returns the related CoGe::Genome::DB::Sequence_type object
 sequence_type    =>  alias for sequence_type_id
 seq_type         =>  alias for sequence_type_id
 type             =>  alias for sequence_type_id

 featrue_id       =>  returns the related CoGe::Genome::DB::Feature object
 feature          =>  alias for feature_id
 feat             =>  alias for feature_id

 sequence_data    =>  sequence data stored in the object
 seq_data         =>  alias for sequence_data
 seq              =>  alias for sequence_data

 new              =>  creates a new object (inherited from Class::Accessor)

=cut


sub sequence_type
  {
    my $self = shift;
    return $self->sequence_type_id();
  }

sub seq_type
  {
    my $self = shift;
    return $self->sequence_type_id();
  }

sub type
  {
    my $self = shift;
    return $self->sequence_type_id();
  }

sub feature
  {
    my $self = shift;
    return $self->feature_id();
  }

sub feat
  {
    my $self = shift;
    return $self->feature_id();
  }

sub seq
  {
    my $self = shift;
    return $self->sequence_data(@_);
  }

sub seq_data
  {
    my $self = shift;
    return $self->sequence_data(@_);
  }

sub id
  {
    my $self = shift;
    return $self->sequence_id();
  }

################################################ subroutine header begin ##

=head2 get_genomic_position

 Usage     : my $pos = $sequence_obj->get_genomic_position(3);
 Purpose   : convert the relative position in a sequence to the appropriate genomic location
 Returns   : an integer that represents a chromosomal position in chromosomal units (usually nucleotides)
 Argument  : an integer that represents the relative position in the sequence contained in the sequence object
 Throws    : returns 0 if there are no locations or the relative position was outside the range of the
           : assoicated feature.
 Comments  : This can be a tricky operation, especially for protein sequences because
           : of the need to convert amino acids into nucleotides and account for 
           : multiple exons.  This routine does check to see if the squence is a
           : protein through the related sequence_type object in order to figure
           : this out.

See Also   : 

=cut

################################################## subroutine header end ##

sub get_genomic_position
  {
    my $self = shift;
#    my %opts = @_;
    my $pos = shift;
    $pos = 1 unless $pos; #assume they want the start?
    #need to figure out which location object the position is within.
    my $npos = $self->type->name =~ /prot/i ? $pos*3-2: $pos; #relative position in nucleotides
    foreach my $loc (sort {$a->start <=> $b->start} $self->feat->locs())
      {
	my $locsize = $loc->stop - $loc->start +1;
	if ( $locsize > $npos) #we have fount the location object in which the position lies
	  {
	    return $loc->start+$npos-1;
	  }
	else
	  {
	    $npos -= $locsize;
	  }
      }
    return 0;
  }

################################################ subroutine header begin ##

=head2 get_genomic_locations

 Usage     : my @locs = $sequence_obj->get_genomic_locations(start=>10, stop=>100);
 Purpose   : This routine is primarily designed to make it easy to go from relative
           : protein coordinates to genomic coordinates especially when the protein
           : location may bridge multiple exons.  
           
 Returns   : CoGe::Genome::DB::Location objects (either a array ref or ref detected by wantarray)
 Argument  : start|begin=> relative start position in the sequence object 
             stop |end  => relative stop position in the sequence object
 Throws    : returns 0 in case of error
 Comments  : 

See Also   : CoGe::Genome::DB::Location

=cut

################################################## subroutine header end ##

sub get_genomic_locations
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{start} || $opts{begin} || 1;
    my $stop = $opts{stop} || $opts{end};
    unless ($start && $stop)
      {
	warn "Need a valid start and stop position\n";
	return 0;
      }
    my @locs;
    #are we dealing with a protein sequence for locaiton conversion?
    my $nstart = $self->type->name =~ /prot/i ? $start*3-2: $start; #relative position in nucleotides
    my $nstop = $self->type->name =~ /prot/i ? $stop*3 : $stop; #relative position in nucleotides
    my $genome_start;
    my $genome_stop;
    my $rev = $self->feat->strand =~ /-/; #are we on the bottom strand?
    my @tlocs = sort {$a->start <=> $b->start} $self->feat->locs();
    if ($rev)
      {
	@tlocs = reverse @tlocs if $rev;
	foreach my $loc (@tlocs)
	  {
	    my $locsize = $loc->stop - $loc->start +1;
	    #need to detect actual start
	    if ( $locsize > $nstart && !$genome_stop) #we have found the location object in which the start lies
	      {
		$genome_stop = $loc->stop-$nstart+1;
	      }
	    elsif ($genome_stop) #we have a start, but the end has yet to be found.  Assign the start to the start of the current loc object
	      {
		$genome_stop = $loc->stop;
	      }
	    #once start is found, need to know if also have the stop in the exon, or if lies in a future exon
	    if ( $locsize > $nstop && $genome_stop) #get have the real stop
	      {
		$genome_start = $loc->stop-$nstop+1;
	      }
	    elsif ( $genome_stop) #we have the start, but the end is in another exon
	      {
		$genome_start = $loc->start; 
	      }
	    if ($genome_start && $genome_stop)
	      {
		my $nloc = CoGe::Genome::DB::Location->new();
		$nloc->start($genome_start);
		$nloc->stop($genome_stop);
		$nloc->strand($loc->strand);
		$nloc->chromosome($loc->chr);
		push @locs, $nloc;
	      }
	    $nstart -= $locsize; 
	    $nstop -= $locsize;
	    last if ($nstop <= 0); #better stop;
	  }
      }

    else
      {
	foreach my $loc (@tlocs)
	  {
	    my $locsize = $loc->stop - $loc->start +1;
	    #need to detect actual start
	    if ( $locsize > $nstart && !$genome_start) #we have found the location object in which the start lies
	      {
		$genome_start = $loc->start+$nstart-1;
	      }
	    elsif ($genome_start) #we have a start, but the end has yet to be found.  Assign the start to the start of the current loc object
	      {
		$genome_start = $loc->start;
	      }
	    #once start is found, need to know if also have the stop in the exon, or if lies in a future exon
	    if ( $locsize > $nstop && $genome_start) #get have the real stop
	      {
		$genome_stop = $loc->start+$nstop-1;
	      }
	    elsif ( $genome_start) #we have the start, but the end is in another exon
	      {
		$genome_stop = $loc->stop; 
	      }
	    if ($genome_start && $genome_stop)
	      {
		my $nloc = CoGe::Genome::DB::Location->new();
		$nloc->start($genome_start);
		$nloc->stop($genome_stop);
		$nloc->strand($loc->strand);
		$nloc->chromosome($loc->chr);
		push @locs, $nloc;
	      }
	    $nstart -= $locsize; 
	    $nstop -= $locsize;
	    last if ($nstop <= 0); #better stop;
	  }
      }
    return wantarray ? @locs : \@locs;
  }


1; #this line is important and will help the module return a true value

