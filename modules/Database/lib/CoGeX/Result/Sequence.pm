package CoGeX::Result::Sequence;

# Created by DBIx::Class::Schema::Loader v0.03009 @ 2006-12-01 18:13:38

use strict;
use warnings;

use base 'DBIx::Class';

=head1 NAME

CoGeX::Sequence

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<sequence> table in the CoGe database.
The C<sequence> table contains an actual genomic sequence, the type of sequence (in the form of a reference to the C<sequence_type> table/the CCoGeX::Result::SequenceType> object, and an ID number that relates it to a record in the C<feature> table/a CCoGeX::Result::Feature> object.

=head1 DESCRIPTION

Has columns:
C<sequence_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11
Primary identification key for table.

C<sequence_type_id>
Type: INT, Default: 0, Nullable: no, Size: 11
Reference key to a record in the C<sequence_type> table.

C<sequence_data>
Type:LONGTEXT, Default: "", Nullable: no, Size: 4294967295
Genomic sequence data.

C<feature_id>
Type: INT, Default: 0, Nullable: no, Size: 11
Reference to a record in the C<feature> table.

Relates to CCoGeX::Result::Feature> via C<feature_id>, one-to-one relationship.
Relates to CCoGeX::Result::SequenceType> via C<sequence_type_id>, one-to-one relationship.


=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->load_components("PK::Auto", "Core");
__PACKAGE__->table("sequence");
__PACKAGE__->add_columns(
  "sequence_id",
  { data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "sequence_type_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "sequence_data",
  {
    data_type => "LONGTEXT",
    default_value => "",
    is_nullable => 0,
    size => 4294967295,
  },
  "feature_id",
  { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);

__PACKAGE__->set_primary_key("sequence_id");

__PACKAGE__->belongs_to("feature" => "CoGeX::Result::Feature", 'feature_id');
__PACKAGE__->belongs_to("sequence_type" => "CoGeX::Result::SequenceType", 'sequence_type_id');



################################################ subroutine header begin ##

=head2 genomic_position

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub genomic_position
  {
    my $self = shift;
    my $pos = shift;
    $pos = 1 unless $pos; #assume they want the start?
    #need to figure out which location object the position is within.
    my $npos = $self->sequence_type->name =~ /prot/i ? $pos*3-2: $pos; #relative position in nucleotides
    foreach my $loc (sort {$a->start <=> $b->start} $self->feature->locs())
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

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

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
    my $nstart = $self->sequence_type->name =~ /prot/i ? $start*3-2: $start; #relative position in nucleotides
    my $nstop = $self->sequence_type->name =~ /prot/i ? $stop*3 : $stop; #relative position in nucleotides
    my $genome_start;
    my $genome_stop;
    my $rev = $self->feature->strand =~ /-/; #are we on the bottom strand?
    my @tlocs = sort {$a->start <=> $b->start} $self->feature->locs();
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
		my $nloc = {
			    start=>$genome_start,
			    stop=>$genome_stop,
			    strand=>$loc->strand,
			    chromosome=>$loc->chromosome,
			    chr=>$loc->chromosome,
			    };
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
		my $nloc = {
			    start=>$genome_start,
			    stop=>$genome_stop,
			    strand=>$loc->strand,
			    chromosome=>$loc->chromosome,
			    chr=>$loc->chromosome,
			    };
		push @locs, $nloc;
	      }
	    $nstart -= $locsize; 
	    $nstop -= $locsize;
	    last if ($nstop <= 0); #better stop;
	  }
      }
    return wantarray ? @locs : \@locs;
  }



1;


=head1 BUGS


=head1 SUPPORT


=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

=cut
