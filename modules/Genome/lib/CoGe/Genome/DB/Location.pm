package CoGe::Genome::DB::Location;
use strict;
use base 'CoGe::Genome::DB';

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->table('location');
    __PACKAGE__->columns(All=>qw{location_id start stop strand chromosome feature_id});
    __PACKAGE__->has_a(feature_id=>'CoGe::Genome::DB::Feature');
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Location - Genome::DB::Location

=head1 SYNOPSIS

  #location objects are usually obtained through a feature object and contain the chromosomal 
  #location information for a feature.  A feature object may be associated with one or more
  #location objects

  use CoGe::Genome;
  #create the master genome database object
  my $db = CoGe::Genome->new;

  #using the master database object, let's find some features of a given name

  foreach my $feat ($db->get_features_by_name("GenePoo"))
   {
     #print out the name(s) of the feature.  Names are stored in a separate table and 
     #and are accessible through the feature object
     print join (", ", map {$_->name} $feat->names)

     #print out the type of the feature as well.  Each feature is of one and only one type.  
     #Featuer types are also stored in a separate table and are accessible through the feature object
     #(Three cheers for Class::DBI!)
     print " ", $feat->type->name,"\n";

     #get and print locations for the gene
     foreach my $loc ($feat->locs)
       {
         #print out something useful
         print "\t", $loc->chromosome," ", $loc->start,"-", $loc->stop," ", $loc->strand,"\n";
       }
   }


=head1 DESCRIPTION

The location table in the genomes database stores the location information associated with a feature.
Since this is a genomic location, the information for a location consists of the start position in
chromosomal units (usually nucleotides), the stop position, the strand, and the chromosome.  This 
table is related to the feature table.

This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.

The columns for this table are:
 location_id
 start
 stop
 strand
 chromosome
 feature_id

Related objects that be accessed through this object are:
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
 CoGe::Genome::DB::Feature
 Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##




=head2 Accessor Functions



new              =>  creates a new object (inherited from Class::Accessor)

location_id      =>  database entry id
id               =>  alias for location_id

start            => returns the start position of the location in chromosomal units (usually nucelotides)
begin            => alias for start

stop             => returns the stop position of the location in chromosomal units (usually nucelotides)
end              => alias for stop

strand           => returns the strand of the location Conventions are "+" or "1" for the top strand,
                    "-" or "-1" for the bottom strand.  Most checks for strandedness in various 
                    applications that need to know if a location is on the top or bottom usuall check
                    for the "-" symbol to know if it is on bottom strand.  Thus, a location of "bottom"
                    may produce some unanticipated results in other parts of the code

chromosome       => returns the chromosome name of the location
chr              => alias for chromosome

feature_id       => returns the related feature object associated with this location
feature          => alias for feature_id
feat             => alias for feature_id

=cut


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
sub id
  {
    my $self = shift;
    return $self->location_id();
  }

################################################ subroutine header begin ##

=head2 combine_overlaps

 Usage     : @locs = $loc_obj->combine_overlaps(\@locs);
 Purpose   : Combines overlapping locations.  This can happen if there are
             multiple splice varients for a given gene.  For example, there 
             are three locations (start-end) 20-200, 175-210, and 10-150.  
             10-210 will be returned.  .
 Returns   : wantarray (an array or a ref to an array) or undef
 Argument  : an array or ref to an array of location objects
 Throws    : undef
 Comments  : This does not check anything other than the start and stop
           : of the location objects.  This routine also modifies Class::DBI
             objects and thus can generate warnings.  Be aware of updating
             those objects in the database! (not recommended)

See Also   : 

=cut

################################################## subroutine header end ##


sub combine_overlaps
  {
    my $self = shift;
    my @locs;
    foreach my $item (@_)
      {
	if (ref ($item) =~ /array/i)
	  {push @locs, @$item;}
	else
	  {push @locs, $item;}
      }
    return unless @locs;
    my $orig_count = scalar @locs;
    my %skip;
    my @new_locs;
    while (my $nloc = shift @locs)
      {
	next if $skip{$nloc->id};
	foreach my $loc (@locs)
	  {
	    if ($nloc->start >= $loc->start && $nloc->start <= $loc->stop)
	      {
		$skip{$loc->id} = 1;
		$nloc->start($loc->start);
		$nloc->discard_changes();
	      }
	    if ($nloc->stop <= $loc->end && $nloc->stop >= $loc->start)
	      {
		$skip{$loc->id} = 1;
		$nloc->stop($loc->stop);
		$nloc->discard_changes();
	      }
	  }
	push @new_locs, $nloc;
      }

    if ($orig_count != scalar @new_locs)
      {

	@new_locs = $self->combine_overlaps(@new_locs);
      }
    return wantarray ? @new_locs : \@new_locs;
  }

1; #this line is important and will help the module return a true value

