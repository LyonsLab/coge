package CoGe::Genome::DB::Sequence_type;
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
    __PACKAGE__->table('sequence_type');
    __PACKAGE__->columns(Primary=>qw{sequence_type_id});
    __PACKAGE__->columns(Others=>qw{name description});
    __PACKAGE__->has_many(sequences=>'CoGe::Genome::DB::Sequence');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature_type 

=head1 SYNOPSIS

  use CoGe::Genome::DB::Feature_type;
  my ($feat_type) = CoGe::Genome::DB::Feature_type->search_like(name=>"protein%");
  print $feat_type->name,"\n";
  print $feat_type->description,"\n";
  print $feat_type->id,"\n";
  my @seq_objs = $feat_type->sequences();
  
  


=head1 DESCRIPTION

The sequence type table in the genomes database stores the name and description 
of the type of sequence that is stored in the sequence table.  The sequence table
in turn stores a specific sequence associated with a "feature" (See the feature
module for more information on the feature table).  Since this database is designed
to store genomic sequences and associated information, the sequences stored for
a feature in the sequence table are usually translated protein sequences.  (The 
actual DNA genomic sequence is stored in the genomic_sequence table of the database.)
For that example, there would be an entry in the sequence_type table such that the 
name would be"protein" and the description would be "translated protein".  

This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.


The columns for this table are:
 sequence_type_id
 name
 description

Related objects that can be accessed through this object are:
 CoGe::Genome::DB::Sequence

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
CoGe::Genome::DB::Sequence
Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

 name             =>  name 

 description      =>  description
 desc             =>  alias for description

 sequence_type_id =>  database entry id
 id               =>  alias for sequence_type_id

 sequences        =>  returns an array of CoGe::Genome::DB::Sequence objects or 
                      a Class::DBI iterator
 seqs             =>  alias for sequences 

 new              =>  creates a new object (inherited from Class::Accessor)

=cut

sub seqs
  {
    my $self = shift;
    return $self->sequences();
  }

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->sequence_type_id();
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

