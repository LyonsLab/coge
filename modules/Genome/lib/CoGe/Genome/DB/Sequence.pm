package CoGe::Genome::DB::Sequence;
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
    return $self->sequence_data();
  }

sub seq_data
  {
    my $self = shift;
    return $self->sequence_data();
  }

sub id
  {
    my $self = shift;
    return $self->sequence_id();
  }

1; #this line is important and will help the module return a true value

