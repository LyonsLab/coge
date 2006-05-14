package CoGe::Genome;
use strict;
use CoGe::Genome::DB::Annotation;
use CoGe::Genome::DB::Annotation_type;
use CoGe::Genome::DB::Annotation_type_group;
use CoGe::Genome::DB::Data_source;
use CoGe::Genome::DB::Data_information;
use CoGe::Genome::DB::Feature;
use CoGe::Genome::DB::Feature_name;
use CoGe::Genome::DB::Feature_type;
use CoGe::Genome::DB::Genomic_sequence;
use CoGe::Genome::DB::Location;
use CoGe::Genome::DB::Organism;
use CoGe::Genome::DB::Sequence;
use CoGe::Genome::DB::Sequence_type;
use Carp qw(cluck);
use Data::Dumper;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();

}


########################################### main pod documentation begin ##



=head1 NAME

Genome - Genome

=head1 SYNOPSIS

  use Genome
  my $genome_obj = Genome->new();


=head1 DESCRIPTION

  Master object for accessing information in the Genome Database from the SBCS
  Group at UC Berkeley


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

perl(1).

=cut

############################################# main pod documentation end ##


################################################ subroutine header begin ##

=head2 new

 Usage     : my $genome_obj = CoGe::Genome->new();
 Purpose   : Creates an CoGe::Genome object
 Returns   : an CoGe::Genome object
 Argument  : none
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub new
{
    my ($class, %parameters) = @_;

    my $self = bless ({}, ref ($class) || $class);

    return ($self);
}

################################################ subroutine header begin ##

=head2 get_annotation_obj

 Usage     : $genome_obj->get_annotation_obj();
 Purpose   : gets an annotation object
 Returns   : an annotation object
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Annotation

=cut

################################################## subroutine header end ##



sub get_annotation_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Annotation->new();
  }

################################################ subroutine header begin ##

=head2 get_anno_obj

 Usage     : $genome_obj->get_anno_obj();
 Purpose   : alias for get_annotation_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_anno_obj
  {
    my $self = shift;
    return $self->get_annotation_obj();
  }

################################################ subroutine header begin ##

=head2 get_annotation_type_obj

 Usage     : $genome_obj->get_annotation_type_obj();
 Purpose   : gets an annotation_type object
 Returns   : an annotation_type object
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Annotation_type

=cut

################################################## subroutine header end ##

sub get_annotation_type_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Annotation_type->new();
  }

################################################ subroutine header begin ##

=head2 get_anno_type_obj

 Usage     : $genome_obj->get_anno_type_obj();
 Purpose   : alias for get_annotation_type_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_anno_type_obj
  {
    my $self = shift;
    return $self->get_annotation_type_obj();
  }

################################################ subroutine header begin ##

=head2 get_annotation_type_group_obj

 Usage     : $genome_obj->get_annotation_type_group_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Annotation_type_group

=cut

################################################## subroutine header end ##

sub get_annotation_type_group_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Annotation_type_group->new();
  }

################################################ subroutine header begin ##

=head2 get_anno_type_obj

 Usage     : $genome_obj->get_anno_type_obj();
 Purpose   : alias
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_anno_type_group_obj
  {
    my $self = shift;
    return $self->get_annotation_type_group_obj();
  }

################################################ subroutine header begin ##

=head2 get_data_source_obj

 Usage     : $genome_obj->get_data_source_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Data_source

=cut

################################################## subroutine header end ##

sub get_data_source_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Data_source->new();
  }

################################################ subroutine header begin ##

=head2 get_source_obj

 Usage     : $genome_obj->get_source_obj();
 Purpose   : alias for get_data_source_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_source_obj
  {
    my $self = shift;
    return $self->get_data_source_obj();
  }

################################################ subroutine header begin ##

=head2 get_data_information_obj

 Usage     : $genome_obj->get_data_information_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Data_information

=cut

################################################## subroutine header end ##

sub get_data_information_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Data_information->new();
  }

################################################ subroutine header begin ##

=head2 get_information_obj

 Usage     : $genome_obj->get_information_obj();
 Purpose   : alias for get_data_information_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_information_obj
  {
    my $self = shift;
    return $self->get_data_information_obj();
  }

################################################ subroutine header begin ##

=head2 get_data_info_obj

 Usage     : $genome_obj->get_data_info_obj();
 Purpose   : alias for get_data_information_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_data_info_obj
  {
    my $self = shift;
    return $self->get_data_information_obj();
  }

################################################ subroutine header begin ##

=head2 get_feature_obj

 Usage     : $genome_obj->get_feature_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##

sub get_feature_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature->new();
  }

################################################ subroutine header begin ##

=head2 get_feat_obj

 Usage     : $genome_obj->get_feat_obj();
 Purpose   : alias for get_feature_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feat_obj
  {
    my $self = shift;
    return $self->get_feature_obj();
  }

################################################ subroutine header begin ##

=head2 get_feature_name_obj

 Usage     : $genome_obj->get_feature_name_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_name

=cut

################################################## subroutine header end ##

sub get_feature_name_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_name->new();
  }

################################################ subroutine header begin ##

=head2 get_feat_name

 Usage     : $genome_obj->get_feat_name_obj();
 Purpose   : alias for get_feature_name_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feat_name_obj
  {
    my $self = shift;
    return $self->get_feature_name_obj();
  }

################################################ subroutine header begin ##

=head2 get_feature_type_obj

 Usage     : $genome_obj->get_feature_type_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature_type

=cut

################################################## subroutine header end ##

sub get_feature_type_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Feature_type->new();
  }

################################################ subroutine header begin ##

=head2 get_feat_type_obj

 Usage     : $genome_obj->get_feat_type_obj();
 Purpose   : alias for get_feature_type_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feat_type_obj
  {
    my $self = shift;
    return $self->get_feature_type_obj();
  }

################################################ subroutine header begin ##

=head2 get_genomic_sequence_obj

 Usage     : $genome_obj->get_genomic_sequence_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Genomic_sequence

=cut

################################################## subroutine header end ##

sub get_genomic_sequence_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Genomic_sequence->new();
  }

################################################ subroutine header begin ##

=head2 get_genomic_seq_obj

 Usage     : $genome_obj->get_genomic_seq_obj();
 Purpose   : aslias for get_genomic_sequence_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_genomic_seq_obj
  {
    my $self = shift;
    return $self->get_genomic_sequence_obj();
  }

################################################ subroutine header begin ##

=head2 get_location_obj

 Usage     : $genome_obj->get_location_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Location

=cut

################################################## subroutine header end ##

sub get_location_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Location->new();
  }

################################################ subroutine header begin ##

=head2 get_loc_obj

 Usage     : $genome_obj->get_loc_obj();
 Purpose   : alias for get_location_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_loc_obj
  {
    my $self = shift;
    return $self->get_location_obj();
  }




=head2 get_organism_obj

 Usage     : get_organism_obj
 Purpose   : gets an organism object
 Returns   : an organism object
 Argument  : none
 Throws    : none
 Comments  : none

See Also   : CoGe::Genome::DB::Organism

=cut


sub get_organism_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Organism->new();
  }

=head2 get_org_obj

 Usage     : get_org_obj
 Purpose   : alias for get_organism_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut


sub get_org_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Organism->new();
  }

=head2 get_species_obj

 Usage     : get_species_obj
 Purpose   : alias for get_organism_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut



sub get_species_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Organism->new();
  }

################################################ subroutine header begin ##

=head2 

 Usage     : $genome_obj->sequence();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Sequence

=cut

################################################## subroutine header end ##

sub get_sequence_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Sequence->new();
  }


################################################ subroutine header begin ##

=head2 get_seq_obj

 Usage     : $genome_obj->get_seq_obj();
 Purpose   : alias for get_sequence_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_seq_obj
  {
    my $self = shift;
    return $self->get_sequence_obj();
  }


################################################ subroutine header begin ##

=head2 get_sequence_type_obj

 Usage     : $genome_obj->get_sequence_type_obj();
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Sequence_type

=cut

################################################## subroutine header end ##

sub get_sequence_type_obj
  {
    my $self = shift;
    return CoGe::Genome::DB::Sequence_type->new();
  }


################################################ subroutine header begin ##

=head2 get_seq_type_obj

 Usage     : $genome_obj->get_seq_type_obj();
 Purpose   : alias for get_sequence_type_obj
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_seq_type_obj
  {
    my $self = shift;
    return $self->get_sequence_type_obj();
  }


################################################ subroutine header begin ##

=head2 get_feature_by_name_and_info_id

 Usage     : my @feats = $genome_obj->get_feature_by_name_and_info_id(name=>$name, id=>$id)
 Purpose   : get feature objects based on a name and a data_information id
 Returns   : wantarray of CoGe::Genome::Feature objects
 Argument  : a hash of key/value pairs
               name=> a string that is the name of a feature to be fetched by $self->get_feature_by_name
               id=>   an integer that is the database id from the data_information table
                      the id can be obtained from a data_information object
 Throws    : none
 Comments  : 

See Also   : sub get_feature_by_name
             CoGe::Data_information for information on that object

=cut

################################################## subroutine header end ##

sub get_feature_by_name_and_info_id
  {
    my $self = shift;
    my %opts = @_;
    my ($name) = $opts{name} || $opts{NAME} || shift;
    my ($id) = $opts{info_id} || $opts{id} || $opts{ID} || shift;
    my @feats;
    foreach my $feat ($self->get_feature_by_name($name))
      {
	push @feats, $feat if $feat->data_information->id eq $id;
      }
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 get_features_by_name_and_version

 Usage     : my @feats = $genome_obj->get_features_by_name_and_version(name=>$name, version=>$version)
 Purpose   : get feature objects based on a name and a data_information version number
 Returns   : wantarray of CoGe::Genome::Feature objects
 Argument  : a hash of key/value pairs
               name=> a string that is the name of a feature to be fetched by $self->get_feature_by_name
               version=>   an integer.  This is obtained from a data_information object.
 Throws    : none
 Comments  : 

See Also   : sub get_feature_by_name
             CoGe::Data_information for information on that object

=cut

################################################## subroutine header end ##

sub get_features_by_name_and_version
  {
    my $self = shift;
    my @feats = $self->get_feature_obj->get_features_by_name_and_version(@_);
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 get_feature_by_name

 Usage     : $genome_obj->get_feature_by_name($name);
 Purpose   : get feature objects using the name of the feature
 Returns   : an array of feature objects (or nothing)
 Argument  : string / the name of the feature(s)
 Throws    : none
 Comments  : same as:
             map {$_->feature} $self->get_feat_name_obj->search(name=>$name)

See Also   : 

=cut

################################################## subroutine header end ##


sub get_feature_by_name
  {
    my $self = shift;
    my $name = shift;
    my @feats = map {$_->feature} $self->get_feat_name_obj->search(name=>$name);
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 get_features_by_name

 Usage     : $genome_obj->get_features_by_name($name);
 Purpose   : alias for get_feature_by_name
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_features_by_name
  {
    my $self = shift;
    my $name = shift;
    return $self->get_feature_by_name($name);
  }

################################################ subroutine header begin ##

=head2 get_feat_by_name

 Usage     : $genome_obj->get_feat_by_name($name);
 Purpose   : alias for get_feature_by_name
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feat_by_name
  {
    my $self = shift;
    my $name = shift;
    return $self->get_feature_by_name($name);
  }

################################################ subroutine header begin ##

=head2 get_feats_by_name

 Usage     : $genome_obj->get_feats_by_name($name);
 Purpose   : alias for get_feature_by_name
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feats_by_name
  {
    my $self = shift;
    my $name = shift;
    return $self->get_feature_by_name($name);
  }


################################################ subroutine header begin ##

=head2 get_protein_seq_by_feat_name

 Usage     : my @protein_seqs = $genome_obj->get_protein_seq_by_feat_name(name=>$name, version=>$version);
 Purpose   : Some features have an associated protein sequence which is stored
             in the sequence table with a sequence type "protein".  This routine 
             finds all features by a specified name and limited by a data_information 
             version and returns the first sequence for that feature it finds of type
             protein.
 Returns   : array of strings (protein sequence)
 Argument  : name => the name you wish to searhc (e.g. At3g010110)
             version => (optional) skips any feature whose version (from data_information)
                        does not equal this specified version.  (e.g. 6)
             hash   => Returns a hash (or hash ref) with keys as name of the feature that
                       matches the submitted name so that version numbers can be used
                       e.g. At3g010110.1=>MAR.. .  
             search_like => flag for "search like".  adds "%" after the accn name for the search
 Throws    : undef
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub get_protein_seq_by_feat_name
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name};
    my $version = $opts{version} || $opts{ver};
    my $hash = $opts{hash};
    my $search_like = $opts{search_like};
    return unless $name;
    my %seqs;
    my @seqs;
    my @no = $search_like ? $self->get_feature_name_obj->search_like(name=>$name."%") : $self->get_feature_name_obj->search(name=>$name);
    
    foreach my $no (@no)
      {
	if ($version)
	  {
	    next unless $no->feat->data_info->version eq $version;
	  }
	
	foreach my $seq ($no->feat->sequences)
	  {
	    push @seqs, $seq->sequence_data if $seq->seq_type->name =~ /prot/i;
	    $seqs{$no->name}=$seq->sequence_data if $seq->seq_type->name =~ /prot/i;
	  }
      }
    if ($hash)
      {
	return wantarray ? %seqs : \%seqs;
      }
    if (@seqs == 1)
      {
	return $seqs[0];
      }
    my @sorted = sort {length($b) <=> length ($a)} @seqs;
    return wantarray ? @sorted : \@sorted;
  }

################################################ subroutine header begin ##

=head2 get_prot_seq_by_feat_name

 Usage     : my @protein_seqs = $genome_obj->get_prot_seq_by_feat_name(name=>$name, version=>$version);
 Purpose   : Alias for get_protein_seq_by_feat_name

See Also   : get_protein_seq_by_feat_name

=cut

################################################## subroutine header end ##

sub get_prot_seq_by_feat_name
  {
    my $self = shift;
    return $self->get_protein_seq_by_feat_name(@_);
  }

################################################ subroutine header begin ##

=head2 get_protein_seq_by_feat

 Usage     : my @protein_seqs = $genome_obj->get_protein_seq_by_feat($feat);
 Purpose   : Some features have an associated protein sequence which is stored
             in the sequence table with a sequence type "protein" (e.g. CDS sequences).  
             This routine finds and returns those sequences.
 Returns   : array of strings (protein sequence) or array ref (based on wantarray)
 Argument  : GoGe::Genomes::DB::Feature object
 Throws    : undef if 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##



sub get_protein_seq_by_feat
  {
    my $self = shift;
    my $feat = shift;
    return unless $feat;
    my @seqs;
    foreach my $seq ($feat->sequences)
      {
	push @seqs, $seq->sequence_data if $seq->seq_type->name =~ /prot/i;
      }
    my @sorted = sort {length($b) <=> length ($a)} @seqs;
    return wantarray ? @sorted : \@sorted;
  }

################################################ subroutine header begin ##

=head2 get_genomic_seq_by_feat_name_and_type_name

 Usage     : my @seqs = $genome_obj->get_genomic_seq_by_feat_name_and_type_name(name=>$name, 
                                                                                type=>$type, 
                                                                                version=>$version);
 Purpose   : Gets the genomic sequence of a feature by the name of the feature and 
             by the type of the feature
 Returns   : array of strings (genomic DNA sequence)
 Argument  : name => the name you wish to searhc (e.g. At3g010110)
             type => the name of the feature type (e.g. CDS)
             version => (optional) skips any feature whose version (from data_information)
                        does not equal this specified version.  (e.g. 6)
             infoid  => (optional) skips any feature whose data information id 
                        does not match
 Throws    : undef, prints error to STDERR
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_genomic_seq_by_feat_name_and_type_name
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name};
    my $version = $opts{version} || $opts{ver};
    my $type = $opts{type};
    my $info_id = $opts{infoid};
    my %locs;
    foreach my $no ($self->get_feature_name_obj->search(name=>$name) )
      {
	foreach my $feat ($no->feat)
	  {
	    if ($version)
	      {
		next unless $feat->data_info->version eq $version;
	      }
	    if ($info_id)
	      {
		next unless $feat->data_info->id eq $info_id;
	      }
	    next unless $feat->feat_type->name =~ /$type/;
	    push @{$locs{$feat->id}}, $feat->locations;
	  }
      }
#    @locs = $self->get_location_obj->combine_overlaps(@locs);
    my $gso = $self->get_genomic_seq_obj;
    my $rev_seq_flag = 0;
    my @seqs;
    foreach my $id (keys %locs)
      {
	my $seq;
	
	foreach my $loc (sort {$a->start <=> $b->start} @{$locs{$id}})
	  {
	    $seq .= $gso->get_sequence(start=>$loc->start,
				       stop =>$loc->stop,
				       chr  =>$loc->chr,
				       strand=>$loc->strand,
				       org_id=>$loc->feat->org->id(),
				       data_info_id=>$loc->feat->data_info->id(),
				      );
	    $rev_seq_flag = 1 if $loc->strand =~ "-";
	  }
	$seq = reverse $seq if $rev_seq_flag;
	push @seqs, $seq;
      }
    @seqs = sort {length($b) <=> length($a)} @seqs;
    return @seqs;
  }

################################################ subroutine header begin ##

=head2 get_genomic_sequence_for_feature

 Usage     : my $seq = $genome_obj->get_genomic_sequence_for_feature($feature_obj)
 Purpose   : Gets the genomic sequence of a feature using a feature object
 Returns   : string (gemomic DNA of feature)
 Argument  : CoGe::Genome::DB::Feature object
 Throws    : none
 Comments  : The sequence is returned in 5'->3' orientation.

See Also   : 

=cut

################################################## subroutine header end ##

sub get_genomic_sequence_for_feature
  {
    my $self = shift;
    my $feat = shift;
    my $seq;
    foreach my $loc ($feat->locs)
      {
	my $tmp_seq = $self->get_genomic_seq_obj->get_sequence(
							 start  => $loc->start,
							 stop   => $loc->stop,
							 chr    => $loc->chr,
							 org_id => $feat->data_info->organism->id,
							 info_id=> $feat->data_info->id,
							 strand => $loc->strand,
							 );
	$tmp_seq = reverse $tmp_seq if $loc->strand =~ /-/;
	$seq .= $tmp_seq if $tmp_seq;
      }
    return $seq;
  }

################################################ subroutine header begin ##

=head2 get_genomic_seq_for_feat

 Usage     : my $seq = $genome_obj->get_genomic_seq_for_feat($feature_obj)
 Purpose   : Alias for get_genomic_sequence_for_feature
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub get_genomic_seq_for_feat
  {
    my $self = shift;
    my $feat = shift;
    return $self->get_genomic_sequence_for_feature($feat);
  }


################################################ subroutine header begin ##
=head2 get_features_for_organism_and_version

 Usage     : my @feats = $db->get_features_for_organism_and_version(org=>$org,
                                                                    version=>6,
								    );
 Purpose   : Get all the features for an organism of a particular version
           : as specified by the data_information table.
 Returns   : an array or array_ref (wantarray) of CoGe::Genome::DB::Feature objects
 Argument  : org => either the name of the organism (which will be looked up in the database) --OR--
                    an CoGe::Genome:DB::Organism object
             version=> an integer that specifies the data inforation version from which to find features
 Throws    : Will return undef if no organism object can be found.

See Also   : 

=cut

################################################## subroutine header end ##


sub get_features_for_organism_and_version
  {
    my $self = shift;
    my %opts = @_;
    my $org = $opts{org} || $opts{organism};
    unless (ref ($org) =~ /organism/i)
      {
         ($org) = $self->get_org_obj->search_like(name=>"$org%");
      }
    return unless ref($org) =~ /organism/i;
    my $version = $opts{version} || $opts{ver} || $opts{v};
    my @feats;
    foreach my $di ($org->data_infos)
      {
         next unless $di->version == $version;
	 push @feats, $di->feats;
      }
    return wantarray ? @feats : \@feats;
  }

1; #this line is important and will help the module return a true value
