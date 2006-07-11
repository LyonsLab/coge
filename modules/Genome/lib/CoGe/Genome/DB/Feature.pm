package CoGe::Genome::DB::Feature;
use strict;
use base 'CoGe::Genome::DB';
use CoGe::Genome::Accessory::Annotation;
use CoGe::Genome;

BEGIN {
    use Exporter ();
    use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = 0.1;
    @ISA         = (@ISA, qw (Exporter));
    #Give a hoot don't pollute, do not export more than needed by default
    @EXPORT      = qw ();
    @EXPORT_OK   = qw ();
    %EXPORT_TAGS = ();
    __PACKAGE__->table('feature');
    __PACKAGE__->columns(All=>qw{feature_id feature_type_id data_information_id});
    __PACKAGE__->has_a('data_information_id'=>'CoGe::Genome::DB::Data_information');
    __PACKAGE__->has_a('feature_type_id'=>'CoGe::Genome::DB::Feature_type');
    __PACKAGE__->has_many('feature_names'=>'CoGe::Genome::DB::Feature_name');
    __PACKAGE__->has_many('locations'=>'CoGe::Genome::DB::Location');
    __PACKAGE__->has_many('sequences'=>'CoGe::Genome::DB::Sequence');
    __PACKAGE__->has_many('annotations'=>'CoGe::Genome::DB::Annotation');
    __PACKAGE__->set_sql('select_features_by_name_and_version' => qq{
SELECT f.feature_id
  FROM feature f
  JOIN data_information di USING (data_information_id)
  JOIN feature_name fn USING (feature_id)
 WHERE fn.name = ?
   AND di.version = ?});
    __PACKAGE__->set_sql ('select_features_in_range' => qq{
SELECT DISTINCT f.feature_id
  FROM feature f
  JOIN location l USING (feature_id)
 WHERE ? <= l.stop
   AND ? >= l.start
   AND chromosome = ?
   AND f.data_information_id = ?;
});
    __PACKAGE__->set_sql ('count_features_in_range' => qq{
SELECT COUNT(DISTINCT f.feature_id)
  FROM feature f
  JOIN location l USING (feature_id)
 WHERE ? <= l.stop
   AND ? >= l.start
   AND chromosome = ?
   AND f.data_information_id = ?;
});
  __PACKAGE__->set_sql ('select_all_feature_ids' => qq{
SELECT feature_id 
  FROM feature
});
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature - CoGe::Genome::DB::Feature

=head1 SYNOPSIS




  #let's write a program that finds features with the name "Kinase" and print some info about them
 
  use CoGe::Genome;
  
  my $db = CoGe::Genome->new;

  foreach my $feat ($db->get_feature_by_name("kinase"))
   {
     print "Feature Type: "  , $feat->type->name,"\n";
     print "Dataset: "       , $feat->data_info->name,"\n";
     print "Dataset version:", $feat->data_info->version,"\n";
     print "Organism:"       , $feat->data_info->organism->name,"\n";
     print "Location:\n";
     foreach my $loc ($feat->locations)
      {
        print "\t", $loc->start,"-", $loc->stop,"\n";
      }
     print "Annotations:\n";
     foreach my $anno ($feat->annotations)
      {
        print $anno->type->name,": " if $anno->type;
        print $anno->annotation,"\n";
      }
     print "Genomic Sequence:\n";
     print $feat->genomic_sequence,"\n";
     print "\n";
   }



=head1 DESCRIPTION

  Feature objects are one of the primary objects in the CoGe::Genome system.
  In a genomic sense, features are any genomic region that have some associated
  set of information or function.  For example, a gene, mRNA, tRNA, SSR, CNS, CRE,
  etc. are all genomic features.  This object provides a central access point 
  for accessing all information for a feature as well as some functions for finding
  particular features.

 This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
 Class::DBI provides the basic methods for creating accessor methods for accessing
 table information.  Please see manual pages for Class::DBI for additional information.


 The columns for this table are:
  feature_id
  feature_type_id
  data_information_id

 Related objects that can be accessed through this object are:
  CoGe::Genome::DB::Feature_type
  CoGe::Genome::DB::Feature_name
  CoGe::Genome::DB::Sequence
  CoGe::Genome::DB::Location
  CoGe::Genome::DB::Annotation
  CoGe::Genome::DB::Data_information


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
 Class::DBI



perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions

new              =>  creates a new object (inherited from Class::Accessor)

feature_id       =>  database entry id
id               =>  alias for location_id

data_information_id => returns a CoGe::Genome::DB::Data_information object
data_information
information
data_info
info
dataset

feature_type_id     => returns a CoGe::Genome::DB::Feature_type object
feature_type
feat_type
type


feature_names       => returns a CoGe::Genome::DB::Feature_name object
feat_names
feat_name
names
name
aliases

locations           => returns a CoGe::Genome::DB::Location object
location
locs
loc

sequences           => returns a CoGe::Genome::DB::Sequence object
seqs

annotations         => returns a CoGe::Genome::DB::Annotation object
annos

=cut

################################################ subroutine header begin ##

sub feature_type
  {
    my $self = shift;
    return $self->feature_type_id();
  }

sub feat_type
  {
    my $self = shift;
    return $self->feature_type_id();
  }

sub type
  {
    my $self = shift;
    return $self->feature_type_id();
  }

sub data_information
  {
    my $self = shift;
    return $self->data_information_id();
  }

sub information
  {
    my $self = shift;
    return $self->data_information_id();
  }

sub data_info
  {
    my $self = shift;
    return $self->data_information_id();
  }

sub info
  {
    my $self = shift;
    return $self->data_information_id();
  }

sub dataset
  {
    my $self = shift;
    return $self->data_information_id();
  }

##legacy table, now stored in data_information
sub organism
  {
    my $self = shift;
    return $self->info->organism_id();
  }

sub org
  {
    my $self = shift;
    return $self->organism();
  }

sub species
  {
    my $self = shift;
    return $self->organism();
  }


sub feat_names
  {
    my $self = shift;
    return $self->feature_names();
  }

sub feat_name
  {
    my $self = shift;
    return $self->feature_names();
  }

sub names
  {
    my $self = shift;
    return $self->feature_names();
  }

sub name
  {
    my $self = shift;
    return $self->feature_names();
  }

sub aliases
  {
    my $self = shift;
    return $self->feature_names();
  }

sub location
  {
    my $self = shift;
    return $self->locations();
  }

sub locs
  {
    my $self = shift;
    return $self->locations();
  }

sub loc
  {
    my $self = shift;
    return $self->locations();
  }

sub seqs
  {
    my $self = shift;
    return $self->sequences();
  }

sub annos
  {
    my $self = shift;
    return $self->annotations();
  }

sub id
  {
    my $self = shift;
    return $self->feature_id();
  }


################################################ subroutine header begin ##

=head2 annotation_pretty_print

 Usage     : my $pretty_annotation = $feat->annotation_pretty_print
 Purpose   : returns a string with information and annotations about a feature
             in a nice format with tabs and new-lines and the like.
 Returns   : returns a string
 Argument  : none
 Throws    : 
 Comments  : uses Coge::Genome::Accessory::Annotation to build the annotations,
           : specifying delimters, and printing to string.   Pretty cool object.

See Also   : CoGe::Genome::Accessory::Annotation

=cut

################################################## subroutine header end ##


sub annotation_pretty_print
  {
    my $self = shift;
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");
    my $start = $self->begin_location;
    my $stop = $self->end_location;
    my $chr = $self->chr;
    my $strand = $self->strand;
    my $info_id = $self->data_info->id;
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} $self->locs);
    $location .="(".$strand.")";
    #my $location = "Chr ".$chr. "".$start."-".$stop.""."(".$strand.")";
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"Location", Values=>[$location], Type_delimit=>": ", Val_delimit=>" "));
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"Name(s)");
    $anno_type->Type_delimit(": ");
    $anno_type->Val_delimit(", ");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot($name->name);
      }
    
    $anno_obj->add_Annot($anno_type);
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>$type->name);
	$anno_type->Val_delimit("\n");

	$anno_type->add_Annot($anno->annotation);
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Genome::Accessory::Annotation(Type=>$group->name);
	    $anno_g->add_Annot($anno_type);
	    $anno_g->Type_delimit(": ");
	    $anno_g->Val_delimit(", ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit(": ");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    return $anno_obj->to_String;
  }


################################################ subroutine header begin ##

=head2 annotation_pretty_print_html

 Usage     : my $pretty_annotation_html = $feat->annotation_pretty_print_html
 Purpose   : returns a string with information and annotations about a feature
             in a nice html format with breaks and class tags (called "annotation")
 Returns   : returns a string
 Argument  : none
 Throws    : 
 Comments  : uses Coge::Genome::Accessory::Annotation to build the annotations,
           : specifying delimters, and printing to string.   Pretty cool object.

See Also   : CoGe::Genome::Accessory::Annotation

=cut

################################################## subroutine header end ##


sub annotation_pretty_print_html
  {
    my $self = shift;
    my %opts = @_;
    my $loc_link = $opts{loc_link};
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n<BR>\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n<BR>\n");
    my $start = $self->begin_location;
    my $stop = $self->end_location;
    my $chr = $self->chr;
    my $strand = $self->strand;
    my $info_id = $self->data_info->id;
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_name = $type->name;
	$anno_name = "<font class=\"annotation\">". $anno_name."</font>" unless ref($group) =~ /group/i;
	
	my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>$anno_name);
	$anno_type->Val_delimit("\n<li>");

	$anno_type->add_Annot($anno->annotation);
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"annotation\">".$group->name."</font>");
	    $anno_g->add_Annot($anno_type);
	    $anno_g->Type_delimit("\n<li>");
	    $anno_g->Val_delimit("\n<li>");
#	    $anno_g->Val_delimit(" ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit("\n<li>");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"annotation\">"."Name(s)"."</font>");
    $anno_type->Type_delimit("\n<BR><li>");
    $anno_type->Val_delimit("\n<li>");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot($name->name);
      }
    
    $anno_obj->add_Annot($anno_type);
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} $self->locs);
    $location .="(".$strand.")";
    $location = qq{<a href="$loc_link?start=$start&stop=$stop&chr=$chr&di=$info_id&strand=$strand">}.$location."</a>\n" if $loc_link;
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"annotation\">Location</font>", Values=>[$location], Type_delimit=>"\n<BR><li>", Val_delimit=>" "));
    return $anno_obj->to_String;
  }




################################################ subroutine header begin ##

=head2 genbank_location_string

 Usage     : my $genbank_loc = $feat->genbank_location_string
 Purpose   : generates a genbank location string for the feature in genomic coordinates or
           : based on a recalibration number that is user specified
           : e.g.: complement(join(10..100,200..400))
 Returns   : a string
 Argument  : hash:  recalibrate => number of positions to subtract from genomic location
 Throws    : none
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##


sub genbank_location_string
  {
    my $self = shift;
    my %opts = @_;
    my $recal = $opts{recalibrate};
    my $string;
    my $count= 0;
    my $comp = 0;
    foreach my $loc (sort {$a->start <=> $b->start}  $self->locs())
      {
	$comp = 1 if $loc->strand =~ "-";
	$string .= "," if $count;
	$string .= $recal ? ($loc->start-$recal+1)."..".($loc->stop-$recal+1): $loc->start."..".$loc->stop;
	$count++;
      }
    $string = "join(".$string.")" if $count > 1;
    $string = "complement(".$string.")" if $comp;
    return $string;
  }


################################################ subroutine header begin ##

=head2 begin_location

 Usage     : my $feat_start = $feat->begin_location
 Purpose   : returns the start of the feature (does not take into account the strand on which
             the feature is located)
 Returns   : a string, number usually
 Argument  : none
 Throws    : 
 Comments  : this simply calles $feat->locs, sorts them based on their starting position, and
           : returns the smallest position

See Also   : 

=cut

################################################## subroutine header end ##


sub begin_location
  {
    my $self = shift;
    my @locs = $self->locs;
    return unless scalar @locs;
    my ($val) = sort {$a->begin <=> $b->begin} @locs;
    return $val->begin;
  }

################################################ subroutine header begin ##

=head2 start_location

 Usage     : my $feat_start = $feat->start_location
 Purpose   : alias for $feat->begin_location

=cut

################################################## subroutine header end ##

sub start_location
  {
    my $self = shift;
    return $self->begin_location;
  }

################################################ subroutine header begin ##

=head2 start

 Usage     : my $feat_start = $feat->start
 Purpose   : alias for $feat->begin_location

=cut

################################################## subroutine header end ##

sub start
  {
    my $self = shift;
    return $self->begin_location;
  }

################################################ subroutine header begin ##

=head2 begin

 Usage     : my $feat_start = $feat->begin
 Purpose   : alias for $feat->begin_location

=cut

################################################## subroutine header end ##

sub begin
  {
    my $self = shift;
    return $self->begin_location;
  }


################################################ subroutine header begin ##

=head2 end_location

 Usage     : my $feat_end = $feat->end_location
 Purpose   : returns the end of the feature (does not take into account the strand on which
             the feature is located)
 Returns   : a string, number usually
 Argument  : none
 Throws    : 
 Comments  : this simply calles $feat->locs, sorts them based on their ending position, and
           : returns the largest position

See Also   : 

=cut

################################################## subroutine header end ##


sub end_location
  {
    my $self = shift;
    my @locs = $self->locs;
    return unless scalar @locs;
    my ($val) = sort {$b->end <=> $a->end} @locs;
    return $val->end;
  }

################################################ subroutine header begin ##

=head2 stop_location

 Usage     : my $feat_end = $feat->stop_location
 Purpose   : alias for $feat->end_location

=cut

################################################## subroutine header end ##

sub stop_location
  {
    my $self = shift;
    return $self->end_location;
  }

################################################ subroutine header begin ##

=head2 stop

 Usage     : my $feat_end = $feat->stop
 Purpose   : alias for $feat->end_location

=cut

################################################## subroutine header end ##

sub stop
  {
    my $self = shift;
    return $self->end_location;
  }

################################################ subroutine header begin ##

=head2 end

 Usage     : my $feat_end = $feat->end
 Purpose   : alias for $feat->end_location

=cut

################################################## subroutine header end ##

sub end
  {
    my $self = shift;
    return $self->end_location;
  }


################################################ subroutine header begin ##

=head2 chromosome

 Usage     : my $chr = $feat->chromosome
 Purpose   : return the chromosome of the feature
 Returns   : a string
 Argument  : none
 Throws    : none
 Comments  : returns $self->locs->next->chr
           : 

See Also   : 

=cut

################################################## subroutine header end ##


sub chromosome
  {
    my $self = shift;
    my $loc = $self->locs->next;
    return unless $loc;
    return $loc->chr();
  }

################################################ subroutine header begin ##

=head2 chr

 Usage     : my $chr = $feat->chr
 Purpose   : alias for $feat->chromosome

=cut

################################################## subroutine header end ##

sub chr
  {
    my $self = shift;
    return $self->chromosome;
  }

################################################ subroutine header begin ##

=head2 strand

 Usage     : my $strand = $feat->strand
 Purpose   : return the chromosome strand of the feature
 Returns   : a string (usally something like 1, -1, +, -, etc)
 Argument  : none
 Throws    : none
 Comments  : returns $self->locs->next->strand
           : 

See Also   : 

=cut

################################################## subroutine header end ##


sub strand
  {
    my $self = shift;
    my $loc = $self->locs->next;
    return unless $loc;
    return $loc->strand();
  }

################################################ subroutine header begin ##

=head2 get_features_in_region

 Usage     : $object->get_features_in_region(start   => $start, 
                                             stop    => $stop, 
                                             chr     => $chr,
                                             info_id => $data_info->id());

 Purpose   : gets all the features in a specified genomic region
 Returns   : an array or an array_ref of feature objects (wantarray)
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             info_id => data_information id in database (obtained from a
                        CoGe::Data_information object)
                        of the dna seq will be returned
             OPTIONAL
             count   => flag to return only the number of features in a region
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_features_in_region
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{'start'} || $opts{'START'} || $opts{begin} || $opts{BEGIN};
    $start = 0 unless $start;
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    $stop = $start unless defined $stop;
    my $chr = $opts{chr} || $opts{CHR} || $opts{chromosome} || $opts{CHROMOSOME};
    my $info_id = $opts{info_id} || $opts{INFO_ID} || $opts{data_info_id} || $opts{DATA_INFO_ID};
    my $count_flag = $opts{count} || $opts{COUNT};
    my $sth = $count_flag ? $self->sql_count_features_in_range : $self->sql_select_features_in_range();
    $sth->execute($start, $stop, $chr, $info_id);
    if ($count_flag)
      {
	my $count = $sth->fetch->[0];
	$sth->finish;
	return $count;
      }
    my @feats;
    while (my $q = $sth->fetch())
      {
	push @feats, $self->search(feature_id=>$q->[0]);
      }
    $sth->finish;
    return wantarray ? @feats : \@feats;
  }

################################################ subroutine header begin ##

=head2 count_features_in_region

 Usage     : $object->count_features_in_region(start   => $start, 
                                             stop    => $stop, 
                                             chr     => $chr,
                                             info_id => $data_info->id());

 Purpose   : counts the features in a specified genomic region
 Returns   : an integer
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             info_id => data_information id in database (obtained from a
                        CoGe::Data_information object)
                        of the dna seq will be returned
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub count_features_in_region
  {
    my $self = shift;
    my %opts = @_;
    return $self->get_features_in_region (%opts, count=>1);
  }

################################################ subroutine header begin ##

=head2 get_features_by_name_and_data_information_version

 Usage     : $object->get_features_by_name_and_data_information_version(name   => $name, 
                                                                        version=> $ver);
 Purpose   : gets all the features based on a name and data information version
 Returns   : an array or an array_ref of feature objects (wantarray)
 Argument  : name   => genomic start position
             version=> data_information version (obtained from a
                        CoGe::Data_information object)
 Throws    : undef if $name or $ver are undefined
 Comments  : the name is obtained from a join to the feature_name table
             the version is obtained from a join to the data_information table

See Also   : CoGe::Genome::DB::Data_information

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_data_information_version
  {
    my $self = shift;
    my %opts = @_;
    my $name = $opts{name} || $opts{NAME};
    my $ver = $opts{version} || $opts{VERSION} || $opts{ver} || $opts{VER};
    my $sth = $self->sql_select_features_by_name_and_version();
    return unless defined $name && defined $ver;
    $sth->execute($name, $ver);
    my @feats;
    while (my $q = $sth->fetch())
      {
	push @feats, $self->search(feature_id=>$q->[0]);
      }
    $sth->finish;
    return wantarray ? @feats : \@feats;
  }


################################################ subroutine header begin ##

=head2 get_features_by_name_and_data_info_version

 Usage     : $object->get_features_by_name_and_data_info_version(name   => $name, 
                                                                 version=> $ver);
 Purpose   : alias for get_features_by_name_and_data_information_version

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_data_info_version
  {
    my $self = shift;
    return ($self->get_features_by_name_and_data_information_version(@_));
  }

################################################ subroutine header begin ##

=head2 get_features_by_name_and_version

 Usage     : $object->get_features_by_name_and_version(name   => $name, 
                                                       version=> $ver);
 Purpose   : alias for get_features_by_name_and_data_information_version

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_version
  {
    my $self = shift;
    return ($self->get_features_by_name_and_data_information_version(@_));
  }

################################################ subroutine header begin ##

=head2 get_all_feature_ids

 Usage     : my @ids = $object->get_all_feature_ids
 Purpose   : get all the feature ids from the database
 Returns   : an array or arrayref based on wantarray
 Argument  : none
 Comments  : a simple method to get all the feature ids without creating 
             CoGe::Genome::DB::Feature objects (which could be potentially
             dangerous)

See Also   : 

=cut

################################################## subroutine header end ##

sub get_all_feature_ids
  {
    my $self = shift;
    my $sth = $self->sql_select_all_feature_ids;
    $sth->execute;
    my @ids = map {$_->[0]} $sth->fetchall;
    $sth->finish;
    return wantarray ? @ids : \@ids;
  }

################################################ subroutine header begin ##

=head2 genomic_sequence

 Usage     : my $genomic_seq = $feat->genomic_sequence
 Purpose   : gets the genomic seqence for a feature
 Returns   : a string
 Argument  : none
 Comments  : This method simply creates a CoGe::Genome object and calls:
             get_genomic_sequence_for_feature($self)
See Also   : CoGe::Genome

=cut

################################################## subroutine header end ##


sub genomic_sequence
  {
    my $self = shift;
    my $db = CoGe::Genome->new;
    return $db->get_genomic_sequence_for_feature($self);
  }

################################################ subroutine header begin ##

=head2 

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


1; #this line is important and will help the module return a true value


