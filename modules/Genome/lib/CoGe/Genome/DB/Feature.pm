package CoGe::Genome::DB::Feature;
use strict;
use base 'CoGe::Genome::DB';
use CoGe::Genome::Accessory::Annotation;
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
    __PACKAGE__->table('feature');
    __PACKAGE__->columns(All=>qw{feature_id feature_type_id dataset_id});
    __PACKAGE__->has_a('dataset_id'=>'CoGe::Genome::DB::Dataset');
    __PACKAGE__->has_a('feature_type_id'=>'CoGe::Genome::DB::Feature_type');
    __PACKAGE__->has_many('feature_names'=>'CoGe::Genome::DB::Feature_name');
    __PACKAGE__->has_many('locations'=>'CoGe::Genome::DB::Location');
    __PACKAGE__->has_many('sequences'=>'CoGe::Genome::DB::Sequence');
    __PACKAGE__->has_many('annotations'=>'CoGe::Genome::DB::Annotation');
    __PACKAGE__->set_sql('select_features_by_name_and_version' => qq{
SELECT f.feature_id
  FROM feature f
  JOIN dataset ds USING (dataset_id)
  JOIN feature_name fn USING (feature_id)
 WHERE fn.name = ?
   AND ds.version = ?});
    __PACKAGE__->set_sql ('select_features_in_range' => qq{
SELECT DISTINCT f.feature_id
  FROM feature f
  JOIN location l USING (feature_id)
 WHERE ? <= l.stop
   AND ? >= l.start
   AND chromosome = ?
   AND f.dataset_id = ?;
});
    __PACKAGE__->set_sql ('count_features_in_range' => qq{
SELECT COUNT(DISTINCT f.feature_id)
  FROM feature f
  JOIN location l USING (feature_id)
 WHERE ? <= l.stop
   AND ? >= l.start
   AND chromosome = ?
   AND f.dataset_id = ?;
});
  __PACKAGE__->set_sql ('select_all_feature_ids' => qq{
SELECT feature_id 
  FROM feature
});

    __PACKAGE__->set_sql ('no_name_features' => qq{
SELECT feature_id from feature
  LEFT OUTER JOIN feature_name using (feature_id)
 WHERE feature_name.name is NULL
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
     print "Dataset: "       , $feat->dataset->name,"\n";
     print "Dataset version:", $feat->dataset->version,"\n";
     print "Organism:"       , $feat->dataset->organism->name,"\n";
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
  dataset_id

 Related objects that can be accessed through this object are:
  CoGe::Genome::DB::Feature_type
  CoGe::Genome::DB::Feature_name
  CoGe::Genome::DB::Sequence
  CoGe::Genome::DB::Location
  CoGe::Genome::DB::Annotation
  CoGe::Genome::DB::Dataset


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

data_information_id => obselete - now returns a CoGe::Genome::DB::Dataset object
data_information
data_info
info


dataset_id => returns a CoGe::Genome::DB::Dataset object
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

aliases             => like the above, but returns the actual name(s) instead of
                       Feature_name_obj

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

sub dataset
  {
    my $self = shift;
    return $self->dataset_id();
  }

sub data_information_id
  {
    my $self = shift;
    carp "data_information_id is obselete. Please use dataset_id";
    return $self->dataset_id();
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

##legacy table, now stored in dataset
sub organism
  {
    my $self = shift;
    return $self->dataset->organism_id();
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
    my @names = map {$_->name} $self->feature_names();
    return wantarray ? @names : \@names;
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
    #look into changing this to set_id
    my $info_id = $self->dataset->id;
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
    my $info_id = $self->dataset->id;
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"title4\">"."Name(s):"."</font>");
    $anno_type->Type_delimit("\n");
    $anno_type->Val_delimit("\n, ");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot("<a class=\"data\" href=FeatView.pl?accn=".$name->name.">".$name->name."</a>");
      }
    
    $anno_obj->add_Annot($anno_type);
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_name = $type->name;
	$anno_name = "<font class=\"title4\">". $anno_name."</font>" unless ref($group) =~ /group/i;
	
	my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>$anno_name);
	$anno_type->Val_delimit(", ");

	$anno_type->add_Annot("<font class=\"data\">".$anno->annotation."</font>");
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"title4\">".$group->name."</font>");
	    $anno_g->add_Annot($anno_type);
	    $anno_g->Type_delimit(": ");
	    $anno_g->Val_delimit(", ");
#	    $anno_g->Val_delimit(" ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit(": ");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} $self->locs);
    $location .="(".$strand.")";
    $location = qq{<a href="$loc_link?start=$start&stop=$stop&chr=$chr&di=$info_id&strand=$strand">}.$location."</a\n" if $loc_link;
    $location = qq{<font class="data">$location</font>};
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"title4\">Location</font>", Values=>[$location], Type_delimit=>": ", Val_delimit=>" "));
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

=head2 version

 Usage     : my $version = $feat->version
 Purpose   : return the dataset version of the feature
 Returns   : an integer
 Argument  : none
 Throws    : none
 Comments  : returns $self->dataset->version
           : 

See Also   : 

=cut

################################################## subroutine header end ##


sub version
  {
    my $self = shift;
    return $self->dataset->version();
  }

################################################ subroutine header begin ##

=head2 get_features_in_region

 Usage     : $object->get_features_in_region(start   => $start, 
                                             stop    => $stop, 
                                             chr     => $chr,
                                             dataset_id => $dataset->id());

 Purpose   : gets all the features in a specified genomic region
 Returns   : an array or an array_ref of feature objects (wantarray)
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             dataset_id => dataset id in database (obtained from a
                        CoGe::Dataset object)
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
    my $dataset_id = $opts{dataset} || $opts{dataset_id} || $opts{info_id} || $opts{INFO_ID} || $opts{data_info_id} || $opts{DATA_INFO_ID} ;
    my $count_flag = $opts{count} || $opts{COUNT};
    my $sth = $count_flag ? $self->sql_count_features_in_range : $self->sql_select_features_in_range();
    $sth->execute($start, $stop, $chr, $dataset_id);
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
                                             dataset_id => $dataset->id());

 Purpose   : counts the features in a specified genomic region
 Returns   : an integer
 Argument  : start   => genomic start position
             stop    => genomic stop position
             chr     => chromosome
             dataset_id => dataset id in database (obtained from a
                        CoGe::Dataset object)
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

=head2 get_features_by_name_and_dataset_version

 Usage     : $object->get_features_by_name_and_dataset_version(name   => $name, 
                                                                        version=> $ver);
 Purpose   : gets all the features based on a name and dataset version
 Returns   : an array or an array_ref of feature objects (wantarray)
 Argument  : name   => genomic start position
             version=> dataset version (obtained from a
                        CoGe::Dataset object)
 Throws    : undef if $name or $ver are undefined
 Comments  : the name is obtained from a join to the feature_name table
             the version is obtained from a join to the dataset table

See Also   : CoGe::Genome::DB::Dataset

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_dataset_version
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

=head2 get_features_by_name_and_data_information_version

 Usage     : $object->get_features_by_name_and_data_information_version(name   => $name, 
                                                                        version=> $ver);
 Purpose   : See get_features_by_name_and_dataset_version
 Returns   :
 Argument  :
 Throws    :
 Comments  : Obselete method. Please use get_features_by_name_and_dataset_version

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_data_information_version
  {
    my $self = shift;
    print STDERR "get_features_by_name_and_data_information_version is obselete. Please use get_features_by_name_and_dataset_version";
    return ($self->get_features_by_name_and_dataset_version(@_));
  }

################################################ subroutine header begin ##

=head2 get_features_by_name_and_data_info_version

 Usage     : $object->get_features_by_name_and_data_info_version(name   => $name, 
                                                                 version=> $ver);
 Purpose   : alias for get_features_by_name_and_dataset_version

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_data_info_version
  {
    my $self = shift;
    print STDERR "get_features_by_name_and_data_info_version is obselete. Please use get_features_by_name_and_dataset_version";
    return ($self->get_features_by_name_and_dataset_version(@_));
  }

################################################ subroutine header begin ##

=head2 get_features_by_name_and_version

 Usage     : $object->get_features_by_name_and_version(name   => $name, 
                                                       version=> $ver);
 Purpose   : alias for get_features_by_name_and_dataset_version

=cut

################################################## subroutine header end ##


sub get_features_by_name_and_version
  {
    my $self = shift;
    return ($self->get_features_by_name_and_dataset_version(@_));
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
    my %opts = @_;
    my $upstream = $opts{upstream} || $opts{up} || $opts{us} || 0;
    my $downstream = $opts{downstream} || $opts{down} || $opts{ds} || 0;
    my $db = CoGe::Genome->new;
    return $db->get_genomic_sequence_for_feature(feat=>$self, up=>$upstream, down=>$downstream);
  }


sub has_genomic_sequence
  {
    my $self = shift;
    return 1 if $self->dataset->has_genomic_sequence;
    return 0;
  }

################################################ subroutine header begin ##

=head2 protein_sequence

 Usage     : my @prots = $feat->protein_sequence
 Purpose   : trys to find protein sequences for this feature
 Returns   : an array or array ref depending on wantarray
 Argument  : none
 Comments  : this routine searches through associated sequence objects for those
             of type protein
See Also   : CoGe::Genome::DB::Sequence

=cut

################################################## subroutine header end ##

sub protein_sequence
  {
    my $self = shift;
    my @seqs;
    foreach my $seq ($self->sequences)
      {
	next unless $seq->type->name =~ /protein/;
	push @seqs, $seq->seq;
      }
    return wantarray ? @seqs : \@seqs;
  }

################################################ subroutine header begin ##

=head2 blast_bit_score

 Usage     : my $bit_score = $feature->blast_bit_score();
 Purpose   : returns the blast bit score for the feature's self-self identical hit
 Returns   : an int -- the blast bit score
 Argument  : optional hash
             match    => the score for a nucleotide match. DEFAULT: 1
             mismatch => the score for a nucleotide mismatch.  DEFAULT: -3
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub blast_bit_score
  {
    my $self = shift;
    my %opts = @_;
    my $match = $opts{match} || 1;
    my $mismatch = $opts{mismatch} || -3;
    my $lambda = $self->_estimate_lambda(match=>$match, mismatch=>$mismatch);
    my $seq = $self->genomic_sequence();
    warn "No genomic sequence could be obtained for this feature object.  Can't calculate a blast bit score.\n" unless $seq;
    my $bs = sprintf("%.0f", $lambda*length($seq)*$match/log(2));
    return $bs;
  }

################################################ subroutine header begin ##

=head2 _estimate_lambda

 Usage     : my $lambda = $feature->_estimate_lambda
 Purpose   : estimates lambda for calculating blast bit scores.  Lambda is
             a matrix-specific constant for normalizing raw blast scores 
 Returns   : a number, lambda
 Argument  : optional hash
             match    => the score for a nucleotide match. DEFAULT: 1
             mismatch => the score for a nucleotide mismatch.  DEFAULT: -3
             precision=> the different between the high and low estimate 
                         of lambda before lambda is returned.  
                         DEFAULT: 0.001
 Throws    : a warning if there is a problem with the calcualted expected_score
             or the match score is less than 0;
 Comments  : Assumes an equal probability for each nucleotide.
           : this routine is based on example 4-1 from 
           : BLAST: An essential guide to the Basic Local Alignment Search Tool 
           : by Korf, Yandell, and Bedell published by O'Reilly press.

See Also   : 

=cut

################################################## subroutine header end ##


sub _estimate_lambda
  {
    #this routine is based on example 4-1 from BLAST: An essential guide to the Basic Local Alignment Search Tool by Korf, Yandell, and Bedell published by O'Reilly press.
    my $self = shift;
    my %opts = @_;
    my $match = $opts{match} || 1;
    my $mismatch = $opts{mismatch} || -3;
    my $precision = $opts{precision} || 0.001;
      
    use constant Pn => 0.25; #prob of any nucleotide
    my $expected_score = $match * 0.25 + $mismatch * 0.75; 
    if ($match <= 0 or $expected_score >= 0)
      {
	warn qq{
Problem with scores.  Match: $match (should be greater than 0).
             Expected score: $expected_score (should be less than 0).
};
	return 0;
      }
    # calculate lambda 
    my ($lambda, $high, $low) = (1, 2, 0); # initial estimates 
    while ($high - $low > $precision) 
      {         # precision 
	# calculate the sum of all normalized scores 
	my $sum = Pn * Pn * exp($lambda * $match) * 4 
	  + Pn * Pn * exp($lambda * $mismatch) * 12; 
	# refine guess at lambda 
	if ($sum > 1) 
	  { 
	    $high = $lambda;
	    $lambda = ($lambda + $low)/2; 
	  } 
	else 
	  {
	  $low = $lambda; 
	  $lambda = ($lambda + $high)/2; 
	}
      }
    # compute target frequency and H 
    my $targetID = Pn * Pn * exp($lambda * $match) * 4; 
    my $H = $lambda * $match    *     $targetID 
      + $lambda * $mismatch * (1 -$targetID); 
    # output 
#    print "expscore: $expected_score\n"; 
#    print "lambda:   $lambda nats (", $lambda/log(2), " bits)\n"; 
#    print "H:        $H nats (", $H/log(2), " bits)\n"; 
#    print "%ID:      ", $targetID * 100, "\n"; 

    return $lambda;
  }


sub get_no_name_features
  {
    my $self = shift;
    my $sth = $self->sql_no_name_features;
    $sth->execute;
    my @ids = map {$self->retrieve($_->[0])} $sth->fetchall;
    $sth->finish;
    return wantarray ? @ids : \@ids;
  }

sub reverse_complement
  {
    my $self = shift;
    my $seq = shift || $self->genomic_sequence;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/; 
    return $rcseq;
  }

sub frame6_trans
  {
    my $self = shift;
    my %opts = @_;
    my(%code) = (
		 'TCA' => 'S',# Serine
		 'TCC' => 'S',# Serine
		 'TCG' => 'S',# Serine
		 'TCT' => 'S',# Serine
		 'TTC' => 'F',# Fenilalanine
		 'TTT' => 'F',# Fenilalanine
		 'TTA' => 'L',# Leucine
		 'TTG' => 'L',# Leucine
		 'TAC' => 'Y',# Tirosine
		 'TAT' => 'Y',# Tirosine
		 'TAA' => '*',# Stop
		 'TAG' => '*',# Stop
		 'TGC' => 'C',# Cysteine
		 'TGT' => 'C',# Cysteine
		 'TGA' => '*',# Stop
		 'TGG' => 'W',# Tryptofane
		 'CTA' => 'L',# Leucine
		 'CTC' => 'L',# Leucine
		 'CTG' => 'L',# Leucine
		 'CTT' => 'L',# Leucine
		 'CCA' => 'P',# Proline
		 'CCC' => 'P',# Proline
		 'CCG' => 'P',# Proline
		 'CCT' => 'P',# Proline
		 'CAC' => 'H',# Hystidine
		 'CAT' => 'H',# Hystidine
		 'CAA' => 'Q',# Glutamine
		 'CAG' => 'Q',# Glutamine
		 'CGA' => 'R',# Arginine
		 'CGC' => 'R',# Arginine
		 'CGG' => 'R',# Arginine
		 'CGT' => 'R',# Arginine
		 'ATA' => 'I',# IsoLeucine
		 'ATC' => 'I',# IsoLeucine
		 'ATT' => 'I',# IsoLeucine
		 'ATG' => 'M',# Methionina
		 'ACA' => 'T',# Treonina
		 'ACC' => 'T',# Treonina
		 'ACG' => 'T',# Treonina
		 'ACT' => 'T',# Treonina
		 'AAC' => 'N',# Asparagina
		 'AAT' => 'N',# Asparagina
		 'AAA' => 'K',# Lisina
		 'AAG' => 'K',# Lisina
		 'AGC' => 'S',# Serine
		 'AGT' => 'S',# Serine
		 'AGA' => 'R',# Arginine
		 'AGG' => 'R',# Arginine
		 'GTA' => 'V',# Valine
		 'GTC' => 'V',# Valine
		 'GTG' => 'V',# Valine
		 'GTT' => 'V',# Valine
		 'GCA' => 'A',# Alanine
		 'GCC' => 'A',# Alanine
		 'GCG' => 'A',# Alanine
		 'GCT' => 'A',# Alanine
		 'GAC' => 'D',# Aspartic Acid
		 'GAT' => 'D',# Aspartic Acid
		 'GAA' => 'E',# Glutamic Acid
		 'GAG' => 'E',# Glutamic Acid
		 'GGA' => 'G',# Glicine
		 'GGC' => 'G',# Glicine
		 'GGG' => 'G',# Glicine
		 'GGT' => 'G',# Glicine
		);
    my $code = $opts{code} || \%code;
    my $seq = $opts{seq} || $self->genomic_sequence;

    my %seqs;
    $seqs{"1"} = $self->_process_seq(seq=>$seq, start=>0, code1=>$code, codonl=>3);
    $seqs{"2"} = $self->_process_seq(seq=>$seq, start=>1, code1=>$code, codonl=>3);
    $seqs{"3"} = $self->_process_seq(seq=>$seq, start=>2, code1=>$code, codonl=>3);
    my $rcseq = $self->reverse_complement;
    $seqs{"-1"} = $self->_process_seq(seq=>$rcseq, start=>0, code1=>$code, codonl=>3);
    $seqs{"-2"} = $self->_process_seq(seq=>$rcseq, start=>1, code1=>$code, codonl=>3);
    $seqs{"-3"} = $self->_process_seq(seq=>$rcseq, start=>2, code1=>$code, codonl=>3);
    return \%seqs;

  }

sub alpha_trans
  {
    my $self = shift;
    my %opts = @_;
    my $alter = $opts{alter};
    $alter = "X" unless defined $alter; #place to store the symbol used for characters not translated
    my $hyb = $opts{hyb};
    my %code1 = (
	     "UU"=>"L",
	     "UA"=>"Y",
	     "UG"=>"C",
	     "AU"=>"I",
	     "AA"=>"K",
	     "AG"=>"R",
	     "CU"=>"L",
#	     "CA"=>"Q",
	     "CG"=>"R",
	     "GU"=>"V",
	     "GA"=>"E",
	    );
    my $code1 = $opts{code1} || \%code1;
    my %code2 = (
	     "UU"=>"F",
	     "UC"=>"S",
#	     "UG"=>"*",
#	     "UA"=>"*",
	     "AC"=>"T",
	     "AA"=>"N",
	     "AG"=>"S",
	     "CC"=>"P",
#	     "CA"=>"H",
	     "GC"=>"A",
	     "GA"=>"D",
	     "GG"=>"G",
	     );
    my $code2 = $opts{code2} || \%code2;

    my $seq = $opts{seq} || $self->genomic_sequence;
    my %seqs;
    $seqs{I1} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>0, code1=>$code1, codonl=>2);
    $seqs{I2} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>1, code1=>$code1, codonl=>2);
    $seqs{I1H} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>0, code1=>$code1, code2=>$code2, codonl=>2) if $hyb;
    $seqs{I2H} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>1, code1=>$code1, code2=>$code2, codonl=>2) if $hyb;;
    $seqs{II1} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>0, code1=>$code2, codonl=>2);
    $seqs{II2} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>1, code1=>$code2, codonl=>2);
    $seqs{II1H} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>0, code1=>$code2, code2=>$code1, codonl=>2) if $hyb;;
    $seqs{II2H} = $self->_process_seq(alter=>$alter, seq=>$seq, start=>1, code1=>$code2, code2=>$code1, codonl=>2) if $hyb;;
    my $rcseq = $self->reverse_complement;
    $seqs{"I-1"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>0, code1=>$code1, codonl=>2);
    $seqs{"I-2"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>1, code1=>$code1, codonl=>2);
    $seqs{"I-1H"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>0, code1=>$code1, code2=>$code2, codonl=>2) if $hyb;;
    $seqs{"I-2H"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>1, code1=>$code1, code2=>$code2, codonl=>2) if $hyb;;
    $seqs{"II-1"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>0, code1=>$code2, codonl=>2);
    $seqs{"II-2"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>1, code1=>$code2, codonl=>2);
    $seqs{"II-1H"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>0, code1=>$code2, code=>$code1, codonl=>2) if $hyb;;
    $seqs{"II-2H"} = $self->_process_seq(alter=>$alter, seq=>$rcseq, start=>1, code1=>$code2, code=>$code1, codonl=>2) if $hyb;;
    return \%seqs;
  }

sub _process_seq
  {
    my $self = shift;
    my %opts = @_;
    my $seq = $opts{seq};
    my $start = $opts{start};
    my $code1 = $opts{code1};
    my $code2 = $opts{code2};
    my $alter = $opts{alter};
    my $codonl = $opts{codonl} || 2;
    my $seq_out;
    for (my $i = $start; $i < length ($seq); $i = $i+$codonl)
      {
	my $codon = uc(substr($seq, $i, $codonl));
	my $chr = $code1->{$codon} || $code2->{$codon};
	unless ($chr)
	  {
	    $chr= $alter if $alter;
	  }
	$seq_out .= $chr if $chr;
      }
    return $seq_out;
  }

sub percent_translation_system
  {
    my $self = shift;
    my %opts = @_;
    
    my %code1 = (
	     "L"=>1,
	     "Y"=>1,
	     "C"=>1,
	     "I"=>1,
	     "K"=>1,
	     "R"=>1,
	     "L"=>1,
	     "Q"=>1,
	     "R"=>1,
	     "V"=>1,
	     "E"=>1,
	    );
    my $code1 = $opts{code1} || \%code1;
    my %code2 = (
	     "F"=>1,
	     "S"=>1,
	     "*"=>1,
	     "*"=>1,
	     "T"=>1,
	     "N"=>1,
	     "S"=>1,
	     "P"=>1,
	     "H"=>1,
	     "A"=>1,
	     "D"=>1,
	     "G"=>1,
	     );
    my $code2 = $opts{code2} || \%code2;
    my ($seq) = $opts{seq} || $self->protein_sequence;
    return (0,0) unless $seq;
    my ($c1, $c2, $total) = (0,0,0);
    foreach (split //, $seq)
      {
	$c1++ if $code1->{$_};
	$c2++ if $code2->{$_};
	$total++;
      }
    return (map {sprintf("%.4f", $_)} $c1/$total, $c2/$total);
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


