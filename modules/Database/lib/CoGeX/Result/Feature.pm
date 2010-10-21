package CoGeX::Result::Feature;

use strict;
use warnings;
use base 'DBIx::Class::Core';
use CoGe::Accessory::genetic_code;
use Text::Wrap;
use Data::Dumper;
use CoGe::Accessory::Annotation;
use base 'Class::Accessor';

=head1 NAME

CoGeX::Feature

=head1 SYNOPSIS

This object uses the DBIx::Class to define an interface to the C<feature> table in the CoGe database.


=head1 DESCRIPTION

Has columns:
C<feature_id> (Primary Key)
Type: INT, Default: undef, Nullable: no, Size: 11

C<feature_type_id>
Type: INT, Default: 0, Nullable: no, Size: 11

C<dataset_id>
Type: INT, Default: 0, Nullable: no, Size: 11

C<start>
Type: INT, Default: 0, Nullable: yes, Size: 11

C<stop>
Type: INT, Default: 0, Nullable: yes, Size: 11

C<strand>
Type: TINYINT, Default: 0, Nullable: yes, Size: 4

C<chromosome>
Type: VARCHAR, Default: 0, Nullable: yes, Size: 255


Belongs to CCoGeX::Result::FeatureType> via C<feature_type_id>
Belongs to CCoGeX::Result::Dataset> via C<dataset_id>

Has many CCoGeX::Result::FeatureName> via C<feature_id>
Has many CCoGeX::Result::Annotation> via C<feature_id>
Has many CCoGeX::Result::Location> via C<feature_id>
Has many CCoGeX::Result::Sequence> via C<feature_id>

=head1 USAGE

  use CoGeX;

=head1 METHODS

=cut

__PACKAGE__->table("feature");
__PACKAGE__->add_columns(
  "feature_id",{ data_type => "INT", default_value => undef, is_nullable => 0, size => 11 },
  "feature_type_id",{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "dataset_id",{ data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
  "start",{ data_type => "INT", default_value => 0, is_nullable => 1, size => 11 },
  "stop",{ data_type => "INT", default_value => 0, is_nullable => 1, size => 11 },
  "strand",{ data_type => "TINYINT", default_value => 0, is_nullable => 1, size => 4 },
  "chromosome",{ data_type => "VARCHAR", default_value => 0, is_nullable => 1, size => 255 },
  "access_count",{ data_type => "INT", default_value => 0, is_nullable => 1, size => 10 },
);

__PACKAGE__->set_primary_key("feature_id");

# feature has many feature_names
__PACKAGE__->has_many( 'feature_names' => "CoGeX::Result::FeatureName", 'feature_id');

# feature has many annotations
__PACKAGE__->has_many( 'annotations' => "CoGeX::Result::Annotation", 'feature_id');

# feature has many locations
__PACKAGE__->has_many( 'locations' => "CoGeX::Result::Location", 'feature_id');

# feature has many sequences - note this is non-genomic sequence (see
# genomic_sequence() and sequence
__PACKAGE__->has_many( 'sequences' => "CoGeX::Result::Sequence", 'feature_id');

# feature_type has many features
__PACKAGE__->belongs_to("feature_type" => "CoGeX::Result::FeatureType", 'feature_type_id');

# dataset has many features
__PACKAGE__->belongs_to("dataset" => "CoGeX::Result::Dataset", 'dataset_id');


__PACKAGE__->mk_accessors('_genomic_sequence', 'gst', 'dsg'); #_genomic_sequence =>place to store the feature's genomic sequence with no up and down stream stuff


################################################ subroutine header begin ##

=head2 type

 Usage     : $returned_featuretype_object = $FeatureObject->type();
 Purpose   : Shortcut to return a FeatureType object from a Feature object.
 Returns   : A FeatureType object.
 Argument  : None
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub type
  {
    my $self = shift;
    return $self->feature_type();
  }


################################################ subroutine header begin ##

=head2 dataset_groups

 Usage     : $returned_dataset_group_objecst = $FeatureObject->dataset_groups();
 Purpose   : Shortcut to return dataset group objects from a Feature object.
 Returns   : An DatasetGroup objects. (array or array ref depending on wantarray
 Argument  : None
 Throws    : 
 Comments  : 
           : 

See Also   : org()

=cut

################################################## subroutine header end ##

sub dataset_groups
  {
    my $self = shift;
    my @dsgs =  $self->dataset->dataset_groups();
    return wantarray ? @dsgs : \@dsgs;
  }

################################################ subroutine header begin ##

=head2 organism

 Usage     : $returned_organism_object = $FeatureObject->organism();
 Purpose   : Shortcut to return an Organism object (name, description, normalized name) from a Feature object.
 Returns   : An Organism object.
 Argument  : None
 Throws    : 
 Comments  : 
           : 

See Also   : org()

=cut

################################################## subroutine header end ##

sub organism
  {
    my $self = shift;
    my ($dsg) =  $self->dataset->dataset_groups();
    return $dsg->organism;
  }


################################################ subroutine header begin ##

=head2 org

 Usage     : 
 Purpose   : Alias to the organism() method.
 Returns   : See organism()
 Argument  : None
 Throws    : 
 Comments  : 
           : 

See Also   : organism()

=cut

################################################## subroutine header end ##

sub org
  {
    my $self = shift;
    return $self->organism();
  }


################################################ subroutine header begin ##

=head2 names

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

sub names
  {
    my $self = shift;
    if ($self->{_names})
      {
	return wantarray ? @{$self->{_names}} : $self->{_names};
      }
#    my @names =  sort $self->feature_names()->get_column('name')->all;
    my @names;
    foreach my $name (sort {$a->name cmp $b->name} $self->feature_names())
      {
	if ($name->primary_name)
	  {
	    unshift @names, $name->name;
	  }
	else
	  {
	    push @names, $name->name;
	  }
      }
    $self->{_names}=\@names;
    return wantarray ? @names : \@names;
  }


################################################ subroutine header begin ##

=head2 primary_name

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

sub primary_name
 {
   my $self = shift;
   my ($nameo) = $self->feature_names({primary_name=>1});
   my ($name) = ref($nameo) =~ /name/i ? $nameo->name : $self->names;
 }


################################################ subroutine header begin ##

=head2 locs

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

sub locs
  {
    my $self = shift;
    return $self->locations();
  }


################################################ subroutine header begin ##

=head2 seqs

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

sub seqs
  {
    my $self = shift;
    return $self->sequences();
  }


################################################ subroutine header begin ##

=head2 eannotations

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

sub eannotations
  {
    my $self = shift;
    return $self->annotations(undef,{prefetch=>['annotation_type']});
  }


################################################ subroutine header begin ##

=head2 annos

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

sub annos
  {
    shift->eannotations(@_);
  }


################################################ subroutine header begin ##

=head2 length

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

sub length
    {
      my $self = shift;
      my $length = 0;
      map {$length+=($_->stop-$_->start+1)} $self->locations;
      return $length;
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
    my %opts = @_;
    my $gstid = $opts{gstid};
    my $anno_obj = new CoGe::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");
    my $start = $self->start;
    my $stop = $self->stop;
    my $chr = $self->chr;
    my $strand = $self->strand;
    #look into changing this to set_id
    my $info_id = $self->dataset->id;
    my $location = "Chr ".$chr." ";
    $location .= join (", ", map {$_->start."-".$_->stop} sort {$a->start <=> $b->start} $self->locs);
    $location .="(".$strand.")";
    #my $location = "Chr ".$chr. "".$start."-".$stop.""."(".$strand.")";
    $anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"Location", Values=>[$location], Type_delimit=>": ", Val_delimit=>" "));
    my $anno_type = new CoGe::Accessory::Annotation(Type=>"Name(s)");
    $anno_type->Type_delimit(": ");
    $anno_type->Val_delimit(" , ");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot($name);
      }
    
    $anno_obj->add_Annot($anno_type);
    foreach my $anno (sort {$b->type->name cmp $a->type->name} $self->annos)
      {
	my $type = $anno->type();
	my $group = $type->group();
	my $anno_type = new CoGe::Accessory::Annotation(Type=>$type->name);
	$anno_type->Val_delimit("\n");

	$anno_type->add_Annot($anno->annotation);
	if (ref ($group) =~ /group/i)
	  {
	    my $anno_g = new CoGe::Accessory::Annotation(Type=>$group->name);
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
    $anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<span class=\"title4\">Length</span>", Values=>[$self->length],Type_delimit=>": ", Val_delimit=>" "));
    my $ds = $self->dataset;
    my $org = $ds->organism->name;
    $org .= ": ".$ds->organism->description if $ds->organism->description;
    
    $anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"Organism", Values=>[$org], Type_delimit=>": ", Val_delimit=>" "));
    my ($gc, $at) = $self->gc_content(gstid=>$gstid);
    $gc*=100;
    $at*=100;
    $anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"DNA content", Values=>["GC: $gc%","AT: $at%"], Type_delimit=>": ", Val_delimit=>" "));
    my ($wgc, $wat) = $self->wobble_content(gstid=>$gstid);
    if ($wgc || $wat)
      {
	$wgc*=100;
	$wat*=100;
	$anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"Wobble content", Values=>["GC: $wgc%","AT: $wat%"], Type_delimit=>": ", Val_delimit=>" "));
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

See Also   : CoGe::Accessory::Annotation

=cut

################################################## subroutine header end ##


sub annotation_pretty_print_html
  {
    my $self = shift;
    my %opts = @_;
    my $loc_link = $opts{loc_link};
    my $minimal = $opts{minimal};
    my $gstid = $opts{gstid};
    $gstid = 1 unless defined $gstid;
    my $skip_GC = $opts{skip_GC};
    $loc_link = "FastaView.pl" unless defined $loc_link;
    my $anno_obj = new CoGe::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("<BR/>");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("<BR/>");
    my $start = $self->start;
    my $stop = $self->stop;
    my $chr = $self->chr;
    my $strand = $self->strand;
    my $dataset_id = $self->dataset->id;
    my $anno_type = new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">"."Name(s)"."</span>");
    $anno_type->Type_delimit(": <td>");
    $anno_type->Val_delimit(" , ");
    my $outname;
    foreach my $name ($self->names)
      {
	$outname = $name unless $outname;
	$anno_type->add_Annot("<span class=\"data5 link\" onclick=\"window.open('FeatView.pl?accn=".$name."');\">".$name."</span>");
      }
    $anno_obj->add_Annot($anno_type);
    unless ($minimal)
      {
	foreach my $anno (sort {uc($a->type->name) cmp uc($b->type->name)} $self->annos({},{prefetch=>{annotation_type=>'annotation_type_group'}}))
	  {
	    my $type = $anno->type();
	    my $group = $type->group();
	    my $anno_name = $type->name;
	    $anno_name .= ", ".$type->description if $type->description;
	    if (ref($group) =~ /group/i && !($type->name eq $group->name) )
	      {
		  {
		    $anno_name .= ":" unless $anno_name =~/:$/;
		    $anno_name = "<span class=\"title5\">". $anno_name."</span>";
		  }
	      }
	    else
	      {
		if ($anno->link)
		  {
		    $anno_name = "<tr><td nowrap='true'><span class=\"coge_link\">".$anno_name."</span>";
		  }
		else
		  {
		    $anno_name = "<tr><td nowrap='true'><span class=\"title5\">". $anno_name."</span>";
		  }
	      }
	    my $anno_type = new CoGe::Accessory::Annotation(Type=>$anno_name);
	    $anno_type->Val_delimit("<br>");
	    $anno_type->Type_delimit(" ");
#	    my $annotation = $anno->annotation;
	    my $annotation = "<span class=\"data5";
	    $annotation .= qq{ link" onclick="window.open('}.$anno->link.qq{')} if $anno->link;
	    $annotation .="\">".$anno->annotation."</span>";
	    $anno_type->add_Annot($annotation) if $anno->annotation;
	    if (ref ($group) =~ /group/i && !($type->name eq $group->name) )
	      {
		my $class = $anno->link ? "coge_link" : "title5";
		my $anno_g = new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"$class\">".$group->name."</span>");
		$anno_g->add_Annot($anno_type);
		$anno_g->Type_delimit(":<td>");
		$anno_g->Val_delimit(", ");
		#	    $anno_g->Val_delimit(" ");
		$anno_obj->add_Annot($anno_g);
	      }
	    else
	      {
		$anno_type->Type_delimit(":<td>");
		$anno_obj->add_Annot($anno_type);
	      }
	  }

	my $location = "Chr ".$chr." ";
	#    $location .= join (", ", map {$_->start."-".$_->stop} sort {$a->start <=> $b->start} $self->locs);
	$location .= $self->commify($self->start)."-".$self->commify($self->stop);
	$location .=" (".$strand.")";
	my $featid = $self->id;
	$anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">Length</span>", Values=>["<span class='data5'>".$self->length." nt</span>"],Type_delimit=>":<td>", Val_delimit=>" ")) unless $minimal;
#	$location = qq{<span class="data5 link" onclick="window.open('$loc_link?featid=$featid&start=$start&stop=$stop&chr=$chr&dsid=$dataset_id&strand=$strand&gstid=$gstid')" >}.$location."</span>" if $loc_link;
	$location = qq{<span class="data5 link" onclick="window.open('$loc_link?featid=$featid&gstid=$gstid')" >}.$location."</span>" if $loc_link;
	$location = qq{<span class="data">$location</span>};
	$anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5 link\"><span onclick=\"window.open('GenomeView.pl?chr=$chr&ds=$dataset_id&x=$start&z=5&gstid=$gstid')\" >Location</span></span>", Values=>[$location], Type_delimit=>":<td>", Val_delimit=>" "));

	my $ds=$self->dataset();
	my $dataset = qq{<span class="data5 link" onclick="window.open('OrganismView.pl?dsid=}.$ds->id."')\">".$ds->name;
#	$dataset .= ": ".$ds->description if $ds->description;
	$dataset .= "</span>";
	$anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">Dataset</span>", Values=>[$dataset], Type_delimit=>":<td>", Val_delimit=>" "));
	my $org = qq{<span class="data5 link" onclick = "window.open('OrganismView.pl?oid=}.$ds->organism->id."')\">".$ds->organism->name;
#	$org .= ": ".$ds->organism->description if $ds->organism->description;
	$org .= "</span>";
	
	$anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">Organism</span>", Values=>[$org], Type_delimit=>":<td>", Val_delimit=>" "));
	my $gst;
	foreach my $dsg ($ds->dataset_groups)
	  {
	    $gst = $dsg->type if $dsg->type->id == $gstid;
	  }
	if ($gst)
	  {
	    $anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">Genomic Sequnce</span>", Values=>["<span class=data5>".$gst->name."</span>"], Type_delimit=>":<td>", Val_delimit=>" "));
	  }
	unless ($skip_GC)
	  {
	    my ($gc, $at) = $self->gc_content(gstid=>$gstid);
	    $gc*=100;
	    $at*=100;
	    $anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">DNA content</span>", Values=>["<span class='data5'>GC: $gc%","AT: $at%</span>"], Type_delimit=>":<td>", Val_delimit=>" "));
	    my ($wgc, $wat) = $self->wobble_content(gstid=>$gstid);
	    if ($wgc || $wat)
	      {
		$wgc*=100;
		$wat*=100;
		$anno_obj->add_Annot(new CoGe::Accessory::Annotation(Type=>"<tr><td nowrap='true'><span class=\"title5\">Wobble content</span>", Values=>["<span class='data5'>GC: $wgc%","AT: $wat%</span>"], Type_delimit=>":<td>", Val_delimit=>" "));
	      }
	  }
      }
    return "<table cellpadding=0 class='ui-widget-content ui-corner-all'>".$anno_obj->to_String."</table>";
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
  #?
	# $comp = 1 if $loc->strand =~ "-";
	$comp = 1 if $loc->strand == "-1";
	$string .= "," if $count;
	$string .= $recal ? ($loc->start-$recal+1)."..".($loc->stop-$recal+1): $loc->start."..".$loc->stop;
	$count++;
      }
    $string = "join(".$string.")" if $count > 1;
    $string = "complement(".$string.")" if $comp;
    return $string;
  }


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

=head2 genomic_sequence

 Usage     : my $genomic_seq = $feat->genomic_sequence
 Purpose   : gets the genomic seqence for a feature
 Returns   : a string
 Argument  : none
 Comments  : This method simply creates a CoGe object and calls:
             get_genomic_sequence_for_feature($self)
See Also   : CoGe

=cut

################################################## subroutine header end ##

sub genomic_sequence {
  my $self = shift;
  my %opts = @_;
  my $up = $opts{up} || $opts{upstream} || $opts{left};
  my $down = $opts{down} || $opts{downstream} || $opts{right};
  my $debug=$opts{debug};
  my $gstid = $opts{gstid};
  my $seq = $opts{seq};
  my $dsgid = $opts{dsgid};
  #have a full sequence? -- pass it in and the locations will be parsed out of it!
  if (!$up && !$down && $self->_genomic_sequence)
    {
      return $self->_genomic_sequence;
    }

  my $dataset = $self->dataset();
  my @sequences;
  my %locs = map {($_->start,$_->stop)}$self->locations() ;#in case a mistake happened when loading locations and there are multiple ones with the same start
  my @locs = map {[$_, $locs{$_}]} sort {$a <=> $b} keys %locs;
  ($up,$down) = ($down, $up) if ($self->strand =~/-/); #must switch these if we are on the - strand;
  if ($up)
    {
      my $start = $locs[0][0]-$up;
      $start = 1 if $start < 1;
      $locs[0][0]=$start;
    }
  if ($down)
    {
      my $stop = $locs[-1][1]+$down;
      $locs[-1][1]=$stop;
    }
  my $chr = $self->chromosome;
  my $start = $locs[0][0];
  my $stop = $locs[-1][1];
  my $full_seq = $seq ? $seq : $dataset->get_genomic_sequence(
							      chromosome=>$chr,
							      start=>$start,
							      stop=>$stop,
							      debug=>$debug,
							      gstid=>$gstid,
							      dsgid=>$dsgid,
							     );
  if ($full_seq)
    {
      foreach my $loc (@locs)
	{
	  if ($loc->[0]-$start+$loc->[1]-$loc->[0]+1 > CORE::length ($full_seq))
	    {
	      print STDERR "#"x20,"\n";
	      print STDERR "Error in feature->genomic_sequence, Sequence retrieved is smaller than the length of the exon being parsed! \n";
	      print STDERR "Organism: ", $self->organism->name,"\n";
	      print STDERR "Dataset: ", $self->dataset->name,"\n";
	      use Data::Dumper;
	      print STDERR "Locations data-structure: ", Dumper \@locs;
	      print STDERR "Retrieved sequence lenght: ", CORE::length ($full_seq),"\n";
	      print STDERR $full_seq,"\n";
	      print STDERR "Feature object information: ",Dumper {
		chromosome=>$chr,
		  skip_length_check=>1,
		    start=>$start,
		      stop=>$stop,
			dataset=>$dataset->id,
			  feature=>$self->id,
			};
	      print STDERR "#"x20,"\n";
	    }
	  
	  my $sub_seq = substr($full_seq, $loc->[0] - $start, $loc->[1] - $loc->[0] + 1);
	  next unless $sub_seq;
	  if ($self->strand == -1){
	    unshift @sequences, $self->reverse_complement($sub_seq);
	  }else{
	    push @sequences, $sub_seq;
	  }
	}      
    }
  my $outseq = join( "", @sequences );
  if (!$up && !$down)
    {
      $self->_genomic_sequence($outseq);
    }
  return $outseq;
}


################################################ subroutine header begin ##

=head2 genome_sequence

 Usage     : 
 Purpose   : See genomic_sequence()
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : Alias for the genomic_sequence() method.

See Also   : genomic_sequence()

=cut

################################################## subroutine header end ##

sub genome_sequence
  {
   shift->genomic_sequence(@_);
  }


################################################ subroutine header begin ##

=head2 has_genomic_sequence

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

sub has_genomic_sequence
  {
    my $self = shift;
    return 1 if $self->dataset->has_genomic_sequence;
    return 0;
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
    my $gstid = $opts{gstid};
    my $match = $opts{match} || 1;
    my $mismatch = $opts{mismatch} || -3;
    my $seq_length = $opts{seq_length} || $opts{length}; #need a sequence length for calculation
    my $lambda = $self->_estimate_lambda(match=>$match, mismatch=>$mismatch);
    my $seq = $self->genomic_sequence(gstid=>$gstid) unless $seq_length; #get sequence for feature if no seq_length has been passed in;
    $seq_length = CORE::length ($seq) if $seq;
    warn "No genomic sequence could be obtained for this feature object.  Can't calculate a blast bit score.\n" unless $seq_length;
    my $bs = sprintf("%.0f", $lambda*$seq_length*$match/log(2));
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


################################################ subroutine header begin ##

=head2 reverse_complement

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

sub reverse_complement
  {
    my $self = shift;
    my $seq = shift;# || $self->genomic_sequence;
    if (ref($self) =~ /Feature/)
      {
	$seq = $self->genomic_sequence unless $seq; #self seq unless we have a seq
      }
    else #we were passed a sequence without invoking self
      {
	$seq = $self unless $seq;
      }
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/; 
    return $rcseq;
  }


################################################ subroutine header begin ##

=head2 reverse_comp

 Usage     : 
 Purpose   : See reverse_complement()
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : Alias for the reverse_complement function.
           : 

See Also   : reverse_complement()

=cut

################################################## subroutine header end ##

sub reverse_comp {
  shift->reverse_complement(@_);
}


################################################ subroutine header begin ##

=head2 protein_sequence

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

sub protein_sequence {
  my $self = shift;
  my %opts = @_;
  my $gstid = $opts{gstid};
  my $siter = $self->sequences();
  my @sequence_objects;
  while ( my $seq = $siter->next() ) {
    push @sequence_objects, $seq;
  }

  if (@sequence_objects == 1) 
    {
      return $sequence_objects[0]->sequence_data();
    } 
  elsif ( @sequence_objects > 1 )  
    {
      return wantarray ? @sequence_objects : \@sequence_objects;
    } 
  else 
    {
      my ($seqs,$type) = $self->frame6_trans(gstid=>$gstid);
      #check to see if we can find the best translation
      my $found=0;
      my @seqs;
      while (my ($k, $v) = each %$seqs)
	{
	  next unless $v;
	  push @seqs, $v;
	  if (($v =~ /^M/ || $v =~ /^L/ )&& $v =~ /\*$/) #need the L for some microbial start positions
	    {
	      next if $v =~ /\*\w/;
	      $found = $k;
	    }
	}
      if ($found)
	{
	  return $seqs->{$found};
	}
      else
	{
	  return wantarray ? @seqs : \@seqs;
	}
    }
}


################################################ subroutine header begin ##

=head2 frame6_trans

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

sub frame6_trans
  {
    my $self = shift;
    my %opts = @_;
    my $trans_type = $opts{trans_type};
    my $gstid = $opts{gstid};
    my $code;
    ($code, $trans_type) = $opts{code} || $self->genetic_code(trans_type=>$trans_type);
    my $seq = $opts{seq} || $self->genomic_sequence(gstid=>$gstid);

    my %seqs;
    $seqs{"1"} = $self->_process_seq(seq=>$seq, start=>0, code1=>$code, codonl=>3);
    $seqs{"2"} = $self->_process_seq(seq=>$seq, start=>1, code1=>$code, codonl=>3);
    $seqs{"3"} = $self->_process_seq(seq=>$seq, start=>2, code1=>$code, codonl=>3);
    my $rcseq = $self->reverse_complement($seq);
    $seqs{"-1"} = $self->_process_seq(seq=>$rcseq, start=>0, code1=>$code, codonl=>3);
    $seqs{"-2"} = $self->_process_seq(seq=>$rcseq, start=>1, code1=>$code, codonl=>3);
    $seqs{"-3"} = $self->_process_seq(seq=>$rcseq, start=>2, code1=>$code, codonl=>3);
    return \%seqs, $trans_type;

  }


################################################ subroutine header begin ##

=head2 genetic_code

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

sub genetic_code
  {
    my $self = shift;
    my %opts = @_;
    my $trans_type = $opts{trans_type};
    unless ($trans_type)
      {
	foreach my $anno ($self->annotations)
	  {
	    next unless $anno->annotation_type;
	    next unless $anno->annotation_type->name eq "transl_table";
	    $trans_type = $anno->annotation;
	  }
      }
    $trans_type = 1 unless $trans_type;
    my $code = code($trans_type);
    return ($code->{code}, $code->{name});
  }
  
  
################################################ subroutine header begin ##

=head2 _process_seq

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
  
sub _process_seq
  {
    my $self = shift;
    my %opts = @_;
    my $seq = $opts{seq};
    my $start = $opts{start};
    my $code1 = $opts{code1};
    my $code2 = $opts{code2};
    my $alter = $opts{alter} || "X";
    my $codonl = $opts{codonl} || 2;
    my $seq_out;
    for (my $i = $start; $i < CORE::length ($seq); $i = $i+$codonl)
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


################################################ subroutine header begin ##

=head2 percent_translation_system

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

sub percent_translation_system
  {
    my $self = shift;
    my %opts = @_;
    my $counts = $opts{counts};
    my %code1 = (
		 "W"=>1,
		 "M"=>1,
		 "L"=>1,
		 "Y"=>1,
		 "C"=>1,
		 "I"=>1,
#		 "K"=>1, #majority in code2
		 "R"=>1,
		 "Q"=>1,
		 "V"=>1,
		 "E"=>1,
		 
		);
    my $code1 = $opts{code1} || \%code1;
    my %code2 = (
		 "F"=>1,
		 "S"=>1,
		 "T"=>1,
		 "N"=>1,
		 "P"=>1,
		 "H"=>1,
		 "A"=>1,
		 "D"=>1,
		 "G"=>1,
		 "K"=>1,
		);
    my $code2 = $opts{code2} || \%code2;
    my ($seq) = $opts{seq} || $self->protein_sequence;
    return (0,0) unless $seq;
    my ($c1, $c2, $total) = (0,0,0);
    foreach (split //, $seq)
      {
	$_ = uc($_);
	$c1++ if $code1->{$_};
	$c2++ if $code2->{$_};
	$total++;
      }
    if ($counts)
      {
	return $c1, $c2, $total;
      }
    else
      {
	return (map {sprintf("%.4f", $_)} $c1/$total, $c2/$total);
      }
  }


################################################ subroutine header begin ##

=head2 aa_frequency

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

sub aa_frequency
  {
    my $self = shift;
    my %opts = @_;
    my $counts = $opts{counts};
    my $gstid = $opts{gstid};
    my $code = $opts{code};
    ($code) = $self->genetic_code unless $code;
    my %data = map {$_=>0} values %$code;
    my ($seq) = $opts{seq} || $self->protein_sequence(gstid=>$gstid);
    return \%data unless $seq;
    foreach (split //,$seq)
      {
	next if $_ eq "*";
	$data{$_}++ if defined $data{$_};
      }
    if ($counts)
      {
	return \%data;
      }
    else
      {
	my $total = 0;
	foreach (values %data)
	  {
	    $total+=$_;
	  }
	foreach my $aa (keys %data)
	  {
	    $data{$aa} = sprintf("%.4f", ($data{$aa}/$total));
	  }
	return \%data;
      }
  }


################################################ subroutine header begin ##

=head2 codon_frequency

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

sub codon_frequency
  {
    my $self = shift;
    my %opts = @_;
    my $counts = $opts{counts};
    my $code = $opts{code};
    my $code_type = $opts{code_type};
    my $gstid = $opts{gstid};
    ($code, $code_type) = $self->genetic_code unless $code;;
    my %codon = map {$_=>0} keys %$code;
    my $seq = $self->genomic_sequence(gstid=>$gstid);
    my $x=0;
    while ($x<CORE::length($seq))
      {
	$codon{uc(substr($seq, $x, 3))}++;
	$x+=3;
      }
    if ($counts)
      {
	return \%codon, $code_type;
      }
    else
      {
	my $total = 0;
	foreach (values %codon)
	  {
	    $total+=$_;
	  }
	foreach my $codon (keys %codon)
	  {
	    $codon{$codon} = sprintf("%.4f", ($codon{$codon}/$total));
	  }
	return (\%codon, $code_type);
      }
  }


################################################ subroutine header begin ##

=head2 gc_content

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

sub gc_content
  {
    my $self = shift;
    my %opts = @_;
    my $counts = $opts{counts};
    my $gstid = $opts{gstid}; #genomic_sequence_type_id
    my $seq = $self->genomic_sequence(gstid=>$gstid);
    my ($gc,$at, $n);
    $gc = $seq =~ tr/gcGC/gcGC/;
    $at = $seq =~ tr/atAT/atAT/;
    $n = $seq =~ tr/nxNX/nxNX/;
    unless ($counts)
      {
	my $total = CORE::length($seq);
	return (0,0,0) unless $total;
	$gc = sprintf("%.4f", ($gc/$total));
	$at = sprintf("%.4f", ($at/$total));
	$n = sprintf("%.4f", ($n/$total));
      }
    return $gc,$at, $n;
  }


################################################ subroutine header begin ##

=head2 percent_gc

 Usage     : 
 Purpose   : See gc_content()
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : Alias for the gc_content() method.
           : 

See Also   : gc_content()

=cut

################################################## subroutine header end ##

sub percent_gc
  {
    shift->gc_content(@_);
  }


################################################ subroutine header begin ##

=head2 wobble_content

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

sub wobble_content
  {
    my $self = shift;
    my %opts = @_;
    my $counts = $opts{counts};
    my $gstid = $opts{gstid}; #genomic_sequence_type_id
    return unless $self->type->name =~ /cds/i;
    my $seq = $self->genomic_sequence(gstid=>$gstid);
    my $codon_count=0;;
    my $at_count=0;
    my $gc_count=0;
    my $n_count=0;
    for (my $i =0; $i < CORE::length($seq); $i+=3)
      {
        my $codon = substr ($seq, $i, 3);
        $codon_count++;
        my ($wobble) = $codon =~ /(.$)/;
        $at_count++ if $wobble =~ /[at]/i;
        $gc_count++ if $wobble =~ /[gc]/i;
        $n_count++ if $wobble =~ /[nx]/i;
      }
    return ($gc_count, $at_count, $n_count) if $counts;
    my $pat = sprintf("%.4f", $at_count/$codon_count) if $codon_count;
    my $pgc = sprintf("%.4f", $gc_count/$codon_count) if $codon_count;
    my $pn = sprintf("%.4f", $n_count/$codon_count) if $codon_count;
    return ($pgc, $pat, $pn);
  }


################################################ subroutine header begin ##

=head2 fasta

 Usage     : 
 Purpose   : returns a fasta formated sequence for the featre
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub fasta
  {
    my $self = shift;
    my %opts = @_;
    my $col;
    $col = $opts{col};
    my $prot = $opts{protein} || $opts{prot};
    #$col can be set to zero so we want to test for defined variable
    $col = $opts{column} unless defined $col;
    $col = $opts{wrap} unless defined $col;
    $col = 100 unless defined $col;
    my $rc = $opts{rc};
    my $upstream = $opts{upstream} || 0;
    my $downstream = $opts{downstream} || 0;
    my $name_only = $opts{name_only};
    my $add_fid = $opts{add_fid};
    my $sep = $opts{sep}; #returns the header and sequence as separate items.
    my $gstid = $opts{gstid};
    my $gst;
    if ($gstid)
      {
	foreach my $dsg ($self->dataset->dataset_groups)
	  {
	    ($gst) = $dsg->type if $dsg->type->id == $gstid;
	  }
      }
    if (!$gst)
      {
	($gst) = $self->sequence_type();
	$gstid = $gst->id if $gst;
      }
    my ($pri_name) = $self->primary_name;
    my $head = ">";
    if ($name_only)
      {
	$head = "fid:".$self->id." ".$head if $add_fid;
	$head .= $pri_name;
      }
    else
      {
	$head .= $self->dataset->organism->name."(v".$self->version;
	$head .= ", ".$gst->name if $gst;
	$head .= ")".", Name: ".(join (", ", $self->names)).", Type: ".$self->type->name.", Feature Location: (Chr: ".$self->chromosome.", ".$self->genbank_location_string.")";
	$head .= "fid:".$self->id if $add_fid;
	$head .= " +up: $upstream" if $upstream;
	$head .= " +down: $downstream" if $downstream;
	$head .= " (reverse complement)" if $rc;
      }
    my ($start, $stop) = ($self->start, $self->stop);
    if ($rc) 
      {
	$start -= $downstream;
	$stop += $upstream;
      }
    else
      {
	$start -= $upstream;
	$stop += $downstream;
      }

    $head .= " Genomic Location: $start-$stop" unless $name_only;
    $Text::Wrap::columns=$col;
    my $fasta;
    if ($prot)
      {
	foreach my $seq ($self->protein_sequence(gstid=>$gstid))
	  {
	    $seq = join ("\n", wrap("","",$seq)) if $col;
	    $fasta .= $head."\n".$seq."\n";
	  }
      }
    else
      {
	my $seq = $self->genomic_sequence(upstream=>$upstream, downstream=>$downstream, gstid=>$gstid);
	$seq = $self->reverse_complement($seq) if $rc;
	$seq = join ("\n", wrap("","",$seq)) if $col;
	$fasta = $head."\n".$seq."\n";
	return $head, $seq if ($sep);
      }
    return $fasta;
  }


################################################ subroutine header begin ##

=head2 sequence_type

 Usage     : 
 Purpose   : returns the genomic_sequence_type object for the sequence
 Returns   : wantarray -- may be more than one genomic_sequence_type sequences associated with this feature
             looked up through dataset->dataset_group_connector->dataset_group->genomic_sequence_type
 Argument  : 
 Throws    : 
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub sequence_type
  {
    my $self = shift;
    my %gst;
    foreach my $dsg ($self->dataset->dataset_groups)
      {
	$gst{$dsg->genomic_sequence_type->id} = $dsg->genomic_sequence_type;
      }
    my @gst = values %gst;
    return shift @gst if scalar @gst == 1; #only one result
    return wantarray ? @gst : \@gst;
  }


################################################ subroutine header begin ##

=head2 commify

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

sub commify 
    {
      my $self = shift;
      my $input = shift;
      $input = reverse $input;
      $input =~ s<(\d\d\d)(?=\d)(?!\d*\.)><$1,>g;
      return scalar reverse $input;
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
