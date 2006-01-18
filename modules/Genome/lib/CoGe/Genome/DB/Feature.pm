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
#    __PACKAGE__->has_a('organism_id'=>'CoGe::Genome::DB::Organism');
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
   AND di.version = ?
});
    __PACKAGE__->set_sql ('select_features_in_range' => qq{
SELECT DISTINCT f.feature_id
  FROM feature f
  JOIN location l USING (feature_id)
 WHERE ? <= l.stop
   AND ? >= l.start
   AND chromosome = ?
   AND f.data_information_id = ?;
});
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Feature - Genome::DB::Feature

=head1 SYNOPSIS

  use Genome::DB::Featur
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

##legacy table, now stored in data_information
sub organism
  {
    my $self = shift;
    return $self->info->organism_id();
  }

sub org
  {
    my $self = shift;
    return $self->info->organism_id();
  }

sub species
  {
    my $self = shift;
    return $self->info->organism_id();
  }


sub feat_names
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

sub locs
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

sub annotation_pretty_print
  {
    my $self = shift;
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Val_delimit("\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n");
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"Location", Values=>["Chr ".$self->chr, "".$self->begin_location."-".$self->end_location.""."(".$self->strand.")"], Type_delimit=>"\n\t", Val_delimit=>" "));
    foreach my $anno ($self->annos)
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
	    $anno_g->Type_delimit("\n\t");
	    $anno_g->Val_delimit(" ");
	    $anno_obj->add_Annot($anno_g);
	  }
	else
	  {
	    $anno_type->Type_delimit("\n\t");
	    $anno_obj->add_Annot($anno_type);
	  }
      }
    my $anno_type = new CoGe::Genome::Accessory::Annotation(Type=>"Name(s)");
    $anno_type->Type_delimit("\n\t");
    $anno_type->Val_delimit("\n\t");
    foreach my $name ($self->names)
      {
	$anno_type->add_Annot($name->name);
      }
    
    $anno_obj->add_Annot($anno_type);
    return $anno_obj->to_String;
  }

sub annotation_pretty_print_html
  {
    my $self = shift;
    my $anno_obj = new CoGe::Genome::Accessory::Annotation(Type=>"anno");
    $anno_obj->Val_delimit("\n<BR>\n");
    $anno_obj->Add_type(0);
    $anno_obj->String_end("\n<BR>\n");
    $anno_obj->add_Annot(new CoGe::Genome::Accessory::Annotation(Type=>"<font class=\"annotation\">Location</font>", Values=>["Chr ".$self->chr, "".$self->begin_location."-".$self->end_location.""."(".$self->strand.")"], Type_delimit=>"\n<BR><li>", Val_delimit=>" "));
    foreach my $anno ($self->annos)
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
    return $anno_obj->to_String;
  }



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

sub begin_location
  {
    my $self = shift;
    my @locs = $self->locs;
    return unless @locs;
    my ($val) = sort {$a->begin <=> $b->begin} @locs;
    return $val->begin;
  }

sub end_location
  {
    my $self = shift;
    my @locs = $self->locs;
    return unless @locs;
    my ($val) = sort {$b->end <=> $a->end} @locs;
    return $val->end;
  }

sub chromosome
  {
    my $self = shift;
    my $loc = $self->locs->next;
    return unless $loc;
    return $loc->chr();
  }

sub chr
  {
    my $self = shift;
    return $self->chromosome;
  }

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
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    my $chr = $opts{chr} || $opts{CHR} || $opts{chromosome} || $opts{CHROMOSOME};
    my $info_id = $opts{info_id} || $opts{INFO_ID} || $opts{data_info_id} || $opts{DATA_INFO_ID};
    my $sth = $self->sql_select_features_in_range();
    $sth->execute($start, $stop, $chr, $info_id);
    my @feats;
    while (my $q = $sth->fetch())
      {
	push @feats, $self->search(feature_id=>$q->[0]);
      }
    return wantarray ? @feats : \@feats;
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

1; #this line is important and will help the module return a true value


