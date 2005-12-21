package CoGe::Genome::DB::Feature;
use strict;
use base 'CoGe::Genome::DB';
use CoGe::Genome::Accessory::Annotation;

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

1; #this line is important and will help the module return a true value

