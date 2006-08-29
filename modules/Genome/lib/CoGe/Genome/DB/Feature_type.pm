package CoGe::Genome::DB::Feature_type;
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
    __PACKAGE__->table('feature_type');
    __PACKAGE__->columns(All=>qw{feature_type_id name description});
    __PACKAGE__->has_many(features=>'CoGe::Genome::DB::Feature');
 }


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Feature_type - CoGe::Genome::DB::Feature_type

=head1 SYNOPSIS

  
  #feature_types are associated with features and store the what kind of feature 
  #it is.  Common examples are Gene, mRNA, CDS, etc"

  use CoGe::Genome;
  #create the master genome database object
  my $db = CoGe::Genome->new;
 
  #using the master database object, let's find some features
  foreach my $feat ($db->get_features_by_name("GenePoo"))
   {
     #print out the name(s) of the feature.  Names are stored in a separate table and 
     #and are accessible through the feature object
     print join (", ", map {$_->name} $feat->names)

     #print out the type of the feature as well.  Each feature is of one and only one type.  
     #Featuer types are also stored in a separate table and are accessible through the feature object
     #(Three cheers for Class::DBI!)
     print " ", $feat->type->name,"\n";
   }


=head1 DESCRIPTION

The feature_type table in the genomes database soters the type of a feature.  A
genomic feature can be defined as a region on a chromosome that has associated 
information (such as a gene).  There are relatively few types of features
(when compared to features themselves) and some common examples would be:
Gene, CDS, mRNA, tRNA, snoRNA, etc.

This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.


The columns for this table are:
 feature_type_id
 name
 description

Related objects that can be accessed through this object are:
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

feature_type_id  =>  database entry id
id               =>  alias for location_id

name             =>  name of feature type

description      =>  description of feature type
desc             =>  alias for description

features         =>  returns an array of CoGe::Genome::DB::Feature objects or
                     a Class::DBI interator
feats            =>  alias for features

=cut


sub feats
  {
    my $self = shift;
    return $self->features();
  }

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->feature_type_id();
  }


1; #this line is important and will help the module return a true value

