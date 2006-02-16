package CoGe::Genome::DB::Feature_name;
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
    __PACKAGE__->table('feature_name');
    __PACKAGE__->columns(All=>qw{feature_name_id name description feature_id});
    __PACKAGE__->has_a(feature_id=>'CoGe::Genome::DB::Feature');
    __PACKAGE__->set_sql(search_name=>qq{
SELECT name
  FROM feature_name
 WHERE name like ?
});
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Genome::DB::Feature_name - Genome::DB::Feature_name

=head1 SYNOPSIS

 #feature name objects are usually obtained through a feature object and contain a single name
 #and description of that name.  A given feature, such as a gene, may have many names and each name
 #is stored as a seperate entry in the feature_name table
 #Note:  In many cases, the description field is empty because a name can say all there is to say.  

 use CoGe::Genome;
 #create the master genome database object
 my $db = CoGe::Genome->new;
 
 #using the master database object, let's find some features.  We'll use the method get_feature_obj
 #to get a Feature object and then use the method get_features_in_region to find all the features 
 #in a particular region for a given data information source.  Data information ids are database ids.
 #

 foreach my $feat ($db->get_feature_obj->get_features_in_region(start=>10000, stop=>50000, chromosome=>3, info_id=>7))
   {
     #let's print some information about the feature
     print "Feature Name:  ", join (", ", map {$_->name} $feat->names),"\n";
     #let's print the location(s) for the feature
     print "\t", join ("\n\t", map {$_->start."-".$_->stop.": ".$_->chr} $feat->locs),"\n";
     print "\n";
   }

=head1 DESCRIPTION

 The feature name table in the genomes database stores the name(s) and description of said names(s)
 for a feature.  A genomic feature is more or less defined as a region on a chromosome that has
 associated information about it.  Some common genomic features are genes, transcripts, CDS (coding
 sequences), tRNAs, etc.  Since any feature may have one or more names, there may be one or more
 names associated with a feature.

This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.

The columns for this table are:
 feature_name_id
 name
 description
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
 CoGe::Genome::Feature
 Class::DBI

perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions



new              =>  creates a new object (inherited from Class::Accessor)

feature_name_id  =>  database entry id
id               =>  alias for location_id

name             =>  name of feature

description      =>  description of feature name
desc             =>  alias for description

feature_id       =>  returns the related feature object associated with this name
feature          =>  alias for feature_id
feat             =>  alias for feature_id

=cut



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

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->feature_name_id();
  }

################################################ subroutine header begin ##

=head2 search_name

 Usage     : my @names = $feat_name_obj->search_name($name);
 Purpose   : fetches all the feature names using a "like" call to the database
           : THIS DOES NOT RETURN feature_name objects!
 Returns   : an a array or arrayref based on wantarray
 Argument  : a string for a database search using a "like" call
           : for example, the wild card for a like call is "%" so
             a search for "at1g01%" will return anything in the databaes
             that begins with "at1g01".  Likewise, if you wish to search for
             somethat that contains "1g" you would submit the string "%1g%"
 Throws    : undef if nothing is passed in.
 Comments  : This routines searches the database for names and returns
           : them without creating feature_name objects.  Thus, this routine
             is much quicker than the usual one supplied by Class::DBI
             E.g.: $feat_name_obj->search_like(name=>'at1g01%')

See Also   : Class::DBI

=cut

################################################## subroutine header end ##


sub search_name
  {
    my $self = shift;
    my $name = shift;
    return unless $name;
    my $sth = $self->sql_search_name();
    $sth->execute($name);
    my @names;
    while (my $q = $sth->fetch)
      {
	push @names, $q->[0];
      }
    return wantarray? @names : \@names;
  }

################################################ subroutine header begin ##

=head2 sample_function

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

