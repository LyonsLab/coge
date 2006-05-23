package CoGe::Genome::DB::Data_information;
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
    __PACKAGE__->table('data_information');
    __PACKAGE__->columns(All=>qw{data_information_id data_source_id organism_id name description link version});
    __PACKAGE__->has_a(data_source_id=>'CoGe::Genome::DB::Data_source');
    __PACKAGE__->has_a('organism_id'=>'CoGe::Genome::DB::Organism');
    __PACKAGE__->has_many(features=>'CoGe::Genome::DB::Feature');
    __PACKAGE__->has_many(genomic_sequences=>'CoGe::Genome::DB::Genomic_sequence');

    __PACKAGE__->set_sql(get_chromosomes=>qq{
SELECT DISTINCT l.chromosome
 FROM location l
 JOIN feature f using (feature_id)
 JOIN data_information di using (data_information_id)
WHERE di.data_information_id = ?

});
}


########################################### main pod documentation begin ##
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

CoGe::Genome::DB::Data_source

=head1 SYNOPSIS

  use Genome::DB;

  my $db = new Genome::DB;

  #let's say we want to make a fasta file of all CDS features for Arabidopsis

  #first, let's find the correct organism object

  my ($org) = $db->get_org_obj->search_like(name=>'Arabidopsis%');

  #let's find all the data_information objects for arabidopsis version 2

  foreach my $di ($db->get_data_info_obj->search(organism_id=>$org->id,version=>2)
    {
      #let's go thru all the features for each data information object
      foreach my $feat ($di->feats)
        {
           #skip any that aren't of type CDS
           next unless $feat->type->name eq "CDS";
           #print the fasta header of names of the feature
           print ">",join (" ", map {$_->name} $feat->names),"\n";
           #print the sequence
           print $db->get_genomic_sequence_for_feature($feat),"\n";
        }
    }

#done!

=head1 DESCRIPTION

Data information objects contain information about the origins of genomic data
such as the file name of the data, a URL link to where it was obtained, the 
version of the data, etc.  The reason that this is separate from the notion of 
a "data source" is that a given data source (e.g. TIGR) can release multiple 
files per genome (e.g. one file per chromosome) , different genomes, and 
different versions of any given genome.  In order to accommodate such variation
in genomic data sources, this table provides a way to track individual data
sets from a given data source.  After writing the database, this module, 
all related modules, many end-user applications, I realize that a better
name for this object would have been data_set.  However, I'm much too busy
with important stuff (tm) to go back and recode all of this with the appropriate
name.  If this is problematic for you, cope.

In terms of functionality, this object is one of the most important for
figuring out which features or other genomic information is needed.  For 
example, let's say that the database has two genomes, Arabidopsis (aka 
Shitty Little Weed, SLW) and Mus musculus (aka mouse) and has two versions 
of each genome (version 1 and version 2).  Each genome also had an individual 
data file for each chromosome (that means 5 for each version of Arabidopsis and 
21 for each version of mouse).  That means that for this example , there are 
52 data information entries in our database!  Understanding if a user is 
interested in chromosome 1 of Arabidopsis version 2 can make our database searches
faster as well as given us ( the programmers ) the power of tracking several
versions of any given genome.



This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.

The columns for this table are:
 data_information_id
 name
 description
 link
 version
 organism_id

Related objects that be accessed through this object are:
 CoGe::Genome::DB::Feature
 CoGe::Genome::DB::Organism
 CoGe::Genome::DB::Genomic_sequence
 CoGe::Genome::DB::Data-source

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
 CoGe::Genome::DB::Data::Source
 CoGe::Genome::DB::Organism
 CoGe::Genome::DB::Genomic_sequence
 Class::DBI


perl(1).

=cut

############################################# main pod documentation end ##


=head2 Accessor Functions



new                 =>  creates a new object (inherited from Class::Accessor)

data_information_id =>  database entry id
id                  =>  alias for location_id

name                =>  name, usually the file name of the data set

description         =>  description of the data set
desc                =>  alias for description

link                =>  web link or similar for where the data was gotten
url                 =>  alias for link

version             =>  version number or similar of the data
ver                 =>  alias for version
v                   =>  alias for version

data_source_id      =>  returns the data source object associated with this data set
data_source         =>  alias for above
source              =>  alias for above

organism_id         =>  returns the organism object for this data set
organism            =>  alias for above
org                 =>  alias for above
species             =>  alias for above

genomic_sequence_id =>  returns all the genomic sequence objects associated with this data set
genomic_seqs        => alias for above
seqs                => alias for above

features            =>  returns the related feature objects associated with this data set
feature             =>  alias for features
feats               =>  alias for features

=cut


sub data_source
  {
    my $self = shift;
    return $self->data_source_id();
  }

sub source
  {
    my $self = shift;
    return $self->data_source_id();
  }

sub feature
  {
    my $self = shift;
    return $self->features();
    }

sub feats
  {
    my $self = shift;
    return $self->features();
    }

sub genomics_seqs
  {
    my $self = shift;
    return $self->genomics_sequences();
  }

sub seqs
  {
    my $self = shift;
    return $self->genomics_sequences();
  }

sub desc
  {
    my $self = shift;
    return $self->description(@_);
  }

sub id
  {
    my $self = shift;
    return $self->data_information_id();
  }

sub organism
  {
    my $self = shift;
    return $self->organism_id();
  }

sub org
  {
    my $self = shift;
    return $self->organism_id();
  }

sub species
  {
    my $self = shift;
    return $self->organism_id();
  }

sub url
  {
    my $self = shift;
    return $self->link();
  }

sub v
  {
    my $self = shift;
    return $self->version();
  }

sub ver
  {
    my $self = shift;
    return $self->version();
  }


################################################ subroutine header begin ##

=head2 get_associated_data_infos

 Usage     : my @di = $data_info_obj->get_associated_data_infos()
 Purpose   : This will find other data information objects that might be important given the current 
             data information object.  For example, let's say you loaded a genome where each chromosome
             was loaded separately so each has its own data_information object.  Then you load another
             set of data for this gene such as annotations from an HMM profile search.  However, this
             data-set covers the entire genome and not any individual chromosome.  So, if you want to
             find if a gene on some chromosome also has an HMM domain annotated, you will need to 
             1. find the gene's data information object
             2. find any other data information objects that might contain relevant information (using this function)
             3. assemble that data together
 Returns   : an array or array ref of CoGe::Genome::DB::Data_information objects depending on value of wantarray.
             This array will include the Data_information objects used in the search
 Argument  : accepts one or more data_information objects or uses self if none are specified.
 Throws    : None
 Comments  : The mechanism by which this routine works is to:
             1. find the organism, version, and chromosome(s) for the give data_information objects(s)
             2. search for other data_information objects with the same organism, version, and chromosome
           : NOTE: Chromosome search has been disabled due to speed constraints.
See Also   : 

=cut

################################################## subroutine header end ##


sub get_associated_data_infos
  {
    my $self = shift;
    my $db = new CoGe::Genome;
    my @dio = @_;
    push @dio, $self unless @dio;
    my %ids;
    foreach my $di (@dio)
      {
	$ids{$di->id} = $di;
#	my (%chr) = map {$_, 1} $di->get_chromosomes;
	foreach my $tdi ($self->search({organism_id=>$di->org->id, version=>$di->version}))
	  {
# 	    my $pass;
# 	    foreach my $chr ($di->get_chromosomes)
# 	      {
# 		if ($chr{$chr})
# 		  {
# 		    $pass = 1;
# 		    last;
# 		  }
# 	      }
# 	    next unless $pass;
	    $ids{$tdi->id}=$tdi;
	  }
      }
    my @new_dios = values %ids;
    return wantarray ? @new_dios : \@new_dios;
  }

################################################ subroutine header begin ##

=head2 get_chromosomes

 Usage     : my @chrs = $data_information_obj->get_chromosomes
 Purpose   : finds all the chromosome names for one or more data information objects
 Returns   : an array or array ref of strings (depending on wantarray)
 Argument  : accepts an array of data information objects or database ids.  
             if none ore specified, it will use itself
 Throws    : none
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_chromosomes
  {
    my $self = shift;
    my @di = @_;
    push @di, $self unless @di;
    my %chrs;
    my $sth = $self->sql_get_chromosomes();
    foreach my $di (@di)
      {
	my $id = ref ($di) =~ /information/i ? $di->id : $di;

	$sth->execute($id);
	my $q = $sth->fetch();
	next unless $q;
	foreach my $i (@$q)
	  {
	    $chrs{$i} = 1;
	  }
      }
    $sth->finish;
    return wantarray ? keys %chrs : [keys %chrs];
  }

1; #this line is important and will help the module return a true value

