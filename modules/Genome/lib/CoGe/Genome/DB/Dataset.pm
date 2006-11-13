package CoGe::Genome::DB::Dataset;
use strict;
use base 'CoGe::Genome::DB';
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
    __PACKAGE__->table('dataset');
    __PACKAGE__->columns(All=>qw{dataset_id data_source_id organism_id name description link version});
    __PACKAGE__->has_a(data_source_id=>'CoGe::Genome::DB::Data_source');
    __PACKAGE__->has_a('organism_id'=>'CoGe::Genome::DB::Organism');
    __PACKAGE__->has_many(features=>'CoGe::Genome::DB::Feature');
    __PACKAGE__->has_many(genomic_sequences=>'CoGe::Genome::DB::Genomic_sequence');

    __PACKAGE__->set_sql(get_chromosomes=>qq{
SELECT DISTINCT l.chromosome
 FROM location l
 JOIN feature f using (feature_id)
 JOIN dataset ds using (dataset_id)
WHERE ds.dataset_id = ?

});
    __PACKAGE__->set_sql(get_feature_type_count=>qq{
SELECT count(distinct(feature_id)), ft.name 
  FROM feature
  JOIN feature_type ft using (feature_type_id)
 WHERE dataset_id = ?
 GROUP BY ft.name
});
    __PACKAGE__->set_sql(get_feature_type_count_chr=>qq{
SELECT count(distinct(feature_id)), ft.name 
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN location l USING (feature_id)
 WHERE dataset_id = ?
   AND l.chromosome = ?
 GROUP BY ft.name
});
    __PACKAGE__->set_sql(get_current_version_for_organsim=>qq{
SELECT DISTINCT(version), chromosome
  FROM dataset
  JOIN genomic_sequence USING (dataset_id)
 WHERE organism_id = ?
 ORDER BY version DESC
});
  
   __PACKAGE__->set_sql(get_features_for_organism => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
 WHERE organism_id = ?
});
  
   __PACKAGE__->set_sql(get_features_for_organism_version => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
 WHERE organism_id = ?
   AND version = ?
});

   __PACKAGE__->set_sql(get_features_for_organism_version_type => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
  JOIN feature_type USING (feature_type_id)
 WHERE organism_id = ?
   AND version = ?
   AND feature_type.name = ?
});

   __PACKAGE__->set_sql(get_features_for_organism_version_chromosome => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
  JOIN location USING (feature_id)
 WHERE organism_id = ?
   AND version = ?
   AND location.chromosome = ?
});

   __PACKAGE__->set_sql(get_features_for_organism_version_type_chromosome => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
  JOIN feature_type USING (feature_type_id)
  JOIN location USING (feature_id)
 WHERE organism_id = ?
   AND version = ?
   AND feature_type.name = ?
   AND location.chromosome = ?
});

   __PACKAGE__->set_sql(get_features_for_organism_type_chromosome => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
  JOIN feature_type USING (feature_type_id)
  JOIN location USING (feature_id)
 WHERE organism_id = ?
   AND feature_type.name = ?
   AND location.chromosome = ?
});

   __PACKAGE__->set_sql(get_features_for_organism_chromosome => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
  JOIN location USING (feature_id)
 WHERE organism_id = ?
   AND location.chromosome = ?
});

   __PACKAGE__->set_sql(get_features_for_organism_type => qq{
SELECT DISTINCT(feature_id) 
  FROM feature
  JOIN dataset USING (dataset_id)
  JOIN feature_type USING (feature_type_id)
 WHERE organism_id = ?
   AND feature_type.name = ?
});

  __PACKAGE__->set_sql(has_genomic_sequence => qq{
SELECT count(*)
 FROM genomic_sequence
 WHERE dataset_id = ?
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

  #let's find all the dataset objects for arabidopsis version 2

  foreach my $ds ($db->get_dataset_obj->search(organism_id=>$org->id,version=>2)
    {
      #let's go thru all the features for each dataset object
      foreach my $feat ($ds->feats)
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

Dataset objects contain information about the origins of genomic data
such as the file name of the data, a URL link to where it was obtained, the 
version of the data, etc.  The reason that this is separate from the notion of 
a "data source" is that a given data source (e.g. TIGR) can release multiple 
files per genome (e.g. one file per chromosome) , different genomes, and 
different versions of any given genome.  In order to accommodate such variation
in genomic data sources, this table provides a way to track individual data
sets from a given data source.  After writing the database, this module, 
all related modules, many end-user applications, I realize that a better
name for this object would have been data_set.  Although it would be great
to update the database, API, and end-user applications to reflect this,
chances are it will not be done anytime soon.

In terms of functionality, this object is one of the most important for
figuring out which features or other genomic information is needed.  For 
example, let's say that the database has two genomes, Arabidopsis and Mus musculus
and has two versions of each genome (version 1 and version 2).  Each genome also 
has an individual data file for each chromosome (that means 5 for each version of 
Arabidopsis and 21 for each version of mouse).  That means that for this example, 
there are 52 data information entries in our database!  Understanding if a user is 
interested in chromosome 1 of Arabidopsis version 2 can make our database searches
faster as well as giving us ( the programmers ) the power of tracking several
versions of any given genome.



This object inherits from CoGe::Genome::DB which in turn inherits from Class::DBI.
Class::DBI provides the basic methods for creating accessor methods for accessing
table information.  Please see manual pages for Class::DBI for additional information.

The columns for this table are:
 dataset_id
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

dataset_id =>  database entry id
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
    return $self->description();
  }

sub id
  {
    my $self = shift;
    return $self->dataset_id();
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

=head2 resolve_dataset

 Usage     : my $ds = resolve_dataset($ds_thing);
 Purpose   : given a dataset name, a dataset database id, or a dataset object,
             this will return the dataset object for you
 Returns   : CoGe::Genome::DB::Dataset object
 Argument  : ds_thing can be a dataset name, a dataset database id, 
             or a dataset object
 Throws    : will throw a warning if a valid object was not created
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub resolve_dataset
  {
    my $self = shift;
    my $dsin = shift;
    my $dsout;
    if (ref ($dsin) =~ /Dataset/i) #we were passed an object
      {
	$dsout = $dsin;
      }
    elsif ($dsin =~ /^\d+$/) #only numbers, probably a database id
      {
	$dsout = $self->retrieve($dsin);
      }
    else #probably a name. . .
      {
	($dsout) = $self->search_like (name=>"%".$dsin."%")
      }
    warn "unable to resolve data $dsin in Dataset->resolve_dataset" unless ref ($dsout) =~ /Dataset/i;
    return $dsout;
  }


################################################ subroutine header begin ##

=head2 get_associated_datasets

 Usage     : my @ds = $dataset_obj->get_associated_datasets()
 Purpose   : This will find other dataset objects that might be important given the current 
             dataset object.  For example, let's say you loaded a genome where each chromosome
             was loaded separately so each has its own dataset object.  Then you load another
             set of data for this gene such as annotations from an HMM profile search.  However, this
             data-set covers the entire genome and not any individual chromosome.  So, if you want to
             find if a gene on some chromosome also has an HMM domain annotated, you will need to 
             1. find the gene's dataset object
             2. find any other dataset objects that might contain relevant information (using this function)
             3. assemble that data together
 Returns   : an array or array ref of CoGe::Genome::DB::Dataset objects depending on value of wantarray.
             This array will include the Dataset objects used in the search
 Argument  : accepts one or more dataset objects or uses self if none are specified.
 Throws    : None
 Comments  : The mechanism by which this routine works is to:
             1. find the organism, version, and chromosome(s) for the give dataset objects(s)
             2. search for other dataset objects with the same organism, version, and chromosome
           : NOTE: Chromosome search has been disabled due to speed constraints.
See Also   : 

=cut

################################################## subroutine header end ##


sub get_associated_datasets
  {
    my $self = shift;
    my $db = new CoGe::Genome;
    my @dso = @_;
    push @dso, $self unless @dso;
    my %ids;
    foreach my $ds (@dso)
      {
	$ids{$ds->id} = $ds;
	foreach my $tds ($self->search({organism_id=>$ds->org->id, version=>$ds->version}))
	  {
	    $ids{$tds->id}=$tds;
	  }
      }
    my @new_dsos = values %ids;
    return wantarray ? @new_dsos : \@new_dsos;
  }

sub get_associated_data_infos
  {
    my $self = shift;
    print STDERR "This method is obsolete, use get_associated_datasets. Also, bug Josh to replace and delete me\n";
    return ($self->get_associated_datasets);
  }

sub get_associated_data_information
  {
    my $self = shift;
    print STDERR "This method is obsolete, use get_associated_datasets. Also, bug Josh to replace and delete me\n";
    return ($self->get_associated_datasets);
  }

sub get_associated_data_informations
  {
    my $self = shift;
    print STDERR "This method is obsolete, use get_associated_datasets. Also, bug Josh to replace and delete me\n";
    return ($self->get_associated_datasets);
  }


################################################ subroutine header begin ##

=head2 get_chromosomes

 Usage     : my @chrs = $dataset_obj->get_chromosomes($id)
 Purpose   : finds all the chromosome names for one or more dataset objects
 Returns   : an array or array ref of strings (depending on wantarray)
 Argument  : accepts an array of dataset objects or database ids.  
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
    my @ds = @_;
    push @ds, $self unless @ds;
    my %chrs;
    my $sth = $self->sql_get_chromosomes();
    foreach my $ds (@ds)
      {
	my $id = ref ($ds) =~ /dataset/i ? $ds->id : $ds;
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

################################################ subroutine header begin ##

=head2 chromosomes

 Usage     : $dataset_obj->chromosomes
 Purpose   : alias for get_chromosomes

See Also   : $self->get_chromosomes

=cut

################################################## subroutine header end ##



sub chromosomes
  {
    my $self = shift;
    return $self->get_chromosomes(@_);
  }

################################################ subroutine header begin ##

=head2 chr

 Usage     : $dataset_obj->chr
 Purpose   : alias for get_chromosomes

See Also   : $self->get_chromosomes

=cut

################################################## subroutine header end ##


sub chr
  {
    my $self = shift;
    return $self->get_chromosomes(@_);
  }

################################################ subroutine header begin ##

=head2 get_feature_type_count

 Usage     : my $hash_ref = $ds->get_feature_type_count()
 Purpose   : get a count of all the feature types contained in a dataset object
             for example, you want to know the number of genes, mRNAs, CDS, etc contained
             in a particular dataset set
 Returns   : a hash ref where keys are the names of a feature type and the values are the 
             counts of that feature type in the dataset
 Argument  : will take a dataset database id but
             if none are specified, it will use itself
 Throws    : none
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_feature_type_count
  {
    my $self = shift;
    my %opts = @_;
    my $ds = $opts{ds};
    my $chr = $opts{chr};
    $ds = $self->id unless $ds;
    return unless $ds;
    my %feats;
    my $sth = $chr ? $self->sql_get_feature_type_count_chr() : $self->sql_get_feature_type_count();
    $chr ? $sth->execute($ds, $chr) : $sth->execute($ds);
    while (my $q = $sth->fetchrow_arrayref)
      {
	$feats{$q->[1]} = $q->[0];
      }
    return \%feats;
  }

################################################ subroutine header begin ##

=head2 get_current_version_for_organism

 Usage     : my $version = $ds->get_current_version_for_organism(org=>$org_obj);
 Purpose   : Find the most current version number for an organism
 Returns   : an integer that is the most current version of data for an organism
 Argument  : a CoGe::Genome::DB::Organism object of the database id of an organism
 Throws    : return undef if
 Comments  : One problem of this method is when an organism has multiple chromosomes
           : and each chromosome is of a different version.  In this case, this
           : method will return 0

See Also   : 

=cut

################################################## subroutine header end ##

sub get_current_version_for_organism
  {
    my $self = shift;
    my %opts = @_;
    my  $db = new CoGe::Genome;
    my $org = $db->get_organism_obj->resolve_organism($opts{org});
    return 0 unless $org;
    my $sth = $self->sql_get_current_version_for_organsim();
    $sth->execute($org->id);
    my %data;
    my $version;
    while (my $row = $sth->fetchrow_arrayref)
      {
	#row[0] version, row[1] chromosome
	$data{$row->[1]} = $row->[0] unless $data{$row->[1]};
	if ($data{$row->[1]} < $row->[0])
	  {
	    $data{$row->[1]} = $row->[0]; #chromosome, version	    
	  }
      }
    foreach my $v (values %data)
      {
	$version = $v unless $version;
	if ($version != $v)
	  {
	    return 0;
	  }
      }
    return $version;
  }

################################################ subroutine header begin ##

=head2 get_features_for_organism

 Usage     : my @features = $ds->get_features_for_organism(org=>$org_obj);
 Purpose   : Returns an array or array ref of CoGe::Genome::DB::Feature objects
             for an organism
 Returns   : an array or array ref of CoGe::Genome::DB::Feature objects
 Argument  : a hash of key->value pairs
             org           => CoGe::Genome::DB::Organism object or organism database 
                              id or organism name
             version | ver => which version of data to use (if not specified this
                            routine will find the current most version of the data
                            and use that)
             type        => the type of feature to return (e.g. "gene").  If this is
                            not specified it will return all features. (optional)
             chr | chromosome => chromosome name to limit search to just one 
                                 chromosome (optional)
 Throws    : none
 Comments  : 
           : 

See Also   : 

=cut

################################################## subroutine header end ##

sub get_features_for_organism
  {
    my $self = shift;
    my %opts = @_;
    my $db = new CoGe::Genome;
    my $org = $db->get_organism_obj->resolve_organism($opts{org});
    return 0 unless $org;
    my $version = $opts{version} || $opts{ver} || $self->get_current_version_for_organism(org=>$org);
    my $type = $opts{type};
    my $chr = $opts{chromosome} || $opts{chr};
    my $sth;
    my @args;
    if ($type && $chr && $version)
      {
	$sth = $self->sql_get_features_for_organism_version_type_chromosome();

	push @args, ($version, $type, $chr);
      }
    elsif ($type && $version)
      {
	$sth = $self->sql_get_features_for_organism_version_type();
	push @args, ($version, $type);
      }
    elsif ($chr && $version)
      {
	$sth = $self->sql_get_features_for_organism_version_chromosome();
	push @args, ($version, $chr);
      }
    elsif ($chr && $type)
      {
	$sth = $self->sql_get_features_for_organism_type_chromosome();
	push @args, ($type, $chr);
      }
    elsif ($type)
      {
	$sth = $self->sql_get_features_for_organism_type();
	push @args, ($type);
      }
    elsif ($chr)
      {
	$sth = $self->sql_get_features_for_organism_chromosome();
	push @args, ($chr);
      }
    elsif ($version)
      {
	$sth = $self->sql_get_features_for_organism_version();
	push @args, $version;
      }
    else
      {
	$sth = $self->sql_get_features_for_organism();
      }
    $sth->execute($org->id, @args);
    my @feats;
    while (my $q = $sth->fetchrow_arrayref)
      {
	push @feats, $db->get_feature_obj->retrieve($q->[0]);
      }
    return wantarray ? @feats : \@feats;
  }

################################################ subroutine header begin ##

=head2 has_genomic_sequence

 Usage     : if ($dataset_obj->has_genomic_sequence) { do stuff . . .}
 Purpose   : find out if a dataset has associated genomic sequence
 Returns   : returns the number of entries in the genomic_sequence table for this dataset_id
 Argument  : none
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub has_genomic_sequence
  {
    my $self = shift;
    my $sth = $self->sql_has_genomic_sequence();
    $sth->execute($self->id);
    my $num = $sth->fetchrow_arrayref->[0];
    return $num;
  }

################################################ subroutine header begin ##

=head2 get_dataset

 Usage     : my $ds = $dataset_obj->get_dataset(org=>"Arabidopsis", version=>5, chr=>1);
 Purpose   : gets a dataset object based on user criteria
 Returns   : CoGe::Genome::DB::Dataset object
 Argument  : org     => organism name or database id
             version => version of the organism's data
             chr     => chromsome
             dataset => dataset name or database id
 Throws    : carps if opts are not correctly specified or a dataset can't
             be found
 Comments  : This method is used to find a dataset object matching
             a variety of criteria.  Usually this will involve either
             an organism name, version, and chromosome OR the name (or
             database id) of a dataset.  If an organism is provided but
             the version is not specified, the most recent version of the 
             organism is used.  If no chromosome is specified, then chr: 1
             is used.

See Also   : CoGe::Genome::DB::Dataset

=cut

################################################## subroutine header end ##



sub get_dataset
  {
    my $self = shift;
    my %opts = @_;
    my $org = $opts{org};
    my $version = $opts{version};
    my $chr = $opts{chr} || 1;
    my $dataset = $opts{dataset};
    if ($dataset)
      {
	return $self->resolve_dataset($dataset);
      }
    my  $db = new CoGe::Genome;
    $org = $db->get_organism_obj->resolve_organism($org);
    unless ($org) 
      {
	carp "No organism or dataset specified.  get_dataset can't proceed.\n";
	return;
      }
    $version = $self->get_current_version_for_organism(org=>$org) unless $version;
    foreach my $ds ($org->datasets)
      {
	next unless $ds->version eq $version; #does version match?
	foreach my $c ($ds->get_chromosomes)
	  {
	    return $ds if $c eq $chr; #does chr match?
	  }
      }
    carp "Unable to find a dataset for org: $org, version: $version, chromsome: $chr\n";
    return;

  }

################################################ subroutine header begin ##

=head2 

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : CoGe::Genome::DB::Feature

=cut

################################################## subroutine header end ##


1; #this line is important and will help the module return a true value

