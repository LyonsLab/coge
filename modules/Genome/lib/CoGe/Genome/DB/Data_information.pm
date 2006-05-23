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

Genome::DB::Data_source - Genome::DB::Data_source

=head1 SYNOPSIS

  use Genome::DB::Data_source
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

