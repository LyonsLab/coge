package CoGeX;

use strict;
use warnings;

use vars qw( $VERSION );

$VERSION = 0.01;

use base 'DBIx::Class::Schema';

__PACKAGE__->load_classes;

=head1 NAME

CoGeX - CoGeX

=head1 SYNOPSIS

  use CoGeX;
  blah blah blah


=head1 DESCRIPTION

Primary object for interacting with CoGe database system.

=head1 USAGE

  use CoGeX;

  my $connstr = 'dbi:mysql:genomes:biocon:3306';
  my $s = CoGeX->connect($connstr, 'cnssys', 'CnS' ); # Biocon's ro user

  my $rs = $s->resultset('Feature')->search(
                {
                  'organism.name' => "Arabidopsis thaliana"
                },
                { join => [ 'dataset', 'organism' ] }
  );

=head1 BUGS


=head1 SUPPORT


=head1 AUTHOR

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut
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

sub get_features_in_region_old
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
    if ($count_flag)
      {
	return $self->resultset('Feature')->count({
						   "locations.chromosome"=>$chr,
						   "dataset_id"=>$dataset_id,
						   -and=>[
							  -or=>[
								-and=>[
								       "locations.stop"=>  {"<=" => $stop},
								       "locations.stop"=> {">=" => $start},
								      ],
								-and=>[
								       "locations.start"=>  {"<=" => $stop},
								       "locations.start"=> {">=" => $start},
								      ],
							       ],
							 ],
						  },
						  {
						   join => ["locations"],
						   distinct=>['feature_id'],
						  }
						 );
      }
    my @feats = $self->resultset('Feature')->search({
						     "locations.chromosome"=>$chr,
						     "me.dataset_id"=>$dataset_id,
						     -and=>[
                   -or=>[
                     -and=>[
                       "locations.stop"=>  {"<=" => $stop},
                       "locations.stop"=> {">=" => $start},
                     ],
                     -and=>[
                       "locations.start"=>  {"<=" => $stop},
                       "locations.start"=> {">=" => $start},
                     ],
                   ],
                 ],
						     },
                  {
                     join => ["locations"],
                     distinct=>['feature_id'],
                     prefetch=>["locations", "feature_type"],
                  }
						   );
    return wantarray ? @feats : \@feats;
  }



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
    if ($count_flag)
      {
	return $self->resultset('Feature')->count(
						  {
						   "me.chromosome" => $chr,
						   "me.dataset_id" => $dataset_id,
						   -and=>[
							  -or=>[
								-and=>[
								       "me.stop"=>  {"<=" => $stop},
								       "me.stop"=> {">=" => $start},
								      ],
								-and=>[
								       "me.start"=>  {"<=" => $stop},
								       "me.start"=> {">=" => $start},
								      ],
							       ],
							 ],
						  },
						  {
#						   prefetch=>["locations", "feature_type"],
						  }
						 );
      }
    my @feats = $self->resultset('Feature')->search({
                 "me.chromosome" => $chr,
                 "me.dataset_id" => $dataset_id,
						     -and=>[
                   -or=>[
                     -and=>[
                       "me.stop"=>  {"<=" => $stop},
                       "me.stop"=> {">=" => $start},
                     ],
                     -and=>[
                       "me.start"=>  {"<=" => $stop},
                       "me.start"=> {">=" => $start},
                     ],
                   ],
                 ],
						     },
                  {
                    prefetch=>["locations", "feature_type"],
                    #prefetch=>["locations", "dataset", "feature_type"],
                    #prefetch=>["locations"],
                    #prefetch=>["dataset"],
                    #prefetch=>["feature_type"],
                  }
						   );
    return wantarray ? @feats : \@feats;
  }

sub get_features_in_region_split
  {
    my $self = shift;
    my %opts = @_;
    my $start = $opts{'start'} || $opts{'START'} || $opts{begin} || $opts{BEGIN};
    $start = 0 unless $start;
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    $stop = $start unless defined $stop;
    my $chr = $opts{chr} || $opts{CHR} || $opts{chromosome} || $opts{CHROMOSOME};
    my $dataset_id = $opts{dataset} || $opts{dataset_id} || $opts{info_id} || $opts{INFO_ID} || $opts{data_info_id} || $opts{DATA_INFO_ID} ;

    my @startfeats = $self->resultset('Feature')->search({
                 "me.chromosome" => $chr,
                 "me.dataset_id" => $dataset_id,
                 -and => [
                   "me.stop"=> {">=" => $start},
                   "me.stop"=>  {"<=" => $stop},
                 ],
						     },
                 {
                   prefetch=>["locations", "feature_type"],
                 }
						   );
    my @stopfeats = $self->resultset('Feature')->search({
                 "me.chromosome" => $chr,
                 "me.dataset_id" => $dataset_id,
                 -and => [
                   "me.start"=>  {">=" => $start},
                   "me.start"=> {"<=" => $stop},
                 ],
						     },
                 {
                   prefetch=>["locations", "feature_type"],
                 }
						   );

    my %seen;
    my @feats;

    foreach my $f ( @startfeats ) {
      if ( not exists $seen{ $f->id() } ) {
        $seen{$f->id()}+=1;
        push( @feats, $f );
      }
    }

    foreach my $f ( @stopfeats ) {
      if ( not exists $seen{ $f->id() } ) {
        $seen{$f->id()}+=1;
        push( @feats, $f );
      }
    }

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


1;
