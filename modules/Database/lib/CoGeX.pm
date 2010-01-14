package CoGeX;

no warnings qw(redefine);
use strict;
use warnings;


use vars qw( $VERSION );

$VERSION = 0.01;

use base 'DBIx::Class::Schema';
use base qw(Class::Accessor);

__PACKAGE__->load_classes;
__PACKAGE__->mk_accessors(qw(db_connection_string db_name db_passwd));


use vars qw($DEFAULT_CONNECTION_STRING $DEFAULT_NAME $DEFAULT_PASSWD);
$DEFAULT_CONNECTION_STRING = 'dbi:mysql:dbname=coge;host=biocon.berkeley.edu;port=3306';

#$DEFAULT_CONNECTION_STRING = 'dbi:mysql:dbname=coge;host=homer.cnr.berkeley.edu;port=3306';
$DEFAULT_NAME = "coge";
$DEFAULT_PASSWD = "123coge321";

=head1 NAME

CoGeX - CoGeX

=head1 SYNOPSIS

  use CoGeX;
  This object is the API to the CoGe genomes database.  It uses DBIx::Class for managing
  relationships and access to the database.  Various other "high-level" functions are
  provided to make getting genomic data easier.


=head1 DESCRIPTION

Primary object for interacting with CoGe database system.

=head1 USAGE

  use CoGeX;

  my $s = CoGeX->connectdb(); # Biocon's ro user

  my $rs = $s->resultset('Feature')->search(
                {
                  'organism.name' => "Arabidopsis thaliana"
                },
                { join => [ 'dataset', 'organism' ] }
  );

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

perl(1).

=cut

################################################ subroutine header begin ##

=head2 dbconnect

 Usage     : use CoGeX;
             my $coge = CoGeX->dbconnect;

 Purpose   : generates the CoGe genomes API object
 Returns   : CoGeX object
 Argument  : db_connection_string=>'dbi:mysql:dbname=genomes;host=dbhostname.edu;port=3306'
             db_name=>'database_user_name'
             db_passwd=>'user_password'
             Arguments should be optional with defaults set in object so that this function can be
             called to generate db API object
 Throws    : none
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##


sub dbconnect
  {
    my ($self, %opts) = self_or_default(@_);
    my $str = $opts{db_connection_string};
    my $name = $opts{db_name};
    my $pwd = $opts{db_passwd};
    $str = $self->db_connection_string unless $str;
    $str = $DEFAULT_CONNECTION_STRING unless $str;
    $name = $self->db_name unless $name;
    $name = $DEFAULT_NAME unless $name;
    $pwd = $self->db_passwd unless $pwd;
    $pwd = $DEFAULT_PASSWD unless $pwd;
    my $cogedb = CoGeX->connect($str, $name, $pwd );
    $cogedb->db_connection_string($str);
    $cogedb->db_name($name);
    $cogedb->db_passwd($pwd);
    return $cogedb;
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
    my $start = 1 unless $start;
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
								-and=>[
								       "me.start"=>  {"<=" => $start},
								       "me.stop"=> {">=" => $stop},
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
			      -and=>[
				     "me.start"=>  {"<=" => $start},
				     "me.stop"=> {">=" => $stop},
				    ],
			     ],
		       ],
						    },
						    {
						     prefetch=>["locations", "feature_type"],
						     order_by=>"me.start",
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

sub get_current_datasets_for_org
  {
    my $self = shift;
    warn 'THIS METHOD (get_current_datasets_for_org) IS OBSELETE.  PLEASE CALL $org_obj->current_datasets()\n';
    my %opts = @_ if @_ >1;
    my $orgid = $opts{org} || $opts{orgid} || $opts{organism};
    $orgid = shift unless $orgid;
    return unless $orgid;
    my ($org) = $self->resultset('Organism')->resolve($orgid);
    return unless $org;
    return $org->current_datasets();
  }

sub log_user
  {
    my $self = shift;
    my %opts = @_;
    my $user = $opts{user};
    my $session = $opts{session};
    my $uid = ref($user) =~ /User/ ? $user->id : $user;
    unless ($uid =~ /^\d+$/)
      {
	warn "Error adding user_id to User_session.  Not a valid uid: $uid\n";
	return;
      }
    #FIRST REMOVE ALL ENTRIES FOR THIS USER
    foreach my $item ($self->resultset('UserSession')->search(session=>$session))
      {
	next unless $item;
	$item->delete;
      }
    #ADD NEW ENTRY
    my $item = $self->resultset('UserSession')->create({user_id=>$uid, date=>\'NOW()', session=>$session});
    return $item;
  }

sub self_or_default 
    { #adapted from CGI.pm
      shift @_ if $_[0] eq "CoGeX";
      my $Q;
      unless (
	      defined($_[0]) && 
	      (ref($_[0]) eq 'CoGeX')# || UNIVERSAL::isa($_[0],'CoGeX')) # slightly optimized for common case
	     ) 
	{
	  $Q = CoGeX->new unless defined($Q);
	  unshift(@_,$Q);
	}
      return wantarray ? @_ : $Q;
    }


1;
