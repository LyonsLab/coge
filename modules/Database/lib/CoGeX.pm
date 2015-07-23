package CoGeX;

no warnings qw(redefine);
use strict;
use warnings;
use Data::Dumper;
use Cwd 'abs_path';

use base 'DBIx::Class::Schema';
use base qw(Class::Accessor);

#__PACKAGE__->load_classes();
__PACKAGE__->load_namespaces();
#__PACKAGE__->mk_accessors();

=head1 NAME

CoGeX - CoGeX

=head1 SYNOPSIS

  use CoGeX;
  This object is the API to the CoGe genomes database.  It uses DBIx::Class for managing
  relationships and access to the database.  Various other "high-level" functions are
  provided to make getting genomic data easier.

=head1 DESCRIPTION

Primary object for interacting with CoGe database system.

=head1 AUTHORS

 Eric Lyons
 Brent Pedersen

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

BEGIN {
    use vars qw( %pool $VERSION ); # persistent db connection pool
    $VERSION = 5.6;
}

################################################ subroutine header begin ##

=head2 dbconnect

 Usage     : use CoGeX;
             my $coge = CoGeX->dbconnect;
 Purpose   : generates the CoGeX API object
 Returns   : CoGeX object
 Argument  : none
 Throws    : none
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub dbconnect {
    my $self      = shift;
    my $conf      = shift;
    my $conn_name = shift;    # optional connection name
    $conn_name = 'default' unless $conn_name;

    # Connect to the database
    unless ( defined $pool{$conn_name} )
        #and $pool{$conn_name}->storage->dbh->ping() )
    {
        my $db      = $conf->{DB};
        my $dbname  = $conf->{DBNAME};
        my $dbhost  = $conf->{DBHOST};
        my $dbport  = $conf->{DBPORT};
        my $dbuser  = $conf->{DBUSER};
        my $dbpass  = $conf->{DBPASS};
        my $connstr = "dbi:$db:dbname=$dbname;host=$dbhost;port=$dbport";
        $pool{$conn_name} = $self->connect( $connstr, $dbuser, $dbpass );

        #print STDERR "CoGeX: new connection '$conn_name'\n";
        #$coge->storage->debugobj(new DBIxProfiler());
        #$coge->storage->debug(1);
    }

    return $pool{$conn_name};
}

################################################ subroutine header begin ##

=head2 get_features_in_region

 Usage     : $object->get_features_in_region(start   => $start,
                                             stop    => $stop,
                                             chr     => $chr,
                                             ftid    => $ftid,
                                             dataset_id => $dataset->id(),);

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
             ftid    => limit features to those with this feature type id
 Throws    : none
 Comments  :

See Also   :

=cut

################################################## subroutine header end ##

sub get_features_in_region {
    my $self = shift;
    my %opts = @_;
    my $start =
         $opts{'start'}
      || $opts{'START'}
      || $opts{begin}
      || $opts{BEGIN};
    $start = 1 unless $start;
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    $stop = $start unless defined $stop;
    my $chr = $opts{chr};
    $chr = $opts{chromosome} unless defined $chr;
    my $dataset_id =
         $opts{dataset}
      || $opts{dataset_id}
      || $opts{info_id}
      || $opts{INFO_ID}
      || $opts{data_info_id}
      || $opts{DATA_INFO_ID};
    my $genome_id  = $opts{gid};
    my $count_flag = $opts{count} || $opts{COUNT};
    my $ftid       = $opts{ftid};

    if ( ref($ftid) =~ /array/i ) {
        $ftid = undef unless @$ftid;
    }
    my @dsids;
    push @dsids, $dataset_id if $dataset_id;
    if ($genome_id) {
        my $genome = $self->resultset('Genome')->find($genome_id);
        push @dsids, map { $_->id } $genome->datasets if $genome;
    }
    if ($count_flag) {
        return $self->resultset('Feature')->count(
            {
                "me.chromosome" => $chr,
                "me.dataset_id" => [@dsids],
                -and            => [
                    "me.start" => { "<=" => $stop },
                    "me.stop"  => { ">=" => $start },
                ],

                #  -and=>[
                # 	  -or=>[
                # 		-and=>[
                # 		       "me.stop"=>  {"<=" => $stop},
                # 		       "me.stop"=> {">=" => $start},
                # 		      ],
                # 		-and=>[
                # 		       "me.start"=>  {"<=" => $stop},
                # 		       "me.start"=> {">=" => $start},
                # 		      ],
                # 		-and=>[
                # 		       "me.start"=>  {"<=" => $start},
                # 		       "me.stop"=> {">=" => $stop},
                # 		      ],
                # 	       ],
                # 	 ],
            },
            {

                #						   prefetch=>["locations", "feature_type"],
            }
        );
    }
    my %search = (
        "me.chromosome" => $chr,
        "me.dataset_id" => [@dsids],
        -and            => [
            "me.start" => { "<=" => $stop },
            "me.stop"  => { ">=" => $start },
        ],

        # -and=>[
        # 	-or=>[
        # 	      -and=>[
        # 		     "me.stop"=>  {"<=" => $stop},
        # 		     "me.stop"=> {">=" => $start},
        # 		    ],
        # 	      -and=>[
        # 		     "me.start"=>  {"<=" => $stop},
        # 		     "me.start"=> {">=" => $start},
        # 		    ],
        # 	      -and=>[
        # 		     "me.start"=>  {"<=" => $start},
        # 		     "me.stop"=> {">=" => $stop},
        # 		    ],
        # 	     ],
        #      ]
    );
    $search{"me.feature_type_id"} = { "IN" => $ftid } if $ftid;
    my @feats = $self->resultset('Feature')->search(
        \%search,
        {

            #					     prefetch=>["locations", "feature_type"],
            #						     order_by=>"me.start",
        }
    );
    return wantarray ? @feats : \@feats;
}

sub get_features_in_region_split {
    my $self = shift;
    my %opts = @_;
    my $start =
         $opts{'start'}
      || $opts{'START'}
      || $opts{begin}
      || $opts{BEGIN};
    $start = 0 unless $start;
    my $stop = $opts{'stop'} || $opts{STOP} || $opts{end} || $opts{END};
    $stop = $start unless defined $stop;
    my $chr = $opts{chr};
    $chr = $opts{chromosome} unless defined $chr;
    my $dataset_id =
         $opts{dataset}
      || $opts{dataset_id}
      || $opts{info_id}
      || $opts{INFO_ID}
      || $opts{data_info_id}
      || $opts{DATA_INFO_ID};

    my @startfeats = $self->resultset('Feature')->search(
        {
            "me.chromosome" => $chr,
            "me.dataset_id" => $dataset_id,
            -and            => [
                "me.stop" => { ">=" => $start },
                "me.stop" => { "<=" => $stop },
            ],
        },
        { prefetch => [ "locations", "feature_type" ], }
    );
    my @stopfeats = $self->resultset('Feature')->search(
        {
            "me.chromosome" => $chr,
            "me.dataset_id" => $dataset_id,
            -and            => [
                "me.start" => { ">=" => $start },
                "me.start" => { "<=" => $stop },
            ],
        },
        { prefetch => [ "locations", "feature_type" ], }
    );

    my %seen;
    my @feats;

    foreach my $f (@startfeats) {
        if ( not exists $seen{ $f->id() } ) {
            $seen{ $f->id() } += 1;
            push( @feats, $f );
        }
    }

    foreach my $f (@stopfeats) {
        if ( not exists $seen{ $f->id() } ) {
            $seen{ $f->id() } += 1;
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

sub count_features_in_region {
    my $self = shift;
    my %opts = @_;
    return $self->get_features_in_region( %opts, count => 1 );
}

sub get_current_datasets_for_org {
    my $self = shift;
    warn
'THIS METHOD (get_current_datasets_for_org) IS OBSELETE.  PLEASE CALL $org_obj->current_datasets()\n';
    my %opts = @_ if @_ > 1;
    my $orgid = $opts{org} || $opts{orgid} || $opts{organism};
    $orgid = shift unless $orgid;
    return unless $orgid;
    my ($org) = $self->resultset('Organism')->resolve($orgid);
    return unless $org;
    return $org->current_datasets();
}

sub log_user {
    my $self    = shift;
    my %opts    = @_;
    my $user    = $opts{user};
    my $session = $opts{session};
    my $uid     = ref($user) =~ /User/ ? $user->id : $user;
    unless ( $uid =~ /^\d+$/ ) {
        warn "Error adding user_id to User_session.  Not a valid uid: $uid\n";
        return;
    }

    #FIRST REMOVE ALL ENTRIES FOR THIS USER
    foreach my $item (
        $self->resultset('UserSession')->search( { session => $session } ) )
    {
        next unless $item;
        $item->delete;
    }

    #ADD NEW ENTRY
    my $item =
      $self->resultset('UserSession')
      ->create( { user_id => $uid, date => \'NOW()', session => $session } );
    return $item;
}

sub self_or_default {    #adapted from CGI.pm
    shift @_ if $_[0] eq "CoGeX";
    my $Q;
    unless (
        defined( $_[0] )
        && (
            ref( $_[0] ) eq 'CoGeX'
        ) # || UNIVERSAL::isa($_[0],'CoGeX')) # slightly optimized for common case
      )
    {
        $Q = CoGeX->new unless defined($Q);
        unshift( @_, $Q );
    }
    return wantarray ? @_ : $Q;
}

sub node_types {
    my $self  = shift;
    my %types = (
        list       => 1,
        notebook   => 1,   # alias for list
        genome     => 2,
        experiment => 3,
        feature    => 4,
        user       => 5,
        user_group => 6,
        group      => 6,   # alias for user_group
        workflow   => 7
    );
    return wantarray ? %types : \%types;
}

sub node_type_name {
    my $self = shift;
    my $type_id = shift;
    my $node_types = $self->node_types;
    my ($type_name) = grep { $node_types->{$_} eq $type_id } keys %$node_types;
    $type_name = 'notebook' if ( $type_name eq 'list' );
    return $type_name;
}

1;
