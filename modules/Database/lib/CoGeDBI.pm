package CoGeDBI;

=head1 NAME

CoGeDBI

=head1 SYNOPSIS

  This module provides a low-level API to the CoGe database.  It was written for
  high-performance cases where the CoGeX-DBIX ORM is too slow.

=head1 DESCRIPTION

Low-level access to CoGe database.

=head1 AUTHORS

 Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;
use DBI;
use Hash::Merge::Simple qw(merge);
use Data::Dumper;
#use Time::HiRes qw(time);

BEGIN {
    use vars qw ($VERSION $MAX_FETCH_ROWS @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    $MAX_FETCH_ROWS = 10_000;
    @ISA     = qw (Exporter);
    @EXPORT = qw( 
        get_table get_user_access_table get_experiments get_distinct_feat_types
        get_genomes_for_user get_experiments_for_user get_lists_for_user
        get_groups_for_user get_group_access_table get_datasets
        get_feature_counts get_feature_count_summary get_features get_feature_types
        get_feature_names get_feature_annotations get_locations 
        get_total_queries get_dataset_ids feature_type_names_to_id
        get_features_by_range get_table_count get_uptime
        get_favorites_for_user
    );
}

# takes: string of comma delimeted, singled quoted feature type names
# returns: string of comma delimited corresponding feature type ids
sub feature_type_names_to_id {
	my $type_names = shift;
	my $dbh = shift;
	return join(',', @{$dbh->selectcol_arrayref('SELECT feature_type_id FROM feature_type WHERE name IN(' . $type_names . ')')});
}

sub get_table {
    my $dbh = shift;         # database connection handle
    my $table = shift;       # table name
    my $hash_fields = shift; # array ref of field names to use as hash keys
    my $conditions = shift;
    my $operator_hash = shift; #hash ref of the comparison operators to be used, with the conditions as keys
    $hash_fields = [$table.'_id'] unless $hash_fields;
     my $operator_string = "=" unless $operator_hash;
    
    # Build query string
    my $query = "SELECT * FROM $table";
    if ($operator_hash) {
    	if ($conditions) {
    		$query .= ' WHERE ' . join(' AND ', map { $_.$operator_hash->{$_}.$conditions->{$_} } keys %$conditions);
    	}
    } else {
	    if ($conditions) {
	        $query .= ' WHERE ' . join(' AND ', map { $_.$operator_string.$conditions->{$_} } keys %$conditions);
	    }
	    #print STDERR $query, "\n";
    }
    
    # Execute query
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    # Fetch results
    my $results = $sth->fetchall_hashref($hash_fields);
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_user_access_table {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    #my $child_type = shift;  # optional child_type
    return unless $user_id;
    
    # Get user connections
    my $query = "SELECT * FROM user_connector WHERE parent_type=5 AND parent_id=$user_id";
    #$query .= " AND child_type=$child_type" if ($child_type);
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results1 = $sth->fetchall_hashref(['child_type', 'child_id']);
    #print STDERR Dumper $results1, "\n";
    
    # Get group connections for user
    my $groups_str = join(',', map { $_->{child_id} } values %{$results1->{6}});
    my $results2 = {};
    if ($groups_str) {
        $query = "SELECT * FROM user_connector WHERE parent_type=6 AND parent_id IN ($groups_str)";
        #$query .= " AND child_type=$child_type" if ($child_type);
        $sth = $dbh->prepare($query);
        $sth->execute();
        $results2 = $sth->fetchall_hashref(['child_type', 'child_id']);
        #print STDERR Dumper $results2, "\n";
    }
    
    # Get list connections for user and user's groups
    my $lists_str = join(',', map { $_->{child_id} } 
        (defined $results1 && defined $results1->{1} ? values %{$results1->{1}} : ()), 
        (defined $results2 && defined $results2->{1} ? values %{$results2->{1}} : ())
    );
    my $results3 = {};
    if ($lists_str) {
        my $query = "SELECT * FROM list_connector WHERE parent_id IN ($lists_str)";
        #$query .= " AND child_type=$child_type" if ($child_type);
        my $sth = $dbh->prepare($query);
        $sth->execute();
        $results3 = $sth->fetchall_hashref(['child_type', 'child_id']);
        #print STDERR Dumper $results3, "\n"; 
    }
    
    Hash::Merge::set_behavior('LEFT_PRECEDENT');
    #my $combined = Hash::Merge::merge($results1, $results2, $results3); # mdb removed 12/11/15
    my $combined = Hash::Merge::merge($results1, $results2); # order is important here, user connections should override group and list connections
    my $combined2 = Hash::Merge::merge($combined, $results3); # mdb added 12/11/15
    #print STDERR Dumper $combined2, "\n";
    
    return $combined2;
}

sub get_group_access_table {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    #my $child_type = shift;  # optional child_type
    return unless $user_id;
    
    # Get user connections
    my $query = "SELECT * FROM user_connector WHERE parent_type=6 AND parent_id=$user_id";
    #$query .= " AND child_type=$child_type" if ($child_type);
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results1 = $sth->fetchall_hashref(['child_type', 'child_id']);
    #print STDERR Dumper $results1, "\n";
    
    # Get list connections for user and user's groups
    my $lists_str = join(',', map { $_->{child_id} } 
        (defined $results1 && defined $results1->{1} ? values %{$results1->{1}} : ())
    );
    my $results2 = {};
    if ($lists_str) {
        my $query = "SELECT * FROM list_connector WHERE parent_id IN ($lists_str)";
        #$query .= " AND child_type=$child_type" if ($child_type);
        my $sth = $dbh->prepare($query);
        $sth->execute();
        $results2 = $sth->fetchall_hashref(['child_type', 'child_id']);
        #print STDERR Dumper $results3, "\n"; 
    }
    
    Hash::Merge::set_behavior('LEFT_PRECEDENT');
    my $combined = Hash::Merge::merge($results1, $results2);
    #print STDERR Dumper \%combined, "\n";
    
    return $combined;
}

sub get_group_str_for_user {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    return unless $user_id;
    
    # Get groups for user
    my $query = qq{
        SELECT child_id
        FROM user_connector AS uc 
        WHERE uc.parent_type=5 AND uc.parent_id=$user_id AND uc.child_type=6
    };
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $groups = $sth->fetchall_arrayref();
    my $group_str = join(', ', map { $_->[0] } @$groups);
    #print STDERR $group_str, "\n";
    return $group_str;
}

sub get_groups_for_user {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    return unless $user_id;
    
    # Get groups for user
    my $query = qq{
        SELECT ug.user_group_id AS id, ug.name AS name, ug.description AS description,
        ug.deleted AS deleted, ug.locked AS locked, ug.role_id AS role_id
        FROM user_connector AS uc 
        JOIN user_group AS ug ON (uc.child_id=ug.user_group_id)
        WHERE uc.parent_type=5 AND uc.parent_id=$user_id AND uc.child_type=6
    };
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $groups = $sth->fetchall_arrayref({});
    return $groups;
}

sub get_genomes_for_user {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    return unless $user_id;
    my $role_ids = shift;    # optional ref to role IDs array
    my $role_id_str;
    $role_id_str = join(', ', @$role_ids) if $role_ids;
    #print STDERR "CoGeDBI::get_genomes_for_user $user_id\n";
    
    # Get groups for user
    my $group_str = get_group_str_for_user($dbh, $user_id) || -1; # default to -1 to prevent empty IN clause in query
        
    # Get user/group genome connections
    # mdb added favorite_connector subquery for COGE-388
    my $query = qq{
        SELECT g.genome_id AS id, g.name AS name, g.description AS description, 
            g.version AS version, g.restricted AS restricted, g.certified AS certified, g.deleted AS deleted, 
            g.date AS date, o.name AS organism, ds.date AS dataset_date, uc.role_id AS role_id, 
            (CASE WHEN EXISTS (SELECT 1 FROM favorite_connector AS fc WHERE fc.user_id=$user_id AND fc.child_id=g.genome_id AND fc.child_type=2)
                THEN 1 ELSE 0
            END) AS favorite
        FROM user_connector AS uc 
        JOIN genome AS g ON (uc.child_id=g.genome_id) 
        JOIN dataset_connector AS dsc ON (dsc.genome_id=g.genome_id) 
        JOIN dataset AS ds ON (ds.dataset_id=dsc.dataset_id) 
        JOIN organism AS o ON (o.organism_id=g.organism_id) 
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=2 
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $query .= " ORDER BY uc.role_id DESC"; # mdb added 4/28/16 -- ensure that roles with greater priveleges have precedence if multiple connections to same item
    #my $t1 = time;
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results1 = $sth->fetchall_hashref(['id']);
    #my $t2 = time;
    #print STDERR 'query1: ', ($t2-$t1)*1000, "\n";
    #print STDERR Dumper $results1, "\n";

    # Get list genome connections
    # mdb added uc2 join for COGE-629
    # mdb added favorite_connector subquery for COGE-388
    $query = qq{
        SELECT g.genome_id AS id, g.name AS name, g.description AS description, 
            g.version AS version, g.restricted AS restricted, g.certified AS certified, g.deleted AS deleted, 
            g.date AS date, o.name AS organism, ds.date AS dataset_date, uc.role_id AS role_id,
            (CASE WHEN EXISTS (SELECT 1 FROM favorite_connector AS fc WHERE fc.user_id=$user_id AND fc.child_id=g.genome_id AND fc.child_type=2)
                THEN 1 ELSE 0
            END) AS favorite
        FROM user_connector AS uc 
        JOIN list AS l ON (uc.child_id=l.list_id) 
        JOIN list_connector AS lc ON (lc.parent_id=l.list_id)
        JOIN genome AS g ON (lc.child_id=g.genome_id)
        JOIN user_connector AS uc2 ON (uc2.child_id=g.genome_id)
        JOIN dataset_connector AS dsc ON (dsc.genome_id=g.genome_id) 
        JOIN dataset AS ds ON (ds.dataset_id=dsc.dataset_id) 
        JOIN organism AS o ON (o.organism_id=g.organism_id) 
        WHERE uc.child_type=1 AND ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc2.child_type=2 AND ((uc2.parent_type=5 AND uc2.parent_id=$user_id) OR (uc2.parent_type=6 AND uc2.parent_id IN ($group_str))) AND uc2.role_id IN (2,3,4)
            AND lc.child_type=2
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $query .= " ORDER BY uc.role_id DESC"; # mdb added 4/28/16 -- ensure that roles with greater priveleges have precedence if multiple connections to same item
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $results2 = $sth->fetchall_hashref(['id']);
    #my $t3 = time;
    #print STDERR 'query2: ', ($t3-$t2)*1000, "\n";
    #print STDERR Dumper $results2, "\n";
    
    #TODO Replace this merge with SQL UNION of queries above?
    Hash::Merge::set_behavior('LEFT_PRECEDENT');
    my $combined = Hash::Merge::merge($results1, $results2); # order is important here, results1 should overwrite results2
    #my $t4 = time;
    #print STDERR 'merge: ', ($t4-$t3)*1000, "\n";
    
    return [ values $combined ];
}

sub get_experiments_for_user {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    return unless $user_id;
    my $role_ids = shift;    # optional ref to role IDs array
    my $role_id_str;
    $role_id_str = join(', ', @$role_ids) if $role_ids;
    #print STDERR "CoGeDBI::get_experiments_for_user $user_id\n";
    
    # Get groups for user
    my $group_str = get_group_str_for_user($dbh, $user_id) || -1; # default to -1 to prevent empty IN clause in query
        
    # Get user/group experiment connections
    my $query = qq{
        SELECT e.experiment_id AS id, e.name AS name, e.description AS description, 
            e.version AS version, e.restricted AS restricted, e.deleted AS deleted, 
            e.date AS date, uc.role_id as role_id,
            (CASE WHEN EXISTS (SELECT 1 FROM favorite_connector AS fc WHERE fc.user_id=$user_id AND fc.child_id=e.experiment_id AND fc.child_type=3)
                THEN 1 ELSE 0
            END) AS favorite
        FROM user_connector AS uc 
        JOIN experiment AS e ON (uc.child_id=e.experiment_id) 
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=3
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $query .= " ORDER BY uc.role_id DESC"; # mdb added 3/31/16 -- ensure that roles with greater priveleges have precedence if multiple connections to same item
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results1 = $sth->fetchall_hashref(['id']);
    #print STDERR Dumper $results1, "\n";

    # Get list experiment connections
    # mdb added uc2 join for COGE-629
    $query = qq{
        SELECT e.experiment_id AS id, e.name AS name, e.description AS description, 
            e.version AS version, e.restricted AS restricted, e.deleted AS deleted, 
            e.date AS date, uc.role_id as role_id,
            (CASE WHEN EXISTS (SELECT 1 FROM favorite_connector AS fc WHERE fc.user_id=$user_id AND fc.child_id=e.experiment_id AND fc.child_type=3)
                THEN 1 ELSE 0
            END) AS favorite
        FROM user_connector AS uc 
        JOIN list AS l ON (uc.child_id=l.list_id) 
        JOIN list_connector AS lc ON (lc.parent_id=l.list_id)
        JOIN experiment AS e ON (lc.child_id=e.experiment_id)
        JOIN user_connector AS uc2 ON (uc2.child_id=e.experiment_id)
        WHERE uc.child_type=1 AND ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc2.child_type=3 AND ((uc2.parent_type=5 AND uc2.parent_id=$user_id) OR (uc2.parent_type=6 AND uc2.parent_id IN ($group_str))) AND uc2.role_id IN (2,3,4)
            AND lc.child_type=3 
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $query .= " ORDER BY uc.role_id DESC"; # mdb added 3/31/16 -- ensure that roles with greater priveleges have precedence if multiple connections to same item
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $results2 = $sth->fetchall_hashref(['id']);
    #print STDERR Dumper $results2, "\n";
    
    Hash::Merge::set_behavior('LEFT_PRECEDENT');
    my $combined = Hash::Merge::merge($results1, $results2); # order is important here, results1 should overwrite results2
    return [ values $combined ];
}

sub get_lists_for_user {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    return unless $user_id;
    my $role_ids = shift;    # optional ref to role IDs array
    my $role_id_str;
    $role_id_str = join(', ', @$role_ids) if $role_ids;
    #print STDERR "CoGeDBI::get_lists_for_user $user_id\n";
    
    # Get groups for user
    my $group_str = get_group_str_for_user($dbh, $user_id) || -1; # default to -1 to prevent empty IN clause in query
        
    # Get user/group list connections
    my $query = qq{
        SELECT l.list_id AS id, l.name AS name, l.description AS description, 
            l.restricted AS restricted, l.deleted AS deleted,
            l.date AS date, uc.role_id as role_id,
            (CASE WHEN EXISTS (SELECT 1 FROM favorite_connector AS fc WHERE fc.user_id=$user_id AND fc.child_id=l.list_id AND fc.child_type=1)
                THEN 1 ELSE 0
            END) AS favorite
        FROM user_connector AS uc 
        JOIN list AS l ON (uc.child_id=l.list_id) 
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str)))
            AND uc.child_type=1 
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $query .= " ORDER BY uc.role_id DESC"; # mdb added 4/28/16 -- ensure that roles with greater priveleges have precedence if multiple connections to same item
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_arrayref({});
    #print STDERR Dumper $results, "\n";

    return $results;
}

sub get_favorites_for_user {
    my $dbh = shift;         # database connection handle
    my $user_id = shift;     # user id
    return unless $user_id;
    my $type_ids = shift;     # optional child type id
    
    # Get
    my $query = qq{
        SELECT lc.child_id AS id
        FROM list AS l
        JOIN list_connector AS lc ON (lc.parent_id=l.list_id) 
        WHERE l.name='Favorites' 
    };
    $query .= ' AND lc.child_type IN (' . join(', ', @$type_ids) . ')' if $type_ids;
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_hashref('id');
    return $results;
}

sub get_experiments {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;   # optional genome id
    
    # Execute query
    my $query = "SELECT * FROM experiment";
    if ($genome_id) {
        $query .= " WHERE genome_id=$genome_id";
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    # Fetch results
    my $results = $sth->fetchall_hashref('experiment_id');
    #print STDERR Dumper $results, "\n";
    
    return wantarray ? values %$results : $results;
}

sub get_feature_counts {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;
    my $dataset_id = shift;
    return unless ($dbh and ($genome_id or $dataset_id));
    my $returnArray = shift; # optional flag to return result as raw array ref

    # Execute query
    my $query = qq{
        SELECT f.chromosome AS chromosome, f.feature_type_id AS type_id, count(*) AS count
        FROM dataset_connector AS dc
        JOIN feature AS f ON (dc.dataset_id=f.dataset_id)
    };
    if ($genome_id) {
        $query .= " WHERE (dc.genome_id=$genome_id)";
    }
    else {
        $query .= " WHERE (dc.dataset_id=$dataset_id)";
    }
    $query .= " GROUP BY f.chromosome, f.feature_type_id";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results;
    if ($returnArray) {
        $results = $sth->fetchall_arrayref({});
    }
    else {
        $results = $sth->fetchall_hashref([ 'chromosome', 'type_id' ]);
    }
    #print STDERR Dumper $results, "\n";

    return $results;
}

sub get_feature_count_summary {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;
    return unless ($dbh and $genome_id);

    # Execute query
    my $query = qq{
        SELECT f.feature_type_id AS type_id, ft.name as type_name, count(*) AS count
        FROM dataset_connector AS dc
        JOIN feature AS f ON (dc.dataset_id=f.dataset_id)
        JOIN feature_type AS ft ON (f.feature_type_id=ft.feature_type_id)
        WHERE (dc.genome_id=$genome_id)
        GROUP BY f.feature_type_id
    };
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_arrayref({});
    #print STDERR Dumper $results, "\n";

    return $results;
}

sub get_features {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;   # genome id (either this or dataset_id required)
    my $dataset_id = shift;  # dataset id
    return unless ($dbh and ($genome_id or $dataset_id));
    my $type_id = shift;     # optional feature type
    my $returnArray = shift; # optional flag to return result as raw array ref

    # Test for presence of primary_name
    my $hasPrimaryName = 0;
    my $query = qq{
        SELECT count(*)
        FROM dataset_connector AS dc
        JOIN feature AS f ON (f.dataset_id=dc.dataset_id)
        JOIN feature_type AS ft ON (ft.feature_type_id=f.feature_type_id)
        JOIN feature_name AS fn ON (fn.feature_id=f.feature_id)
        WHERE fn.primary_name=1
    };
    if ($genome_id) {
        $query .= " AND dc.genome_id=$genome_id";
    }
    else {
        $query .= " AND dc.dataset_id=$dataset_id";
    }
    if ($type_id) {
        $query .= " AND f.feature_type_id=$type_id";
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    ($hasPrimaryName) = $sth->fetchrow_array();

    # Execute query
    $query = qq{
        SELECT f.feature_id AS id, f.feature_type_id AS type_id, ft.name AS type_name,
            f.start AS start, f.stop AS stop, f.strand AS strand, f.chromosome AS chromosome,
            fn.name as name
        FROM dataset_connector AS dc 
        JOIN feature AS f ON (f.dataset_id=dc.dataset_id) 
        JOIN feature_type AS ft ON (ft.feature_type_id=f.feature_type_id)
        JOIN feature_name AS fn ON (fn.feature_id=f.feature_id)
        WHERE 1=1
    };
    if ($genome_id) {
        $query .= " AND dc.genome_id=$genome_id";
    }
    else {
        $query .= " AND dc.dataset_id=$dataset_id";
    }
    if ($type_id) {
        $query .= " AND f.feature_type_id=$type_id";
    }
    if ($hasPrimaryName) {
        $query .= " AND fn.primary_name=1";
    }

    warn $query;
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $results;
    if ($returnArray) {
        $results = $sth->fetchall_arrayref({});
    }
    else {
        $results = $sth->fetchall_hashref([ 'chromosome', 'id' ]);
    }
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_features_by_range { # for JBrowse::Annotation
    my %opts = @_;
    my $dbh       = $opts{dbh}; # database connection handle
    my $gid       = $opts{gid} || $opts{genome_id};
    my $chr       = $opts{chr};
    my $start     = $opts{start};
    my $stop      = $opts{stop} || $opts{end};
    my $feat_type = $opts{feat_type};
    my $dsid      = $opts{dsid} || $opts{dataset_id};
    die unless ($gid && defined $chr && defined $start && defined $stop);
    
    # mdb 4/24/14 - added fn.primary_name=1 constraint to keep from returning arbitrary name
    my $query = qq{
        SELECT l.start as locstart, l.stop as locstop, l.strand as locstrand, 
            ft.name as type, fn.name, l.location_id, f.start, f.stop, 
            f.feature_id, fn.primary_name, l.chromosome, f.chromosome, f.strand
            FROM genome g
            JOIN dataset_connector dc ON dc.genome_id = g.genome_id
            JOIN dataset d on dc.dataset_id = d.dataset_id
            JOIN feature f ON d.dataset_id = f.dataset_id
            JOIN location l ON f.feature_id = l.feature_id
            JOIN feature_name fn ON f.feature_id = fn.feature_id
            JOIN feature_type ft ON f.feature_type_id = ft.feature_type_id
            WHERE g.genome_id = $gid
                AND f.chromosome = '$chr'
                AND f.stop > $start AND f.start <= $stop
                AND ft.feature_type_id != 4
    };
    
    if ($feat_type) {
        $query .=  " AND ft.name = '$feat_type'";
    }
    # mdb added 4/21/14 issue 363
    if ($dsid) {
        $query .= " AND dc.dataset_id = $dsid";
    }
    
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_arrayref({});
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_feature_types {
    my $dbh = shift; # database connection handle
    return unless $dbh;
    
    # Execute query
    my $query = qq{
        SELECT ft.feature_type_id AS ftid, ft.name AS name, ft.description AS description
        FROM feature_type AS ft
    };
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_hashref(['ftid']);
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_feature_names {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;   # genome id
    my $dataset_id = shift;  # dataset id
    return unless ($dbh and ($genome_id or $dataset_id));
    
    # Execute query
    my $query = qq{
        SELECT f.feature_id AS fid, f.chromosome AS chr, f.feature_type_id AS ftid,
            fn.feature_name_id AS fnid, fn.name AS name, fn.primary_name AS primary_name,
            dc.dataset_id AS dsid
        FROM dataset_connector AS dc 
        JOIN feature AS f ON (f.dataset_id=dc.dataset_id) 
        JOIN feature_name AS fn ON (fn.feature_id=f.feature_id)
    };
    if ($genome_id) {
        $query .= " WHERE (dc.genome_id=$genome_id)";
    }
    else {
        $query .= " WHERE (dc.dataset_id=$dataset_id)";
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_hashref(['fid', 'fnid']);
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_feature_annotations {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;   # genome id
    my $dataset_id = shift;  # dataset id
    return unless ($dbh and ($genome_id or $dataset_id));
    
    # Execute query
    my $query = qq{
        SELECT f.feature_id AS fid, fa.feature_annotation_id AS faid,
            fa.annotation AS annotation, at.name AS atname, atg.name AS atgname
        FROM dataset_connector AS dc 
        JOIN feature AS f ON (f.dataset_id=dc.dataset_id) 
        JOIN feature_annotation AS fa ON (fa.feature_id=f.feature_id)
        LEFT JOIN annotation_type AS at ON (fa.annotation_type_id=at.annotation_type_id)
        LEFT JOIN annotation_type_group AS atg ON (at.annotation_type_group_id=atg.annotation_type_group_id)
    };
    if ($genome_id) {
        $query .= " WHERE (dc.genome_id=$genome_id)";
    }
    else {
        $query .= " WHERE (dc.dataset_id=$dataset_id)";
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_hashref(['fid', 'faid']);
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_locations {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;   # genome id
    my $dataset_id = shift;  # dataset id
    return unless ($dbh and ($genome_id or $dataset_id));
    
    # Execute query
    my $query = qq{
        SELECT f.feature_id AS fid, l.location_id AS lid, l.start AS start,
            l.stop AS stop, l.strand AS strand, l.chromosome AS chr
        FROM dataset_connector AS dc 
        JOIN feature AS f ON (f.dataset_id=dc.dataset_id) 
        JOIN location AS l ON (l.feature_id=f.feature_id)
    };
    if ($genome_id) {
        $query .= " WHERE (dc.genome_id=$genome_id)";
    }
    else {
        $query .= " WHERE (dc.dataset_id=$dataset_id)";
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_hashref(['fid', 'lid']);
    #print STDERR Dumper $results, "\n";
    
    return $results;
}

sub get_datasets {
    my $dbh = shift;         # database connection handle
    my $genome_id = shift;   # genome id
    
    # Execute query
    my $query = qq{
        SELECT d.dataset_id AS dataset_id, d.name AS name, d.description AS description,
            d.link AS link, d.version AS version, d.date AS date
        FROM dataset_connector AS dc
        JOIN dataset AS d ON (dc.dataset_id=d.dataset_id)
        WHERE genome_id=$genome_id
    };
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    # Fetch results
    my $results = $sth->fetchall_hashref(['dataset_id']);
    #print STDERR Dumper $results, "\n";
    
    return wantarray ? values %$results : $results;
}

sub get_dataset_ids {
	my $gid = shift;
	my $dbh = shift;
	$dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $gid);
}

# doesn't seem to be used anywhere
#sub get_chromosomes {
#    my $dbh = shift;         # database connection handle
#    my $genome_id = shift;   # genome id
#    my $dataset_id = shift;  # dataset id
#    return unless ($dbh and ($genome_id or $dataset_id));
#    
#    # Execute query
#    my $query = qq{
#        SELECT f.feature_id AS fid, f.start AS start, f.stop AS stop, 
#            f.chromosome AS chr
#        FROM dataset_connector AS dc 
#        JOIN feature AS f ON (f.dataset_id=dc.dataset_id) 
#        WHERE feature_type_id=4
#    };
#    if ($genome_id) {
#        $query .= " AND dc.genome_id=$genome_id";
#    }
#    else { # dataset_id
#        $query .= " AND dc.dataset_id=$dataset_id";
#    }
#    my $sth = $dbh->prepare($query);
#    $sth->execute();
#    my $results = $sth->fetchall_arrayref({});
#    #print STDERR Dumper $results, "\n";
#    
#    return wantarray ? @$results : $results;
#}

# doesn't seem to be used anywhere
#sub get_chromosomes_from_features {
#    my $dbh = shift;         # database connection handle
#    my $genome_id = shift;   # genome id
#    my $dataset_id = shift;  # dataset id
#    return unless ($dbh and ($genome_id or $dataset_id));
#    
#    # Execute query
#    my $query = qq{
#        SELECT f.chromosome AS chr
#        FROM dataset_connector AS dc 
#        JOIN feature AS f ON (f.dataset_id=dc.dataset_id) 
#    };
#    if ($genome_id) {
#        $query .= " WHERE dc.genome_id=$genome_id";
#    }
#    else { # dataset_id
#        $query .= " WHERE dc.dataset_id=$dataset_id";
#    }
#    my $sth = $dbh->prepare($query);
#    $sth->execute();
#    my $results = $sth->fetchall_hashref(['chr']);
#    #print STDERR Dumper $results, "\n";
#    
#    return wantarray ? keys %$results : [ keys %$results ];
#}

# Estimate table count (for large InnoDB tables)
sub get_table_count {
    my $dbh = shift;   # database connection handle
    my $table = shift; # table name
    my $where = shift; # optional where clause
    return unless ($dbh and $table);
    my $id = $table . '_id';
    
    my $query = "SELECT $id FROM $table " . ($where ? "WHERE $where" : '') . " ORDER BY $id DESC LIMIT 1";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my ($count) = $sth->fetchrow_array();
    
    return $count;
}

sub get_total_queries {
	my $dbh = shift;
	
	my $query = 'SHOW STATUS WHERE Variable_name="Queries"';
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $results = $sth->fetchall_hashref("Variable_name");
	
	return $results;
}

sub get_uptime {
	my $dbh = shift;

	my $query = 'SHOW STATUS WHERE Variable_name="Uptime"';
	my $sth = $dbh->prepare($query);
	$sth->execute();
	my $results = $sth->fetchall_hashref("Variable_name");

	return $results;
}

1;
