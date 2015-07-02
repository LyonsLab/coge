package CoGeDBI;

=head1 NAME

CoGeDBI

=head1 SYNOPSIS

  This object is a low-level API to the CoGe database.  It was written for
  high-performance cases where the CoGeX ORM is too slow.

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
use Benchmark;

BEGIN {
    use vars qw ($VERSION $MAX_FETCH_ROWS @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    $MAX_FETCH_ROWS = 10*1000;
    @ISA     = qw (Exporter);
    @EXPORT = qw( 
        get_table get_user_access_table get_experiments get_distinct_feat_types
        get_genomes_for_user get_experiments_for_user get_lists_for_user
        get_groups_for_user get_group_access_table
    );
}

sub get_table {
    my $dbh = shift;         # database connection handle
    my $table = shift;       # table name
    my $hash_fields = shift; # array ref of field names to use as hash keys
    my $conditions = shift;
    my $search_string = shift;
    $hash_fields = [$table.'_id'] unless $hash_fields;
    $search_string = "=" unless $search_string;
    
    # Build query string
    my $query = "SELECT * FROM $table";
    if ($conditions) {
        $query .= ' WHERE ' . join(' AND ', map { $_.$search_string.$conditions->{$_} } keys %$conditions);
    }
    #print STDERR $query, "\n";
    
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
    
    my $combined = merge($results1, $results2, $results3);
    #print STDERR Dumper \%combined, "\n";
    
    return $combined;
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
    
    my $combined = merge($results1, $results2);
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
    my $query = qq{
        SELECT g.genome_id AS id, g.name AS name, g.description AS description, 
            g.version AS version, g.restricted AS restricted, g.deleted AS deleted, 
            g.date AS date, o.name AS organism, ds.date AS dataset_date, uc.role_id as role_id
        FROM user_connector AS uc 
        JOIN genome AS g ON (uc.child_id=g.genome_id) 
        JOIN dataset_connector AS dsc ON (dsc.genome_id=g.genome_id) 
        JOIN dataset AS ds ON (ds.dataset_id=dsc.dataset_id) 
        JOIN organism AS o ON (o.organism_id=g.organism_id) 
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=2 
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    #my $t1    = new Benchmark;
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results1 = $sth->fetchall_hashref(['id']);
    #my $t2    = new Benchmark;
    #print STDERR 'query1: ', timestr( timediff( $t2, $t1 ) ), "\n";
    #print STDERR Dumper $results1, "\n";

    # Get list genome connections
    $query = qq{
        SELECT g.genome_id AS id, g.name AS name, g.description AS description, 
            g.version AS version, g.restricted AS restricted, g.deleted AS deleted, 
            g.date AS date, o.name AS organism, ds.date AS dataset_date, uc.role_id as role_id
        FROM user_connector AS uc 
        JOIN list AS l ON (uc.child_id=l.list_id) 
        JOIN list_connector AS lc ON (lc.parent_id=l.list_id)
        JOIN genome AS g ON (lc.child_id=g.genome_id)
        JOIN dataset_connector AS dsc ON (dsc.genome_id=g.genome_id) 
        JOIN dataset AS ds ON (ds.dataset_id=dsc.dataset_id) 
        JOIN organism AS o ON (o.organism_id=g.organism_id) 
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=1 AND lc.child_type=2
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $results2 = $sth->fetchall_hashref(['id']);
    #print STDERR Dumper $results2, "\n";
    
    my $combined = merge($results1, $results2);
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
            e.date AS date, uc.role_id as role_id
        FROM user_connector AS uc 
        JOIN experiment AS e ON (uc.child_id=e.experiment_id) 
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=3
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results1 = $sth->fetchall_hashref(['id']);
    #print STDERR Dumper $results1, "\n";

    # Get list experiment connections
    $query = qq{
        SELECT e.experiment_id AS id, e.name AS name, e.description AS description, 
            e.version AS version, e.restricted AS restricted, e.deleted AS deleted, 
            e.date AS date, uc.role_id as role_id
        FROM user_connector AS uc 
        JOIN list AS l ON (uc.child_id=l.list_id) 
        JOIN list_connector AS lc ON (lc.parent_id=l.list_id)
        JOIN experiment AS e ON (lc.child_id=e.experiment_id)
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=1 AND lc.child_type=3 
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    $sth = $dbh->prepare($query);
    $sth->execute();
    my $results2 = $sth->fetchall_hashref(['id']);
    #print STDERR Dumper $results2, "\n";
        
    my $combined = merge($results1, $results2);
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
            l.restricted AS restricted, l.deleted AS deleted, lt.name AS type_name,
            l.date AS date, uc.role_id as role_id
        FROM user_connector AS uc 
        JOIN list AS l ON (uc.child_id=l.list_id) 
        JOIN list_type AS lt ON (lt.list_type_id=l.list_type_id)
        WHERE ((uc.parent_type=5 AND uc.parent_id=$user_id) OR (uc.parent_type=6 AND uc.parent_id IN ($group_str))) 
            AND uc.child_type=1 
    };
    $query .= " AND uc.role_id in ($role_id_str)" if $role_id_str;
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my $results = $sth->fetchall_arrayref({});
    #print STDERR Dumper $results, "\n";

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

1;