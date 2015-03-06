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
use Data::Dumper;

BEGIN {
    use vars qw ($VERSION $MAX_FETCH_ROWS @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    $MAX_FETCH_ROWS = 10*1000;
    @ISA     = qw (Exporter);
    @EXPORT = qw( 
        get_table get_user_access_table get_experiments get_distinct_feat_types
    );
}

sub get_table {
    my $dbh = shift;         # database connection handle
    my $table = shift;       # table name
    my $hash_fields = shift; # array ref of field names to use as hash keys
    my $conditions = shift;  
    $hash_fields = [$table.'_id'] unless $hash_fields;
    
    # Build query string
    my $query = "SELECT * FROM $table";
    if ($conditions) {
        $query .= ' WHERE ' . join(' AND ', map { $_.'='.$conditions->{$_} } keys %$conditions);
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
    
    my %combined = (%$results1, %$results2, %$results3);
    #print STDERR Dumper \%combined, "\n";
    
    return \%combined;
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