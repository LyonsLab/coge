#------------------------------------------------------------------------------
# Copy missing entries from source sqlite db to destination sqlite db.
# For when coge crashed (and we didn't have a backup of the JEX db) and we
# switched to geco, then needed to switch back to coge.
#------------------------------------------------------------------------------
use DBI;
use Data::Dumper;

if (@ARGV < 2) {
    print STDERR "Usage: perl sync_jex_db.pl <src_sqlite> <dest_sqlite>\n";
    exit;
}

my $src_file = $ARGV[0];
my $dest_file = $ARGV[1];

my $src_dbh = DBI->connect("dbi:SQLite:dbname=$src_file","","");
my $dest_dbh = DBI->connect("dbi:SQLite:dbname=$dest_file","","");

my $src_workflows = load_workflows($src_dbh);
print STDERR "Source workflows: ", scalar(keys $src_workflows), "\n";
#print STDERR Dumper $src_workflows, "\n";

my $dest_workflows = load_workflows($dest_dbh);
print STDERR "Destination workflows: ", scalar(keys $dest_workflows), "\n";
#print STDERR Dumper $dest_workflows, "\n";

my $added = 0;
foreach my $wid (keys $src_workflows) {
    next if ($dest_workflows->{$wid}); # workflow already exists in both
    add_workflow($dest_dbh, $src_workflows->{$wid});
    $added++;
}

print "$added workflows added to dest db\n";

exit;
#------------------------------------------------------------------------------

sub load_workflows {
    my $dbh = shift;
    
    my $query = "SELECT * FROM workflows";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    
    my $workflows = $sth->fetchall_hashref('id');
    return $workflows;
}

sub add_workflow {
    my $dbh = shift;
    my $workflow = shift;
    
    my $query = 'INSERT INTO workflows (id, name, log, jobs, submitted, completed, priority, status) '
                . 'VALUES (?, ?, ?, ?, ?, ?, ?, ?)';
    my $sth = $dbh->prepare($query);
    $sth->execute($workflow->{id}, $workflow->{name}, $workflow->{log}, $workflow->{jobs}, $workflow->{submitted}, $workflow->{completed}, $workflow->{priority}, $workflow->{status});
}
