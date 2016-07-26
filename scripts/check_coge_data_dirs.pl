#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./check_coge_data_dirs.pl -db XXXXXXX -u XXXXXXX -p XXXXXXX
#
# Scan for coge_data directories for all users and report missing directories.
#
#-------------------------------------------------------------------------------

use strict;
use Getopt::Long;
use IPC::System::Simple qw(capture system $EXITVAL EXIT_ANY);
use CoGeX;

my ($db, $user, $pass);
GetOptions (
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw|p=s"	=> \$pass
);
die "Missing DB params\n" unless ($db and $user and $pass);

$| = 1;

print STDERR "Running $0\n";
print STDERR "PERL5LIB=" . $ENV{PERL5LIB} . "\n";

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
my $DB = CoGeX->connect($connstr, $user, $pass);

# Verify users
print STDERR "Verifying users' coge_data directories ...\n";
foreach my $user ($DB->resultset('User')->all) {
    my $data_path = '/iplant/home/'.$user->name.'/coge_data';
    my $cmd = "ils -l '$data_path' 2>&1";
    my @ils = capture( EXIT_ANY, $cmd );
    if ($EXITVAL) {
        print STDERR $user->name, "\n";
    }
}

print STDERR "All done!\n";
exit;
