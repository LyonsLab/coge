#!/usr/bin/perl
#-------------------------------------------------------------------------------
# Usage:
#    ./audit_datasets.pl -database XXXXXXX -user XXXXXXX -password XXXXXXX
#
#	perl -I /home/mbomhoff/perl5 scripts/audit_datasets.pl -db coge_matt -u XXXXXXX -p XXXXXXX
#
# Scan the datasets and report inconsistencies.
#
#-------------------------------------------------------------------------------

use strict;
use DBI;
use CoGeX;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $coge);

my ($db, $user, $pass, $cleanup);
GetOptions (
	"debug=s"			=> \$DEBUG,
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw|p=s"	=> \$pass,
);
die "Missing DB params\n" unless ($db and $user and $pass);

$| = 1;
$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on

print STDERR "Running $0\n";
print STDERR "PERL5LIB=" . $ENV{PERL5LIB} . "\n";

#-------------------------------------------------------------------------------
# Connect to database
#-------------------------------------------------------------------------------

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $dbh = $coge->storage->dbh;

#-------------------------------------------------------------------------------
# Audit
#-------------------------------------------------------------------------------

print STDERR "Verifying datasets ---------------------------------------------\n";
foreach my $ds (sort {$b->id <=> $a->id} ($coge->resultset('Dataset')->all)) {
	my $count = 0;
	print STDERR "Testing dsid ", $ds->id, " ", $ds->date, "\n";
	foreach my $f ($ds->features) {
	    if ($f->start =~ /[a-zA-Z]/ || $f->stop =~ /[a-zA-Z]/ || 
	       ($f->strand < -1 and $f->strand > 1) || $f->start > $f->stop) 
	    {
	        print STDERR "   Bad feature ", $f->id, "\n";
	        goto NEXT;
	    }
	    
	    foreach my $loc ($f->locations) {
	        if ($loc->start =~ /[a-zA-Z]/ || $loc->stop =~ /[a-zA-Z]/ || 
               ($loc->strand < -1 and $loc->strand > 1) || $loc->start > $loc->stop) 
            {
                print STDERR "   Bad location ", $loc->id, " in fid ", $f->id, "\n";
                goto NEXT;
            }
	    }
	    
	    #last if ($count++ > 1000);
	}
	NEXT:
}

#-------------------------------------------------------------------------------
# All done!
#-------------------------------------------------------------------------------

print STDERR "All done!\n";
exit;
