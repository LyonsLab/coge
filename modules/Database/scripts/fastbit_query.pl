#!/usr/bin/perl

use warnings;
use strict;

use lib '/home/mbomhoff/CoGeX/lib';

use Data::Dumper;
use Getopt::Long;
use DBI;
#use DBIxProfiler;
use CoGeX;
use CoGe::Accessory::Web;

use vars qw($DEBUG $coge);

my $cmdPath = '/home/mbomhoff/fastbit-ibis1.3.0/install/bin';

my ($exp_id, $db, $user, $pass, $host, $port, $coge_pref, $chr, $start, $stop);

GetOptions ( 
			"debug=s" 			=> \$DEBUG,
	    
	    "exp_id=i" 			=> \$exp_id,

	    "database|db=s" => \$db,
	    "user|u=s" 			=> \$user,
			"password|pw=s" => \$pass,
			"host|h=s" 				=> \$host,
			"port|p=i" 			=> \$port,

			"chr=s" 				=> \$chr,
			"start=i" 			=> \$start,
			"stop=i" 				=> \$stop,
			
			"coge_pref=s" 	=> \$coge_pref   
	   );

$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on
#print STDERR "Running $0\n";
print help() unless ($exp_id and $chr and $start and $stop);

# Load config file
my ($P, $DBNAME, $DBHOST, $DBPORT, $DBUSER, $DBPASS);
if ($coge_pref) {
	$P = CoGe::Accessory::Web::get_defaults($coge_pref);
	$db = $P->{DBNAME};
	$host = $P->{DBHOST};
	$port = $P->{DBPORT};
	$user = $P->{DBUSER};
	$pass = $P->{DBPASS};
}
$host = 'localhost' unless ($host);
$port = 3307 unless ($port);

# Connect to the database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port";
$coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

# Get experiment's storage path
my ($exp) = $coge->resultset("Experiment")->find($exp_id);
print help("Couldn't find experiment") unless ($exp);

my $exp_storage_path = $exp->storage_path;

# Call FastBit to do query
my $cmd = "$cmdPath/ibis -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where chr=$chr and start > $start and stop < $stop\" 2>&1";
my $cmdOut = qx{$cmd};
my $cmdStatus = $?;
die "Error executing command $cmdStatus" if ($cmdStatus != 0);
print $cmdOut;

exit;

#-------------------------------------------------------------------------------

sub help # FIXME
{
	my $msg = shift;
	
	print $msg if ($msg);
	
	print qq
	{
    Options:
        exp_id        experiment ID
        database|db   database name
        user|u        database username
        password|pw   database password
        chr           chromosome
        start         start position
        stop          stop position
		    coge_prefs    configuration file
};
	exit;
}

__END__
