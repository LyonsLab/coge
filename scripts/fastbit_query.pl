#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;
use Getopt::Long;
use DBI;
#use DBIxProfiler;
use CoGeX;
use CoGe::Accessory::Web;

use vars qw($DEBUG $coge);

my ($exp_id, $coge_conf, $chr, $start, $stop, $P, $DBNAME, $DBHOST, $DBPORT, $DBUSER, $DBPASS, $CMDPATH);

GetOptions (
			"debug=s" 			=> \$DEBUG,

	    "exp_id=i" 			=> \$exp_id,

	    "database|db=s" => \$DBNAME,
	    "user|u=s" 			=> \$DBUSER,
			"password|pw=s" => \$DBPASS,
			"host|h=s" 			=> \$DBHOST,
			"port|p=i" 			=> \$DBPORT,

			"chr=s" 				=> \$chr,
			"start=i" 			=> \$start,
			"stop=i" 				=> \$stop,

      "coge_conf=s" 	=> \$coge_conf
	   );

$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on
#print STDERR "Running $0\n";
print help() unless (defined $exp_id and defined $chr and defined $start and defined $stop);

# Load config file
if ($coge_conf) {
	$P = CoGe::Accessory::Web::get_defaults($coge_conf);
	$DBNAME = $P->{DBNAME};
	$DBHOST = $P->{DBHOST};
	$DBPORT = $P->{DBPORT};
	$DBUSER = $P->{DBUSER};
	$DBPASS = $P->{DBPASS};
	$CMDPATH = $P->{FASTBIT_QUERY};
}
$DBHOST = 'localhost' unless ($DBHOST);
$DBPORT = 3307 unless (defined $DBPORT);
$CMDPATH = '/usr/local/bin/ibis' unless ($CMDPATH);

# Connect to the database
my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

# Get experiment's storage path
my ($exp) = $coge->resultset("Experiment")->find($exp_id);
print help("Couldn't find experiment") unless ($exp);

my $exp_storage_path = $exp->storage_path;

# Call FastBit to do query
my $cmd = "$CMDPATH -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where chr=$chr and start > $start and stop < $stop\" 2>&1";
my $cmdOut = qx{$cmd};
my $cmdStatus = $?;
die "Error executing command $cmdStatus" if ($cmdStatus != 0);
print $cmdOut;

exit;

#-------------------------------------------------------------------------------

sub help
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
        coge_conf     configuration file
};
	exit;
}

__END__
