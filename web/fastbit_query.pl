#!/usr/bin/perl

use warnings; # FIXME remove
use strict; # FIXME remove

use lib '/home/mbomhoff/CoGeX/lib'; # FIXME remove

#use Roman;
use Data::Dumper;
use Getopt::Long;
use DBI;
#use DBIxProfiler;
use CoGeX;
use CoGe::Accessory::Web;
use CGI;

my $cmdPath = '/home/mbomhoff/fastbit-ibis1.3.0/install/bin'; # FIXME get from conf file
my ($exp_id, $coge_pref, $chr, $start, $stop);

GetOptions ( # FIXME remove
	    "exp_id=i" 			=> \$exp_id,
			"chr=s" 				=> \$chr,
			"start=i" 			=> \$start,
			"stop=i" 				=> \$stop,
			"coge_pref=s" 	=> \$coge_pref   
	   );

# Load config file
my ($DBNAME, $DBHOST, $DBPORT, $DBUSER, $DBPASS);
my $P = CoGe::Accessory::Web::get_defaults($coge_pref);
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};

# Parse URL
my $FORM = CGI->new;
$exp_id = $FORM->param('exp_id')	unless defined $exp_id;
$chr 		= $FORM->param('chr') 		unless defined $chr;
$start 	= $FORM->param('start') 	unless defined $start;
$stop 	= $FORM->param('stop') 		unless defined $stop;

# Connect to the database
my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

# Get experiment's storage path
my ($exp) = $coge->resultset("Experiment")->find($exp_id);
# FIXME handle experiment not found case
my $exp_storage_path = $exp->storage_path;

# Call FastBit to do query
my $cmd = "$cmdPath/ibis -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where chr=$chr and start > $start and stop < $stop\" 2>&1";
my $cmdOut = qx{$cmd};
my $cmdStatus = $?;
die "Error executing command $cmdStatus" if ($cmdStatus != 0);

# Convert FastBit output into JSON
my ($numHits, $queryStr, $header, $results);

my @lines = split("\n", $cmdOut);
foreach (@lines) 
	{
		if (/(.+) \(with counts\)/) { # first line
			$header = $1;
			$header =~ s/\"?([\w\.\+\-]+)\"?/\"$1\"/g;
		}
		elsif (/^doQuery:: evaluate\((.+)\) produced (\d+) hit/) { # last line
			$queryStr = $1;
			$numHits = $2;
		}
		elsif (my @fields = /\"?([\w\.\+\-]+)\"?\,/g) { # result line
			$results .= ',' if ($results);
			$results .= '[' . join(',', map {"\"$_\""} @fields) . ']';
		}
	}

if (defined $numHits and defined $queryStr) 
	{
		$header = '' unless ($header);
		$results = '' unless ($results);
		print qq{{"query" : "$queryStr","result_count" : "$numHits","header" : [$header],"results" : [$results]}};
	}

exit;

__END__
