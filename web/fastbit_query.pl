#!/usr/bin/perl

use strict;

use Data::Dumper;
use Getopt::Long;
use DBI;
#use DBIxProfiler;
use CoGeX;
use CoGe::Accessory::Web;
use CGI;
use Benchmark;

my ($exp_id, $coge_conf, $chr, $start, $stop);
my $t1    = new Benchmark;

GetOptions ( # FIXME for testing, remove someday
	    "exp_id=i" 			=> \$exp_id,
			"chr=s" 				=> \$chr,
			"start=i" 			=> \$start,
			"stop=i" 				=> \$stop,
			"coge_conf=s" 	=> \$coge_conf   
	   );

# Load config file
my $P = CoGe::Accessory::Web::get_defaults($coge_conf);
my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $CMDPATH = $P->{FASTBIT_QUERY};

# Parse URL
my $FORM = new CGI;
$exp_id = $FORM->param('exp_id')	unless defined $exp_id;
$chr 		= $FORM->param('chr') 		unless defined $chr;
$start 	= $FORM->param('start') 	unless defined $start;
$stop 	= $FORM->param('stop') 		unless defined $stop;
# Connect to the database
my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);
my $t2    = new Benchmark;
# Get experiment's storage path
my ($exp) = $coge->resultset("Experiment")->find($exp_id);
# FIXME handle experiment not found case
my $exp_storage_path = $exp->storage_path;
my $t3    = new Benchmark;
# Call FastBit to do query
my $cmd = "$CMDPATH -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where chr=$chr and start > $start and stop < $stop\" 2>&1";
print STDERR $cmd;
my $cmdOut = qx{$cmd};
my $cmdStatus = $?;
die "Error executing command $CMDPATH ($cmdStatus)" if ($cmdStatus != 0);
my $t4    = new Benchmark;
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
		elsif (my @fields = /\"?([\w\.\+\-]*)\"?\,/g) { # result line
		  s/"//g;
		  s/, /,/g;
		  my @items = split/,/;
		  for (my $i =0; $i<@items; $i++)
		    {
		      $items[$i] = 1 if $items[$i] !~ /\w/;
		    }
		  $results .= ',' if ($results);
		  $results .= '[' . join(',', map {"\"$_\""} @items) . ']';
		}
	}
my $t5    = new Benchmark;
print "Content-type: text/json\n\n";
if (defined $numHits and defined $queryStr) 
	{
		$header = '' unless ($header);
		$results = '' unless ($results);
		print qq{{"query" : "$queryStr","result_count" : "$numHits","header" : [$header],"results" : [$results]}};
	}
print STDERR "Initialize CoGe: ". timestr( timediff( $t2, $t1 ) ),"\n";
print STDERR "Get experiment info from CoGe: ".timestr( timediff( $t3, $t2 ) ),"\n";
print STDERR "Get data from fastbit: ".timestr( timediff( $t4, $t3 ) ),"\n";
print STDERR "Assemble results: ". timestr( timediff( $t5, $t4 ) ),"\n";


exit;

__END__
