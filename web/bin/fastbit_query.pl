#!/usr/bin/perl

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!! AS OF AUGUST 2013 THIS FILE IS DEPRECATED, USE THIS MODULE INSTEAD:
#!!!! Accessory::Storage
#!!!! This file will be kept to preserve legacy Open Layers GenomeView
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use strict;

use Data::Dumper;
use Getopt::Long;
use DBI;
#use DBIxProfiler;
use CoGeX;
use CoGe::Accessory::Web;
use CGI;
use Benchmark;
use Cwd 'abs_path';

my $NUM_COL = 6;

my $t1 = new Benchmark;

my ( $exp_id, $coge_conf, $chr, $start, $stop );
GetOptions(    # FIXME for testing, remove someday
    "exp_id=i"    => \$exp_id,
    "chr=s"       => \$chr,
    "start=i"     => \$start,
    "stop=i"      => \$stop,
    "coge_conf=s" => \$coge_conf
);

# Load config file
unless ($coge_conf) {    # mdb added 3/26/13, tempfix for issue 61
    $coge_conf = abs_path($0);
    $coge_conf =~ s/bin\/fastbit_query\.pl//;
    $coge_conf .= 'coge.conf';
}
my $P       = CoGe::Accessory::Web::get_defaults($coge_conf);
my $DBNAME  = $P->{DBNAME};
my $DBHOST  = $P->{DBHOST};
my $DBPORT  = $P->{DBPORT};
my $DBUSER  = $P->{DBUSER};
my $DBPASS  = $P->{DBPASS};
my $CMDPATH = $P->{FASTBIT_QUERY};

# Parse URL
my $FORM = new CGI;
$exp_id = $FORM->param('exp_id') unless defined $exp_id;
$chr    = $FORM->param('chr')    unless defined $chr;
$chr    = $FORM->param('ref')    unless defined $chr;
$start  = $FORM->param('start')  unless defined $start;
$stop   = $FORM->param('stop')   unless defined $stop;
$stop   = $FORM->param('end')    unless defined $stop;

#TODO need some param checking here

# Connect to the database
my $connstr = "dbi:$P->{DB}:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
my $coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);
my $t2 = new Benchmark;

# Get experiment's storage path
my ($exp) = $coge->resultset("Experiment")->find($exp_id);
unless ($exp) {
    print STDERR "fastbit_query: experiment $exp_id not found\n";
    exit;
}

my $exp_storage_path = $exp->storage_path;
my $t3               = new Benchmark;

# Call FastBit to do query
# Issue 61: query string must contain a "." for fastbit to use consistent
# format (see jira thread for explanation) so "0.0=0.0" was added along with
# -v option.  Output parsing was modified accordingly for new output format.
#my $cmd = "$CMDPATH -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where chr='$chr' and start > $start and stop < $stop\" 2>&1"; # mdb removed 3/27/13 issue 61
my $cmd =
"$CMDPATH -v 1 -d $exp_storage_path -q \"select chr,start,stop,strand,value1,value2 where 0.0=0.0 and chr='$chr' and start > $start and stop < $stop order by start limit 999999999\" 2>&1"
  ;    # mdb added 3/27/13 issue 61
#print STDERR "$cmd\n";
my @cmdOut = qx{$cmd};
#print STDERR @cmdOut;
my $cmdStatus = $?;
die "Error executing command $CMDPATH ($cmdStatus)" if ( $cmdStatus != 0 );
my $t4 = new Benchmark;

# Convert FastBit output into JSON
my $results = '';
foreach (@cmdOut) {    # mdb rewritten 3/27/13 issue 61
    chomp;
    if (/^\"/) {       #if (/^\"$chr\"/) { # potential result line
        s/"//g;
        my @items = split(/,\s*/);
        next
          if ( @items != $NUM_COL )
          ;    # || $items[0] !~ /^\"?$chr/); # make sure it's a row output line
        for ( my $i = 0 ; $i < @items ; $i++ ) {
            $items[$i] = 1 if $items[$i] !~ /\w/;    # what's this for?
        }
        $results .= ',' if ($results);
        $results .= '[' . join( ',', map { "\"$_\"" } @items ) . ']';
    }
}

my $t5 = new Benchmark;
print "Content-type: text/json\n\n", qq{{"results" : [$results]}}, "\n";
#print STDERR qq{{"results" : [$results]}}."\n";

#print STDERR "Initialize CoGe: ". timestr( timediff( $t2, $t1 ) ),"\n",
#			 "Get experiment info from CoGe: ".timestr( timediff( $t3, $t2 ) ),"\n",
#			 "Get data from fastbit: ".timestr( timediff( $t4, $t3 ) ),"\n",
#			 "Assemble results: ". timestr( timediff( $t5, $t4 ) ),"\n";

exit;

__END__
