#!/usr/bin/perl -w

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!! AS OF AUGUST 2013 THIS FILE IS DEPRECATED, USE THESE MODULES INSTEAD:
#!!!! Accessory::Storage or Services::Sequence (Web Service)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use strict;
use warnings;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGeX;
use CoGe::Accessory::Web;
no warnings 'redefine';

delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

my $P    = CoGe::Accessory::Web::get_defaults();
my $FORM = new CGI;

my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};
my $connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
my $coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

my $chr    = $FORM->param('chr');
my $start  = $FORM->param('start');
my $stop   = $FORM->param('stop');
my $strand = $FORM->param('strand');
my $dsgid  = $FORM->param('dsgid');
$dsgid = $FORM->param('gid') unless $dsgid;

my $usage = qq{
USAGE:
  $0?chr=1;dsgid=1;start=1;stop=10;strand=1

chr  :   chromosome
dsgid:   database id for dataset group
start:   start position of sequence (DEFAULT = 1)
stop:    stop position of sequence (DEFAULT = end of chromosome)
strand:  strand '1' or '-1' (DEFAULT = 1)

};

print $FORM->header, "\n";

$dsgid =~ s/\s+//g if $dsgid;
unless ( $dsgid && $dsgid =~ /^\d+$/ ) {
    print
"No Dataset Group Database ID specified.  Can't retrieve sequence without this!\n";
    print $usage;
    exit;
}

unless ($chr) {
    print qq{
Can't retrieve a sequence without a chromosome!  Please specify a chromosome!
  };
    print $usage;
    exit;
}

my ($dsg) = $coge->resultset('Genome')->find($dsgid);
unless ($dsg) {
    print "Unable to retrieve a valid dataset group object for id $dsgid\n";
    print $usage;
    exit;
}

unless ( $dsg->sequence_length($chr) ) {
    print "Chromosome '$chr' does not exist for this dataset group\n";
    print $usage;
    exit;
}

$start = 1 unless $start && $start =~ /^\d+$/;
$start = 1 if $start < 1;
$strand = 1 unless $strand;
$stop = $dsg->last_chromosome_position($chr) unless $stop && $stop =~ /^\d+$/;

my $seq = $dsg->get_seq(
    chr          => $chr,
    start        => $start,
    stop         => $stop,
    strand       => $strand,
    storage_path => $P->{SEQDIR}
);
if ($seq) {
    print $seq;
}
else {
    print "Unable to retrieve sequence.  Please check parameters\n";
    print $usage;
}
