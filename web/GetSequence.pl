#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGeX;
$ENV{PATH} = "/opt/apache/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

my $FORM = new CGI;
my  $coge = CoGeX->dbconnect();

my $chr = $FORM->param('chr');
my $start = $FORM->param('start');
my $stop = $FORM->param('stop');
my $strand = $FORM->param('strand');
my $dsgid = $FORM->param('dsgid');


my $usage = qq{
USAGE:
  $0?chr=1;dsgid=1;start=1;stop=10;strand=1

chr  :   chromosome
dsgid:   database id for dataset group
start:   start position of sequence (DEFAULT = 1)
stop:    stop position of sequence (DEFAULT = end of chromosome)
strand:  strand '1' or '-1' (DEFAULT = 1)

};
$dsgid =~ s/\s+//g if $dsgid;
print $FORM->header,"\n";
unless ($dsgid && $dsgid =~ /^\d+$/)
  {
    print "No Dataset Group Database ID specified.  Can't retrieve sequence without this!\n";
    print $usage;
    exit;
  }

unless ($chr)
  {
    print qq{
Can't retrieve a sequence without a chromosome!  Please specify a chromosome!
  };
    print $usage;
    exit;
  }
my ($dsg) = $coge->resultset('DatasetGroup')->find($dsgid);
unless ($dsg)
  {
    print "Unable to retrieve a valid dataset group object for id $dsgid\n";
    print $usage;
    exit;
  }


$start = 1 unless $start;
$start = 1 if $start < 1;
$strand = 1 unless $strand;
$stop = $dsg->last_chromosome_position($chr) unless $stop;

my $seq = $dsg->get_seq(chr=>$chr, start=>$start, stop=>$stop, strand=>$strand);
if ($seq)
  {
    print $seq;
  }
else
  {
    print "Unable to retrieve sequence.  Please check parameters\n";
    print $usage;
  }
