#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use CoGeX;

my ($help, $dsgid, $version, $type, $chr, $DEBUG, $user, $pass, $db);

GetOptions ("h|help" =>  \$help,
	    "dsgid=i" => \$dsgid,
	    "t|type=i"    => \$type,
	    "d|debug"     => \$DEBUG,
            "c|chr=s"     => \$chr,
	    "u=s" => \$user,
	    "p=s"=>\$pass,
	    "db=s"=>\$db
	    );

$help = 1 unless ($dsgid);
help() if $help;

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
my $coge = CoGeX->connect($connstr, $user, $pass );

my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
my $ft = $coge->resultset('FeatureType')->find($type) if $type;
my $search = {};
$search->{feature_type_id} = $ft->id if $ft;

print join ("\t", qw(TYPE CHR START STOP STRAND NAME)),"\n";

foreach my $ds ($dsg->datasets)
  {
    print STDERR "working on ",$ds->name,"\n";
    foreach my $feat ($ds->features($search))
      {
	print join ("\t", $feat->type->name, $feat->chr, $feat->start, $feat->stop, $feat->strand, join (", ", $feat->names)),"\n";
      }
  }

sub help
  {
    print qq{
Welcome to $0

This program generates a table of all the features in a genome

Options:

  -h | -help               prints this message

  -dsgid                   CoGe dataset group database id

  -t | -type               specify the type_id of feature from which to extract the sequence.
                           e.g.: 3 (for CDS)

  -u                        database user name

  -p                        database user password

  -db                       database database name

  -d | -debug              print debugging messages

};
exit(0);
  }
