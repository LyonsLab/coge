#!/usr/bin/perl -w

use strict;
use CoGeX;
use Getopt::Long;

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

my ($name, $desc, $version, $org_id, $gst_id, $file_path);

GetOptions (
	    "name|n=s"=>\$name,
	    "desc|d=s"=>\$desc,
	    "org_id|o=i"=>\$org_id,
	    "gst_id|t=i"=>\$gst_id,
	    "file|f=s"=>\$file_path,
	    "version|v=s"=>\$version,
	   );

$version = 1 unless defined $version;
$gst_id = 1 unless defined $gst_id;

unless ($org_id && $file_path)
  {
    usage();
    exit;
  }

my $dsg = $coge->resultset('DatasetGroup')->create({
							    name=>$name,
							    description=>$desc,
							    organism_id=>$org_id,
							    version=>$version,
							    genomic_sequence_type_id=>$gst_id,
							    file_path=>$file_path
							   });

if ($dsg)
  {
    print "Created dataset group with database id ",$dsg->id,"\n";
  }
else
  {
    print "Error occurred, did not create dataset group object or database entry\n";
    exit;
  }

$/ = "\n>";
open (IN, $file_path);
while (<IN>)
  {
    s/>//g;
    my ($name, $seq) = split /\n/,$_,2;
    $seq =~ s/>//g;
    $seq =~ s/\n//g;
    my @name = split/\|/, $name if $name =~ /\|/;
    $name = $name[1] if $name[1];
    $coge->resultset('GenomicSequence')->find_or_create({
							 sequence_length => length($seq),
							 chromosome => $name,
							 dataset_group_id => $dsg->id,
							});
  }
close IN;

sub usage
  {
    print qq{
Welcome to $0.
This program creates a dataset group entry in CoGe's database.

Options:

-name      | -n            name of dataset group (optional)

-desc      | -d            description of dataset group (optional)

-org_id    | -o            database id of the organism to which this dataset group belongs
                           MANDATORY!  You can get this value from CoGe's web app OrganismView

-gst_id    | -t            genomic sequence type is of the dataset group.  1 is usualy unmasked
                           sequence, 2 is 50x repeat masked, and then there are others.
                           DEFAULT: 1

-file      | -f            path to the genomic sequence file in fasta format.  Fasta headers must
                           be the names of the chromosomes (or contigs)
                           MANDATORY!  Example: /opt/apache/CoGe/data/genomic_sequence/0/0/0/1/1.faa

-version   | -v            version of dataset group.  DEFAULT: 1

};
  }
