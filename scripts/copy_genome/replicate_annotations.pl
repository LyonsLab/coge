#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;

use vars qw($DEBUG $db $user $pass $coge $dsgid1 $dsgid2);

$DEBUG = 0;

GetOptions(
	"database|db=s"   => \$db,
	"user|u=s"        => \$user,
	"password|pw|p=s" => \$pass,
	"dsgid1=i"        => \$dsgid1,
	"dsgid2=i"        => \$dsgid2,
);

print STDOUT "Starting $0 (pid $$)\n";

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect( $connstr, $user, $pass );

unless ( $coge && $dsgid1 && $dsgid2 )
{
	print qq{
Welcome to $0

Purpose:  clone gene models and annotations between two dataset groups.  Does check to see if there is a matching entry in the database to avoid duplicates.  This permits adding new annotations, updating annotations, etc.

Usage:  $0  -db <coge database name> -u <coge database user> -pw <coge database password> -dsgid1 <coge database id for dataset_group to be master> -dsgid2 <coge database id for dataset_group to be copied into>

Copyright Eric Lyons 2012

};
	exit;
}

my ($dsg1) = $coge->resultset('Genome')->find($dsgid1);
my ($dsg2) = $coge->resultset('Genome')->find($dsgid2);

unless ( $dsg1 && $dsg2 ) {
	print STDERR "Error:  unable to create a datasetgroup object for one of the IDs\n" if $DEBUG;
	exit;
}

my $ds1 = map_chr_to_ds($dsg1);
my $ds2 = map_chr_to_ds($dsg2);

unless ( scalar keys %$ds1 == scalar keys %$ds2 ) {
	print STDERR "Chromosome count mismatch!  Continuing.  Check errors below.\n" if $DEBUG;
}

foreach my $chr1 ( sort $dsg1->chromosomes ) {
	unless ( $ds2->{$chr1} ) {
		print STDERR "Chromosome $chr1 does not exist in dataset_group_2.  Please check the chromosome names if you think this is in error.  Skipping to next chromosome\n" if $DEBUG;
		next;
	}
	foreach my $ds_1 ($dsg1->datasets)
	  {
	    #print STDERR $ds_1->id,"\t", $ds_1->name,"\n";
	    foreach my $f1 ( $ds_1->features->search( { chromosome => $chr1 } ) )
	      {
	        #print STDERR $f1->type->name,"\n";
		my ($f2) = $ds2->{$chr1}->features->search(
							   {
							    chromosome      => $chr1,
							    feature_type_id => $f1->feature_type_id,
							    start           => $f1->start,
							    stop            => $f1->stop,
							    strand          => $f1->strand,
							   }
							  );
		if ($f2) {
		  print STDERR "Feature already exists.\n" if $DEBUG;
		}
		else {
		  ($f2) = replicate_feature( f1 => $f1, ds2 => $ds2->{$chr1} );
		}
		next unless $f2;    #can't do much without the annotation
		replicate_names( f1 => $f1, f2 => $f2 );
		replicate_annotations( f1 => $f1, f2 => $f2 );
	      }
	  }
}

sub replicate_annotations {
	my %opts = @_;
	my $f1   = $opts{f1};
	my $f2   = $opts{f2};
	print STDERR "Processing annotations\n" if $DEBUG;
	foreach my $anno ( $f1->annotations ) {
		my (@annos) = $f2->annotations->search(
			{
				annotation         => $anno->annotation,
				annotation_type_id => $anno->annotation_type_id
			}
		);
		if ( scalar @annos ) {
			print STDERR "\tAnnotation exists\n" if $DEBUG;
			next;
		}
		print STDERR "\tAdding annotation\n" if $DEBUG;
		$f2->add_to_feature_annotations(
			{
				annotation         => $anno->annotation,
				annotation_type_id => $anno->annotation_type_id
			}
		  );
	}
}

sub replicate_names {
	my %opts = @_;
	my $f1   = $opts{f1};
	my $f2   = $opts{f2};
	print STDERR "Processing names\n" if $DEBUG;
	foreach my $name ( $f1->feature_names ) {
		my @names =
		  $f2->feature_names->search(
			{ name => $name->name, description => $name->description } );
		if ( scalar @names ) {
			print STDERR "\tName exists\n" if $DEBUG;
			next;
		}
		print STDERR "\tAdding name\n" if $DEBUG;
		$f2->add_to_feature_names(
			{
				name         => $name->name,
				description  => $name->description,
				primary_name => $name->primary_name
			}
		  );
	}
}

sub replicate_feature {
	my %opts = @_;
	my $f1   = $opts{f1};
	my $ds2  = $opts{ds2};
	print STDERR "Creating feature\n" if $DEBUG;
	my ($f2) = $ds2->add_to_features(
		{
			chromosome      => $f1->chromosome,
			feature_type_id => $f1->feature_type_id,
			start           => $f1->start,
			stop            => $f1->stop,
			strand          => $f1->strand,
		}
	  );
	foreach my $loc ( $f1->locations ) {
		print STDERR "\tadding location\n" if $DEBUG;
		$f2->add_to_locations(
			{
				start      => $loc->start,
				stop       => $loc->stop,
				strand     => $loc->strand,
				chromosome => $loc->chromosome,
			}
		  );
	}
	return $f2;
}

sub map_chr_to_ds {
	my $dsg = shift;
	my %ds;
	foreach my $ds ( $dsg->datasets ) {
		foreach my $chr ( $ds->chromosomes ) {
		  $ds{$chr}=$ds;
		}
	}
	return \%ds;
}
