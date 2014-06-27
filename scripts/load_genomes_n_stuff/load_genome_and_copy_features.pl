#!/usr/bin/perl -w

use strict;
use CoGeX;
use Data::Dumper;
use Getopt::Long;

use vars qw($DEBUG $GO $db $user $pass $coge $dsgid $faa $fasta_genome_loader $replicate_annotations $seq_type_id $seq_type_name $seq_type_desc $restricted);

#../delete_dataset_group.pl -u xxx -db coge -pw xxx -delete_seqs -dsgid 19478
#./load_and_genome.pl  -u xxx -p xxx -db coge -dsgid 782 -faa ~/tmp/autog-test.faa -seq_type_id 3 -debug 1 -go 1

GetOptions ( "debug=s" => \$DEBUG,
             "go=s"    => \$GO,
             "database|db=s"=>\$db,
	     "user|u=s"=>\$user,
	     "password|pw|p=s"=>\$pass,
	     "dsgid=i"=>\$dsgid,
	     "faa|nt|fasta=s"=>\$faa,
	     "fasta_genome_loader=s"=>\$fasta_genome_loader,
	     "replicate_annotations=s"=>\$replicate_annotations,
	     "seq_type_name=s" => \$seq_type_name,
             "seq_type_desc=s" => \$seq_type_desc,
             "seq_type_id=i"=>\$seq_type_id, # masked50 == id 2
	     "restricted=i"=> \$restricted,
	   );

my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

$fasta_genome_loader = "/home/elyons/projects/CoGeX/scripts/load_genomes_n_stuff/fasta_genome_loader.pl" unless $fasta_genome_loader;
$replicate_annotations = "/home/elyons/projects/CoGeX/scripts/load_genomes_n_stuff/replicate_annotations.pl" unless $replicate_annotations;
$restricted = 0 unless $restricted;

unless ($coge && $dsgid && $faa)
  {
    print qq{
Welcome to $0

Purpose:  load a new sequence for an existing genome.  For example, you have created masked genome from an existing genome in CoGe.  Now you want to load that genome and map over the features from the existing genome in CoGe.

Usage:  $0  -db <coge database name> -u <coge database user> -pw <coge database password> -dsgid <coge database id for dataset_group to be master> -faa <fasta formated genome sequence> -seq_type_id <coge database id for genomic sequence type.  e.g. 1 is for unmasked>

Options:
   -go 1        |      Make the database calls.  Default 0

   -fasta_genome_loader  |   Path to CoGe's script fasta_genome_loader.pl for loading the fasta sequence.  Default: /home/elyons/projects/CoGeX/scripts/load_genomes_n_stuff/fasta_genome_loader.pl

   -replicate_annotations | Path to CoGe's script replicate_annotations.pl for copying features and annotations between two dataset groups in CoGe.  Default: /home/elyons/projects/CoGeX/scripts/load_genomes_n_stuff/replicate_annotations.pl

   -seq_type_id       |    coge data genomic_sequence_type id (id is required:  defaults to 1 (unmasked)

   -seq_type_name     |    coge data genomic_sequence_type name (this or id is required)

   -seq_type_desc     |    coge data genomic_sequence_type desc (optional)

   -restricted        |    mark genome as restricted (set to 1; default 0)
Copyright Eric Lyons 2012

};
    exit;
  }

unless (-r $faa)
  {
    print qq{
Error:  unable able to read fasta file: $faa
};
    exit;
  }

unless (-r $fasta_genome_loader)
  {
    print qq{
Error:  unable to read fasta_genome_loader: $fasta_genome_loader
};
    exit;
  }

my ($dsg) = $coge->resultset('DatasetGroup')->find($dsgid);

unless ($dsg)
  {
    print qq{
Error: Unable to create a datasetgroup object for $dsgid
};
    exit;
  }

$seq_type_id = get_seq_type(seq_type_id => $seq_type_id, seq_type_name => $seq_type_name, seq_type_desc => $seq_type_desc);

print "Loading genome\n";
my $new_dsgid = load_genome (faa=>$faa, dsg=>$dsg, stid=>$seq_type_id);
print "NEW DSGID: $new_dsgid\n" if $DEBUG;

add_annotations(dsgid1=>$dsgid, dsgid2=>$new_dsgid);

sub add_annotations
  {
    my %opts = @_;
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    my $cmd = $replicate_annotations;
    $cmd .= " -db " . $db;
    $cmd .= " -u " . $user;
    $cmd .= " -pw " . $pass;
    $cmd .= " -dsgid1 " . $dsgid1;
    $cmd .= " -dsgid2 " . $dsgid2;
    $cmd .= " -go " . $GO;
    $cmd .= " -debug ". $DEBUG;
    print "Running: $cmd\n";
    open (CMD, "$cmd |");
    while (<CMD>)
      {
	print $_ if $DEBUG;
      }
    close CMD;
  }

sub load_genome
  {
    my %opts = @_;
    my $faa = $opts{faa};
    my $dsg = $opts{dsg};
    my $stid = $opts{stid};
    my $cmd = $fasta_genome_loader;
    $cmd .= " -org_id ".$dsg->organism->id;
    my ($ds) = $dsg->datasets;
    $cmd .= " -source_id ".$ds->data_source->id;
    $cmd .= " -ds_version ". $dsg->version;
    $cmd .= " -seq_type_id ".$stid;
    $cmd .= " -u ". $user;
    $cmd .= " -pw ".$pass;
    $cmd .= " -db ".$db;
    $cmd .= " -restricted ". $restricted if $restricted;
    my $seq_dir = $dsg->file_path;
    $seq_dir =~ s/[^\/]*$//;
    $seq_dir =~ s/\d+\///g;
    $cmd .= " -sd ". $seq_dir;
    $cmd .= " -nt ". $faa;
    print "Running: ",$cmd,"\n";
    my $dsgid;
    if ($GO)
      {
	open (CMD, "$cmd |");
	while (<CMD>)
	  {
	    if (/dataset_group_id:\s+(\d+)/)
	      {
		$dsgid = $1;
		print "Captured dsgid: $dsgid\n" if $DEBUG;
	      }
	    print $_ if $DEBUG;;
	  }
	close CMD;
      }
    return $dsgid;
  }

sub get_seq_type
  {
    my %opts = @_;
    my $seq_type_id = $opts{seq_type_id};
    my $seq_type_name = $opts{seq_type_name};
    my $seq_type_desc = $opts{seq_type_desc};
    if ($seq_type_name && $GO)
      {
	print "Generating seq_type: $seq_type_name\n" if $DEBUG;
	my $gst = $coge->resultset('GenomicSequenceType')->find_or_create({name=>$seq_type_name,
									   description=>$seq_type_desc,
									  });
	$seq_type_id = $gst->id;
      }
    $seq_type_id = 1 unless $seq_type_id;  #default to unmasked sequence data
    return $seq_type_id;
  }
