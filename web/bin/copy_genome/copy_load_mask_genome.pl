#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web;
use Data::Dumper;
use Getopt::Long;
use File::Path;

use vars qw($DEBUG $GO $conf_file $coge $gid $restricted $P $mask $uid);


GetOptions ( "debug=s" => \$DEBUG,
             "go=s"    => \$GO,
	     "conf_file|cf=s" => \$conf_file,
	     "gid=i"=>\$gid,
	     "restricted|r=i"=> \$restricted,
	     "mask|m"=> \$mask,
	     "uid=i"=> \$uid, #coge user id to whom the genome will be assigned
	   );

$P = CoGe::Accessory::Web::get_defaults($conf_file);

unless ($P && $P->{DBNAME}) {usage();}

my $TEMPDIR = $P->{TEMPDIR}."copy_genome";
mkpath( $TEMPDIR,    1, 0777 );

my $DBNAME = $P->{DBNAME};
my $DBHOST = $P->{DBHOST};
my $DBPORT = $P->{DBPORT};
my $DBUSER = $P->{DBUSER};
my $DBPASS = $P->{DBPASS};

my $connstr = "dbi:mysql:dbname=".$DBNAME.";host=".$DBHOST.";port=".$DBPORT;
$coge = CoGeX->dbconnect(db_connection_string=>$connstr, db_name=>$DBUSER, db_passwd=>$DBPASS );

#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my $fasta_genome_loader = $P->{FASTA_GENOME_LOADER};
my $replicate_annotations = $P->{REPLICATE_ANNOTATIONS};
my $windowmasker = $P->{WINDOWMASKER};
my $hard_mask = $P->{HARD_MASK};


unless ($gid)
  {
    usage();
  }

unless (-r $fasta_genome_loader)
  {
    print STDERR qq{
Error:  unable to read fasta_genome_loader: $fasta_genome_loader
};
    exit;
  }



my ($dsg) = $coge->resultset('Genome')->find($gid);
$restricted = $dsg->restricted unless defined $restricted; #use whatever is specified by the original genome unless user specifies otherwise
unless ($dsg)
  {
    print STDERR qq{
Error: Unable to create a datasetgroup object for $gid
};
    exit;
  }


my $new_faa;
my $stid;
if ($mask)
  {
    $new_faa = mask_genome(dsg=>$dsg);
    $stid = get_sequence_type_id_for_windowmasker();
  }
else
  {
    $new_faa = get_and_clean_sequence($dsg->file_path);
    $stid = $dsg->genomic_sequence_type_id;
  }


print STDERR "Loading genome\n";
my $new_gid = load_genome (faa=>$new_faa, dsg=>$dsg, stid=>$stid);
print STDERR "NEW DSGID: $new_gid\n" if $DEBUG;
unless ($new_gid)
  {
    print STDERR "Error: unable to capture a new genome id for the masked genome.  Exiting.\n";
    exit;
  }

#assign genome to user
if ($uid)
  {
    my $node_types = CoGeX::node_types();
    my $uc = $coge->resultset('UserConnector')->find_or_create({
								parent_id=>$uid,
								parent_type => $node_types->{user},
								child_id=>$new_gid,
								child_type=>$node_types->{genome},
								role_id=>2
								});
  }

#copy_permissions (gid1=>$gid, gid2=>$new_gid);
add_annotations(gid1=>$gid, gid2=>$new_gid);

sub copy_permissions 
  {
    my %opts = @_;
    my $gid1 = $opts{gid1};
    my $gid2 = $opts{gid2};
    #permissions will be copied from g1 to g2;
    my $g1 = $coge->resultset('Genome')->find($gid1);
    my $g2 = $coge->resultset('Genome')->find($gid2);
  }

sub add_annotations
  {
    my %opts = @_;
    my $dsgid1 = $opts{gid1};
    my $dsgid2 = $opts{gid2};
    my $cmd = $replicate_annotations;
    $cmd .= " -db " . $DBNAME;
    $cmd .= " -u " . $DBUSER;
    $cmd .= " -pw " . $DBPASS;
    $cmd .= " -dsgid1 " . $dsgid1;
    $cmd .= " -dsgid2 " . $dsgid2;
    $cmd .= " -go " . $GO;
    $cmd .= " -debug ". $DEBUG if $DEBUG;
    print STDERR "Running: $cmd\n";
    open (CMD, "$cmd |");
    while (<CMD>)
      {
	print STDERR $_ if $DEBUG;
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
    $cmd .= " -genome_name '". $dsg->name."'" if $dsg->name;
    $cmd .= " -genome_message 'Copy of genome gid:".$dsg->id."'";
    $cmd .= " -genome_link 'OrganismView.pl?gid=".$dsg->id."'";

    my ($ds) = $dsg->datasets;
    $cmd .= " -source_id ".$ds->data_source->id;
    $cmd .= " -ds_version '". $dsg->version."'";
    $cmd .= " -ds_name '". $ds->name."'";
    $cmd .= " -ds_desc '". $ds->description."'" if $ds->description;
    $cmd .= " -seq_type_id ".$stid;
    $cmd .= " -u ". $DBUSER;
    $cmd .= " -pw ".$DBPASS;
    $cmd .= " -db ".$DBNAME;
    $cmd .= " -restricted ". $restricted if $restricted;
    my $seq_dir = $dsg->file_path;
    $seq_dir =~ s/[^\/]*$//;
    $seq_dir =~ s/\d+\///g;
    $cmd .= " -sd ". $seq_dir;
    $cmd .= " -nt ". $faa;
    print STDERR "Running: ",$cmd,"\n";
    my $dsgid;
    if ($GO)
      {
	open (CMD, "$cmd |");
	while (<CMD>)
	  {
	    if (/genome_id:\s+(\d+)/)
	      {
		$dsgid = $1;
		print STDERR "Captured dsgid: $dsgid\n" if $DEBUG;
	      }
	    print STDERR $_ if $DEBUG;;
	  }
	close CMD;
      }
    return $dsgid;
  }

sub get_and_clean_sequence
  {
    my $faa = shift;
    my $rnd = int(rand(1000000000000));
    my $new_faa = $TEMPDIR."/".$rnd.".clean.faa";
    open (OUT, ">$new_faa");
    open (IN, $faa);
    while (<IN>)
      {
	if (/^>/)
	  {
	    s/lcl\|//g;
	    s/gi\|//g;
	    print OUT $_;
	  }
	else
	  {
	    print OUT $_;
	  }
      }
    close OUT;
    close IN;
    return $new_faa;
  }

sub mask_genome
  {
    my %opts = @_;
    my $dsg = $opts{dsg};
    my $rnd = int(rand(1000000000000));
    my $counts = $TEMPDIR."/".$rnd. ".counts";
    my $masked = $TEMPDIR."/".$rnd. ".masked";
    my $hard = $masked.".hard";
    my $cmd = $windowmasker ." -in ". $dsg->file_path." -mk_counts -out ". $counts;
    print STDERR "running $cmd\n";
    `$cmd` if $GO;
    $cmd = $windowmasker ." -in ". $dsg->file_path." -ustat ". $counts ." -outfmt fasta -dust T -out ".$masked;
    print STDERR "running $cmd\n";
    `$cmd` if $GO;
    $cmd = $hard_mask. " < $masked > $hard";
    print STDERR "running $cmd\n";
    `$cmd` if $GO;
    return $hard;
  }

sub get_sequence_type_id_for_windowmasker
  {
    my $st = $coge->resultset('GenomicSequenceType')->find_or_create({name=>"NCBI WindowMasker (Hard)"});
    return $st->id;
  }

sub usage
  {
    print STDERR qq{
Welcome to $0

Purpose:  take a genome ID from coge, generate a masked version of the genome, load it into CoGe, and map over the annotations

Usage:  $0  -conf_file <coge configuration file> -dsgid <coge database id for genome to be masked and copies> -uid <coge user id> -go 1

Options:
   -go 1                  |      Make the database calls.  Default 0

   -conf_file | cf        |      CoGe conf file

   -gid                   |      CoGe genome id

   -mask                  |      Mask the genome with NCBI's windowmasker

   -restricted            |      mark genome as restricted (Defaults to whatever is set for gid's genome)

   -uid                   |      CoGe user id to whom the new genome will be assigned


Copyright Eric Lyons 2013

};
    exit;
  }
