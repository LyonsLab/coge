#!/usr/bin/perl -w

use DBI;
use strict;
use CoGeX;
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;
use URI::Escape::JavaScript qw(escape unescape);

use vars qw($staging_dir $install_dir $fasta_files 
			$name $description $version $type_id $restricted 
			$org_name $source_name $user_name
			$host $port $db $user $pass);

GetOptions(
	"staging_dir=s"	=> \$staging_dir,
	"install_dir=s" => \$install_dir,
	"fasta_files=s" => \$fasta_files,	# comma-separated list (JS escaped)
	"name=s"		=> \$name,			# genome name (JS escaped)
	"desc=s"		=> \$description,	# genome description (JS escaped)
	"version=s"		=> \$version,		# genome version (JS escaped)
	"type_id=i"		=> \$type_id,		# genomic_sequence_type_id
	"restricted=i"	=> \$restricted,	# genome restricted flag
	"org_name=s"	=> \$org_name,		# organism name (JS escaped)
	"source_name=s"	=> \$source_name,	# data source name
	"user_name=s"	=> \$user_name,		# user name
	
	# Database params
	"host|h=s"			=> \$host,
	"port|p=s"			=> \$port,
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw=s"		=> \$pass,
);

$fasta_files = unescape($fasta_files);
$name = unescape($name);
$description = unescape($description);
$version = unescape($version);
$org_name = unescape($org_name);

# Open log file
$| = 1;
my $logfile = "$staging_dir/log.txt";
open(my $log, ">>$logfile") or die "Error opening log file";
$log->autoflush(1);

# Process each file into staging area
my @files = split(',', $fasta_files);
my %sequences;
foreach my $file (@files) {
	process_fasta_file(\%sequences, $file, $staging_dir);
}

# If we've made it this far without error then we can feel confident about our
# ability to parse all of the input files.  Now we can go ahead and 
# create the db entities and install the files.

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge    = CoGeX->connect( $connstr, $user, $pass );
unless ($coge) {
	print $log "Couldn't connect to database\n";
	die;	
}

# Retrieve organism
my $organism = $coge->resultset('Organism')->find( { name => $org_name } );
die "Error finding organism" unless ($organism);
print STDERR "organism id: " . $organism->id . "\n";
	
# Create datasource
my $datasource = $coge->resultset('DataSource')->find_or_create( { name => $source_name, description => "Loaded into CoGe via LoadGenome" } );
die "Error creating data source" unless $datasource;
print STDERR "datasource id: " . $datasource->id . "\n";
	
# Create genome
my $genome = $coge->resultset('Genome')->create( 
  { name => $name, 
  	description => $description, 
  	version => $version, 
  	organism_id => $organism->id,
  	genomic_sequence_type_id => $type_id,
  	restricted => $restricted
  } );
unless ($genome) {
	print $log "Error creating genome\n";
	die;	
}
print STDERR "genome id: " . $genome->id . "\n";
print $log "genome id: " . $genome->id . "\n"; # don't change, gets parsed by calling code

$install_dir = "$install_dir/" . $genome->get_path . "/";# . $genome->id . ".faa";
$genome->file_path($install_dir);
$genome->update;

# Add new genome to user's owner list
my $user = $coge->resultset('User')->find( { user_name => $user_name } );
unless ($user) {
	print $log "Error finding user '$user_name'\n";
	die;
}
my $child_types = CoGeX::list_child_types();
my $listconn = $coge->resultset('ListConnector')->create(
	{ parent_id => $user->owner_list->id,
	  child_id => $genome->id,
	  child_type => $child_types->{genome}
	} );
unless ($listconn) {
	print $log "Error creating list connector\n";
	die;
}

# Create datasets
my %datasets;
foreach my $file (@files) {
	my ($filename) = $file =~ /^.+\/([^\/]+)$/;
	my $dataset = $coge->resultset('Dataset')->create( 
	  { data_source_id => $datasource->id,
	  	name => $filename,
	  	description => $description,
		version => $version,
		restricted => $restricted,
	  } );
	unless ($dataset) {
		print $log "Error creating dataset\n";
		die;
	}
	#TODO set link field if from FTP
	print STDERR "dataset id: " . $dataset->id . "\n";
	$datasets{$file} = $dataset->id;
		
	my $dsconn = $coge->resultset('DatasetConnector')->create( { dataset_id => $dataset->id, genome_id => $genome->id } );
	unless ($dsconn) {
		print $log "Error creating dataset connector\n";
		die;
	}
}

# Create genomic_sequence/feature/location for ea. chromosome
foreach my $chr (sort keys %sequences) {
	my $seqlen = $sequences{$chr}{size};
	my $dsid = $datasets{ $sequences{$chr}{file} };
	
	$genome->add_to_genomic_sequences( { sequence_length => $seqlen, chromosome => $chr } );

	#must add a feature of type chromosome to the dataset so the dataset "knows" its chromosomes
	my $feat_type = $coge->resultset('FeatureType')->find_or_create( { name => 'chromosome' } );
	my $feat = $coge->resultset('Feature')->find_or_create(
	  { dataset_id => $dsid, 
	  	feature_type_id => $feat_type->id, 
	  	start => 1, 
	  	stop => $seqlen, 
	  	chromosome => $chr, 
	  	strand => 1 
	  } );
	my $feat_name = $coge->resultset('FeatureName')->find_or_create( { name => "chromosome $chr", feature_id => $feat->id } );	

	my $loc = $coge->resultset('Location')->find_or_create( 
	  { feature_id => $feat->id, 
	  	start => 1, 
	  	stop => $seqlen, 
	  	strand => 1, 
	  	chromosome => $chr 
	  } );
}

# Copy files from staging directory to installation directory
print STDERR "install_dir=$install_dir\n";
mkpath($install_dir);
mkpath( $install_dir . "/chr" );
my $cmd = "cp -r $staging_dir/chr $install_dir";
print STDERR "$cmd\n";
`$cmd`;
my $genome_filename = $genome->id . ".faa";
$cmd = "cp $staging_dir/genome.faa $install_dir/$genome_filename";
print STDERR "$cmd\n";
`$cmd`;

# Yay!
print $log "All done!";
print STDERR "All done!";
close($log);
exit;

#-------------------------------------------------------------------------------
sub process_fasta_file {
	my $pSeq = shift;
	my $filepath = shift;
	my $target_dir = shift;
	
	print STDERR "process_fasta_file: $filepath\n";
	$/ = "\n>";
	open( my $in, $filepath ) || die "can't open $filepath for reading: $!";

	while (<$in>) {
		s/>//g;
		my ( $name, $seq ) = split /\n/, $_, 2;
		$seq =~ s/\n//g;

		my ($chr) = split(/\s+/, $name);
		$chr =~ s/^lcl\|//;
		$chr =~ s/chromosome//i;
		$chr =~ s/^chr//i;
		$chr =~ s/^0+//;
		$chr =~ s/^_+//;
		$chr =~ s/\s+/ /;
		$chr =~ s/^\s//;
		$chr =~ s/\s$//;
		$chr = 0 unless $chr;

		#TODO add more checks on chr name and sequence here
		die 'Error parsing section header' if (not $chr);
		die "Duplicate section name '$chr'" if (defined $pSeq->{$chr});
		
		my ($filename) = $filepath =~ /^.+\/([^\/]+)$/;
		print $log "log: Processed '$chr' in $filename\n";

		# Append sequence to master file
		open( my $out, ">>$target_dir/genome.faa" );
		my $head = $chr =~ /^\d+$/ ? ">gi" : ">lcl";
		$head .= "|" . $chr;
		print $out "$head\n$seq\n";
		close($out);
	
		# Create individual file for chromosome
		mkpath("$target_dir/chr");
		open( $out, ">$target_dir/chr/$chr" );
		print $out $seq;
		close($out);
		
		$pSeq->{$chr} = { size => length $seq, file => $filepath };
	}
	close($in);
}
