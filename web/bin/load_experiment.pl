#!/usr/bin/perl -w

use DBI;
use strict;
use CoGeX;
use Roman;
use Data::Dumper;
use Getopt::Long;
use File::Path;
use URI::Escape::JavaScript qw(escape unescape);
use CoGe::Accessory::Web qw(get_defaults);

use vars qw($staging_dir $install_dir $data_file 
			$name $description $version $restricted 
			$gid $source_name $user_name
			$host $port $db $user $pass $P);

my $MIN_COLUMNS = 5;
my $MAX_COLUMNS = 6;

GetOptions(
	"staging_dir=s"	=> \$staging_dir,
	"install_dir=s" => \$install_dir,
	"data_file=s" 	=> \$data_file,		# data file (JS escape)
	"name=s"		=> \$name,			# experiment name (JS escaped)
	"desc=s"		=> \$description,	# experiment description (JS escaped)
	"version=s"		=> \$version,		# experiment version (JS escaped)
	"restricted=i"	=> \$restricted,	# experiment restricted flag
	"gid=s"			=> \$gid,			# genome id
	"source_name=s"	=> \$source_name,	# data source name (JS escaped)
	"user_name=s"	=> \$user_name,		# user name
	
	# Database params
	"host|h=s"			=> \$host,
	"port|p=s"			=> \$port,
	"database|db=s"		=> \$db,
	"user|u=s"			=> \$user,
	"password|pw=s"		=> \$pass,
);

$data_file = unescape($data_file);
$name = unescape($name);
$description = unescape($description);
$version = unescape($version);
$source_name = unescape($source_name);

# Open log file
$| = 1;
my $logfile = "$staging_dir/log.txt";
open(my $log, ">>$logfile") or die "Error opening log file";
$log->autoflush(1);

# Open system config file
$P = CoGe::Accessory::Web::get_defaults();
my $FASTBIT_LOAD = $P->{FASTBIT_LOAD};
my $FASTBIT_QUERY = $P->{FASTBIT_QUERY};
if (not $FASTBIT_LOAD or not $FASTBIT_QUERY or not -e $FASTBIT_LOAD or not -e $FASTBIT_QUERY) {
	print STDERR "FASTBIT_LOAD: $FASTBIT_LOAD\n";
	print STDERR "FASTBIT_QUERY: $FASTBIT_QUERY\n";
	print $log "log: error: can't find fastbit commands\n";
	exit(-1);
}

# Validate the data file
my ($filename) = $data_file =~ /^.+\/([^\/]+)$/;
print $log "log: Validating $filename\n";
my ($count, $pChromosomes) = validate_data_file($data_file);
if (not $count) {
	exit(-1);	
}
print $log "log: Successfully read $count lines\n";

# Copy to staging area and generate fastbit database/index
my $cmd = "cp -f $data_file $staging_dir";
`$cmd`;

my $staged_data_file = $staging_dir . '/' . $filename;

print $log "log: Generating database\n";
$cmd = "$FASTBIT_LOAD -d $staging_dir -m \"chr:key, start:unsigned long, stop:unsigned long, strand:byte, value1:double, value2:double\" -t $staged_data_file";
print STDERR $cmd, "\n";
my $rc = system($cmd);
if ($rc != 0) {
	print $log "log: error executing ardea command: $rc\n";
	exit(-1);
}

print $log "log: Indexing database (may take a few minutes)\n";
$cmd = "$FASTBIT_QUERY -d $staging_dir -v -b \"<binning precision=2/><encoding equality/>\"";
print STDERR $cmd, "\n";
$rc = system($cmd);
if ($rc != 0) {
	print $log "log: error executing ibis command: $rc\n";
	exit(-1);
}

# If we've made it this far without error then we can feel confident about
# the input data.  Now we can go ahead and create the db entities and install the files.

# Connect to database
my $connstr = "dbi:mysql:dbname=$db;host=$host;port=$port;";
my $coge    = CoGeX->connect( $connstr, $user, $pass );
unless ($coge) {
	print $log "log: couldn't connect to database\n";
	exit(-1);
}

# Retrieve genome
my $genome = $coge->resultset('Genome')->find( { genome_id => $gid } );
unless ($genome) {
	print $log "log: error finding genome\n";
	exit(-1);
}

# Verify that chromosome names in input file match those for genome
my %genome_chr = map { $_ => 1 } $genome->chromosomes;
foreach (keys %genome_chr) {
	print $log "genome chromosome $_\n";
}
foreach (keys %$pChromosomes) {
	print $log "input chromosome $_\n";
}
foreach (keys %$pChromosomes) {
	if (not defined $genome_chr{$_}) {
		print $log "log: error: chromosome '$_' not found in genome\n";
		exit(-1);
	}
}

# Create datasource
my $datasource = $coge->resultset('DataSource')->find_or_create( { name => $source_name, description => "Loaded into CoGe via LoadExperiment" } );
die "Error creating data source" unless $datasource;

# Create experiment
my $experiment = $coge->resultset('Experiment')->create(
	{ name				=> $name,
	  description		=> $description,
	  version			=> $version,
	  #link				=> $link, #FIXME 
	  data_source_id	=> $datasource->id,
	  genome_id			=> $gid,
	  restricted		=> $restricted
	});
my $storage_path = "$install_dir/".$experiment->get_path;
print STDERR 'Storage path: ', $storage_path, "\n";
$experiment->storage_path($storage_path);			
$experiment->update;

print STDERR "experiment id: " . $experiment->id . "\n";
print $log "experiment id: " . $experiment->id . "\n"; # don't change, gets parsed by calling code

#TODO create experiment type & connector

# Add new experiment to user's owner list
my $user = $coge->resultset('User')->find( { user_name => $user_name } );
unless ($user) {
	print $log "Error finding user '$user_name'\n";
	exit(-1);
}
my $child_types = CoGeX::list_child_types();
my $listconn = $coge->resultset('ListConnector')->create(
	{ parent_id => $user->owner_list->id,
	  child_id => $experiment->id,
	  child_type => $child_types->{experiment}
	} );
unless ($listconn) {
	print $log "Error creating list connector\n";
	exit(-1);
}

## Copy files from staging directory to installation directory
mkpath($install_dir);
$cmd = "cp -r $staging_dir $storage_path";
print STDERR "$cmd\n";
`$cmd`;

# Yay!
print $log "log: All done!";
print STDERR "All done!";
close($log);
exit;

#-------------------------------------------------------------------------------
sub validate_data_file {
	my $filepath = shift;
	
	my %chromosomes;
	my $line_num = 1;
	
	print STDERR "validate_data_file: $filepath\n";
	open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
	while (my $line = <$in>) {
		$line_num++;
		next if ($line =~ /^\s*#/);
		chomp $line;
		my @tok = split(',', $line);
		
		# Validate format
		if (@tok < $MIN_COLUMNS) {
			print $log "log: error at line $line_num: more columns expected\n";
			return;
		}
		elsif (@tok > $MAX_COLUMNS) {
			print $log "log: error at line $line_num: fewer columns expected\n";
			return;
		}
		
		# Validate values
		my ($chr, $start, $stop, $strand, $val1, $val2) = @tok;
		if ($val1 < 0 or $val1 > 1) {
			print $log "log: error at line $line_num: value 1 not between 0 and 1\n";
			return;
		}
		
		# Fix chromosome identifier
		$chr =~ s/^lcl\|//;
		$chr =~ s/chromosome//i;
		$chr =~ s/^chr//i;
		$chr =~ s/^0+//;
		$chr =~ s/^_+//;
		$chr =~ s/\s+/ /;
		$chr =~ s/^\s//;
		$chr =~ s/\s$//;
		
		if (not $chr) {
			print $log "log: error at line $line_num: trouble parsing chromosome\n";
			return;			
		}
		
		$chromosomes{$chr}++;
	}
	close($in);
	
	return ($line_num, \%chromosomes);
}
