#!/usr/bin/perl -w


use strict;
use Data::Dumper;
use CoGeX;
use Getopt::Long;
use File::Path;
use URI::Escape::JavaScript qw(escape unescape);
use CoGe::Accessory::Web qw(get_defaults);

use vars qw($staging_dir $install_dir $data_file 
			$name $description $version $restricted 
			$gid $source_name $user_name $config
			$host $port $db $user $pass $P);

my $DATA_TYPE_QUANT = 1;
my $DATA_TYPE_VCF	= 2;

my $MIN_QUANT_COLUMNS = 5;
my $MAX_QUANT_COLUMNS = 6;

my $MIN_VCF_COLUMNS = 8;
#my $MAX_VCF_COLUMNS = 10;

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
	
	# Or use config file
	"config=s"			=> \$config
);

if ($config) {
	my $P = CoGe::Accessory::Web::get_defaults($config);
	$db   = $P->{DBNAME};
	$host = $P->{DBHOST};
	$port = $P->{DBPORT};
	$user = $P->{DBUSER};
	$pass = $P->{DBPASS};	
}

$data_file = unescape($data_file);
$name = unescape($name);
$description = unescape($description);
$version = unescape($version);
$source_name = unescape($source_name);

# Open log file
$| = 1;
my $logfile = "$staging_dir/log.txt";
open(my $log, ">>$logfile") or die "Error opening log file $logfile";
$log->autoflush(1);

# Open system config file
$P = CoGe::Accessory::Web::get_defaults();
my $FASTBIT_LOAD = $P->{FASTBIT_LOAD};
my $FASTBIT_QUERY = $P->{FASTBIT_QUERY};
if (not $FASTBIT_LOAD or not $FASTBIT_QUERY or not -e $FASTBIT_LOAD or not -e $FASTBIT_QUERY) {
	print $log "FASTBIT_LOAD: $FASTBIT_LOAD\n";
	print $log "FASTBIT_QUERY: $FASTBIT_QUERY\n";
	print $log "log: error: can't find fastbit commands\n";
	exit(-1);
}

# Determine file type (.csv quant values or .vcf variant values)
print $log "log: Detecting file type\n";
my $data_type = detect_data_type($data_file);
if (!$data_type) {
	print $log "log: error: unable to detect data type from file\n";
	exit(-1);
}
elsif ($data_type == $DATA_TYPE_QUANT) {
	print $log "log: Detected a CSV quantitative file\n";
}
elsif ($data_type == $DATA_TYPE_VCF) {
	print $log "log: Detected a VCF variant file\n";
}

# Validate the data file
print $log "log: Validating data file\n";
my $count = 0;
my $pChromosomes;
my $data_spec;
if ($data_type == $DATA_TYPE_QUANT) {
	($data_file, $data_spec, $count, $pChromosomes) = validate_quant_data_file($data_file);
}
elsif ($data_type == $DATA_TYPE_VCF) {
	($data_file, $data_spec, $count, $pChromosomes) = validate_vcf_data_file($data_file);
}
if (not $count) {
	print $log "log: error: file contains no data\n";
	exit(-1);
}

my ($filename) = $data_file =~ /^.+\/([^\/]+)$/;
print $log "log: Successfully read " . commify($count) . " lines\n";

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
	print $log "log: error finding genome id$gid\n";
	exit(-1);
}

# Verify that chromosome names in input file match those for genome
my %genome_chr = map { $_ => 1 } $genome->chromosomes;
foreach (sort keys %genome_chr) {
	print $log "genome chromosome $_\n";
}
foreach (sort keys %$pChromosomes) {
	print $log "input chromosome $_\n";
}
my $error = 0;
foreach (sort keys %$pChromosomes) {
	if (not defined $genome_chr{$_}) {
		print $log "log: chromosome '$_' not found in genome\n";
		$error++;
	}
}
if ($error) {
	print $log "log: error: input chromosome names don't match genome\n";
	exit(-1);
}

# Copy input file to staging area and generate fastbit database/index
my $cmd = "cp -f $data_file $staging_dir";
`$cmd`;

my $staged_data_file = $staging_dir . '/' . $filename;

#TODO redirect fastbit output to log file instead of stderr

print $log "log: Generating database\n";
$cmd = "$FASTBIT_LOAD -d $staging_dir -m \"$data_spec\" -t $staged_data_file";
print $log $cmd, "\n";
my $rc = system($cmd);
if ($rc != 0) {
	print $log "log: error executing ardea command: $rc\n";
	exit(-1);
}

print $log "log: Indexing database (may take a few minutes)\n";
$cmd = "$FASTBIT_QUERY -d $staging_dir -v -b \"<binning precision=2/><encoding equality/>\"";
print $log $cmd, "\n";
$rc = system($cmd);
if ($rc != 0) {
	print $log "log: error executing ibis command: $rc\n";
	exit(-1);
}

################################################################################
# If we've made it this far without error then we can feel confident about
# the input data.  Now we can go ahead and create the db entities and 
# install the files.
################################################################################

# Create data source
my $data_source = $coge->resultset('DataSource')->find_or_create( { name => $source_name, description => "" } );#, description => "Loaded into CoGe via LoadExperiment" } );
unless ($data_source) {
	print $log "log: error creating data source\n";
	exit(-1);
}

# Create experiment
my $experiment = $coge->resultset('Experiment')->create(
	{ name				=> $name,
	  description		=> $description,
	  version			=> $version,
	  #link				=> $link, #FIXME 
	  data_source_id	=> $data_source->id,
	  data_type			=> $data_type,
	  genome_id			=> $gid,
	  restricted		=> $restricted
	});
my $storage_path = "$install_dir/".$experiment->get_path;
print $log 'Storage path: ', $storage_path, "\n";
$experiment->storage_path($storage_path);			
$experiment->update;

print $log "experiment id: " . $experiment->id . "\n"; # don't change, gets parsed by calling code

#TODO create experiment type & connector

# Make user owner of new experiment
my $user = $coge->resultset('User')->find( { user_name => $user_name } );
unless ($user) {
	print $log "log: error finding user '$user_name'\n";
	exit(-1);
}
my $node_types = CoGeX::node_types();
my $conn = $coge->resultset('UserConnector')->create(
  { parent_id => $user->id,
	parent_type => $node_types->{user},
	child_id => $experiment->id,
	child_type => $node_types->{experiment},
	role_id => 2 # FIXME hardcoded
  } );
unless ($conn) {
	print $log "log: error creating user connector\n";
	exit(-1);
}

## Copy files from staging directory to installation directory
mkpath($storage_path);
$cmd = "cp -r $staging_dir/* $storage_path";
print $log "$cmd\n";
`$cmd`;

# Yay!
CoGe::Accessory::Web::log_history( db => $coge, user_id => $user->id, page => "LoadExperiment", description => 'load experiment id' . $experiment->id, link => 'ExperimentView.pl?eid=' . $experiment->id );
#print STDERR "experiment id: ".$experiment->id,"\n";
print $log "log: All done!";
close($log);
exit;

#-------------------------------------------------------------------------------
sub detect_data_type {
	my $filepath = shift;

	print $log "detect_data_type: $filepath\n";
	
	my $type;
	
	my ($fileext) = $filepath =~ /\.([^\.]+)$/;
	if ($fileext eq 'csv') {
		return $DATA_TYPE_QUANT;
	}
	elsif ($fileext eq 'vcf') {
		return $DATA_TYPE_VCF;
	}
	else {
		print $log "detect_data_type: unknown file ext '$fileext'\n";
	}
	
	open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
	my $line_num = 1;
	while (my $line = <$in>) {
		last if ($line_num++ > 100); # only check the beginning of the file
		chomp $line;
		next if $line =~ /^\s*#/;
		
		if ($line =~ /fileformat\=VCF/) {
			return $DATA_TYPE_VCF;
		}
		
		my @tok = split(/\s+/, $line);
		if (@tok >= $MIN_QUANT_COLUMNS && @tok <= $MAX_QUANT_COLUMNS) {
			return $DATA_TYPE_QUANT;
		}
		elsif (@tok >= $MIN_VCF_COLUMNS) {
			return $DATA_TYPE_VCF;
		}
	}
	close($in);
	
	return;
}

sub validate_quant_data_file {
	my $filepath = shift;
	
	my %chromosomes;
	my $line_num = 1;
	
	print $log "validate_quant_data_file: $filepath\n";
	open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
	my $outfile = $filepath.".csv";
	open (my $out, ">$outfile");
	while (my $line = <$in>) {
		$line_num++;
		next if ($line =~ /^\s*#/);
		chomp $line;
		my @tok = split(/,/, $line);
		
		# If comma-delimited parsing didn't work try white-space delimited
		unless (@tok >= $MIN_QUANT_COLUMNS && @tok <= $MAX_QUANT_COLUMNS) {
			@tok = split (/\s+/,$line);
		}
		# Validate format
		if (@tok < $MIN_QUANT_COLUMNS) {
			print $log "log: error at line $line_num: more columns expected (" . @tok . " < $MIN_QUANT_COLUMNS)\n";
			return;
		}
		elsif (@tok > $MAX_QUANT_COLUMNS) {
			print $log "log: error at line $line_num: fewer columns expected (" . @tok . " > $MAX_QUANT_COLUMNS)\n";
			return;
		}
		
		# Validate values
		my ($chr, $start, $stop, $strand, $val1, $val2) = @tok;
		if (not defined $chr or not defined $start or not defined $stop or not defined $strand) {
			print $log "log: error at line $line_num: missing value in a column\n";
			return;
		}
		if (not defined $val1 or $val1 < 0 or $val1 > 1) {
			print $log "log: error at line $line_num: value 1 not between 0 and 1\n";
			return;
		}
		if (not defined $val2) {
			$val2 = 0;
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
		$strand = $strand =~ /-/ ? -1 : 1;
		print $out join (",", $chr, $start, $stop, $strand, $val1, $val2),"\n";
		$chromosomes{$chr}++;
	}
	close($in);
	close($out);
	my $format = "chr:key, start:unsigned long, stop:unsigned long, strand:byte, value1:double, value2:double";
	return ($outfile, $format, $line_num, \%chromosomes);
}

# For VCF format specification v4.1, see http://www.1000genomes.org/node/101
sub validate_vcf_data_file {
	my $filepath = shift;
	
	my %chromosomes;
	my $line_num = 1;
	
	print $log "validate_vcf_data_file: $filepath\n";
	open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
	my $outfile = $filepath.".vcf";
	open (my $out, ">$outfile");
	while (my $line = <$in>) {
		$line_num++;
		#TODO load VCF metadata for storage as experiment annotations in DB (lines that begin with '##')
		next if ($line =~ /^#/);
		chomp $line;
		next unless $line;
		my @tok = split(/\s+/, $line);
		
		# Validate format
		if (@tok < $MIN_VCF_COLUMNS) {
			print $log "log: error at line $line_num: more columns expected (" . @tok . " < $MIN_VCF_COLUMNS)\n";
			return;
		}
		
		# Validate values and set defaults
		my ($chr, $pos, $id, $ref, $alt, $qual, undef, $info) = @tok;
		if (not defined $chr ||
			not defined $pos ||
			not defined $ref || 
			not defined $alt) 
		{
			print $log "log: error at line $line_num: missing required value in a column\n";
			return;
		}
		next if ($alt eq '.'); # skip monomorphic sites
		$id = '' if (not defined $id or $id eq '.');
		$qual = 0 if (not defined $qual);
		$info = '' if (not defined $info);

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
		
		# Each line could encode multiple alleles
		my @alleles = split(',', $alt);
		foreach my $a (@alleles) {
			# Determine site type
			my $type = detect_type($ref, $a);
			# Save to file
			print $out join (",", $chr, $pos, $pos+length($ref)-1, $type, $id, $ref, $a, $qual, $info),"\n";
		}
	}
	close($in);
	close($out);
	my $format = "chr:key, start:unsigned long, stop:unsigned long, type:key, id:text, ref:key, alt:key, qual:double, info:text";
	return ($outfile, $format, $line_num, \%chromosomes);
}

sub detect_type {
	my $ref = shift;
	my $alt = shift;
	
	return 'snp' if (length $ref == 1 and length($alt) == 1);
	return 'deletion' if (length $ref > length $alt);
	return 'insertion' if (length $ref < length $alt);
	return 'unknown';
}

sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}
