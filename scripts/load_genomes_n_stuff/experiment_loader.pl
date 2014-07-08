#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;
use Getopt::Long;
use File::Path;
use DBI;
#use DBIxProfiler;
use CoGeX;

use vars qw($DEBUG $coge $GO $ERASE);

my $cmdPath = '/usr/local/bin';

my ($exp_file, $exp_dir, $exp_name, $exp_version, $exp_desc, $exp_storage_path, $exp_link,
		$exp_type_name, $exp_type_desc, $exp_type_id,
		$source_name, $source_desc, $source_link, $source_id,
		$dsg_id, $restricted, $db, $user, $pass);

GetOptions (
			"debug=s" => \$DEBUG,
	    "go=s"    => \$GO,
	    "erase|e" => \$ERASE,

	    # Input file
	    "exp_file=s" => \$exp_file, 			# input data file
	    "exp_dir=s" => \$exp_dir,					# target directory in CoGe installation

	    # 'experiment' table values
	    "exp_name=s" => \$exp_name,
	    "exp_version=s" => \$exp_version,
	    "exp_desc=s" => \$exp_desc,
	    "exp_link=s" => \$exp_link,

	    # 'experiment_type' table values - UNUSED
	    "exp_type_name=s" => \$exp_type_name,
	    "exp_type_desc=s" => \$exp_type_desc,
	    "exp_type_id=i" => \$exp_type_id,

	    "source_name=s" => \$source_name, # datasource
	    "source_desc=s" => \$source_desc,
	    "source_link=s" => \$source_link,
	    "source_id=s"   => \$source_id,

	    "dsg_id=s"=>\$dsg_id,
	    "restricted=i"=>\$restricted,
	    "database|db=s" => \$db,
	    "user|u=s" => \$user,
			"password|pw=s" => \$pass,
	   );

$DEBUG = 1 unless defined $DEBUG; # set to '1' to get updates on what's going on
$GO = 1 unless defined $GO; #set to 1 to actually make db calls.

$restricted = 0 unless defined $restricted;
print STDERR "Running $0\n";

print help("Missing experiment name/file/dir\n") unless ($exp_name and $exp_file and $exp_dir);
$exp_version = 1 unless defined $exp_version;

# Connect to the database
my $connstr = "dbi:mysql:dbname=$db;host=localhost;port=3307";
$coge = CoGeX->connect($connstr, $user, $pass);
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my ($dsg) = $coge->resultset("DatasetGroup")->find($dsg_id);
print help("Couldn't find dataset group") unless ($dsg);

$exp_type_name = "unclassified" unless ($exp_type_name);
unless ($exp_type_id)
	{
    my $et = $coge->resultset("ExperimentType")->find_or_create({name=>$exp_type_name,description=>$exp_type_desc});
    $exp_type_id = $et->id;
	}

if ($source_name)
  {
    my $source = $coge->resultset("DataSource")->find_or_create({name=>$source_name,description=>$source_desc,link=>$source_link});
    $source_id = $source->id;
  }
print help("Couldn't find data source\n") unless ($source_id);

my $exp = $coge->resultset('Experiment')->create({
					   	name              => $exp_name,
					   	description       => $exp_desc,
					   	version						=> $exp_version,
					   	link              => $exp_link,
					   	data_source_id    => $source_id,
					   	dataset_group_id	=> $dsg_id,
							restricted        => $restricted
					  }) if $GO;
$exp_storage_path = "$exp_dir/".$exp->get_path;
print 'Storage path: ', $exp_storage_path, "\n";
$exp->storage_path($exp_storage_path);
$exp->update;

my $etc = $coge->resultset('ExperimentTypeConnector')->create({
							experiment_type_id 	=> $exp_type_id,
							experiment_id 			=> $exp->id
						}) if $GO;

# Generate FastBit database and indexes
die "Database directory already exists, please rename or remove it. $!\n" if (-d $exp_storage_path);

mkpath($exp_storage_path);
my $cmd = "$cmdPath/ardea -d $exp_storage_path -m \"chr:key, start:unsigned long, stop:unsigned long, strand:byte, value1:double, value2:double\" -t $exp_file";
print "Executing: $cmd\n";
my $rc = system($cmd);
die "Error executing command $rc" if ($rc != 0);

$cmd = "$cmdPath/ibis -d $exp_storage_path -v -b \"<binning none/><encoding equality/>\"";
print "Executing: $cmd\n";
$rc = system($cmd);
die "Error executing command $rc" if ($rc != 0);

print "all done!\n";

exit;

#-------------------------------------------------------------------------------

sub help # FIXME
{
	my $msg = shift;

	print $msg if $msg;

	print qq
	{
    Options (required):
        database|db   database name
        user|u        database username
        password|pw   database password
        exp_file      raw data file in .csv format
        exp_dir       base directory where experiment databases are stored
        exp_name      name of experiment
        dsg_id        data set group id (genome)
        source_name   data source name

};
	exit;
}

__END__
