#!/usr/bin/perl -w

use strict;

use CoGeX;
use CoGe::Core::Metadata qw(export_annotations);

use File::Basename;
use File::Path;
use File::Spec;
use Getopt::Long;
use Data::Dumper;

our ( $id, $type, $config, $workdir, $output_filename, $resdir, $files);

GetOptions(
    "id=i"             => \$id,       # exeriment/genome id
    "type=s"           => \$type,     # 'experiment' or 'genome'
    "output=s"         => \$output_filename, # output tarball filename
    "directory|dir=s"  => \$workdir,  # output directory,
    "files=s"          => \$files,    # optional additional input files
    "config=s"         => \$config    # CoGe configuration file
);

$| = 1; # enable autoflushing

# Open log file
mkpath($workdir, 0, 0755) unless -r $workdir;
my $logfile = File::Spec->catdir($workdir, "$output_filename.log");
open (my $logh, ">", $logfile) or die "Error opening log file";

# Verify required options
unless ($id) {
    say $logh "log: error: experiment/genome not specified, use 'id' option";
    exit(-1);
}
$type = lc($type);
unless ($type && ($type eq 'experiment' || $type eq 'genome')) {
    say $logh "log: error: type not specified, use 'type' option";
    exit(-1);
}
unless ($output_filename) {
    say $logh "log: error: output file not specified, use 'output' option";
    exit(-1);
}
unless ($config) {
    say $logh "log: error: config file not specified, use 'config' option";
    exit(-1);
}

# Open config file
my $P = CoGe::Accessory::Web::get_defaults($config);

$resdir = $P->{RESOURCEDIR};
unless ($resdir) {
    say $logh "missing RESOURCEDIR setting in config file";
    exit(-1);
}

# Connect to DB
my $connstr = "dbi:mysql:dbname=".$P->{DBNAME}.";host=".$P->{DBHOST}.";port=".$P->{DBPORT}.";";
my $coge = CoGeX->connect( $connstr, $P->{DBUSER}, $P->{DBPASS} );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);
unless ($coge) {
    say $logh "log: error: couldn't connect to database";
    exit(-1);
}

# Generate output tarball
my $archive = File::Spec->catdir($workdir, $output_filename);

my $sourceObj;
if ($type eq 'experiment') {
    $sourceObj = $coge->resultset('Experiment')->find($id);
}
elsif ($type eq 'genome') {
    $sourceObj = $coge->resultset('Genome')->find($id);
}

my @annotations = $sourceObj->annotations;
@annotations = () unless @annotations;

unless (-r $archive and -r "$archive.finished") {
    my @file_list = export_annotations( annotations => \@annotations, export_path => $workdir );
    push @file_list, split(',', $files);
    copy_readme();
    my $info = export_info( $sourceObj->info_file );

    my $input_dir = $sourceObj->storage_path;
    my $cmd = "tar -czf $archive --exclude=log.txt --directory $input_dir .";
    $cmd .= " --transform 's,^,/".$type."_".$id."/,'";
    $cmd .= " --directory $workdir $info";
    $cmd .= " " . join " ", @file_list if @file_list;
    $cmd .= " README.txt";

    system($cmd);
    system("touch $archive.finished");
}

exit;
#-------------------------------------------------------------------------------

sub copy_readme {
    my $readme = File::Spec->catdir($workdir, "README.txt");
    system("cp $resdir/$type.README.txt $readme");
}

sub export_info {
    my $info = shift;
    
    my $info_file = File::Spec->catdir($workdir, "info.csv");
    unless (-r $info_file) {
        open(my $fh, ">", $info_file);
        say $fh $info;
        close($fh);
    }

    return basename($info_file);
}
