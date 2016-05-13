#!/usr/bin/env perl
use v5.14;
use strict;
use warnings;

use Data::Dumper qw(Dumper);
use File::Path qw(mkpath);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Accessory::Utils qw(to_pathname);
use CoGe::Core::Storage qw(create_image);
use CoGeX;

my ($image_file, $config_file, $log_file) = @ARGV;

$| = 1;
print STDOUT "Starting $0 (pid $$)\n", qx/ps -o args $$/;

# Check required parameters
die "ERROR: missing argument" unless ($image_file && $config_file);
$log_file = "$image_file.log" unless $log_file;

# Connect to DB
my $conf = get_defaults($config_file);
my $db = CoGeX->dbconnect($conf);
die "ERROR: couldn't connect to the database" unless $db;

# Create the image
my $image = create_image(filename => $image_file, db => $db);
die "Unabled to create image" unless $image;

# Save image ID in log file -- signals task completion to JEX
my $log_path = to_pathname($log_file);
mkpath($log_path);
open(my $fh, ">>", $log_file);
say $fh "image id: " . $image->id;
close($fh);

exit;
