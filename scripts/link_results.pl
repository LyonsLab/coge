#!/usr/bin/perl -w

use v5.10;
use strict;
use warnings;

use CoGeX;
use Getopt::Long;
use File::Path;
use File::Basename qw( basename );
use File::Spec::Functions qw( catdir catfile );
use URI::Escape::JavaScript qw(unescape);

our ($result_dir, $input_files);

GetOptions(
    "result_dir=s"  => \$result_dir,     # results path
    "input_files=s" => \$input_files,    # comma-separated list (JS escaped) of files
);

mkpath($result_dir);
$input_files = unescape($input_files) if $input_files;
my @files = split( ',', $input_files );

# We require at least one result file
exit 1 unless @files;

foreach my $file (@files) {
    my $filename = basename($file);

    # Exit failure if symlink fails
    exit 1 unless symlink $file, catfile($result_dir, $filename);
}
