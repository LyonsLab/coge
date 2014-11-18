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

our ($result_dir, $input_files, $output_files);

GetOptions(
    "result_dir=s"   => \$result_dir,     # results path
    "input_files=s"  => \$input_files,    # comma-separated list (JS escaped) of files
    "output_files=s" => \$output_files,   # optional comma-separated list (JS escaped) of link names
);

mkpath($result_dir);
$input_files  = unescape($input_files)  if $input_files;
$output_files = unescape($output_files) if $output_files;
my @infiles  = split(',', $input_files);
my @outfiles = split(',', $output_files);

# We require at least one result file
exit 1 unless @infiles;

for (my $i = 0;  $i < scalar(@infiles);  $i++) {
    my $inf  = $infiles[$i];
    my $outf = $inf;
    $outf = $outfiles[$i] if (@outfiles && $outfiles[$i]);

    # Exit failure if symlink fails
    exit 1 unless symlink( $inf, catfile($result_dir, basename($outf)) );
}
