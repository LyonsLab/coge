#!/usr/bin/perl -w

use v5.10;
use strict;
use warnings;

use CoGeX;
use Getopt::Long;
use File::Path;
use File::Basename qw( basename dirname );
use File::Spec::Functions qw( catdir catfile );
use URI::Escape::JavaScript qw(unescape);

our ($result_dir, $input_files, $dest_type, @DEST_TYPES);

@DEST_TYPES = qw(irods http);

GetOptions(
    "result_dir=s"  => \$result_dir,     # results path
    "input_files=s" => \$input_files,    # comma-separated list (JS escaped) of files
    "type=s"        => \$dest_type
);

$input_files = unescape($input_files) if ($input_files);

my @files = split( ',', $input_files );
my @result;

# Check if the dest_type is supported
exit 1 unless grep { $_ eq $dest_type } @DEST_TYPES;

# We require at least one result file
exit 1 unless @files;

foreach my $file (@files) {
    my $name = basename($file);
    say "add file=$name to the results";

    push @result, {
        type => $dest_type,
        path => dirname($file),
        name => $name
    }
}

mkpath($result_dir);
my $output = catfile($result_dir, '1');
CoGe::Accessory::TDS::write($output, @result);
