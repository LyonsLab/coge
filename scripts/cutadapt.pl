#!/usr/bin/env perl
#-------------------------------------------------------------------------------
# Purpose:  Wrapper script for Cutadapt for use with JEX.  Cutadapt's command-line
#     usage with paired-end reads will not work with JEX's workflow interface.
# Author:   mdb
# Created:  5/12/15
#------------------------------------------------------------------------------

use strict;
use warnings;

use File::Path qw(make_path);
use File::Copy qw(move);
use File::Touch;
use File::Spec::Functions qw(catfile);
use CoGe::Accessory::Utils qw(execute is_fastq_file add_fastq_ext to_filename detect_paired_end);

my $read_type = shift;   # 'single' or 'paired'
my $output_path = shift; # path for output files
my $cmd_args = shift;    # cutadapt executable and all command line arguments except for file lists
my @files = @ARGV;       # fastq files

make_path($output_path); # prevent cutadapt from erroring on missing output path

# Rename files if missing FASTQ extension required by cutadapt - mdb added 11/6/15 COGE-673
my %verified_files;
foreach my $orig_file (@files) {
    my $new_file = $orig_file;
    unless (is_fastq_file($orig_file)) {
        $new_file = add_fastq_ext($orig_file);
        move($orig_file, $new_file);
        print STDERR "Renaming $orig_file to $new_file\n";
        die "error: couldn't rename file '$orig_file'" unless (-r $new_file);
    }
    $verified_files{$new_file} = $orig_file;
}

# Split files into pair groups if paired-end
if ($read_type eq 'paired') {
    my ($m1, $m2) = detect_paired_end([keys %verified_files]);
    die "error: invalid paired-end files" unless (@$m1 and @$m2);
    $cmd_args .= ' ' . join(' ', sort @$m1) . ' ' . join(' ', sort @$m2);
}
else {
    $cmd_args .= ' ' . join(' ', keys %verified_files);
}

# Execute cutadapt
print STDOUT "cutadapt.pl: $cmd_args\n";
my $rc = execute($cmd_args);

# Restore original filenames - mdb added 11/6/15 COGE-673
foreach my $new_file (keys %verified_files) {
    my $orig_file = $verified_files{$new_file};
    if ($orig_file ne $new_file) {
        move($new_file, $orig_file);
    }
}

exit($rc);
