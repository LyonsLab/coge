#!/usr/bin/perl -w

use strict;
use File::Basename qw( basename );

my $data_file = shift;

if (!validate_fastq_data_file($data_file)) {
    print STDOUT "log: error: file contains no data\n";
    exit(-1);
}
#print STDOUT "log: file has been verified successfully\n";

system("touch $data_file.validated"); # signal completion to JEX
exit;

#-------------------------------------------------------------------------------
sub validate_fastq_data_file {
    my $filepath = shift;
    my $filename = basename($filepath);
    my $count = 0;
    
    # Ensure size
    unless ( -s $filepath ) {
        print STDOUT "log: error: '$filename' is zero-length\n";
        exit(-1);
    }
    
    # Ensure text file
    if ( -B $filepath ) {
        print STDOUT "log: error: '$filename' is a binary file\n";
        exit(-1);
    }
    
    # Validate records
    my $in;
    unless (open($in, $filepath)) {
        print STDOUT "Error: can't open $filepath for reading: $!";
    };

    my $line;
    while( $line = <$in> ) {
        if ($line =~ /^\@/) { # header
            #my $seqID = $line;   chomp $seqID;
            my $seq   = <$in>;   chomp $seq;
            my $line3 = <$in>;   #chomp $line3;
            my $qual  = <$in>;   chomp $qual;

            $count++;

# mdb removed 1/25/16 COGE-700
#            if (length $seq != length $qual) {
#                print STDOUT "log: error: invalid record \#$count: sequence line length differs from quality line length\n",
#                             "$line",
#                             "$seq\n",
#                             "$line3",
#                             "$qual\n";
#                exit(-1);
#            }
        }
    }
    close($in);

    return $count;
}
