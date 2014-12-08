#!/usr/bin/perl -w

use strict;

my $data_file = shift;

if (!validate_fastq_data_file($data_file)) {
    print STDOUT "log: error: file contains no data\n";
    exit(-1);
}

system("touch $data_file.validated");
exit;

#-------------------------------------------------------------------------------
sub validate_fastq_data_file {
    my $filepath = shift;
    my $count = 0;

    # Validate file
    open( my $in, $filepath ) || die "can't open $filepath for reading: $!";
    my $line;
    while( $line = <$in> ) {
        if ($line =~ /^\@/) { # header
            #my $seqID = $line;   chomp $seqID;
            my $seq   = <$in>;   chomp $seq;
            my $line3 = <$in>;   #chomp $line3;
            my $qual  = <$in>;   chomp $qual;

            $count++;

            if (length $seq != length $qual) {
                print STDOUT "log: error: invalid record (length seq != length qual) line=$count\n";
                exit(-1);
            }
        }
    }
    close($in);

    return $count;
}
