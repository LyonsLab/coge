use strict;

my $QUERY_FASTA      = "arabidopsis.fasta";
my $SUBJECT_DATABASE = "somedb.fasta";
my $NCPUS            = 8;
my $OUTFILE          = "report.%i.blast"
  ;    # becomes report.1.blast, report.2.blast, report.3.blast, etc.
my $BLAST_CMD =
  'blastall -p blastn -i %s -d %s -e 1 -G 2 -E 5 -r1 -q -2 -K 5 -m 8 -W 7 -o '
  . $OUTFILE;

#========================================================================================================
#========================================================================================================
#========================================================================================================
#========================================================================================================

my $lines = `grep '>' $QUERY_FASTA | wc -l`;
chomp $lines;
my $lines_per_file = int( $lines / $NCPUS + 1 );

# something.fasta becomes something.1.fasta, something.2.fasta, etc.
# something becomes something.1, something.2 etc.
my @dots = split( /\./, $QUERY_FASTA );
my $filename = join(
    ".",
    (
        @dots[ 0 .. $#dots - $#dots >= 1 ? 0 : 1 ],
        "%i",
        $#dots >= 1 ? $dots[$#dots] : ""
    )
);
$filename =~ s/\.$//;

open( QUERY_FASTA, $QUERY_FASTA );
my $filei = 1;
my $linei = 0;

open( OUT, ">", sprintf( $filename, $filei ) );

while ( my $line = <QUERY_FASTA> ) {
    chomp $line;
    next unless $line;

    print OUT $line . "\n";

    if ( $linei > $lines_per_file ) {
        close(OUT);
        printf( $BLAST_CMD . " &\n",
            sprintf( $filename, $filei ),
            $SUBJECT_DATABASE, $filei
        );

        $linei = 0;
        ++$filei;
        open( OUT, ">", sprintf( $filename, $filei ) );
    }
    if ( $line =~ /^>/ ) { ++$linei; }
}
printf( $BLAST_CMD . " &\n",
    sprintf( $filename, $filei ),
    $SUBJECT_DATABASE, $filei
);
