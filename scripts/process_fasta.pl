#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Path qw(mkpath);
use File::Spec::Functions qw( catdir catfile );
use File::Basename qw(basename);
use URI::Escape::JavaScript qw(unescape);
use CoGe::Core::Genome qw(fix_chromosome_id);
use CoGe::Accessory::Utils qw(commify units print_fasta);

use vars qw( 
    $input_fasta_file $output_fasta_file $ignore_chr_limit 
    $keep_headers $keep_all_contigs
);
  
use constant MAX_CHROMOSOMES     => 200_000; # max number of chromosomes/contigs
use constant MIN_CHROMOSOME_SIZE => 1_000;   # min chromosome size (bp)
use constant MAX_PRINT           => 50;
use constant MAX_SEQUENCE_SIZE   => 100_000_000_000; # 100 Gbp total
use constant MAX_CHR_NAME_LENGTH => 255;

GetOptions(
    # Required arguments
    "input_fasta_file=s"  => \$input_fasta_file,
    "output_fasta_file=s" => \$output_fasta_file,
    
    # Optional arguments
    "ignore_chr_limit=i"  => \$ignore_chr_limit, # flag to ignore chromosome/contig limit # mdb added 3/9/15 COGE-595
    "keep_headers=i"      => \$keep_headers,     # flag to keep original headers (no parsing)
    "keep_all_contigs=i"  => \$keep_all_contigs  # disable min length check
);
$keep_all_contigs = 1; # mdb debug 8/9/16

# Process and verify parameters
unless ($input_fasta_file) {
    print STDOUT "log: error: input_fasta_file not specified\n";
    exit(-1);
}
$input_fasta_file = unescape($input_fasta_file);

unless ($output_fasta_file) {
    print STDOUT "log: error: output_fasta_file not specified\n";
    exit(-1);
}
$output_fasta_file = unescape($output_fasta_file);

# Ensure non-zero-length file
unless ( -s $input_fasta_file > 0 ) {
    print STDOUT "log: error: '$input_fasta_file' is zero-length\n";
    exit(-1);
}

# Ensure text file
if ( -B $input_fasta_file ) {
    print STDOUT "log: error: '$input_fasta_file' is a binary file\n";
    exit(-1);
}

# Process/validate the file
my ($seqLength, $numSequences, $numRemoved, $removedLength) = process_fasta_file( $input_fasta_file, $output_fasta_file );
if ($numRemoved && $removedLength) {
    print STDOUT "log: " . commify($numRemoved) . " sequences (", commify($removedLength) , " bp) removed due to size less than ", commify(MIN_CHROMOSOME_SIZE), " bp \n";
}
if ($numRemoved == $numSequences) {
    print STDOUT "log: error: all sequences (contigs) were removed due to being too small.\n";
    exit(-1);
}

if ( $seqLength > MAX_SEQUENCE_SIZE ) {
    print STDOUT "log: error: total sequence size exceeds limit of " . units(MAX_SEQUENCE_SIZE) . "\n";
    exit(-1);
}
#if ( !$ignore_chr_limit && $numSequences > MAX_CHROMOSOMES ) {
#    print STDOUT "log: error: too many sequences, limit is MAX_CHROMOSOMES\n";
#    exit(-1);
#}

if ( $numSequences == 0 or $seqLength == 0 ) {
    print STDOUT "log: error: couldn't parse sequences\n";
    exit(-1);
}

print STDOUT "log: " . commify($numSequences) . " sequences allowed totaling " . commify($seqLength) . " bp\n";

exit(0);

#-------------------------------------------------------------------------------

sub process_fasta_file {
    my $input_file  = shift;
    my $output_file = shift;
    print STDOUT "process_fasta_file: $input_file\n";
    
    my %seq;
    my $fileSize    = -s $input_file;
    my $lineNum     = 0;
    my $totalLength = 0;
    my $chrRemoved  = 0;
    my $removedLength = 0;
    
    # Open fasta file
    my $in; # file handle
    unless (open( $in, $input_file )) {
        print STDOUT "log: error: Error opening file for reading: $!\n";
        exit(-1);
    }
    
    # Process fasta file by sections
    $/ = '>'; # set file parsing delimiter
    while (my $section = <$in>) {
        $lineNum++;
        
        # Process the section in chunks.  There is a known problem where
        # Perl substitions fail on strings larger than 1GB.
        my $sectionName;
        my $sectionLen = length($section);
        my $filteredSeq;
        my $CHUNK_LEN = 10_000; # 8/19/15 COGE-647 mdb increased from 1000 for section headers >1000 chars
        my $ofs = 0;
        while ($ofs < $sectionLen) {
            my $chunk = substr($section, $ofs, $CHUNK_LEN);
            $ofs += $CHUNK_LEN;
            
            $chunk =~ s/>//g;
            $chunk =~ s/^\n+//m;
            $chunk =~ s/\n+$//m;
            next unless $chunk;
            
            if ($ofs == $CHUNK_LEN) { # first chunk
                ( $sectionName, $chunk ) = split(/\n/, $chunk, 2);
            }
            elsif ($sectionLen < $ofs) { # last chunk
                $chunk =~ s/\s+$//; # trim trailing whitespace
            }
            $chunk =~ s/\n//g;
            $chunk =~ s/\r//g;
            $filteredSeq .= $chunk;
        }
        next unless $filteredSeq;
    
        # Convert refseq (chromosome) name
        my $chr;
        if ($keep_headers) { # Don't modify refseq name
            $chr = $sectionName;
        }
        else {
            ($chr) = split(/\s+/, $sectionName);
            $chr = fix_chromosome_id($chr);
        }

        # Check validity of chr name and sequence
        if ( not defined $chr ) {
            print STDOUT "log: error parsing section header, line $lineNum, name='$sectionName'\n";
            exit(-1);
        }
        if ( length $filteredSeq == 0 ) {
            print STDOUT "log: warning: skipping zero-length section '$chr'\n";
            next;
        }
        if ( !$keep_all_contigs && length $filteredSeq < MIN_CHROMOSOME_SIZE ) { # mdb added 5/19/16
            #print STDOUT "Skipping section '$chr' less than ", MIN_CHROMOSOME_SIZE, " bp in size\n"; # removed, prints too much text and JEX truncates
            $chrRemoved++;
            $removedLength += length $filteredSeq;
            next;
        }
        if ( length($chr) > MAX_CHR_NAME_LENGTH ) {
            print STDOUT "log: error: section header name '$chr' is too long (>", MAX_CHR_NAME_LENGTH, " characters)\n";
            exit(-1);
        }
        if ( defined $seq{$chr} ) {
            print STDOUT "log: error: Duplicate section name '$chr'\n";
            exit(-1);
        }
        if ( $filteredSeq =~ /[^ACGTURYSWKMBDHVNX]/i ) { #/\W/ ) { # mdb modified 4/28/16 COGE-715, see http://www.bioinformatics.org/sms/iupac.html
            print STDOUT "log: error: sequence in section '$chr' contains invalid characters, perhaps this is not a nucleotide FASTA file?\n";
            #print STDOUT "log: error: chromosome name $sectionName\n";
            #print STDOUT "log: error: $filteredSeq\n";
            exit(-1);
        }

        # Append sequence to master file
        my $out;
        unless (open( $out, ">>$output_fasta_file" )) {
            print STDOUT "log: error: Couldn't open genome.faa\n";
            exit(-1);
        }
        my $head = $chr =~ /^\d+$/ ? "gi" : "lcl";
        $head .= "|" . $chr;
        print_fasta($out, $head, \$filteredSeq);
        close($out);

        $seq{$chr} = { size => length $filteredSeq, file => $input_file };

        # Print log message
        my $count = keys %seq;
        $totalLength += length $filteredSeq;
        if ( !$ignore_chr_limit && ($count >= MAX_CHROMOSOMES or $totalLength > MAX_SEQUENCE_SIZE) ) {
            print STDOUT "log: warning: skipping remaining sections due to number exceeding ", MAX_CHROMOSOMES, "\n" if ($count >= MAX_CHROMOSOMES);
            goto DONE;
        }
        if ( $count <= MAX_PRINT ) {
            my $filename = basename($input_file);
            print STDOUT "log: Processed chr '$chr' (", commify(length($filteredSeq)), " bp)\n";
        }
        elsif ( $count == MAX_PRINT + 1 ) {
            print STDOUT "log: (only showing first ", MAX_PRINT, " chromosomes)\n";
        }
        elsif ( ( $count % 10000 ) == 0 ) {
            print STDOUT "log: Processed "
              . commify($count) . " ("
              . units($totalLength) . ", "
              . int( 100 * $totalLength / $fileSize )
              . "%) sequences so far ...\n";
        }
    }

    DONE:    
    close($in);
    return $totalLength, scalar(keys %seq), $chrRemoved, $removedLength;
}


