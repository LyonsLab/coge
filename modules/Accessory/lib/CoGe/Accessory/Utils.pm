package CoGe::Accessory::Utils;

=head1 NAME

CoGe::Accessory::Utils

=head1 SYNOPSIS

Miscellaneous utility functions.

=head1 DESCRIPTION

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use POSIX qw( ceil );
use Data::GUID;
use Data::Dumper;
use File::Basename qw(fileparse basename);
use File::Find;

BEGIN {
    use vars qw ($VERSION $FASTA_LINE_LEN @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    $FASTA_LINE_LEN = 80;
    @ISA     = qw (Exporter);
    @EXPORT = qw( 
        units commify print_fasta get_unique_id get_link_coords 
        format_time_diff sanitize_name execute directory_size
        trim js_escape html_escape to_filename to_pathname
        is_fastq_file add_fastq_ext detect_paired_end to_number
        to_filename_without_extension is_gzipped is_bzipped2
        to_filename_base to_compressed_ext remove_fastq_ext
        fastq_description
    );
}

sub units {
    my $val = shift;
    
    # mdb added 5/2/16 for genomic values
    my $divisor = shift; # optional divisor
       $divisor = 1024 unless $divisor;
    my $type = shift; # optional unit type (e.g. 'bp' or 'B')
       $type = '' unless $type;

    if ( $val < $divisor ) {
        return $val . $type;
    }
    elsif ( $val < $divisor * $divisor ) {
        return ceil( $val / $divisor ) . 'K' . $type;
    }
    elsif ( $val < $divisor * $divisor * $divisor ) {
        return ceil( $val / ( $divisor * $divisor ) ) . 'M' . $type;
    }
    else {
        return ceil( $val / ( $divisor * $divisor * $divisor ) ) . 'G' . $type;
    }
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub trim {
    my $s = shift;
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;
    $s =~ s/^\"//;
    $s =~ s/\"$//;
    return $s;
}

sub js_escape {
    my $s = shift;
    $s =~ s/[\x00-\x1f]/ /g; # remove non-printable ascii chars
    $s =~ s/\'/\\'/g;
    $s =~ s/\"/\\"/g;
    return $s;
}

sub html_escape {
    my $s = shift;
    $s =~ s/[\x00-\x1f]/ /g; # remove non-printable ascii chars
    $s =~ s/\'/\&\#8216\;/g; # convert apostrophe char
    return $s;
}

sub print_fasta {
	my $fh = shift;
	my $name = shift;	# fasta section name
	my $pIn = shift; 	# reference to section data
    my $lineLength = shift; # optional
       $lineLength = $FASTA_LINE_LEN unless $lineLength;

	my $len = length $$pIn;
	my $ofs = 0;

	print {$fh} ">$name\n";
    while ($ofs < $len) {
    	print {$fh} substr($$pIn, $ofs, $lineLength) . "\n";
    	$ofs += $lineLength;
    }
}

sub get_unique_id {
	my $id = Data::GUID->new->as_hex;
	$id =~ s/^0x//;
	return $id;
}

sub get_link_coords { # mdb added 11/20/13 issue 254
	my ($start, $stop) = @_;
	return ($start, $stop) unless (defined $start and defined $stop);

	my $offset = 500;#int( abs($stop-$start+1) / 4 );
	($start, $stop) = ($stop, $start) if ($start > $stop);
	$start -= $offset;
	$stop  += $offset;
	return ($start, $stop);
}

# Convert a string to filename-friendly version
sub sanitize_name { 
    my $name = shift;
    
    return unless (defined $name);

    $name =~ s/\///g;   # remove /
    $name =~ s/\s+/_/g; # replace whitespace with _
    $name =~ s/\(//g;   # remove (
    $name =~ s/\)//g;   # remove )
    $name =~ s/://g;    # remove :
    $name =~ s/;//g;    # remove ;
    $name =~ s/#/_/g;   # replace # with _
    $name =~ s/'//g;    # remove '
    $name =~ s/"//g;    # remove "

    return $name;
}

sub format_time_diff {
    my $diff = shift;

    my $d = int($diff / (60*60*24));
    $diff -= $d * (60*60*24);
    my $h = int($diff / (60*60));
    $diff -= $h * (60*60);
    my $m = int($diff / 60);
    $diff -= $m * 60;
    my $s = $diff % 60;

    my $elapsed = '';
    $elapsed .= "${d}d " if $d > 0;
    $elapsed .= "${h}h " if $h > 0;
    $elapsed .= "${m}m " if $m > 0 && $d <= 0;
    $elapsed .= "${s}s" if $s > 0 && $d <= 0;

    return $elapsed;
}

# Return filename without path, e.g. /home/user/test/abc.xyz ==> abc.xyz
sub to_filename {
    my ($name, undef, undef) = fileparse(shift);
    return $name;
}

# Return filename without path and last extension, e.g. /home/user/test/abc.xyz.123 ==> abc.xyz
sub to_filename_without_extension {
    my ($name, undef, undef) = fileparse(shift, qr/\.[^.]*/);
    return $name;
}

# Return filename without path and everything before first period, e.g. /home/user/test/abc.xyz.123 ==> abc
sub to_filename_base {
    my ($base) = fileparse(shift, '\..*');
    return $base;
}

# Return path, e.g. /home/user/test/abc.xyz ==> /home/user/test
sub to_pathname {
    my (undef, $path, undef) = fileparse(shift);#, qr/\.[^.]*/);
    return $path;
}

sub execute {
    my $cmd = shift;
    my $error_msg = shift; # optional

    my @cmdOut = qx{$cmd};
    my $cmdStatus = $?;

    if ($cmdStatus != 0) {
        if ($error_msg) {
            print STDERR $error_msg;
        }
        else {
            say STDERR "error: command failed with rc=$cmdStatus: $cmd";
        }
    }

    return $cmdStatus;
}

sub is_fastq_file {
    my $filename = shift;
    return (
        $filename =~ /\.fastq$/      ||
        $filename =~ /\.fastq\.gz$/  ||
        $filename =~ /\.fastq\.bz2$/ ||
        $filename =~ /\.fq$/         ||
        $filename =~ /\.fq\.gz$/     ||
        $filename =~ /\.fq\.bz2$/
    );
}

sub add_fastq_ext {
    my $filename = shift;
    
    if ($filename =~ /\.gz/) {
        $filename =~ s/\.gz/\.fastq\.gz/;
    }
    elsif ($filename =~ /\.bz2/) {
        $filename =~ s/\.bz2/\.fastq\.bz2/;
    }
    else {
        $filename .= '.fastq';
    }
    
    return $filename;
}

sub remove_fastq_ext {
    my $filename = shift;

    foreach my $ext ('.fastq', '.fastq.gz', '.fastq.bz2', '.fq', '.fq.gz', '.fq.bz2') {
        if ($filename =~ /$ext$/) {
            $filename =~ s/$ext$//;
            return $filename;
        }
    }

    return $filename;
}

sub is_gzipped {
    my $filename = shift;
    return ( $filename =~ /\.gz$/ );
}

sub is_bzipped2 {
    my $filename = shift;
    return ( $filename =~ /\.bz2$/ );
}

sub to_compressed_ext {
    my $filename = shift;
    return '.gz' if ($filename =~ /\.gz$/);
    return '.bz2' if ($filename =~ /\.bz2$/);
    return '';
}

sub fastq_description {
    my ($fastq, $read_type) = @_;
    return (
        @$fastq > 2 ?
            join(' ', scalar(@$fastq), ($read_type eq 'paired' ? 'paired-end' : 'single-ended'), 'files') :
            join(', ', map { basename($_) } @$fastq)
    );
}

# Separate files based on last occurrence of _R1 or _R2 in filename
sub detect_paired_end {
    my $files = shift;
    my (@p1, @p2);
    foreach my $file (@$files) {
        my ($pair_id) = $file =~ /.+\_R?([12])/;
        if (defined $pair_id && $pair_id eq '1') { push @p1, $file; }
        else { push @p2, $file; }
    }
    return (\@p1, \@p2);
}

# Returns the total size of all files in the directory (recursively)
sub directory_size {
    my $path = shift;
    my $size = 0;
    find(sub { $size += -s if -f $_ }, $path);
    return $size;
}

# Force string to int/float
sub to_number {
    my $s = shift;
    return $s + 0;    
}

1;
