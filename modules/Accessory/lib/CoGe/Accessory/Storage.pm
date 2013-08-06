package CoGe::Accessory::Storage;

=head1 NAME

CoGe::Accessory::Storage

=head1 SYNOPSIS

Abstraction layer on top of genome storage sub-system created for issues
77 and 157.  All accesses to genome FASTA sequences should go through this
module.

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

use CoGe::Accessory::Web qw(get_defaults);
use File::Basename qw(fileparse);
use POSIX qw(floor);
use File::Spec::Functions;
use List::Util qw[min max];
use Data::Dumper;

BEGIN {
    use vars qw ($VERSION @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw (Exporter);
    @EXPORT =
      qw( get_genome_path get_genome_file index_genome_file get_genome_seq reverse_complement);

    #__PACKAGE__->mk_accessors();
}

################################################ subroutine header begin ##

=head2 get_genome_path

 Usage     : 
 Purpose   : This method determines the correct directory structure for storing
 			 the sequence files for a dataset group.
 Returns   : 
 Argument  : 
 Throws    : none
 Comments  : The idea is to build a dir structure that holds large amounts of
 			 files, and is easy to lookup based on genome ID number.
			 The strucuture is three levels of directorys, and each dir holds
			 1000 files and/or directorys.
			 Thus:
			 ./0/0/0/ will hold files 0-999
			 ./0/0/1/ will hold files 1000-1999
			 ./0/0/2/ will hold files 2000-2999
			 ./0/1/0 will hold files 1000000-1000999
			 ./level0/level1/level2/

See Also   : 

=cut

################################################## subroutine header end ##

sub get_genome_path {
    my $gid = shift;
    return unless $gid;

    my $level0 = floor( $gid / 1000000000 ) % 1000;
    my $level1 = floor( $gid / 1000000 ) % 1000;
    my $level2 = floor( $gid / 1000 ) % 1000;
    my $path   = catdir( $level0, $level1, $level2, $gid );

    return $path;
}

sub get_genome_file {
    my $gid = shift;
    return unless $gid;
    my $base_path = shift;    # optional

    my $seqdir = CoGe::Accessory::Web::get_defaults()->{'SEQDIR'};
    unless ($seqdir) {
        print STDERR
"Storage::get_genome_file: WARNING, conf file parameter SEQDIR is blank!\n";
    }

    unless ($base_path) {
        $base_path = $seqdir . '/' . get_genome_path($gid) . '/';
    }

    my $file_path = $base_path . 'genome.faa';
    return $file_path if ( -r $file_path );

    $file_path = $base_path . "$gid.faa";
    return $file_path if ( -r $file_path );

    return;
}

sub index_genome_file {
    my %opts      = @_;
    my $file_path = $opts{file_path};
    return unless $file_path;
    my $compress = $opts{compress};

    # Index fasta file
    my $samtools = CoGe::Accessory::Web::get_defaults()->{'SAMTOOLS'};
    unless ($samtools) {
        print STDERR
"Storage::index_genome_file: WARNING, conf file parameter SAMTOOLS is blank!\n";
    }
    my $cmd = "$samtools faidx $file_path";
    qx{ $cmd };
    if ( $? != 0 ) {
        print STDERR
          "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
        return $?;
    }

    # Optionally generate compressed version of fasta/index files
    if ($compress) {
        my $razip = CoGe::Accessory::Web::get_defaults()->{'RAZIP'};
        unless ($razip) {
            print STDERR
"Storage::index_genome_file: WARNING, conf file parameter RAZIP is blank!\n";
        }

        # Compress fasta file into RAZF using razip
        $cmd = "$razip -c $file_path > $file_path.razf";
        qx{ $cmd };
        if ( $? != 0 ) {
            print STDERR
              "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
            return $?;
        }

        # Index compressed fasta file
        $cmd = "$samtools faidx $file_path.razf";
        qx{ $cmd };
        if ( $? != 0 ) {
            print STDERR
              "Storage::index_genome_file: command failed with rc=$?: $cmd\n";
            return $?;
        }
    }

    return 0;
}

sub get_genome_seq {
    my %opts = @_;
    my $gid  = $opts{gid};    # required
    unless ($gid) {
        print STDERR "Storage::get_genome_seq: genome id not specified!\n";
        return;
    }
    my $chr   = $opts{chr};
    my $start = $opts{start};
    $start = 1 unless ( $start and $start > 0 );
    my $stop = $opts{stop};
    $stop = $opts{end} if ( not defined $stop );
    $stop = $start unless $stop;
    my $strand            = $opts{strand};
    my $format            = $opts{format};
    my $fasta             = ( defined $format and $format eq 'fasta' );
    my $FASTA_LINE_LENGTH = 60;
    my $seq;

#print STDERR "Storage::get_genome_seq gid=$gid chr=" . ($chr ? $chr : '') . " start=" . (defined $start ? $start : '') . " stop=" . (defined $stop ? $stop : '') . "\n";

    # Validate params
    my $len;
    if ( defined $start and defined $stop ) {
        ( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
        $len = abs( $stop - $start ) + 1;
    }

    # No chromosome specified, return whole genome fasta file
    unless ($chr) {
        my $file_path = get_genome_file($gid);
        open( my $fh, $file_path ) or die;
        read( $fh, $seq, -s $file_path );
        close($fh);
        return $seq;
    }

    # Determine file path - first try new indexed method, otherwise
    # revert to old method
    my $file_path = get_genome_file($gid);

    #print STDERR "file_path=$file_path\n";

    if ( -r "$file_path.fai" ) {    # new indexed method
                                    # Kluge chr/contig name
        if   ( $chr =~ /\D+/ ) { $chr = 'lcl|' . $chr; }
        else                   { $chr = 'gi|' . $chr; }

        # Extract requested piece of sequence file
        my $region =
          $chr . ( defined $start && defined $stop ? ":$start-$stop" : '' );
        my $samtools = CoGe::Accessory::Web::get_defaults()->{'SAMTOOLS'};
        unless ($samtools) {
            print STDERR
"Storage::get_genome_seq: WARNING, conf file parameter SAMTOOLS is blank!\n";
        }
        my $cmd = "$samtools faidx $file_path '$region'";

        #print STDERR "$cmd\n";
        $seq = qx{$cmd};    #my @cmdOut = qx{$cmd};
        unless ($fasta) {

            # remove header line
            $seq =~ s/^(?:.*\n)//;
            # remove end-of-lines
            $seq =~ s/\n//;
        }

        #print STDERR "$seq\n";
        my $cmdStatus = $?;
        if ( $cmdStatus != 0 ) {
            print STDERR "Error: command failed with rc=$cmdStatus: $cmd\n";
            return;
        }
    }
    else {    # old method
        $file_path =
            CoGe::Accessory::Web::get_defaults()->{'SEQDIR'} . '/'
          . get_genome_path($gid)
          . "/chr/$chr";

        # Extract requested piece of sequence file
        open( my $fh, $file_path ) or die "File not found: '$file_path'";
        seek( $fh, $start - 1, 0 ) if ( $start > 1 );

        $len = -s $file_path unless ( defined $len );
        if ($fasta) {
            $seq = ">$chr\n";
            while ( $len > 0 ) {
                my $count = read( $fh, my $s, min( $len, $FASTA_LINE_LENGTH ) );
                $seq .= $s . "\n";
                $len -= $count;
            }
        }
        else {
            read( $fh, $seq, $len );
        }

        close($fh);
    }

    $seq = reverse_complement($seq)
      if ( defined $strand and $strand =~ /-/ );  #FIXME broken for fasta format

    return $seq;
}

################################################ subroutine header begin ##

=head2 reverse_complement

 Usage     : 
 Purpose   : 
 Returns   : 
 Argument  : 
 Throws    : 
 Comments  : 

See Also   : 

=cut

################################################## subroutine header end ##

sub reverse_complement {
    my $seq   = shift;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCGatcg/TAGCtagc/;
    return $rcseq;
}

1;
