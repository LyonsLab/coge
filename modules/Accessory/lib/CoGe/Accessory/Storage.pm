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

use CoGeX;
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
    @EXPORT  = qw( get_genome_file get_genome_seq reverse_complement);

    #__PACKAGE__->mk_accessors();
}

################################################ subroutine header begin ##

=head2 _get_genome_path

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

sub _get_genome_path {    # private function, not exported
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

    my $base_path =
      CoGeX->get_conf('SEQDIR') . '/' . _get_genome_path($gid) . '/';

    my $file_path = $base_path . 'genome.faa';
    return $file_path if ( -r $file_path );

    $file_path = $base_path . "$gid.faa";
    return $file_path if ( -r $file_path );

    return;
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
    $start = 1 unless $start;
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
        my $cmd = CoGeX->get_conf('SAMTOOLS') . " faidx $file_path '$region'";

        #print STDERR "$cmd\n";
        $seq = qx{$cmd};            #my @cmdOut = qx{$cmd};
        unless ($fasta) {

            # remove header line
            $seq =~ s/^(?:.*\n)//;
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
            CoGeX->get_conf('SEQDIR') . '/'
          . _get_genome_path($gid)
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
