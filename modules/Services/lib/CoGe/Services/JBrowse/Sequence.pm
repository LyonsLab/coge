package CoGe::Services::JBrowse::Sequence;
use base 'CGI::Application';

use CoGe::Accessory::Web;
use CoGe::Accessory::Storage qw( get_genome_seq );

#use File::Basename 'fileparse';

sub setup {
    my $self = shift;
    $self->run_modes(
        'stats_global' => 'stats_global',
        'features'     => 'features',
    );
    $self->mode_param('rm');
}

sub stats_global {
    print STDERR "Sequence::stats_global\n";
    return qq{{}};
}

sub features {
    my $self  = shift;
    my $gid   = $self->param('gid');
    my $chr   = $self->param('chr');
    my $size  = $self->query->param('seqChunkSize');
    my $start = $self->query->param('start');
    my $end   = $self->query->param('end');
    my $len   = $end - $start;
    print STDERR
      "Sequence::features gid=$gid chr=$chr size=$size start=$start end=$end\n";

    # Check params
    return qq{{"features" : []}} if ( $end < 1 );
    $start = 1 if ( $start < 1 );

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        return qq{{"features" : []}};
    }

    # Extract requested piece of sequence file
    #    my ( undef, $storagepath ) = fileparse( $genome->file_path );
    #    my $seqfile = $storagepath . '/chr/' . $chr;
    #    open( my $fh, $seqfile ) or die;
    #    seek( $fh, $start, 0 );
    #    read( $fh, my $seq, $len );
    #    close($fh);
    my $seq = get_genome_seq(
        gid   => $gid,
        chr   => $chr,
        start => $start,
        stop  => $end - 1
    );

    return qq{{"features" : [{"start": $start, "end": $end, "seq": "$seq"}]}};
}

1;
