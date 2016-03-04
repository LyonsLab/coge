package CoGe::Services::Data::JBrowse::Sequence;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Core::Storage qw( get_genome_seq );
use URI::Escape qw(uri_unescape);
use List::Util qw( min );

#sub stats_global {
#    print STDERR "JBrowse::Sequence::stats_global\n";
#    return qq{{}};
#}

sub features {
    my $self = shift;
    my $gid  = $self->stash('id');
    my $chr  = $self->stash('chr');
    $chr = uri_unescape($chr) if (defined $chr);
    my $size  = $self->param('seqChunkSize');
    my $start = $self->param('start');
    my $end   = $self->param('end');
    my $len   = $end - $start;
    print STDERR "JBrowse::Sequence::features gid=$gid chr=$chr size=$size start=$start end=$end\n";

    # Check params
    my $null_response = $self->render(json => { "features" => [] });
    if ( $end < 0 ) {
        return $null_response;
    }

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init; #TODO switch to CoGe::Services::Auth

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Adjust location - note: incoming coordinates are interbase!
    $start = 0 if ( $start < 0 );
    my $chrLen = $genome->get_chromosome_length($chr);
    $start = min( $start, $chrLen );
    $end   = min( $end,   $chrLen );
    if ( $start == $end ) {
        return $null_response;
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        return $null_response;
    }

    # Extract requested piece of sequence file
    my $seq = get_genome_seq(
        gid   => $gid,
        chr   => $chr,
        start => $start + 1, # convert from interbase to base
        stop  => $end
    );

    $self->render(json => {
        "features" => [ 
            { "start" => $start, "end" => $end, "seq" => $seq }
        ]
    });
}

1;
