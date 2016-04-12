package CoGe::Services::API::JBrowse::Sequence;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Core::Storage qw( get_genome_seq );
use CoGe::Services::Auth qw(init);
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
    print STDERR "matt $chr\n";
    $chr = uri_unescape($chr) if (defined $chr);
#    my $size  = $self->param('seqChunkSize'); # mdb removed 4/11/16 -- unused
    my $start = $self->param('start');
    my $end   = $self->param('end');
    my $len   = $end - $start;
    print STDERR "JBrowse::Sequence::features gid=$gid chr=$chr start=$start end=$end\n";

    # Check params
    my $null_response = $self->render(json => { "features" => [] });
    if ( $end < 0 ) {
        print STDERR "JBrowse::Sequence::features ERROR, invalid 'end' value\n";
        return $null_response;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    unless ($genome) {
        print STDERR "JBrowse::Sequence::features ERROR, genome '$gid' not found\n";
        return $null_response;
    }

    # Adjust location - note: incoming coordinates are interbase!
    $start = 0 if ( $start < 0 );
    my $chrLen = $genome->get_chromosome_length($chr);
    $start = min( $start, $chrLen );
    $end   = min( $end,   $chrLen );
    if ( $start == $end ) {
        print STDERR "JBrowse::Sequence::features ERROR, start==end: start=$start end=$end chrLen=$chrLen\n";
        return $null_response;
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        print STDERR "JBrowse::Sequence::features ERROR, access denied\n";
        return $null_response;
    }

    # Extract requested piece of sequence file
    my $seq = get_genome_seq(
        gid   => $gid,
        chr   => $chr,
        start => $start + 1, # convert from interbase to base
        stop  => $end
    );
    unless ($seq && length($seq)) {
        print STDERR "JBrowse::Sequence::features ERROR, no sequence returned\n";
        return $null_response;
    }

    $self->render(json => {
        "features" => [ 
            { "start" => $start, "end" => $end, "seq" => $seq }
        ]
    });
}

1;
