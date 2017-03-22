package CoGe::Services::API::Sequence;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Core::Storage qw( get_genome_seq );
use Data::Dumper;

sub setup {
    my $self = shift;
    $self->run_modes( 'get' => 'get' );
    $self->mode_param('rm');
}

sub get {
    my $self = shift;
    my $gid  = $self->param('gid');
    return unless $gid;
    my $chr   = $self->param('chr');
    my $start = $self->query->param('start');
    my $stop  = $self->query->param('stop');
    $stop = $self->query->param('end') if ( not defined $stop );
    my $strand = $self->query->param('strand');
    print STDERR "Data::Sequence::get gid=$gid chr=$chr start=$start stop=$stop\n";

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init();

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    unless ($genome) {
    	print STDERR "Data::Sequence::get genome $gid not found in db\n";
    	return; # empty response (no data)
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        print STDERR "Data::Sequence::get access denied to genome $gid\n";
        return; # empty response (no data)
    }

    # Force browser to download as attachment
    if ( not defined $chr or $chr eq '' ) {
        my $genome_name = sanitize_name($genome->organism->name);
        $genome_name = 'genome_'.$gid unless $genome_name;
        $self->res->headers->content_disposition("attachment; filename=$genome_name.faa;");
    }

    # Get sequence from file
    $self->render( text => get_genome_seq(
        gid   => $gid,
        chr   => $chr,
        start => $start,
        stop  => $stop
    ) );
}

1;
