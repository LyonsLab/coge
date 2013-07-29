package CoGe::Services::JBrowse::Sequence;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON qq{encode_json};
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename 'fileparse';

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

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    #FIXME need permissions check here

    # Extract requested piece of sequence file
    my ( undef, $storagepath ) = fileparse( $genome->file_path );
    my $seqfile = $storagepath . '/chr/' . $chr;
    open( my $fh, $seqfile ) or die;
    seek( $fh, $start, 0 );
    read( $fh, my $seq, $len );
    close($fh);

    #print STDERR "$seqfile $len $seq\n";

    return qq{
		{ "features" : [
			{ "start": $start, "end": $end, "seq": "$seq" }
			]
		}
	};
}

1;
