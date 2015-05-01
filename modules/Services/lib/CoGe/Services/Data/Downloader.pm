package CoGe::Services::Data::Downloader;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use File::Spec;
use File::Slurp;

sub setup {
    my $self = shift;
    $self->run_modes( 'get' => 'get' );
    $self->mode_param('rm');
}

sub get {
    my $self = shift;
    my $filename = $self->query->param('filename');
    my $gid = $self->query->param('gid');
    my $eid = $self->query->param('eid');
    my $uuid = $self->query->param('uuid') || ''; # optional
    
    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init();

    # Determine path to file
    my @path;
    if ($gid) {
        my $genome = $db->resultset('Genome')->find($gid);
        if ( $genome->restricted
            and ( not defined $user or not $user->has_access_to_genome($genome) ) )
        {
            print STDERR "Data::Sequence::get access denied to genome $gid\n";
            return;
        }

        @path = ($conf->{SECTEMPDIR}, "downloads", "genome", $gid, $uuid, $filename);
    } 
    elsif ($eid) {
        my $exp = $db->resultset('Experiment')->find($eid);
        if ($exp->restricted and
            (not defined $user or not $user->has_access_to_experiment($exp))) {
            print STDERR "Data::Sequence::get access denied to experiment $eid\n";
            return;
        }

        @path = ($conf->{SECTEMPDIR}, "downloads", "experiment", $eid, $uuid, $filename);
    }
    
    my $file_path = File::Spec->catdir(@path);
    say STDERR "CoGe::Services::Data::Downloader file: $file_path";
    return unless (@path);

    # Send file as attachment
    $self->header_add( -attachment => $filename );
    my $content;
    eval {
        $content = read_file($file_path) if -r $file_path;
    };

    return $content;
}

1;
