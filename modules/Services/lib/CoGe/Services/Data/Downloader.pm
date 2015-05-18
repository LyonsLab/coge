package CoGe::Services::Data::Downloader;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Storage qw(get_download_path get_workflow_log_file);
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
    my $wid = $self->query->param('wid');
    my $username = $self->query->param('username');
    my $uuid = $self->query->param('uuid') || ''; # optional
    my $attachment = $self->query->param('attachment') || 1;
    
    # Validate inputs
    unless ($gid || $eid || ($wid && $username)) {
        print STDERR "CoGe::Services::Data::Downloader invalid request\n";
        return;
    }
    
    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init();

    # Determine path to file
    my $file_path = '';
    if ($gid) { # genome sequence export
        my $genome = $db->resultset('Genome')->find($gid);
        if ( $genome->restricted
            and ( not defined $user or not $user->has_access_to_genome($genome) ) )
        {
            print STDERR "CoGe::Services::Data::Downloader access denied to genome $gid\n";
            return;
        }

        my $dl_path = get_download_path('genome', $gid, $uuid);
        $file_path = File::Spec->catdir($dl_path, $filename);
    } 
    elsif ($eid) { # experiment tarball export
        my $exp = $db->resultset('Experiment')->find($eid);
        if ($exp->restricted and
            (not defined $user or not $user->has_access_to_experiment($exp))) {
            print STDERR "CoGe::Services::Data::Downloader access denied to experiment $eid\n";
            return;
        }

        my $dl_path = get_download_path('experiment', $eid, $uuid);
        $file_path = File::Spec->catdir($dl_path, $filename);
    }
    elsif ($wid && $username) { # workflow debug log file
        $file_path = get_workflow_log_file($username, $wid);
    }
    
    say STDERR "CoGe::Services::Data::Downloader file=$file_path";
    return unless ($file_path);

    # Send file
    $self->header_add( -attachment => $filename ) if $attachment; # tell browser to download file
    my $content;
    eval {
        $content = read_file($file_path) if -r $file_path;
    };

    return $content;
}

1;
