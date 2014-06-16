package CoGe::Services::Data::Downloader;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Core::Storage qw( get_genome_seq );
use File::Spec;
use File::Slurp;
use Data::Dumper;

sub setup {
    my $self = shift;
    $self->run_modes( 'get' => 'get' );
    $self->mode_param('rm');
}

sub get {
    my $self = shift;
    my $page = $self->param('page');
    my $file = $self->query->param('file');
    my $gid = $self->query->param('gid');
    my $eid = $self->query->param('eid');
    my $dir = $self->query->param('dir');
    $dir = "" unless $dir;

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init();

    if ($gid) {
        my $genome = $db->resultset('Genome')->find($gid);

        if ( $genome->restricted
            and ( not defined $user or not $user->has_access_to_genome($genome) ) )
        {
            print STDERR "Data::Sequence::get access denied to genome $gid\n";
            return;
        }

        @path = ($conf->{SECTEMPDIR}, $page, "downloads", $gid, $dir, $file);
    } elsif ($eid) {
        my $exp = $db->resultset('Experiment')->find($eid);
        if ($exp->restricted and
            (not defined $user or not $user->has_access_to_experiment($exp))) {
            print STDERR "Data::Sequence::get access denied to experiment $eid\n";
            return;
        }

        @path = ($conf->{SECTEMPDIR}, $page, "downloads", $eid, $dir, $file);
    }

    $self->header_add( -attachment => $file );

    my $file_path = File::Spec->catdir(@path);
    my $content;

    say STDERR "FILE: $file_path";
    eval {
        $content = read_file($file_path) if -r $file_path;
    };

    return $content;
}

1;
