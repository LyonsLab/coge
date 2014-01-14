package CoGe::Services::Data::Downloader;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Storage qw( get_genome_seq );
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

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init();
    my $genome = $db->resultset('Genome')->find($gid);

    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        print STDERR "Data::Sequence::get access denied to genome $gid\n";
        return;
    }

    $self->header_add( -attachment => $file );

    my @path = ($conf->{SECTEMPDIR}, $page, "downloads", $gid, $file);
    my $file_path = File::Spec->catdir(@path);
    my $content;

    say STDERR "FILE: $file_path";
    eval {
        $content = read_file($file_path) if -r $file_path;
    };

    return $content;
}

1;
