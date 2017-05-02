package CoGe::Accessory::BisQue;

use v5.14;
use strict;
use warnings;

use Data::Dumper;
use File::Basename qw(basename dirname);
use File::Spec::Functions qw(catfile);
use HTTP::Request;
use HTTP::Request::Common;

use CoGe::Accessory::IRODS qw(irods_get_base_path irods_imeta_ls irods_imkdir irods_iput irods_irm);
use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Accessory::Web qw(get_defaults);
use CoGe::Core::Storage qw(get_upload_path);

BEGIN {
    our (@ISA, $VERSION, @EXPORT);
    require Exporter;

    $VERSION = 0.0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw( create_bisque_image delete_bisque_image get_bisque_data_url );
}

sub create_bisque_image {
    my ($object, $upload, $user) = @_;

    my $dest = _get_dir($object->type, $object->id, $user);
    irods_imkdir($dest);
    $dest = catfile($dest, basename($upload->filename));
    my $source;
    if ($upload->asset->is_file) {
         $source = $upload->asset->path;
    }
    else {
        $source = get_upload_path('coge', get_unique_id());
        system('mkdir', '-p', $source);
        $source = catfile($source, $upload->filename);
        $upload->asset->move_to($source);
    }
    irods_iput($source, $dest);
    my $bisque_id = _register_image($dest, basename($upload->filename));
    _init_image($bisque_id, $object, $user);
    return $bisque_id, $upload->filename;
}

# object_type must be experiment, genome, or notebook
sub delete_bisque_image {
    my ($object_type, $object_id, $bisque_file, $bisque_id, $user) = @_;
    my $path = catfile(_get_dir($object_type, $object_id, $user), $bisque_file);
    my $res = irods_irm($path);
    my $ua = LWP::UserAgent->new();
    my $req = HTTP::Request->new(DELETE => 'https://bisque.cyverse.org/data_service/' . $bisque_id);
    $req->authorization_basic('coge', CoGe::Accessory::Web::get_defaults()->{BISQUE_PASS});
    $res = $ua->request($req);
}

sub get_bisque_data_url {
    my $id = shift;
    return 'https://bisque.cyverse.org/data_service/' . $id;
}

# sub _get_bisque_id {
#     my $dest = shift;
#     for my $i (0..9) {
#         sleep 5;
#         my $result = irods_imeta_ls($dest, 'ipc-bisque-id');
#         if (@$result == 4 && substr($result->[2], 0, 6) eq 'value:') {
#             my $bisque_id = substr($result->[2], 7);
#             chomp $bisque_id;
#             _init_image($bisque_id);
#             return $bisque_id;
#         }
#     }
#     warn 'unable to get bisque id';
# }

sub _get_dir {
    my ($object_type, $object_id, $user) = @_;
    return catfile(dirname(irods_get_base_path('coge')), 'bisque', $object_type, $object_id);
}

sub _init_image {
    my ($bisque_id, $object, $user) = @_;
    my $ua = LWP::UserAgent->new();
    my $req = HTTP::Request->new(GET => 'https://bisque.cyverse.org/data_service/' . $bisque_id);
    my $conf = CoGe::Accessory::Web::get_defaults();
    my $bisque_pass = $conf->{BISQUE_PASS};
    $req->authorization_basic('coge', $bisque_pass);
    my $res = $ua->request($req);
    my $content = $res->{_content};
    my $start = index($content, 'permission="') + 12;
    my $end = index($content, '"', $start);
    $content = substr($content, 0, $start) . 'published" hidden="true' . substr($content, $end);
    $content = substr($content, 0, length($content) - 2) . '>';
    $content .= '<tag name="Added by CoGe" type="link" value="' . $conf->{SERVER} . '" />';
    $content .= '<tag name="uploaded by" value="' . $user->info . '" />';
    $content .= '<tag name="' . $object->info . '" type="link" value="' . $conf->{SERVER} . $object->page . '?' . substr($object->type, 0, 1) . 'id=' . $object->id . '" />';
    if ($object->type eq 'genome') {
        $content .= '<tag name="' . $object->organism->info . '" type="link" value="' . $conf->{SERVER} . 'OrganismView.pl?gid=' . $object->id . '" />';
    } elsif ($object->type eq 'experiment') {
        my $genome = $object->genome;
        my $organism = $genome->organism;
        $content .= '<tag name="' . $organism->info . '" type="link" value="' . $conf->{SERVER} . 'OrganismView.pl?gid=' . $genome->id . '" />';
        $content .= '<tag name="' . $genome->info . '" type="link" value="' . $conf->{SERVER} . 'GenomeInfo.pl?gid=' . $genome->id . '" />';
    }
    $content .= '</image>';
    warn $content;
    $req = HTTP::Request->new(POST => 'https://bisque.cyverse.org/data_service/' . $bisque_id, ['Content-Type' => 'application/xml']);
    $req->authorization_basic('coge', $bisque_pass);
    $req->content($content);
    $res = $ua->request($req);
}

sub _register_image {
    my ($path, $name) = @_;
    my $ua = LWP::UserAgent->new();
    my $req = POST 'https://bisque.cyverse.org/import/transfer/', 'Content-Type' => 'form-data', 'Content' => [ file1_resource => '<image name="' . $name . '" value="irods://data.iplantcollaborative.org' . $path . '" owner="https://bisque.cyverse.org/data_service/00-h8Xeaz4KAexEt7L8kgEzwe"/>'];
    $req->authorization_basic('coge', CoGe::Accessory::Web::get_defaults()->{BISQUE_PASS});
    my $res = $ua->request($req);
    my $content = $res->{_content};
    my $start = index($content, 'resource_uniq="') + 15;
    my $end = index($content, '"', $start);
    my $bisque_id = substr($content, $start, $end - $start);
    return $bisque_id;
}

# sub set_bisque_visiblity {
#     my ($bisque_id, $public) = @_;
#     my $ua = LWP::UserAgent->new();
#     my $req = HTTP::Request->new(GET => 'https://bisque.cyverse.org/data_service/' . $bisque_id);
#     $req->authorization_basic('coge', CoGe::Accessory::Web::get_defaults()->{BISQUE_PASS});
#     my $res = $ua->request($req);
#     my $content = $res->{_content};
#     my $start = index($content, 'permission="') + 12;
#     my $end = index($content, '"', $start);
#     $content = substr($content, 0, $start) . ($public ? 'published' : 'private') . substr($content, $end);
#     $req = HTTP::Request->new(POST => 'https://bisque.cyverse.org/data_service/' . $bisque_id, ['Content-Type' => 'application/xml']);
#     $req->authorization_basic('coge', CoGe::Accessory::Web::get_defaults()->{BISQUE_PASS});
#     $req->content($content);
#     $res = $ua->request($req);
# }

# sub _share_bisque_image {
#     my ($bisque_id, $user) = @_;
#     my $ua = LWP::UserAgent->new();
#     my $req = HTTP::Request->new(GET => 'https://bisque.cyverse.org/data_service/user?resource_name=' . $user->name . '&wpublic=1');
#     my $res = $ua->request($req);
#     my $content = $res->{_content};
#     my $index = index($content, 'resource_uniq="') + 15;
#     if ($index != -1) {
#         my $coge_user_uniq = substr($content, $index, index($content, '"', $index) - $index);
#         $req = HTTP::Request->new(POST => 'https://bisque.cyverse.org/data_service/' . $bisque_id . '/auth?notify=false', ['Content-Type' => 'application/xml']);
#         $req->authorization_basic('coge', CoGe::Accessory::Web::get_defaults()->{BISQUE_PASS});
#         $req->content('<auth user="' . $coge_user_uniq . '" permission="edit" />');
#         $res = $ua->request($req);
#     }
# }

1;
