package CoGe::Services::Data::User;

use Mojo::Base 'Mojolicious::Controller';
#use IO::Compress::Gzip 'gzip';
use CoGeX;
use CoGe::Accessory::Web;

sub search {
    my $self = shift;
    $self->render(json => { success => Mojo::JSON->true });
}

sub fetch {
    my $self = shift;
    $self->render(json => { success => Mojo::JSON->true });
}

1;
