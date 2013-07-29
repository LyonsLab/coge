package CoGe::Services::Data::Notebook;
use base 'CGI::Application';

use CoGe::Accessory::Web;

sub setup {
    my $self = shift;
    $self->run_modes( 'create' => 'create', );
    $self->mode_param('rm');
}

sub create {
    my $self         = shift;
    my $name         = $self->query->param('name');
    my $description  = $self->query->param('description');
    my $list_type_id = $self->query->param('list_type_id');
    my $restricted   = $self->query->param('restricted');
    print STDERR "Data::Notebook::create\n";

    my ( $db, $user ) = CoGe::Accessory::Web->init;
    return $user->display_name;

    # Create the new list
    #    my $list = $coge->resultset('List')->create({
    #    	name => $name,
    #        description => $desc,
    #        list_type_id => $list_type_id,
    #        restricted => 1
    #    });
    #    return unless $list;

    # Set user as owner
    #    my $conn = $coge->resultset('UserConnector')->create({
    #    	parent_id => $USER->id,
    #    	parent_type => 5, #FIXME hardcoded to "user"
    #    	child_id => $list->id,
    #    	child_type => 1, #FIXME hardcoded to "list"
    #    	role_id => 2 #FIXME hardcoded to "owner"
    #    });
    #    return unless $conn;

    return qq{
		"blah"
	};
}

1;
