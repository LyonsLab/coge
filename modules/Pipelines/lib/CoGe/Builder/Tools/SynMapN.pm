package CoGe::Builder::Tools::SynMap3D;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::JEX::Jex;
use CoGe::Builder::Tools::SynMap qw( defaults gen_org_name );
use File::Spec::Functions;

sub pre_build { # override superclass method
	my $self = shift;

	# Connect to JEX
    my $jex = CoGe::JEX::Jex->new( host => $self->conf->{JOBSERVER}, port => $self->conf->{JOBPORT} );
    unless ($jex) {
        CoGe::Exception::Generic->throw(message => "Couldn't connect to JEX");
    }
    $self->jex($jex);

	# Initialize workflow -- NOTE: init => 0 means that a previous identical workflow will be reused when submitted
    $self->workflow( $jex->create_workflow(name => $self->get_name, init => 0 ) );
    return unless $self->workflow;

	# Set site_url attribute
	my %opts = ( %{ defaults() }, %{ $self->params } );
	$self->site_url( $opts{tinylink} || get_query_link( $self->conf, $self->db, %opts ) );
	my $requester = $self->request->requester;
    if ($requester) { # request is from internal web page - external API requests will not have a 'requester' field
        my $page = $requester->{page}; # page name used for logging
        $self->page($page) if $page;
    }

	return;
}

sub build {
	my $self = shift;

    my $synmap = CoGe::Builder::Tools::SynMap->new($self);
    $synmap->build();
    $self->add_to_all($synmap);
}

sub get_name {
    my $self = shift;
    
    my $description;
    foreach my $key (keys $self->params) {
        next unless $key =~ /^genome_id/;
        
        my ($genome_id) = $key =~ /(\d+)$/;
        my ( $org_name ) = gen_org_name(
            db        => $self->db,
            genome_id => $genome_id
        );
        
        $description .= ($description ? 'v. ' : '') . $org_name . ' ';
    }
    
	return $description;
}

__PACKAGE__->meta->make_immutable;

1;
