package CoGe::Services::Data::Genome;
use base 'CGI::Application';

################################################################################
# This is a legacy endpoint in support of genome loads from the DE.  Need to 
# migrate the DE to the new API and deprecate this at some point.
################################################################################

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Accessory::IRODS;
use CoGe::Core::Storage;
use CoGe::Builder::CommonTasks;
use Data::Dumper;
use JSON qw(decode_json encode_json);
use URI::Escape::JavaScript qw(escape);
use File::Path qw(mkpath);
use File::Spec::Functions qw(catfile);

sub setup {
    my $self = shift;
    $self->run_modes( 'load' => 'load' );
    $self->mode_param('rm');
}

sub load {
    my $self = shift;
    my $q = $self->query();
	#my $ticket = $q->param('ticket'); # ugh: when POSTDATA is present, it doesn't have the ticket from URL
	my $post = $q->param('POSTDATA');
    print STDERR "Data::Genome::load ticket=$ticket\n",
        Dumper $q, "\n",
	    "POSTDATA: $post", "\n",
	    Dumper \%ENV, "\n";

	# Problem: when status is anything but 200 it adds on an html doc
	# to the reponse even if type is 'application/json'
	# Need to migrate to Apache::REST

	# Validate POST data
	if (!$post) {
		#$self->header_props( -status => 500, -type => 'application/json' );
		return encode_json({
			success => JSON::false,
			error => 'Missing data'
		});
	}
	my $data = decode_json($post);
	if (!$data) {
		#$self->header_props( -status => 500, -type => 'application/json' );
		return encode_json({
			success => JSON::false,
			error => 'Invalid request format'
		});
	}
	if (!$data->{items} or !@{$data->{items}}) {
		#$self->header_props( -status => 500, -type => 'application/json' );
		return encode_json({
			success => JSON::false,
			error => 'No data items'
		});
	}

 	# Connect to the database
 	my $ticket = $data->{ticket};
 	my $url = 'https://genomevolution.org'.$ENV{REQUEST_URI}; #FIXME hardcoded !!!! # mdb changed 6/26/14 from http://coge.iplantcollaborative.org per request from DE
 	$url =~ s/\/$//;
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init(ticket => $ticket, ticket_type => 'proxy', url => $url);
	print STDERR "Data::Genome::load user=", $user->name, "\n";

	# Must be authenticated to load a genome
	if (not $user or $user->is_public) {
		$self->header_props(-status => '401 Unauthorized');
		return;
	}

	# Set defaults #TODO review these values
	$data->{name} = '' unless (defined $data->{name});
	$data->{description} = '' unless (defined $data->{desription});
	$data->{link} = '' unless (defined $data->{link});
	$data->{version} = '1' unless (defined $data->{version});
	$data->{restricted} = (defined $data->{restricted} && ($data->{restricted} eq 'true' || $data->{restricted} eq '1') ? 1 : 0 );
	my $type_id = 1; # hardcode to "unmasked sequence"
	my $source_name = $user->name; #FIXME change to user's display_name
	my $organism_id = 38378; # FIXME hardcoded to "test" organism, need to change to "unknown" or something

    my ($workflow_id, $error_msg) = create_genome_from_file(
        user => $user,
        metadata => {
            name => escape($data->{name}),
            description => escape($data->{description}),
            version => escape($data->{version}),
            source_name => escape($source_name),
            restricted => $data->{restricted},
            type_id => $type_id
        },
        organism_id => $organism_id,
        irods => $data->{items}
    );
    unless ($workflow_id) {
        return encode_json({
            success => JSON::false,
            error => "Workflow submission failed: " . $error_msg
        });
    }

    # Get tiny link
    my $tiny_link = CoGe::Accessory::Web::get_tiny_link(
        url => $conf->{SERVER} . "LoadGenome.pl?job_id=$workflow_id"
    );
    unless ($tiny_link) {
    	#$self->header_props( -status => 500 );
		return encode_json({
			success => JSON::false,
			error => 'Link generation failed'
		});
    }

	print STDERR "link: $tiny_link\n";
	return encode_json({
		success => JSON::true,
		link => $tiny_link
	});
}

sub create_genome_from_file { #TODO use CoGe::Builder::Load::Genome pipeline instead
    my %opts = @_;
    my $user = $opts{user};
    my $irods = $opts{irods};
    my $metadata = $opts{metadata};
    my $organism_id = $opts{organism_id};
    #print STDERR "Storage::create_genome_from_file ", Dumper $metadata, " ", Dumper $files, "\n";

    # Connect to workflow engine
    my $conf = CoGe::Accessory::Web::get_defaults();
    my $jex = CoGe::Accessory::Jex->new( host => $conf->{JOBSERVER}, port => $conf->{JOBPORT} );
    unless (defined $jex) {
        return (undef, "Could not connect to JEX");
    }

    # Create the workflow
    my $workflow = $jex->create_workflow( name => 'Create Genome', init => 1 );
    unless ($workflow and $workflow->id) {
        return (undef, 'Could not create workflow');
    }

    # Setup log file, staging, and results paths
    my ($staging_dir, $result_dir) = get_workflow_paths($user->name, $workflow->id);
    $workflow->logfile( catfile($result_dir, 'debug.log') );

    # Create jobs to retrieve irods files
    my @staged_files;
    foreach my $item (@$irods) {
        next unless ($item->{type} eq 'irods');
        my $iget_task = create_iget_job(irods_path => $item->{path}, local_path => $staging_dir);
        unless ( $iget_task ) {
            return (undef, "Could not create iget task");
        }
        $workflow->add_job($iget_task);
        push @staged_files, $iget_task->{outputs}[0];
    }

    # Create load job
    my $load_task = create_load_genome_job(
        user => $user,
        staging_dir => $staging_dir,
        wid => $workflow->id,
        organism_id => $organism_id,
        input_files => \@staged_files,
        irods_files => $irods, # for metadata markup
        metadata => $metadata,
    );
    unless ( $load_task ) {
        return (undef, "Could not create load task");
    }
    $workflow->add_job($load_task);

    # Submit the workflow
    my $result = $jex->submit_workflow($workflow);
    if ($result->{status} =~ /error/i) {
        return (undef, "Could not submit workflow");
    }

    return ($result->{id}, undef);
}

1;
