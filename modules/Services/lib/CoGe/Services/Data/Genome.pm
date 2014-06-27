package CoGe::Services::Data::Genome;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw(get_unique_id);
use CoGe::Accessory::IRODS;
use CoGe::Core::Storage qw(create_genome_from_file);
use Data::Dumper;
use JSON qw(decode_json encode_json);
use URI::Escape::JavaScript qw(escape);
use File::Path qw(mkpath);
#use CGI::Application::Plugin::JSON 'to_json';

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
    print STDERR "Data::Genome::load ticket=$ticket\n";
    print STDERR Dumper $q, "\n";
	print STDERR "POSTDATA: $post", "\n";
	print STDERR Dumper \%ENV, "\n";

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

    # Setup staging area
    my $LOAD_ID = get_unique_id();
	my $TEMPDIR = $conf->{SECTEMPDIR} . '/LoadGenome/' . $user->name . '/' . $LOAD_ID . '/';
    my $stagepath = $TEMPDIR . '/staging/';
    mkpath $stagepath;
	if (not -r $stagepath) {
		#$self->header_props( -status => 500 );
		return encode_json({
			success => JSON::false,
			error => 'Cannot create staging path'
		});
	}

	# Open log file
    my $logfile = $stagepath . '/log.txt';
    open( my $log, ">$logfile" ) or die "Error creating log file";
    print $log "Starting load genome $stagepath\n"; #. "name=$name description=$description version=$version type_id=$type_id restricted=$restricted org_id=$organism_id\n";

	# Copy files from iPlant Data Store to staging path
#	foreach my $item (@{$data->{items}}) {
#		print STDERR 'item: ', $item->{type}, ' ', $item->{path}, "\n";
#		next unless ($item->{type} and $item->{type} eq 'irods');
#		my $local_path = irods_get_file($item->{path}, $TEMPDIR);
#		print STDERR 'local_path=', $local_path, "\n";
#		unless (-r $local_path) {
#			#$self->header_props( -status => 500 );
#			return encode_json({
#				success => JSON::false,
#				error => 'Failed to retrieve file: ' . $item->{path}
#			});
#		}
#		push @files, $local_path;
#	}
#	print STDERR Dumper \@files, "\n";

	#FIXME: what happens when files take a long time to transfer? !!!!!!!!!!!!!!

    # Verify and decompress files #TODO move this into scripts/load_genome.pl
#    my @files;
#    foreach my $item (@$items) {
#        my $fullpath = $TEMPDIR . $item->{path};
#		 $self->header_props( -status => 500 );
#        return encode_json({ error => "File doesn't exist: $fullpath" }) if ( not -e $fullpath );
#        my ( $path, $filename ) = $item->{path} =~ /^(.+)\/([^\/]+)$/;
#        my ($fileext) = $filename =~ /\.([^\.]+)$/;
#
#        #print STDERR "$path $filename $fileext\n";
#        if ( $fileext eq 'gz' ) {
#            my $cmd = $P->{GUNZIP} . ' ' . $fullpath;
#            #print STDERR "$cmd\n";
#            `$cmd`;
#            $fullpath =~ s/\.$fileext$//;
#        }
#
#        #TODO support detecting/untarring tar files also
#        push @files, $fullpath;
#    }

#    my $CONFIGFILE = $ENV{COGE_HOME} . '/coge.conf';
#    my $cmd = $conf->{SCRIPTDIR} . "/load_genome.pl ";
#    $cmd .= '-user_name "' . $user->name . '" ';
#    $cmd .= '-name "' . escape($data->{name}) . '" '        if $data->{name};
#    $cmd .= '-desc "' . escape($data->{description}) . '" ' if $data->{description};
#    $cmd .= '-link "' . escape($data->{link}) . '" '        if $data->{link};
#    $cmd .= '-version "' . escape($data->{version}) . '" '  if $data->{version};
#    $cmd .= '-type_id ' . $type_id . ' ';
#    $cmd .= "-restricted " . $data->{restricted} . ' '      if $data->{restricted};
#    $cmd .= '-organism_id ' . $organism_id . ' ';
#    $cmd .= '-source_name "' . escape($source_name) . '" ';
#    $cmd .= "-staging_dir $stagepath ";
#    $cmd .= '-fasta_files "' . escape( join( ',', @files ) ) . '" ';
#    $cmd .= '-irods_files "' . escape( join( ',', map { $_->{path} } @{$data->{items}} ) ) . '" ';
#    $cmd .= "-config $CONFIGFILE"; #"-host $DBHOST -port $DBPORT -database $DBNAME -user $DBUSER -password $DBPASS";

    ($workflow_id, $error_msg) = create_genome_from_file(
        user => $user,
        metadata => {
            name => escape($data->{name}),
            description => escape($data->{description}),
            version => escape($data->{version}),
            source_name => escape($source_name),
            restricted => $data->{restricted},
            organism_id => $organism_id,
            type_id => $type_id
        },
        #files => \@files,
        irods => $data->{items}
    );
    unless ($workflow_id) {
        return encode_json({
            success => JSON::false,
            error => "Workflow submission failed: " . $error_msg
        });
    }

    # Call load script
#    print STDERR "$cmd\n";
#    print $log "$cmd\n";
#    close($log);
#    if ( !defined( my $child_pid = fork() ) ) {
#    	#$self->header_props( -status => 500 );
#        return encode_json({
#			success => JSON::false,
#			error => "Cannot fork process: $!"
#		});
#    }
#    elsif ( $child_pid == 0 ) {
#        print STDERR "child running: $cmd\n";
#        `$cmd`;
#        return;
#    }

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

sub irods_get_file {
    my ($src_path, $dest_path) = @_;

    my ($filename)   = $src_path =~ /([^\/]+)\s*$/;
    my ($remotepath) = $src_path =~ /(.*)$filename$/;

    my $localpath     = 'irods/' . $remotepath;
    my $localfullpath = $dest_path . $localpath;
    $localpath .= '/' . $filename;
    my $localfilepath = $localfullpath . '/' . $filename;
    print STDERR "irods_get_file $src_path $filename $localfilepath\n";

    my $do_get = 1;

	# Verify checksum -- slow for large files
    #	if (-e $localfilepath) {
    #		my $remote_chksum = irods_chksum($path);
    #		my $local_chksum = md5sum($localfilepath);
    #		$do_get = 0 if ($remote_chksum ne $local_chksum);
    #		print STDERR "$remote_chksum $local_chksum\n";
    #	}

    if ($do_get) {
        mkpath($localfullpath);
        CoGe::Accessory::IRODS::irods_iget( $src_path, $localfullpath );
        return $localfilepath;
    }
}

1;
