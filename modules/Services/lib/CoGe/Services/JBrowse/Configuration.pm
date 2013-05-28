package CoGe::Services::JBrowse::Configuration;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON;
use Data::Dumper;
use Sort::Versions;
use Cwd 'abs_path';

my $coge_conf;

sub setup {
	my $self = shift;

	#FIXME - move this into service.pl
	$coge_conf = abs_path($0);
	$coge_conf =~ s/services\/JBrowse\/service\.pl/coge\.conf/;

	$self->run_modes(
		'refseq_config'	=> 'refseq_config',
		'track_config' 	=> 'track_config',
	);
	$self->start_mode('track_config');
	$self->mode_param('rm');
}

sub refseq_config {
	my $self = shift;
	my $gid = $self->query->param('gid');
	my $SEQ_CHUNK_SIZE = 20000;
	print STDERR "Configuration::refseq_config gid=$gid\n";

	# Load config file
	my $P = CoGe::Accessory::Web::get_defaults($coge_conf);
	my $DBNAME = $P->{DBNAME};
	my $DBHOST = $P->{DBHOST};
	my $DBPORT = $P->{DBPORT};
	my $DBUSER = $P->{DBUSER};
	my $DBPASS = $P->{DBPASS};

	# Connect to the database
	my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
	my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
	#$coge->storage->debugobj(new DBIxProfiler());
	#$coge->storage->debug(1);

	my $genome = $coge->resultset('Genome')->find($gid);
	return unless $genome;

	my @chromosomes;
	foreach my $chr ($genome->genomic_sequences) {
		push @chromosomes, {
			name => $chr->chromosome,
			length => $chr->sequence_length,
			seqChunkSize => $SEQ_CHUNK_SIZE,
			start => 0,
			end => $chr->sequence_length - 1
		};
	}

	return encode_json(\@chromosomes);
}

sub track_config {
	my $self = shift;
	my $gid = $self->query->param('gid');
	my $cas_ticket = $self->query->param('ticket');
	print STDERR "Configuration::track_config gid=$gid\n";

	my $P = CoGe::Accessory::Web::get_defaults($coge_conf);
	my $DBNAME = $P->{DBNAME};
	my $DBHOST = $P->{DBHOST};
	my $DBPORT = $P->{DBPORT};
	my $DBUSER = $P->{DBUSER};
	my $DBPASS = $P->{DBPASS};

	# Connect to the database
	my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
	my $coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);
	#$coge->storage->debugobj(new DBIxProfiler());
	#$coge->storage->debug(1);

	my $COOKIE_NAME = $P->{COOKIE_NAME};
	my $USER = undef;
	($USER) = CoGe::Accessory::Web->login_cas( ticket => $cas_ticket, coge => $coge, this_url => $FORM->url() ) if ($cas_ticket);
	($USER) = CoGe::Accessory::LogUser->get_user( cookie_name => $COOKIE_NAME, coge => $coge ) unless $USER;

	my $genome = $coge->resultset('Genome')->find($gid);
	return unless $genome;

	my @tracks;

	# Add reference sequence track
	push @tracks,
	{
		chunkSize => 20000,
		baseUrl => "services/JBrowse/service.pl/sequence/$gid/",
		type => "SequenceTrack",
		storeClass => "JBrowse/Store/SeqFeature/REST",
		label => "sequence",
		key => "Sequence",
		formatVersion => 1,
		# CoGe-specific stuff
		coge => {
			id => $gid,
			type => 'sequence'
		}
	};

	# Add GC content track
	push @tracks,
	{
		#baseUrl => "http://geco.iplantcollaborative.org/rchasman/gc/$gid/",
		baseUrl => "services/JBrowse/track/gc/$gid/",
		type => "JBrowse/View/Track/Wiggle/Density",
		storeClass => "JBrowse/Store/SeqFeature/REST",
		label => "gc_content",
		key => "GC Content",
		style => {
			height => 50,
			pos_color => 'rgb(0, 135, 0)',
			neg_color => '#f00',
			bg_color  => 'rgba(232, 255, 220, 0.4)'
		},
		# CoGe-specific stuff
		coge => {
			id => $gid,
			type => ''
		}
	};

	# Add gene annotation tracks
	my $num;
	foreach $ds ($genome->datasets) {
		next if (not $ds->has_gene_annotation);
		my $dsid = $ds->id;
		push @tracks,
		{
			#baseUrl => "services/JBrowse/service.pl/annotation/$dsid/",
			#baseUrl => "http://geco.iplantcollaborative.org/rchasman/annotation/$dsid/",
			baseUrl => "services/JBrowse/track/annotation/$gid/",
            autocomplete => "all",
            track => "genes",
            label => "genes",
            key => "Genes" . ($num++ ? " $num" : ''),
            #type => "FeatureTrack",
            type => "CoGe/View/Track/CoGeFeatures",
            storeClass => "JBrowse/Store/SeqFeature/REST",
            onClick => "FeatAnno.pl?dsg=$gid;chr={chr};start={start};stop={end}",
            maxFeatureScreenDensity => 2,
         	style => {
            	className => "cds",
            	featureScale => 0,
            	histScale => 0
         	},
         	# CoGe-specific stuff
         	coge => {
         		id => $dsid,
         		type => 'annotation'
         	}
        };
	}

	# Add experiment tracks
	my %all_notebooks;
	foreach $e (sort experimentcmp $genome->experiments) {
		#FIXME need permission check here
		next if ($e->deleted);
		my $eid = $e->id;

		# Make a list of notebook id's
		my @notebooks = map {$_->id} $e->notebooks;
		map { $all_notebooks{$_->id} = $_ } $e->notebooks;

		# Make a list of annotations
		my @annotations;
		foreach my $a ($e->annotations) {
			push @annotations,
			{
				type  => $a->type->name,
				text  => $a->annotation,
				image => ($a->image_id ? 'image.pl?id='.$a->image_id : undef),
				link  => $a->link
			}
		}

		push @tracks,
		{
			baseUrl => "services/JBrowse/service.pl/experiment/$eid/",
		    autocomplete => "all",
		    track => "exp$eid",
		    label => "exp$eid",
		    key => $e->name,
		    type => "CoGe/View/Track/Wiggle/MultiXYPlot",#"JBrowse/View/Track/Wiggle/XYPlot",
		    storeClass => "JBrowse/Store/SeqFeature/REST",
#		    style => {
#		    	pos_color => $color,
#		    	neg_color => $color
#		    },
		    # CoGe-specific stuff
		    coge => {
		    	id => $eid,
		    	type => 'experiment',
		    	name => $e->name,
		    	description => $e->description,
		    	notebooks => (@notebooks ? \@notebooks : undef),
		    	annotations => (@annotations ? \@annotations : undef),
		    	menuOptions => [
					{ label => 'ExperimentView',
					  action => "function() { window.open( 'ExperimentView.pl?eid=$eid' ); }"
					  # url => ... will open link in a dialog window
					}
				]
		    }
		};
	}

	# Add notebook tracks
	foreach my $n (sort {$a->name cmp $b->name} values %all_notebooks) {
		next if ($n->restricted and not $USER->has_access_to_list($n));

		my $nid = $n->id;

		# Make a list of experiments
		my @experiments;
		foreach my $e ($n->experiments) {
			push @experiments,
			{
				id => $e->id,
				name => $e->name
			}
		}

		push @tracks,
		{
			key => $n->name,
			baseUrl => "services/JBrowse/service.pl/experiment/notebook/$nid/",
		    autocomplete => "all",
		    track => "notebook$nid",
		    label => "notebook$nid",
		    type => "CoGe/View/Track/Wiggle/MultiXYPlot",
		    storeClass => "JBrowse/Store/SeqFeature/REST",
			# CoGe-specific stuff
			show_average => 0,
			coge => {
				id => $nid,
				type => 'notebook',
				name => $n->name,
				description => $n->description,
				experiments => (@experiments ? \@experiments : undef),
				menuOptions => [
					{ label => 'NotebookView',
					  action => "function() { window.open( 'NotebookView.pl?nid=$nid' ); }"
					  # url => ... will open link in a dialog window
					}
				]
			}
		};
	}

	return encode_json({
		formatVersion => 1,
		dataset_id => 'coge',
		plugins => [ 'CoGe' ],
		trackSelector => {
			type => 'CoGe/View/TrackList/CoGe',
			# plugin-specific stuff
			# <none>
		},
		tracks => \@tracks,
	});
}

# FIXME this comparison routine is duplicated elsewhere
sub experimentcmp {
	no warnings 'uninitialized'; # disable warnings for undef values in sort
	versioncmp($b->version, $a->version) || $a->name cmp $b->name || $b->id cmp $a->id
}

1;
