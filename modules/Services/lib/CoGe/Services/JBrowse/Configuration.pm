package CoGe::Services::JBrowse::Configuration;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON;
use Data::Dumper;
use Sort::Versions;
use Cwd 'abs_path';
use Time::HiRes qw ( time );

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
	foreach my $chr (sort {$b->sequence_length <=> $a->sequence_length} $genome->genomic_sequences) {
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

	my $start_time = time;

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
		type => "CoGe/View/Track/GC_Content",
		storeClass => "JBrowse/Store/SeqFeature/REST",
    	track => "gc_content",
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
	
	print STDERR 'time1: ' . (time - $start_time) . "\n";

	# Add gene annotation tracks
	my $num;
	foreach $ds ($genome->datasets) {
		next if (not $ds->has_gene_annotation);
		my $dsid = $ds->id;
		push @tracks,
		{
			#baseUrl => "services/JBrowse/service.pl/annotation/$dsid/",
			#baseUrl => "http://geco.iplantcollaborative.org/rchasman/annotation/$dsid/",
			baseUrl => "services/JBrowse/track/annotation/$dsid/",
            autocomplete => "all",
            track => "genes",
            label => "genes",
            key => "Genes" . ($num++ ? " $num" : ''),
            #type => "FeatureTrack",
            type => "CoGe/View/Track/CoGeFeatures",
            storeClass => "JBrowse/Store/SeqFeature/REST",
            onClick => "FeatAnno.pl?dsg=$gid;chr={chr};start={start};stop={end}",
            maxFeatureScreenDensity => 100,
            maxHeight => 100000,
         	style => {
            	className => "cds",
            	histScale => 0,
            	labelScale => 0.02,
            	featureScale => 0.00000005,
         	},
         	# CoGe-specific stuff
         	coge => {
         		id => $dsid,
         		type => 'annotation'
         	}
        };
	}

	print STDERR 'time2: ' . (time - $start_time) . "\n";
	
	# Create a fake "all experiments" notebook for all genome's experiments
#	push @tracks,
#		{
#			key => 'All Experiments',
#			baseUrl => "services/JBrowse/service.pl/experiment/genome/$gid/",
#		    autocomplete => "all",
#		    track => "notebook0",
#		    label => "notebook0",
#		    type => "CoGe/View/Track/Wiggle/MultiXYPlot",
#		    storeClass => "JBrowse/Store/SeqFeature/REST",
#			# CoGe-specific stuff
#			showAverage => 0,
#			coge => {
#				id => 0, # use id of 0 to represent all experiments
#				type => 'notebook',
#				name => 'All Experiments',
#				description => '',
#				experiments => ,
##				menuOptions => [
##					{ label => 'NotebookView',
##					  action => "function() { window.open( 'NotebookView.pl?nid=$nid' ); }"
##					  # url => ... will open link in a dialog window
##					}
##				]
#			}
#		};

	# Add experiment tracks
	my %all_notebooks;
	my %expByNotebook;
	foreach $e (sort experimentcmp $genome->experiments) { #{}, {prefetch => { 'experiment_annotations' => 'annotation_type' }})) {
		#FIXME need permission check here
		next if ($e->deleted);
		my $eid = $e->id;

		# Make a list of notebook id's
		my @notebooks;
		foreach my $n ($e->notebooks) {
			push @notebooks, $n->id;
			$all_notebooks{$n->id} = $n;
			push @{ $expByNotebook{$n->id} }, { id => $n->id, name => $n->name };
		}
#		push @notebooks, 0; # add fake "all experiments" notebook

		# Make a list of annotations
#		my @annotations;
#		foreach my $a ($e->experiment_annotations) {
#			push @annotations,
#			{
#				type  => $a->annotation_type->name,
#				text  => $a->annotation,
#				image => ($a->image_id ? 'image.pl?id='.$a->image_id : undef),
#				link  => $a->link
#			}
#		}

		push @tracks,
		{
			baseUrl => "services/JBrowse/service.pl/experiment/$eid/",
		    autocomplete => "all",
		    track => "experiment$eid",
		    label => "experiment$eid",
		    key => $e->name,
		    type => "CoGe/View/Track/Wiggle/MultiXYPlot",#"JBrowse/View/Track/Wiggle/XYPlot",
		    storeClass => "JBrowse/Store/SeqFeature/REST",
#		    style => {
#		    	pos_color => $color,
#		    	neg_color => $color
#		    },
		    # CoGe-specific stuff
		    onClick => "ExperimentView.pl?embed=1&eid=$eid",
		    showHoverScores => 1,
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

	print STDERR 'time3: ' . (time - $start_time) . "\n";
	
	# Add notebook tracks
	foreach my $n (sort {$a->name cmp $b->name} values %all_notebooks) {
		next if ($n->restricted and not $USER->has_access_to_list($n));
		my $nid = $n->id;
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
			onClick => "NotebookView.pl?embed=1&lid=$nid",
			showAverage => 0,
			showHoverScores => 1,
			coge => {
				id => $nid,
				type => 'notebook',
				name => $n->name,
				description => $n->description,
				experiments => (@{$expByNotebook{$nid}} ? $expByNotebook{$nid} : undef),
				menuOptions => [
					{ label => 'NotebookView',
					  action => "function() { window.open( 'NotebookView.pl?lid=$nid' ); }"
					  # url => ... will open link in a dialog window
					}
				]
			}
		};
	}

	print STDERR 'time4: ' . (time - $start_time) . "\n";print STDERR 'time1: ' . (time - $start_time) . "\n";

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
