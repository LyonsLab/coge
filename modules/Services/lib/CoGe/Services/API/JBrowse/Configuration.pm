package CoGe::Services::API::JBrowse::Configuration;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use URI::Escape qw(uri_escape);
use Data::Dumper;
use Sort::Versions;
use Time::HiRes qw(time);

use CoGe::Services::Auth qw(init);
use CoGeDBI qw(get_table get_user_access_table get_experiments get_distinct_feat_types);
use CoGe::Core::Chromosomes;
use CoGe::Core::Experiment qw(experimentcmp);

my %expTypeToName = (
    1 => 'quant',
    2 => 'snp'
);

my $DEBUG_PERFORMANCE = 0;

sub refseq_config {
    my $self           = shift;
    my $gid            = $self->param('gid');
    my $SEQ_CHUNK_SIZE = 20000;

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
      	print STDERR "JBrowse::Configuration::refseq_config access denied to genome $gid\n";
       	return '{}';
    }

    my @chromosomes;
	my $c = CoGe::Core::Chromosomes->new($genome->id);
	while ($c->next) {
		push @chromosomes, {
            name         => uri_escape($c->name),
            length       => $c->length,
            seqChunkSize => $SEQ_CHUNK_SIZE,
            start        => 0,
            end          => $c->length - 1
		};
	}

    $self->render(json => \@chromosomes);
}

sub _annotations {
    my ($type, $eid, $db) = @_;
    my $sth = $db->storage->dbh->prepare('SELECT name,annotation FROM ' . $type . '_annotation JOIN annotation_type ON annotation_type.annotation_type_id=' . $type . '_annotation.annotation_type_id WHERE ' . $type . '_id=' . $eid . ' ORDER BY name');
    $sth->execute();
    my $annotations;
    while (my $row = $sth->fetch) {
        $annotations .= "\n" . $row->[0] . ': ' . $row->[1];
    }
    return $annotations;
}

sub track_config {
    my $self = shift;
    my $gid  = $self->param('gid');
    my $experiment_ids = $self->param('experiments');
    my $start_time = time; # for performance testing

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

	# Admins have ability to simulate other users using the "user_id" query parameter
	my $user_id = $self->param('user_id');
	if (defined $user_id && $user->is_admin && $user_id != $user->id) {
	    my $u = $db->resultset('User')->find($user_id);
	    if (defined $u) {
	        warn "Switching to user '", $u->name, "'";
	        $user = $u;
	    }
	}

    my $SERVER_NAME = $conf->{SERVER};
    my $JBROWSE_API = $SERVER_NAME . 'api/v1/jbrowse'; #TODO move to config file

    # Get genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Check permissions
    if ($genome->restricted and (not defined $user or not $user->has_access_to_genome($genome))) {
      	$self->render({error => "JBrowse::Configuration::track_config access denied to genome $gid"});
       	return;
    }

    my @tracks;

    #
    # Add reference sequence track
    #
    push @tracks, {
        chunkSize     => 20000,
        baseUrl       => "$JBROWSE_API/sequence/$gid/",
        type          => "SequenceTrack",
        storeClass    => "JBrowse/Store/SeqFeature/REST",
        track         => 'sequence',
        label         => "sequence",
        key           => "Sequence",
        formatVersion => 1,
        coge => {
            id   => $gid,
            type => 'sequence'
        }
    };

    #
    # Add GC content track
    #
    push @tracks, {
        baseUrl    => "$JBROWSE_API/track/gc/$gid/",
        type       => "CoGe/View/Track/GC_Content",
        storeClass => "JBrowse/Store/SeqFeature/REST",
        track      => "gc_content",
        label      => "gc_content",
        key        => "GC Content",
        style      => {
            height    => 50,
            pos_color => 'rgb(0, 135, 0)',
            neg_color => '#f00',
            bg_color  => 'rgba(232, 255, 220, 0.4)'
        },
        coge => {
            id   => $gid,
            type => 'gc_content'
        }
    };

    print STDERR 'time1: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

    #
    # Add overall feature group for all datasets
    #
    my @feat_type_names = grep( !/chromosome/, $genome->distinct_feature_type_names ); # exclude "chromosome" features
    if (@feat_type_names) {
        # Add main overall composite track
        push @tracks, {
            baseUrl      => "$JBROWSE_API/track/annotation/$gid/",
            autocomplete => "all",
            track        => "feature_group0",
            label        => "feature_group0",
            key          => "Features: all",
            type         => "CoGe/View/Track/CoGeFeatures",
            description  => "note, description",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            onClick      => $SERVER_NAME . 'FeatAnno.pl?dsg=' . $gid . ';chr={chr};start={start};stop={end}',
            maxFeatureScreenDensity => 20,
            maxHeight               => 100000,
            minSubfeatureWidth      => 4,
            style                   => {
#                labelScale               => 0.02,
                arrowheadClass           => "arrowhead",
                className                => "generic_parent",
                maxDescriptionLength     => 70,
                showLabels               => JSON::true,
                centerChildrenVertically => JSON::true,
                subfeatureClasses        => { match_part => "match_part7" }
            },
            coge => {
                id        => 0,
                type      => 'feature_group',
                collapsible => 1
            }
        };

        print STDERR 'time1a: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

        # Add a track for each feature type
        foreach my $type_name ( sort @feat_type_names ) {
            push @tracks, {
                baseUrl => "$JBROWSE_API/track/annotation/$gid/types/$type_name/",
                autocomplete => "all",
                track        => "features$type_name",
                label        => "features$type_name",
                key          => $type_name,
                type         => "JBrowse/View/Track/HTMLFeatures",
                storeClass   => "JBrowse/Store/SeqFeature/REST",
                region_stats => 1, # see HTMLFeatures.js, force calls to stats/region instead of stats/global
                onClick      => $SERVER_NAME . 'FeatAnno.pl?dsg=' . $gid . ';chr={chr};start={start};stop={end};type=' . $type_name,
                maxFeatureScreenDensity => 1000,     #50,
                maxHeight               => 100000,
                style                   => {
                    arrowheadClass           => "arrowhead",
                    className                => "generic_parent",
                    histScale                => 0.002,
                    minSubfeatureWidth       => 6,
                    maxDescriptionLength     => 70,
                    showLabels               => JSON::true,
                    description              => "note, description",
                    centerChildrenVertically => JSON::true,
                    subfeatureClasses        => { match_part => "match_part7" }
                },
                coge => {
                    id      => "$type_name",
                    type    => 'features',
                    dataset_id => 0
                }
            };
        }
    }

    print STDERR 'time1b: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

    #
    # Add a feature group for each dataset
    #
    if ( @{$genome->datasets} > 1) {
        foreach my $ds ( sort { $a->name cmp $b->name } $genome->datasets ) {
            @feat_type_names = grep( !/chromosome/, $ds->distinct_feature_type_names ); # exclude "chromosome" features
            if (@feat_type_names) {
                my $dsid = $ds->id;
                my $dsname = $ds->name || 'ds'.$dsid;

                # Add overall feature track
                push @tracks, {
                    baseUrl      => "$JBROWSE_API/track/annotation/$gid/datasets/$dsid",
                    autocomplete => "all",
                    track        => "feature_group".$dsid,
                    label        => "feature_group".$dsid,
                    key          => "Features: ".$dsname,
                    type         => "CoGe/View/Track/CoGeFeatures",
                    description  => "note, description",
                    storeClass   => "JBrowse/Store/SeqFeature/REST",
                    onClick      => $SERVER_NAME . 'FeatAnno.pl?ds=' . $dsid . ';chr={chr};start={start};stop={end}',
                    maxFeatureScreenDensity => 20,
                    maxHeight               => 100000,
                    minSubfeatureWidth      => 4,
                    style                   => {
                        labelScale               => 0.02,
                        arrowheadClass           => "arrowhead",
                        className                => "generic_parent",
                        maxDescriptionLength     => 70,
                        showLabels               => JSON::true,
                        centerChildrenVertically => JSON::true,
                        subfeatureClasses        => { match_part => "match_part7" }
                    },
                    coge => {
                        id         => $dsid,
                        type       => 'feature_group',
                        collapsible => 1
                    }
                };

                # Add a track for each feature type
                foreach my $type_name ( sort @feat_type_names ) {
                    push @tracks, {
                        baseUrl => "$JBROWSE_API/track/annotation/$gid/types/$type_name/",
                        autocomplete => "all",
                        track        => 'features'.$dsid.'_'.$type_name,
                        label        => 'features'.$dsid.'_'.$type_name,
                        key          => $type_name,
                        type         => "JBrowse/View/Track/HTMLFeatures",
                        storeClass   => "JBrowse/Store/SeqFeature/REST",
                        region_stats => 1, # see HTMLFeatures.js, force calls to stats/region instead of stats/global
                        onClick      => $SERVER_NAME . 'FeatAnno.pl?ds=' . $dsid . ';chr={chr};start={start};stop={end};type=' . $type_name,
                        maxFeatureScreenDensity => 1000,     #50,
                        maxHeight               => 100000,
                        style                   => {
                            arrowheadClass           => "arrowhead",
                            className                => "generic_parent",
                            histScale                => 0.002,
                            minSubfeatureWidth       => 6,
                            maxDescriptionLength     => 70,
                            showLabels               => JSON::true,
                            description              => "note, description",
                            centerChildrenVertically => JSON::true,
                            subfeatureClasses        => { match_part => "match_part7" }
                        },
                        coge => {
                            id         => $dsid.'_'.$type_name,
                            type       => 'features',
                            dataset_id => $dsid
                        }
                    };
                }
            }
        }
    }

    print STDERR 'time2: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

    #
    # Add experiment tracks
    #
    my %experiments;    # all experiments hashed by id -- used later for creating "All Experiments" section
    my %notebooks;      # all notebokos hashed by id -- used later for creating individual notebooks
    my %expByNotebook;  # all experiments hashed by notebook id -- used later for creating individual notebooks

    my $connectors;
    $connectors = get_user_access_table($db->storage->dbh, $user->id) if $user;
    my $allNotebooks = get_table($db->storage->dbh, 'list');
    my $allNotebookConn = get_table($db->storage->dbh, 'list_connector', ['child_id', 'list_connector_id'], {child_type => 3});
    my @genome_experiments = get_experiments($db->storage->dbh, $genome->id);
    if ($experiment_ids) {
        my @ids = split(/,/, $experiment_ids);
        @genome_experiments = grep { my $id = $_->{experiment_id}; grep { $id == $_ } @ids; } @genome_experiments;
    }
    foreach my $e ( sort experimentcmp @genome_experiments ) { # sort experimentcmp $genome->experiments
        next if ( $e->{deleted} );
        my $eid = $e->{experiment_id};
        my $role = $connectors->{3}{$eid} if $connectors;
        $role = $role->{role_id} if $role;
        next if $e->{restricted} && !($user && $user->is_admin) && !($connectors && $connectors->{3}{$eid});
        $experiments{$eid} = $e;

        # Build a list of notebook id's
        my @notebooks;
        foreach my $conn (values %{$allNotebookConn->{$eid}}) {
            my $nid = $conn->{parent_id};
            my $n = $allNotebooks->{$nid};
            next if $n->{deleted};
            next if $n->{restricted} && !($user && $user->is_admin) && !($connectors && $connectors->{1}{$nid});
            push @notebooks, $nid;
            $notebooks{$nid} = $n;
            push @{ $expByNotebook{$nid} },
              {
                value => $eid,
                label => $e->{name}
              };
        }
        push @notebooks, 0;    # add fake "all experiments" notebook

        my ($type, $featureScale, $histScale, $labelScale, $ext);
        if (!$e->{data_type} or $e->{data_type} == 1) { #FIXME hardcoded data_type 'quantitative'
			$type = 'CoGe/View/Track/Wiggle/XYPlot';
			$featureScale = 0.001;
			#$histScale = 0.05;
			#$labelScale = 0.1;
		}
		elsif ($e->{data_type} == 2) { #FIXME hardcoded data_type 'polymorphism'
			$type = 'CoGe/View/Track/CoGeVariants';
			$featureScale = 0.0001;
			$histScale = 0.01;
			$labelScale = 0.5;
		}
		elsif ($e->{data_type} == 3) { #FIXME hardcoded data_type 'alignment'
			$type = 'CoGe/View/Track/CoGeAlignment';
			$featureScale = 0.005;
			$histScale = 0.01;
			$labelScale = 0.5;
		}
        elsif ($e->{data_type} == 4) { #FIXME hardcoded data_type 'marker'
            $type = 'CoGe/View/Track/Markers'; #"JBrowse/View/Track/HTMLFeatures";
            $histScale = 0.002;
            $labelScale = 0.5;
        }

        my $style = {
            featureScale => $featureScale,
            histScale    => $histScale,
            labelScale   => $labelScale,
            className    => '{type}',
            histCss      => 'background-color:' . _getFeatureColor($eid),
            featureCss   => 'background-color:' . _getFeatureColor($eid)
        };
        $style->{showMismatches} = JSON::true if $e->{data_type} == 3;

        my $coge = {
            id          => $eid,
            type        => 'experiment',
            data_type   => $e->{data_type} || 1,
            editable    => (($user && $user->is_admin) || ($role && ($role == 2 || $role == 3))) ? 1 : 0, # mdb added 2/6/15 #TODO move this obscure code into an API
            name        => $e->{name},
            description => $e->{description},
            notebooks   => ( @notebooks ? \@notebooks : undef ),
            onClick     => "ExperimentView.pl?embed=1&eid=$eid",
            annotations => _annotations('experiment', $eid, $db)
        };
        $coge->{editable} = 1 if (($user && $user->is_admin) || ($role && ($role == 2 || $role == 3)));
        $coge->{restricted} = 1 if $e->{restricted};
        $coge->{role} = $role if $role;
        my $config = {
            baseUrl      => "$JBROWSE_API/experiment/$eid/",
            autocomplete => "all",
            track        => "experiment$eid",
            label        => "experiment$eid",
            key          => ( $e->{restricted} ? '🔒 ' : '' ) . $e->{name},
            type         => $type,
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            region_feature_densities => 1, # enable histograms in store
            showLabels   => JSON::false,
            style => $style,
            coge => $coge
        };
        push @tracks, $config;
    }

    print STDERR 'time3: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

    #
    # Create a fake "All Experiments" notebook track
    #
    if ( keys %experiments ) {
        push @tracks, {
            key     => 'All Experiments',
            baseUrl => "$JBROWSE_API/experiment/genome/$gid/",
            autocomplete => "all",
            track        => "notebook0",
            label        => "notebook0",
            type         => "CoGe/View/Track/Wiggle/XYPlot",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            style        => { featureScale => 0.001 },
            coge => {
                id          => 0, # use id of 0 to represent all experiments
                type        => 'notebook',
                collapsible => 1,
                name        => 'All Experiments',
                description => ''
            }
        };
    }

    #
    # Add notebook tracks
    #
    foreach my $n ( sort { $a->{name} cmp $b->{name} } values %notebooks ) {
        my $nid = $n->{list_id};
        my $role = $connectors->{1}{$nid} if $connectors;
        $role = $role->{role_id} if $role;
        my $coge = {
            id      => $nid,
            type    => 'notebook',
            collapsible => 1,
            name        => $n->{name},
            description => $n->{description},
            experiments => ( @{ $expByNotebook{$nid} } ? $expByNotebook{$nid} : undef ),
            onClick     => "NotebookView.pl?embed=1&lid=$nid",
            annotations => _annotations('list', $nid, $db)
        };
        $coge->{editable} = 1 if (($user && $user->is_admin) || ($role && ($role == 2 || $role == 3)));
        $coge->{restricted} = 1 if $n->{restricted};
        $coge->{role} = $role if $role;
        push @tracks, {
            key     => ( $n->{restricted} ? '🔒 ' : '' ) . $n->{name},
            baseUrl => "$JBROWSE_API/experiment/genome/$gid/notebook/$nid",
            autocomplete => "all",
            track        => "notebook$nid",
            label        => "notebook$nid",
            type         => "CoGe/View/Track/Wiggle/XYPlot",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            style        => { featureScale => 0.001 },
            showHoverScores => 1,
            coge            => $coge
        };
    }

    print STDERR 'time4: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

    $self->render(json => {
        formatVersion => 1,
        dataset_id    => 'coge',
        plugins       => ['CoGe','HideTrackLabels'],
        trackSelector => {
            type => 'CoGe/View/TrackList/CoGe'
        },
        tracks => \@tracks,
    });
}

# FIXME this is duplicated in XYPlot.js, need to either make totally client side or totally server side
sub _getFeatureColor {
    my $id = shift;
    return '#'
      . sprintf( "%06X",
        ( ( ( ( $id * 1234321 ) % 0x1000000 ) | 0x444444 ) & 0xe7e7e7 ) );
}

1;
