package CoGe::Services::JBrowse::Configuration;
use base 'CGI::Application';

use Switch;
use JSON;
use URI::Escape qw(uri_escape);
use Data::Dumper;
use Sort::Versions;
use Cwd qw(abs_path);
use Time::HiRes qw(time);
use File::Spec::Functions qw(catdir);

use CoGeX;
use CoGeDBI qw(get_table get_user_access_table get_experiments get_distinct_feat_types);
use CoGe::Accessory::Web;
use CoGe::Core::Chromosomes;
use CoGe::Core::Experiment qw(experimentcmp);

my %expTypeToName = (
    1 => 'quant',
    2 => 'snp'
);

my $DEBUG_PERFORMANCE = 0;

sub setup {
    my $self = shift;
    $self->run_modes(
        'refseq_config' => 'refseq_config',
        'track_config'  => 'track_config',
    );
    $self->mode_param('rm');
}

sub refseq_config {
    my $self           = shift;
    my $gid            = $self->query->param('gid');
    my $SEQ_CHUNK_SIZE = 20000;
    print STDERR "JBrowse::Configuration::refseq_config gid=$gid\n";

    # Connect to the database
    my ( $db, $user ) = CoGe::Accessory::Web->init;

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
#    foreach my $chr ( sort { $b->sequence_length <=> $a->sequence_length }
#        $genome->genomic_sequences )
#    {
#        push @chromosomes,
#          {
#            name         => uri_escape($chr->chromosome), # mdb changed 12/17/13 issue 266
#            length       => $chr->sequence_length,
#            seqChunkSize => $SEQ_CHUNK_SIZE,
#            start        => 0,
#            end          => $chr->sequence_length - 1
#          };
#    }
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

    return encode_json( \@chromosomes );
}

sub track_config {
    my $self = shift;
    my $gid  = $self->query->param('gid');
    print STDERR "JBrowse::Configuration::track_config gid=$gid\n";
    my $start_time = time; # for performance testing
    
    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;
    
    # Admins have ability to simulate other users using the "user_id" query parameter
#    my $user_id = $self->query->param('user_id');
#    if (defined $user_id && $user->is_admin && $user_id != $user->id) {
#        my $simulated_user = $db->resultset('User')->find($user_id);
#        if (defined $simulated_user) {
#            print STDERR "Switching to user '", $simulated_user->name, "'\n";
#            $user = $simulated_user;
#        }
#    }

    # Get server name for constructing URLs
    my $SERVER_NAME = $conf->{SERVER};#$ENV{SERVER_NAME}; # mdb added 12/11/14 for COGE-568
    my $JBROWSE_API = $SERVER_NAME . 'api/v1/jbrowse'; #TODO move to config file
    #print STDERR "SERVER_NAME = $SERVER_NAME\n", "JBROWSE_API = $JBROWSE_API\n";

    # Get genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
      	print STDERR "JBrowse::Configuration::track_config access denied to genome $gid\n";
       	return '{}';
    }

    my @tracks;

    #
    # Add reference sequence track
    #
    push @tracks, {
        chunkSize     => 20000,
        baseUrl       => "$JBROWSE_API/sequence/$gid/", #"https://$SERVER_NAME/services/JBrowse/service.pl/sequence/$gid/",
        type          => "SequenceTrack",
        storeClass    => "JBrowse/Store/SeqFeature/REST",
        label         => "sequence",
        key           => "Sequence",
        formatVersion => 1,

        # CoGe-specific stuff
        coge => {
            id   => $gid,
            type => 'sequence'
        }
    };

    #
    # Add GC content track
    #
    push @tracks, {
        baseUrl    => "$JBROWSE_API/track/gc/$gid/", #"$SERVER_NAME/coge/services/JBrowse/track/gc/$gid/",
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

        # CoGe-specific stuff
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
            track        => "features",
            label        => "features",
            key          => "Features: all",
            type         => "CoGe/View/Track/CoGeFeatures",
            description  => "note, description",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            onClick      => "$SERVER_NAME/FeatAnno.pl?dsg=$gid;chr={chr};start={start};stop={end}",
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

            # CoGe-specific stuff
            coge => {
                id        => 0,
                type      => 'feature_group',
                classes   => ['coge-tracklist-collapsible'],
                collapsed => 1 #FIXME move into CSS
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
                onClick      => "$SERVER_NAME/FeatAnno.pl?dsg=$gid;chr={chr};start={start};stop={end};type=$type_name",
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

                # CoGe-specific stuff
                coge => {
                    id      => "$type_name",
                    type    => 'features',
                    classes => ['coge-tracklist-indented'],
                    collapsed => 1, #FIXME move into CSS
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
                    track        => "features_ds".$dsid,
                    label        => "features_ds".$dsid,
                    key          => "Features: ".$dsname,
                    type         => "CoGe/View/Track/CoGeFeatures",
                    description  => "note, description",
                    storeClass   => "JBrowse/Store/SeqFeature/REST",
                    onClick      => "$SERVER_NAME/FeatAnno.pl?ds=$dsid;chr={chr};start={start};stop={end}",
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

                    # CoGe-specific stuff
                    coge => {
                        id         => $dsid,
                        type       => 'feature_group',
                        classes    => ['coge-tracklist-collapsible'],
                        collapsed  => 1, #FIXME move into CSS
                    }
                };

                # Add a track for each feature type
                foreach my $type_name ( sort @feat_type_names ) {
                    push @tracks, {
                        baseUrl => "$JBROWSE_API/track/annotation/$gid/types/$type_name/",
                        autocomplete => "all",
                        track        => 'features_ds'.$dsid.'_'.$type_name,
                        label        => 'features_ds'.$dsid.'_'.$type_name,
                        key          => $type_name,
                        type         => "JBrowse/View/Track/HTMLFeatures",
                        storeClass   => "JBrowse/Store/SeqFeature/REST",
                        region_stats => 1, # see HTMLFeatures.js, force calls to stats/region instead of stats/global
                        onClick      => "$SERVER_NAME/FeatAnno.pl?ds=$dsid;chr={chr};start={start};stop={end};type=$type_name",
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

                        # CoGe-specific stuff
                        coge => {
                            id         => $dsid.'_'.$type_name,
                            type       => 'features',
                            classes    => ['coge-tracklist-indented'],
                            collapsed  => 1, #FIXME move into CSS
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
   
    # mdb added 2/9/15 for performance improvement, COGE-166
    my $connectors = get_user_access_table($db->storage->dbh, $user->id);
    my $allNotebooks = get_table($db->storage->dbh, 'list');
    my $allNotebookConn = get_table($db->storage->dbh, 'list_connector', ['child_id', 'list_connector_id'], {child_type => 3});

    foreach $e ( sort experimentcmp get_experiments($db->storage->dbh, $genome->id) ) { # sort experimentcmp $genome->experiments
        next if ( $e->{deleted} );
        my $eid = $e->{experiment_id};
        next if ($e->{restricted} && !$user->admin && defined $connectors && !$connectors->{3}{$eid}); #next unless $user->has_access_to_experiment($e); #TODO move into an API
        $experiments{$eid} = $e;

        # Build a list of notebook id's
        my @notebooks;
        foreach my $conn (values %{$allNotebookConn->{$eid}}) {
            my $nid = $conn->{parent_id};
            my $n = $allNotebooks->{$nid};
            next if $n->{deleted};
            next if ($n->{restricted} && !$user->admin && !$connectors->{1}{$nid});
            push @notebooks, $nid;
            $notebooks{$nid} = $n;
            push @{ $expByNotebook{$nid} },
              {
                id   => $eid,
                name => $e->{name},
                type => $expTypeToName{ $e->{data_type} }
              };
        }
        push @notebooks, 0;    # add fake "all experiments" notebook

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

        my ($type, $featureScale, $histScale, $labelScale);
        if (!$e->{data_type} or $e->{data_type} == 1) { #FIXME hardcoded data_type 'quantitative'
			$type = 'CoGe/View/Track/Wiggle/MultiXYPlot';
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
			$type = 'CoGe/View/Track/CoGeAlignment';#"JBrowse/View/Track/Alignments2";
			$featureScale = 0.005;
			$histScale = 0.01;
			$labelScale = 0.5;
		}
        elsif ($e->{data_type} == 4) { #FIXME hardcoded data_type 'marker'
            $type = "JBrowse/View/Track/HTMLFeatures";
            $histScale = 0.002;
            $labelScale = 0.5;
        }

        push @tracks, {
            baseUrl      => "$JBROWSE_API/experiment/$eid/",
            autocomplete => "all",
            track        => "experiment$eid",
            label        => "experiment$eid",
            key          => ( $e->{restricted} ? '&reg; ' : '' ) . $e->{name},
            type         => $type,
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            region_feature_densities => 1, # enable histograms in store
            style => {
                featureScale => $featureScale,
                histScale    => $histScale,
                labelScale   => $labelScale,
                showLabels   => JSON::true,
                className    => '{type}',
                histCss      => 'background-color:' . getFeatureColor($eid),
                featureCss   => 'background-color:' . getFeatureColor($eid)
            },

            histograms => {
            	store => "JBrowse/Store/SeqFeature/REST"
            },

            # CoGe-specific stuff
            coge => {
                id      => $eid,
                type    => 'experiment',
                classes => [
                    'coge-tracklist-indented',
                    'coge-tracklist-deletable',
                    'coge-tracklist-info'
                ],
                collapsed   => 1, #FIXME move into CSS
                name        => $e->{name},
                description => $e->{description},
                notebooks   => ( @notebooks ? \@notebooks : undef ),
                annotations => ( @annotations ? \@annotations : undef ),
                onClick     => "ExperimentView.pl?embed=1&eid=$eid",
                menuOptions => [
                    {
                        label => 'ExperimentView',
                        action => "function() { window.open( 'ExperimentView.pl?eid=$eid' ); }"
                    }
                ]
            }
        };
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
            type         => "CoGe/View/Track/Wiggle/MultiXYPlot",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            style        => { featureScale => 0.001 },

            # CoGe-specific stuff
            coge => {
                id          => 0, # use id of 0 to represent all experiments
                type        => 'notebook',
                classes     => ['coge-tracklist-collapsible'],
                collapsed   => 1, #FIXME move into CSS
                name        => 'All Experiments',
                description => '',
                count       => keys %experiments,
                #experiments => [ values %experiments ],
            }
        };
    }

    #
    # Add notebook tracks
    #
    foreach my $n ( sort { $a->{name} cmp $b->{name} } values %notebooks ) {
        my $nid = $n->{list_id};
        push @tracks, {
            key     => ( $n->{restricted} ? '&reg; ' : '' ) . $n->{name},
            baseUrl => "$JBROWSE_API/experiment/notebook/$nid/",
            autocomplete => "all",
            track        => "notebook$nid",
            label        => "notebook$nid",
            type         => "CoGe/View/Track/Wiggle/MultiXYPlot",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            style        => { featureScale => 0.001 },

            # CoGe-specific stuff
            showHoverScores => 1,
            coge            => {
                id      => $nid,
                type    => 'notebook',
                classes => [
                    'coge-tracklist-collapsible',
                    'coge-tracklist-deletable',
                    'coge-tracklist-info'
                ],
                collapsed   => 1, #FIXME move into CSS
                name        => $n->{name},
                description => $n->{description},
                #editable    => $user->is_admin || $user->is_owner_editor( list => $n ) || undef, # mdb removed 2/6/15
                editable    => $user->is_admin || (defined $connectors && ($connectors->{1}{$nid} == 2 || $connectors->{1}{$nid} == 3)) || undef, # mdb added 2/6/15 #TODO move this obscure code into an API
                experiments => ( @{ $expByNotebook{$nid} } ? $expByNotebook{$nid} : undef ),
                count       => scalar @{ $expByNotebook{$nid} },
                onClick     => "NotebookView.pl?embed=1&lid=$nid",
                menuOptions => [
                    {
                        label => 'NotebookView',
                        action => "function() { window.open( 'NotebookView.pl?lid=$nid' ); }"
                    }
                ]
            }
        };
    }

    print STDERR 'time4: ' . (time - $start_time) . "\n" if $DEBUG_PERFORMANCE;

    return encode_json(
        {
            formatVersion => 1,
            dataset_id    => 'coge',
            plugins       => ['CoGe'],
            trackSelector => {
                type => 'CoGe/View/TrackList/CoGe',
            },
            tracks => \@tracks,
        }
    );
}

# FIXME this is duplicated in JBrowse MultiXYPlot.js, need to either make totally client side or totally server side
sub getFeatureColor {
    my $id = shift;
    return '#'
      . sprintf( "%06X",
        ( ( ( ( $id * 1234321 ) % 0x1000000 ) | 0x444444 ) & 0xe7e7e7 ) );
}

1;
