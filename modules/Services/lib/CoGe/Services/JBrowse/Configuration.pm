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
        'refseq_config' => 'refseq_config',
        'track_config'  => 'track_config',
    );
    $self->start_mode('track_config');
    $self->mode_param('rm');
}

sub refseq_config {
    my $self           = shift;
    my $gid            = $self->query->param('gid');
    my $SEQ_CHUNK_SIZE = 20000;
    print STDERR "Configuration::refseq_config gid=$gid\n";

    # Load config file
    my $P      = CoGe::Accessory::Web::get_defaults($coge_conf);
    my $DBNAME = $P->{DBNAME};
    my $DBHOST = $P->{DBHOST};
    my $DBPORT = $P->{DBPORT};
    my $DBUSER = $P->{DBUSER};
    my $DBPASS = $P->{DBPASS};

    # Connect to the database
    my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
    my $coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

    #$coge->storage->debugobj(new DBIxProfiler());
    #$coge->storage->debug(1);

    # Get genome
    my $genome = $coge->resultset('Genome')->find($gid);
    return unless $genome;

    my @chromosomes;
    foreach my $chr ( sort { $b->sequence_length <=> $a->sequence_length }
        $genome->genomic_sequences )
    {
        push @chromosomes,
          {
            name         => $chr->chromosome,
            length       => $chr->sequence_length,
            seqChunkSize => $SEQ_CHUNK_SIZE,
            start        => 0,
            end          => $chr->sequence_length - 1
          };
    }

    return encode_json( \@chromosomes );
}

sub track_config {
    my $self       = shift;
    my $gid        = $self->query->param('gid');
    my $cas_ticket = $self->query->param('ticket');
    print STDERR "Configuration::track_config gid=$gid\n";

    # Get config
    my $P      = CoGe::Accessory::Web::get_defaults($coge_conf);
    my $DBNAME = $P->{DBNAME};
    my $DBHOST = $P->{DBHOST};
    my $DBPORT = $P->{DBPORT};
    my $DBUSER = $P->{DBUSER};
    my $DBPASS = $P->{DBPASS};

    # Connect to the database
    my $connstr = "dbi:mysql:dbname=$DBNAME;host=$DBHOST;port=$DBPORT";
    my $coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

    #$coge->storage->debugobj(new DBIxProfiler());
    #$coge->storage->debug(1);

    # Get user
    my $COOKIE_NAME = $P->{COOKIE_NAME};
    my $USER        = undef;
    ($USER) = CoGe::Accessory::Web->login_cas(
        ticket   => $cas_ticket,
        coge     => $coge,
        this_url => $FORM->url()
    ) if ($cas_ticket);
    ($USER) = CoGe::Accessory::LogUser->get_user(
        cookie_name => $COOKIE_NAME,
        coge        => $coge
    ) unless $USER;

    #	my $start_time = time;

    # Get genome
    my $genome = $coge->resultset('Genome')->find($gid);
    return unless $genome;

    my @tracks;

    # Add reference sequence track
    push @tracks, {
        chunkSize     => 20000,
        baseUrl       => "services/JBrowse/service.pl/sequence/$gid/",
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

    # Add GC content track
    push @tracks, {

        #baseUrl => "http://geco.iplantcollaborative.org/rchasman/gc/$gid/",
        baseUrl    => "services/JBrowse/track/gc/$gid/",
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

    #	print STDERR 'time1: ' . (time - $start_time) . "\n";

    # Add feature tracks
    my @feat_type_names =
      grep( !/chromosome/, $genome->distinct_feature_type_names );
    if (@feat_type_names) {

        # Add main "Features" track
        push @tracks, {
            baseUrl      => "services/JBrowse/track/annotation/$gid/",
            autocomplete => "all",
            track        => "features",
            label        => "features",
            key          => "Features",
            type         => "CoGe/View/Track/CoGeFeatures",
            description  => "note, description",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            onClick =>
              "FeatAnno.pl?dsg=$gid;chr={chr};start={start};stop={end}",
            maxFeatureScreenDensity => 50,
            maxHeight               => 100000,
            minSubfeatureWidth      => 4,
            style                   => {
                labelScale               => 0.02,
                arrowheadClass           => "arrowhead",
                className                => "generic_parent",
                maxDescriptionLength     => 70,
                showLabels               => 'true',
                centerChildrenVertically => 'true',
                subfeatureClasses        => { match_part => "match_part7" }
            },

            # CoGe-specific stuff
            coge => {
                id        => 0,
                type      => 'features',
                classes   => ['coge-tracklist-collapsible'],
                collapsed => 1                              #FIXME move into CSS
            }
        };

        # Add a track for each feature type
        foreach my $type_name ( sort @feat_type_names ) {
            push @tracks, {
                baseUrl =>
                  "services/JBrowse/track/annotation/$gid/types/$type_name/",
                autocomplete => "all",
                track        => "features$type_name",
                label        => "features$type_name",
                key          => $type_name,
                type         => "JBrowse/View/Track/HTMLFeatures",
                storeClass   => "JBrowse/Store/SeqFeature/REST",
                onClick =>
"FeatAnno.pl?dsg=$gid;chr={chr};start={start};stop={end};type=$type_name",
                maxFeatureScreenDensity => 1000,     #50,
                maxHeight               => 100000,
                style                   => {
                    arrowheadClass => "arrowhead",
                    className      => "generic_parent",
                    histScale      => 0.05,

#_defaultHistScale        => 4,   #FIXME: not supposed to be set by user?  HTMLFeatures.js
#_defaultLabelScale       => 30,  #FIXME: not supposed to be set by user?  HTMLFeatures.js
#_defaultDescriptionScale => 120, #FIXME: not supposed to be set by user?  HTMLFeatures.js
                    minSubfeatureWidth       => 6,
                    maxDescriptionLength     => 70,
                    showLabels               => 'true',
                    description              => "note, description",
                    centerChildrenVertically => 'true',
                    subfeatureClasses        => { match_part => "match_part7" }
                },

                # CoGe-specific stuff
                coge => {
                    id      => "features_$type_name",
                    type    => 'features',
                    classes => ['coge-tracklist-indented'],
                    collapsed => 1    #FIXME move into CSS
                }
            };
        }
    }

    #	print STDERR 'time2: ' . (time - $start_time) . "\n";

    # Add experiment tracks
    my %all_notebooks;
    my %all_experiments;
    my %expByNotebook;
    foreach $e ( sort experimentcmp $genome->experiments ) {
        next if ( $e->deleted );
        next if ( $e->restricted and not $USER->has_access_to_experiment($e) );
        my $eid = $e->id;
        $all_experiments{$eid} = $e;

        # Make a list of notebook id's
        my @notebooks;
        foreach my $n ( $e->notebooks ) {
            push @notebooks, $n->id;
            $all_notebooks{ $n->id } = $n;
            push @{ $expByNotebook{ $n->id } },
              { id => $e->id, name => $e->name };
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

        my $isSNP = ( $e->data_type == 2 );    #FIXME hardcoded data_type
        push @tracks, {
            baseUrl      => "services/JBrowse/service.pl/experiment/$eid/",
            autocomplete => "all",
            track        => "experiment$eid",
            label        => "experiment$eid",
            key          => ( $e->restricted ? '&reg; ' : '' ) . $e->name,
            type         => (
                $isSNP
                ? "JBrowse/View/Track/HTMLVariants"
                : "CoGe/View/Track/Wiggle/MultiXYPlot"
            ),
            storeClass => "JBrowse/Store/SeqFeature/REST",
            region_stats => 1,    # see HTMLFeatures.js
            style        => {
                featureScale => 0.0001,
                histScale    => 0.005,
                labelScale   => 0.05,
                showLabels   => 'true',
                className    => '{type}',
                histCss      => 'background-color:' . getFeatureColor($eid)
            },

            # CoGe-specific stuff
            #onClick         => "ExperimentView.pl?embed=1&eid=$eid",
            showHoverScores => 1,
            coge            => {
                id      => $eid,
                type    => 'experiment',
                classes => [
                    'coge-tracklist-indented', 'coge-tracklist-deletable',
                    'coge-tracklist-info'
                ],
                name        => $e->name,
                description => $e->description,
                notebooks   => ( @notebooks ? \@notebooks : undef ),
                annotations => ( @annotations ? \@annotations : undef ),
                menuOptions => [
                    {
                        label => 'ExperimentView',
                        action =>
"function() { window.open( 'ExperimentView.pl?eid=$eid' ); }"

                          # url => ... will open link in a dialog window
                    }
                ]
            }
        };
    }

    #	print STDERR 'time3: ' . (time - $start_time) . "\n";

    # Create a fake "All Experiments" notebook track
    if ( keys %all_experiments ) {
        push @tracks, {
            key     => 'All Experiments',
            baseUrl => "services/JBrowse/service.pl/experiment/genome/$gid/",
            autocomplete => "all",
            track        => "notebook0",
            label        => "notebook0",
            type         => "CoGe/View/Track/Wiggle/MultiXYPlot",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            style        => { featureScale => 0.001 },

            # CoGe-specific stuff
            showAverage => 0,
            coge        => {
                id   => 0,            # use id of 0 to represent all experiments
                type => 'notebook',
                classes     => ['coge-tracklist-collapsible'],
                name        => 'All Experiments',
                description => '',
                count       => keys %all_experiments,

                #experiments => [ values %all_experiments ],
            }
        };
    }

    # Add notebook tracks
    foreach my $n ( sort { $a->name cmp $b->name } values %all_notebooks ) {
        next if ( $n->restricted and not $USER->has_access_to_list($n) );
        my $nid = $n->id;
        push @tracks, {
            key => ( $n->restricted ? '&reg; ' : '' ) . $n->name,
            baseUrl => "services/JBrowse/service.pl/experiment/notebook/$nid/",
            autocomplete => "all",
            track        => "notebook$nid",
            label        => "notebook$nid",
            type         => "CoGe/View/Track/Wiggle/MultiXYPlot",
            storeClass   => "JBrowse/Store/SeqFeature/REST",
            style        => { featureScale => 0.001 },

            # CoGe-specific stuff
            onClick         => "NotebookView.pl?embed=1&lid=$nid",
            showAverage     => 0,
            showHoverScores => 1,
            coge            => {
                id      => $nid,
                type    => 'notebook',
                classes => [
                    'coge-tracklist-collapsible', 'coge-tracklist-deletable',
                    'coge-tracklist-info'
                ],
                name        => $n->name,
                description => $n->description,
                editable    => $USER->is_admin
                  || $USER->is_owner_editor( list => $n ),
                experiments =>
                  ( @{ $expByNotebook{$nid} } ? $expByNotebook{$nid} : undef ),
                count       => scalar @{ $expByNotebook{$nid} },
                menuOptions => [
                    {
                        label => 'NotebookView',
                        action =>
"function() { window.open( 'NotebookView.pl?lid=$nid' ); }"

                          # url => ... will open link in a dialog window
                    }
                ]
            }
        };
    }

#	print STDERR 'time4: ' . (time - $start_time) . "\n";print STDERR 'time1: ' . (time - $start_time) . "\n";

    return encode_json(
        {
            formatVersion => 1,
            dataset_id    => 'coge',
            plugins       => ['CoGe'],
            trackSelector => {
                type => 'CoGe/View/TrackList/CoGe',

                # plugin-specific stuff
                # <none>
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

# FIXME this comparison routine is duplicated elsewhere
sub experimentcmp {
    no warnings 'uninitialized';    # disable warnings for undef values in sort
    versioncmp( $b->version, $a->version )
      || $a->name cmp $b->name
      || $b->id cmp $a->id;
}

1;
