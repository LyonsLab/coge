package CoGe::Services::Data::JBrowse::Annotation;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Core::Storage qw( get_genome_seq );
use CoGeDBI qw(get_features_by_range);
use URI::Escape qw(uri_unescape);
use List::Util qw( min max );

sub stats_global {
    my $self = shift;
    print STDERR "JBrowse::Annotation::stats_global\n";
    $self->render(json => { 
        "featureDensity" => 0.02,
        "scoreMin" => 0,
        "scoreMax" => 1,
        "featureDensityByRegion" => 50000,
    });
}

sub features {
    my $self = shift;
    my $gid  = $self->stash('gid');
    my $chr  = $self->stash('chr');
    $chr = uri_unescape($chr) if (defined $chr);
    my $start = $self->param('start');
    my $end   = $self->param('end');
    my $scale = $self->param('scale');
    my $basesPerSpan = $self->param('basesPerSpan');
    my $len   = $end - $start;
    print STDERR "JBrowse::Annotation::features gid=$gid chr=$chr start=$start end=$end\n";

    # Check params
    my $null_response = $self->render(json => { "features" => [] });
    if ( $end <= 0 ) {
        return $null_response;
    }

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init; #TODO switch to CoGe::Services::Auth

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Adjust location - note: incoming coordinates are interbase!
    $start = 0 if ( $start < 0 );
    my $chrLen = $genome->get_chromosome_length($chr);
    $start = min( $start, $chrLen );
    $end   = min( $end,   $chrLen );
    if ( $start == $end ) {
        return $null_response;
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        return $null_response;
    }

    # Get features
    my $features = get_features_by_range(dbh => $db->storage->dbh, gid => $gid, chr => $chr, start => $start, end => $end);

    my @results;
    if ($feat_type) { # feature type was specified
        if results:
            response_body["features"].append({"subfeatures" : []})
        my $i = 0;
        my $lastID = 0;
        my $lastEnd = 0;
        my $lastStart = 0;
        my $result;
        foreach my $f (@$features) {
            if ($f->{feature_id} != lastID and lastID != 0) {
                $result = {
                    "chr" => $f->{chromosome},
                    "start" => $lastStart,
                    "end" => $lastEnd,
                    "uniqueID" => $lastID,
                    "name" => $lastName,
                    "strand" => $lastStrand,
                    "type" => $lastType,
                    "subfeatures" : []
                };
                $i += 1;
            }
            elsif (lastID == 0) {
                $result = {
                    "chr" => $f->{chromosome},
                    "start" => $f->{start},
                    "end" => $f->{stop},
                    "uniqueID" => $f->{feature_id},
                    "name" => $f->{name},
                    "strand" => $f->{strand},
                    "type" => $f->{type}
                }
            }
        }

        if ($f->{strand} == $f->{locstrand} && $f->{locstart} >= $f->{start} && $f->{locstart} <= $f->{stop} && $f->{locstop} >= $f->{start} && $f->{locstop} <= $f->{stop}) {
            push @{$result->{subfeatures}}, {
                "chr" => $f->{chromosome}, # mdb added 11/18/13 issue 240 - add chr to FeatAnno.pl onclick url
                "start" => $f->{locstart},
                "end" => $f->{locstop},
                "strand" => $f->{locstrand},
                "type" => $f->{type},
                "name" => $f->{name}
            };
        }

        $lastStrand = $f->{locstrand};
        $lastType = $f->{type};
        $lastName = $f->{name};
        $lastStart = $f->{start};
        $lastEnd = $f->{stop};
        $lastID = $f->{feature_id};
    }
    else: # process all feature types
        # Filter out redundant results due to multiple feature names
        feats = {}
        min_start = None
        max_stop = None
        for row in results:
    if row[2] == row[12] and row[10] == row[11] and row[0] >= row[6] and row[0] <= row[7] and row[1] >= row[6] and row[1] <= row[7]:
                # Find bounds for sequence retrieval later (show_wobble == 1)
                if min_start is None or min_start > row[0]:
                    min_start = row[0]
                if max_stop is None or max_stop < row[1]:
                    max_stop = row[1]
            
                # Hash unique features
                location_id = row[5]
                try:
                    feat = feats[location_id]
                    if row[9] == 1: # is feature name the primary name?
                        feats[location_id] = row
                except KeyError:
                    feats[location_id] = row

        # Calculate wobble GC for CDS features
        wcount = {}
        wtotal = {}
        if show_wobble:
            seq = '';
            rcseq = '';
            for row in feats.values():
                if row[3] == 'CDS' and is_overlapping(int(row[0]), int(row[1]), int(start), int(end)):
                    if not seq:
                        # Get chromosome subsequence using interbase coordinates
                        cookie = get_cookie(environ)
                        seq = fetch_sequence(genome_id, chr_id, int(min_start), int(max_stop), cookie)
                        
                    if (int(row[2]) == 1): # plus strand
                        count, total = calc_wobble(seq, int(row[0])-int(min_start), int(row[1])-int(min_start))
                    else: # minus strand
                        # Reverse complement the seq and coords
                        if not rcseq:
                            rcseq = reverse_complement(seq)
                        rcstart = int(max_stop)-int(row[1])
                        rcstop  = int(max_stop)-int(row[0])
                        count, total = calc_wobble(rcseq, rcstart, rcstop)
                        
                    id = row[8]
                    if id in wcount:
                        wcount[id] += count
                        wtotal[id] += total
                    else:
                        wcount[id] = count
                        wtotal[id] = total
        
        # Build response
        for row in feats.values():
            response_body["features"].append({
                "chr"    : chr_id, # mdb added 11/18/13 issue 240 - add chr to FeatAnno.pl onclick url
                "start"  : row[0],
                "end"    : row[1],
                "strand" : row[2],
                "type"   : row[3],
                "name"   : row[4] if row[3] == 'gene' else None,
                "wobble" : round(wcount[row[8]] / float(wtotal[row[8]]), 3) if (show_wobble and row[8] in wcount) else 0
            })

    $self->render(json => { "features" => \@results });
}

1;
