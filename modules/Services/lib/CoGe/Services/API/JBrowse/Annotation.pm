package CoGe::Services::API::JBrowse::Annotation;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;
use CoGe::Core::Storage qw( get_genome_seq reverse_complement );
use CoGe::Services::Auth;
use CoGeDBI qw(get_features_by_range);
use URI::Escape qw(uri_unescape);
use List::Util qw( min max );

sub stats_global {
    my $self = shift;
    #print STDERR "JBrowse::Annotation::stats_global\n";
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
    my $feat_type    = $self->stash('type');
    my $dsid         = $self->stash('dsid');
    my $start        = $self->param('start');
    my $end          = $self->param('end');
    my $show_wobble  = $self->param('showWobble');
    my $scale        = $self->param('scale');
    my $basesPerSpan = $self->param('basesPerSpan');
    my $len = $end - $start;
#    print STDERR "JBrowse::Annotation::features gid=$gid ", ($dsid ? "dsid=$dsid " : ''), ($feat_type ? "type=$feat_type " : ''), "chr=$chr start=$start end=$end\n";

    # Check params
    if ( $end <= 0 ) {
        warn 'end <= 0';
        $self->render(json => { "features" => [] });
        return;
    }

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Retrieve genome
    my $genome = $db->resultset('Genome')->find($gid);
    return unless $genome;

    # Adjust location - note: incoming coordinates are interbase!
    $start = 0 if ( $start < 0 );
    my $chrLen = $genome->get_chromosome_length($chr);
    $start = min( $start, $chrLen );
    $end   = min( $end,   $chrLen );
    if ( $start == $end ) {
        warn 'start == end';
        $self->render(json => { "features" => [] });
        return;
    }

    # Check permissions
    if ( $genome->restricted
        and ( not defined $user or not $user->has_access_to_genome($genome) ) )
    {
        warn 'user does not have access';
        $self->render(json => { "features" => [] });
        return;
    }

    # Get features
    my $features = get_features_by_range(
        dbh => $db->storage->dbh, 
        gid => $gid, 
        dsid => $dsid,
        chr => $chr, 
        start => $start, 
        end => $end,
        feat_type => $feat_type
    );

    my @results;
    if ($feat_type) { # feature type was specified
        my $i = 0;
        my $lastID = 0;
        my $lastEnd = 0;
        my $lastStart = 0;
        my ($lastName, $lastStrand, $lastType);
        my $result;
        foreach my $f (@$features) {
            if ($f->{feature_id} != $lastID && $lastID != 0) {
                $result = {
                    "chr" => $f->{chromosome},
                    "start" => $lastStart,
                    "end" => $lastEnd,
                    "uniqueID" => $lastID,
                    "name" => $lastName,
                    "strand" => $lastStrand,
                    "type" => $lastType,
                    "subfeatures" => []
                };
                $i += 1;
            }
            elsif ($lastID == 0) {
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

            if ($f->{strand} == $f->{locstrand} && 
                $f->{locstart} >= $f->{start} && $f->{locstart} <= $f->{stop} && 
                $f->{locstop} >= $f->{start} && $f->{locstop} <= $f->{stop}) 
            {
                push @{$result->{subfeatures}}, {
                    "chr" => $f->{chromosome}, # mdb added 11/18/13 issue 240 - add chr to FeatAnno.pl onclick url
                    "start" => $f->{locstart},
                    "end" => $f->{locstop},
                    "strand" => $f->{locstrand},
                    "type" => $f->{type},
                    "name" => $f->{name}
                };
            }
    
            push @results, $result;
    
            $lastStrand = $f->{locstrand};
            $lastType = $f->{type};
            $lastName = $f->{name};
            $lastStart = $f->{start};
            $lastEnd = $f->{stop};
            $lastID = $f->{feature_id};
        }
    }
    else { # process all feature types
        # Filter out redundant results due to multiple feature names
        my $feats = {};
        my $min_start;
        my $max_stop;
        foreach my $f (@$features) {
            if ($f->{locstrand} == $f->{strand} && 
                $f->{locstart} >= $f->{start} && $f->{locstart} <= $f->{stop} && 
                $f->{locstop} >= $f->{start} && $f->{locstop} <= $f->{stop}) 
            {
                # Find bounds for sequence retrieval later (show_wobble == 1)
                if (!defined($min_start) || $min_start > $f->{locstart}) {
                    $min_start = $f->{locstart};
                }
                if (!defined($max_stop) || $max_stop < $f->{locstop}) {
                    $max_stop = $f->{locstop};
                }
            
                # Hash unique features
                my $location_id = $f->{location_id};
                if (defined $feats->{$location_id}) {
                    my $feat = $feats->{$location_id};
                    if ($f->{primary_name} == 1) { # is feature name the primary name?
                        $feats->{$location_id} = $f;
                    }
                }
                else {
                    $feats->{$location_id} = $f;
                }
            }
        }

        # Calculate wobble GC for CDS features
        my (%wcount, %wtotal);
        if ($show_wobble) {
            my $seq = '';
            my $rcseq = '';
            foreach my $f (values %$feats) {
                if ($f->{type} eq 'CDS' && $f->{locstart} <= $end && $f->{locstop} <= $start) {
                    if (!$seq) {
                        # Get chromosome subsequence using interbase coordinates
                        $seq = get_genome_seq(
                            gid   => $gid,
                            chr   => $chr,
                            start => $min_start,
                            stop  => $max_stop
                        );
                    }
                    
                    my ($count, $total);
                    if ($f->{locstrand} eq '1') { # plus strand
                        ($count, $total) = _calc_wobble($seq, $f->{locstart} - $min_start, $f->{locstop} - $min_start);
                    }
                    else { # minus strand
                        # Reverse complement the seq and coords
                        if (!$rcseq) {
                            $rcseq = reverse_complement($seq);
                        }
                        my $rcstart = $max_stop - $f->{locstop};
                        my $rcstop  = $max_stop - $f->{locstart};
                        ($count, $total) = _calc_wobble($rcseq, $rcstart, $rcstop);
                    }
                    
                    my $id = $f->{feature_id};
                    $wcount{$id} += $count;
                    $wtotal{$id} += $total;
                }
            }
        }
        
        # Build response
        foreach my $f (values %$feats) {
            my $fid = $f->{feature_id};
            my $wobble_gc = ($show_wobble && $wtotal{$fid}) ? ($wcount{$fid} / $wtotal{$fid}) : 0;
            push @results, {
                "chr"    => $chr, # mdb added 11/18/13 issue 240 - add chr to FeatAnno.pl onclick url
                "start"  => $f->{locstart},
                "end"    => $f->{locstop},
                "strand" => $f->{locstrand},
                "type"   => $f->{type},
                "name"   => ($f->{type} eq 'gene' ? $f->{name} : undef),
                "wobble" => $wobble_gc
            };
        }
    }

    $self->render(json => { "features" => \@results });
}

# mdb added 11/4/13 issue 246 - add wobble shading
sub _calc_wobble {
    my ($seq, $start_offset, $stop_offset) = @_;
    #print STDERR "_calc_wobble: $start_offset-$stop_offset, $seq\n";
    
    if ($start_offset < 0) {
        $start_offset = 0;
    }
    if ($stop_offset < 0 or $start_offset == $stop_offset) {
        return 0;
    }

    # Calculate GC percent in wobble position
    my ($count, $total) = (0, 0);
    for (my $i = $start_offset+2;  $i < $stop_offset; $i += 3) { #for i in xrange(start_offset+2, stop_offset, 3):
        $total += 1;
        my $c = uc(substr($seq, $i, 1));
        if ($c eq 'G' or $c eq 'C') {
            $count += 1;
        }
    }

    return ($count, $total);
}

1;
