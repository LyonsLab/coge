package CoGe::Services::API::JBrowse::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;

use CoGeX;
use CoGe::Accessory::histogram;
use CoGe::Accessory::IRODS qw( irods_iput );
use CoGe::Accessory::Utils qw( to_number );
use CoGe::Services::Auth qw( init );
use CoGe::Core::Experiment qw( get_data );
use CoGe::Core::Storage qw( $DATA_TYPE_QUANT $DATA_TYPE_POLY $DATA_TYPE_ALIGN $DATA_TYPE_MARKER get_experiment_path get_experiment_metadata );
use Data::Dumper;
use JSON::XS;

my $NUM_QUANT_COL   = 6;
my $NUM_VCF_COL     = 9;
my $NUM_GFF_COL     = 7;
my $MAX_EXPERIMENTS = 20;
my $MAX_WINDOW_SIZE = 1000000;#500000;

my $QUAL_ENCODING_OFFSET = 33;

sub stats_global { # mdb rewritten 8/26/16 COGE-270
    my $self = shift;
    my $eid  = $self->stash('eid');
    my $nid  = $self->stash('nid');
    my $gid  = $self->stash('gid');

    my ($db, $user) = CoGe::Services::Auth::init($self);

    # Get experiments
	my $experiments = _get_experiments($db, $user, $eid, $gid, $nid);

	# Determine min/max values from all experiments' metadata
	my ($globalMax, $globalMin) = (0, 0);
	foreach my $experiment (@$experiments) {
	    # Get range of values for this experiment
	    my ($min, $max);
	    my $md = get_experiment_metadata($experiment->id);
	    if (defined $md && defined $md->{max} && defined $md->{min}) {
	        $min = $md->{min};
	        $max = $md->{max};
	    }
	    else { # legacy default values
	        $min = -1;
	        $max = 1;
	    }

	    # Keep track of widest range
	    if ($max > $globalMax) {
	        $globalMax = $max;
	    }
	    if ($min < $globalMin) {
	        $globalMin = $min;
	    }
	}

    $self->render(json => {
		"scoreMin" => to_number($globalMin),
		"scoreMax" => to_number($globalMax)
	});
}

sub stats_regionFeatureDensities { #FIXME lots of code in common with features()
    my $self     = shift;
    my $eid      = $self->stash('eid');
    my $nid      = $self->stash('nid');
    my $gid      = $self->stash('gid');
    my $chr      = $self->stash('chr');
    my $start    = $self->param('start');
    my $end      = $self->param('end');
    my $bpPerBin = $self->param('basesPerBin');

	# print STDERR "JBrowse::Experiment::stats_regionFeatureDensities eid="
    #  . ( $eid ? $eid : '' ) . " nid="
    #  . ( $nid ? $nid : '' ) . " gid="
    #  . ( $gid ? $gid : '' )
    #  . " $chr:$start:$end ("
    #  . ( $end - $start + 1 )
    #  . ") bpPerBin=$bpPerBin\n";

    my ($db, $user) = CoGe::Services::Auth::init($self);

	my $experiments = _get_experiments($db, $user, $eid, $gid, $nid);

	# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my @bins;
    foreach my $exp (@$experiments) {
        my $data_type = $exp->data_type;
		if ( !$data_type or $data_type == $DATA_TYPE_QUANT ) {
			next;
		}
        my $pData = CoGe::Core::Experiment::get_data(
            eid   => $eid,
            data_type  => $exp->data_type,
            chr   => $chr,
            start => $start,
            end   => $end
        );
        if ( $data_type == $DATA_TYPE_POLY || $data_type == $DATA_TYPE_MARKER ) {
            # Bin and convert to JSON
            foreach my $d (@$pData) {
                my ($s, $e) = ($d->{start}, $d->{stop});
                my $mid = int( ( $s + $e ) / 2 );
                my $b   = int( ( $mid - $start ) / $bpPerBin );
                $b = 0 if ( $b < 0 );
                $bins[$b]++;
            }
        }
        elsif ( $data_type == $DATA_TYPE_ALIGN ) {
	        # Convert SAMTools output into JSON
	        foreach (@$pData) {
	        	chomp;
	        	my (undef, undef, undef, $pos, undef, undef, undef, undef, undef, $seq, undef, undef) = split(/\t/);

	        	my $len = length($seq);
	        	my $s = $pos;
	        	my $e = $pos + $len;

				# Bin results
                my $mid = int( ( $s + $e ) / 2 );
                my $b   = int( ( $mid - $start ) / $bpPerBin );
                $b = 0 if ( $b < 0 );
                $bins[$b]++;
	        }
        }
        else {
        	print STDERR "JBrowse::Experiment::stats_regionFeatureDensities unknown data type for experiment $eid\n";
        	next;
        }
    }

    my ( $sum, $count );
    foreach my $x (@bins) {
        $sum += $x;
         #$max = $x if ( not defined $max or $x > $max ); # mdb removed 1/14/14 for BAM histograms
        $count++;
    }
    my $mean = 0;
    $mean = $sum / $count if ($count);
    my $max = $bpPerBin; # mdb added 1/14/14 for BAM histograms

    $self->render(json => {
        bins => \@bins,
        stats => { basesPerBin => $bpPerBin, max => $max, mean => $mean }
    });
}

sub _get_experiments {
	my $db = shift;
	my $user = shift;
	my $eid = shift;
	my $gid = shift;
	my $nid = shift;

    my @_experiments;
    if ($eid) {
        foreach (split(/,/, $eid)) {
            my $experiment = $db->resultset('Experiment')->find($_);
            push @_experiments, $experiment if $experiment;
        }
    }
    elsif ($nid) {
        my $notebook = $db->resultset('List')->find($nid);
        push @_experiments, $notebook->experiments if $notebook;
    }
    elsif ($gid) {
        my $genome = $db->resultset('Genome')->find($gid);
        push @_experiments, $genome->experiments if $genome;
    }

    # Filter experiments based on permissions
    my @experiments;
    foreach my $e (@_experiments) {
        if (!$e->restricted() || ($user && $user->has_access_to_experiment($e))) {
            if (!$nid || $e->genome_id == $gid) {
            	push @experiments, $e;
            }
        }
        else {
            warn 'JBrowse::Experiment::_get_experiments access denied to experiment ' . $e->id . ( ($user && $user->name) ? ' for user ' . $user->name : '');
        }
    }
    splice(@experiments, $MAX_EXPERIMENTS);
    return \@experiments;
}

sub histogram {
    my $self = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->stash('chr');
    my $storage_path = get_experiment_path($eid);
    my $hist_file = "$storage_path/value1_$chr.hist";
    if (!-e $hist_file) {
        $chr = undef if $chr eq 'Any';
		my $result = CoGe::Core::Experiment::query_data(
			eid => $eid,
			col => 'value1',
			chr => $chr
		);
	    my $bins = CoGe::Accessory::histogram::_histogram_bins($result, 20);
	    my $counts = CoGe::Accessory::histogram::_histogram_frequency($result, $bins);
	    if (open my $fh, ">", $hist_file) {
    	    print {$fh} encode_json({
    	    	first => 0 + $bins->[0][0],
    	    	gap => $bins->[0][1] - $bins->[0][0],
    	    	counts => $counts
    	    });
    	    close $fh;
        }
        else {
            warn "error opening $hist_file";
        }
    }
    open(my $fh, $hist_file);
    my $hist = <$fh>;
    close($fh);
	$self->render(json => decode_json($hist));
}

sub features {
    my $self  = shift;
    my $eid   = $self->stash('eid');
    my $nid   = $self->stash('nid');
    my $gid   = $self->stash('gid');
    my $chr   = $self->stash('chr');
    my $start = $self->param('start');
    my $end   = $self->param('end');
    return unless (($eid or $nid or $gid) and defined $chr and defined $start and defined $end);

    my ($db, $user) = CoGe::Services::Auth::init($self);

	my $experiments = _get_experiments($db, $user, $eid, $gid, $nid);

    if (!@$experiments) {
        return $self->render(json => { features => [] });
    }

	# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my @results;
    foreach my $exp (@$experiments) {
        my $eid       = $exp->id;
        my $data_type = $exp->data_type;

        my $pData = CoGe::Core::Experiment::get_data(
            eid   => $eid,
            data_type  => $data_type,
            chr   => $chr,
            start => $start,
            end   => $end
        );

        if ( !$data_type || $data_type == $DATA_TYPE_QUANT ) {
            foreach my $d (@$pData) {
                #next if ($d->{value1} == 0 and $d->{value2} == 0); # mdb removed 1/15/15 for user Un-Sa - send zero values (for when "Show background" is enabled in JBrowse)
                my %result = (
                    id     => int($eid),
                    start  => $d->{start},
                    end    => $d->{stop},
                );
                $d->{strand} = -1 if ($d->{strand} == 0);
                $result{end} = $result{start} + 1 if ( $result{end} == $result{start} ); #FIXME revisit this
                $result{score} = $d->{strand} * $d->{value1};
                $result{score2} = $d->{value2} if (defined $d->{value2});
                $result{label} = $d->{label} if (defined $d->{label});
                push(@results, \%result);
            }
        }
        elsif ( $data_type == $DATA_TYPE_POLY ) {
            foreach my $d (@$pData) {
                my %result = (
                    id    => $eid,
                    name  => $d->{name},
                    type  => $d->{type},
                    start => $d->{start},
                    end   => $d->{stop},
                    ref   => $d->{ref},
                    alt   => $d->{alt},
                    score => $d->{qual},
                    info  => $d->{info}
                );
                $result{score2} = $d->{value2} if (defined $d->{value2});
                $result{label} = $d->{label} if (defined $d->{label});
                $result{name} = (( $d->{id} && $d->{id} ne '.' ) ? $d->{id} . ' ' : '')
                    . $d->{type} . ' ' . $d->{ref} . " > " . $d->{alt};
                $result{type} = $d->{type} . $d->{ref} . 'to' . $d->{alt}
                    if ( lc($d->{type}) eq 'snp' );
                push(@results, \%result);
            }
        }
        elsif ( $data_type == $DATA_TYPE_MARKER ) {
            foreach my $d (@$pData) {
                warn $d->{attr};
                my ($name) = $d->{attr} =~ /ID\=([\w\.\+\-\,\:\%]+)\;/; # mdb added 4/24/14 for Amanda
                my %result = (
                    uniqueID => $d->{start} . '_' . $d->{stop},
                    type     => $d->{type},
                    start    => $d->{start},
                    end      => $d->{stop},
                    strand   => $d->{strand},
                    score    => $d->{value1},
                    attr     => $d->{attr},
                    name     => $name
                );
                $result{'strand'} = -1 if ($result{'strand'} == 0);
                $result{'end'} = $result{'start'} + 1 if ( $result{'end'} == $result{'start'} ); #FIXME revisit this
                $result{'score'} = $result{'strand'} * $result{'score'} if $result{'score'};
                push(@results, \%result);
            }
        }
        elsif ( $data_type == $DATA_TYPE_ALIGN ) {
	        # Convert SAMTools output into JSON
	        if ($nid) { # return read count if being displayed in a notebook
	           # Bin the read counts by position
                my %bins;
                foreach (@$pData) {
                    chomp;
                    my (undef, undef, undef, $pos, undef, undef, undef, undef, undef, $seq) = split(/\t/);
                    for (my $i = $pos; $i < $pos+length($seq); $i++) {
                        $bins{$i}++;
                    }
                }
                # Convert to intervals
                my ($start, $stop, $lastCount);
                foreach my $pos (sort {$a <=> $b} keys %bins) {
                    my $count = $bins{$pos};
                    if (!defined $start || ($pos-$stop > 1) || $lastCount != $count) {
                        if (defined $start) {
                            $stop++;
                            push(@results, {
                                "id"    => $eid,
                                "start" => $start,
                                "end"   => $stop,
                                "score" => $lastCount
                            });
                        }
                        $start = $pos;
                    }
                    $lastCount = $count;
                    $stop = $pos;
                }
                if (defined $start and $start != $stop) { # do last interval
                    $stop++;
                    push(@results, {
                        "id"    => $eid,
                        "start" => $start,
                        "end"   => $stop,
                        "score" => $lastCount
                    });
                }
	        }
	        else { # else return list reads with qual and seq
    	        foreach (@$pData) {
    	        	chomp;
                    my @fields = split(/\t/);
    	        	my $qname = $fields[0];
                    my $flag = $fields[1];
                    my $pos = $fields[3];
                    my $mapq = $fields[4];
                    my $cigar = $fields[5];
                    my $seq = $fields[9];
                    my $qual = $fields[10];
                    my $md;
                    for (my $i=11; $i<@fields; $i++) {
                        if (substr($fields[$i], 0, 3) eq 'MD:') {
                            $md = substr($fields[$i], 5);
                            last;
                        }
                    }

    	        	my $len = length($seq);
    	        	my $start = $pos;
    	        	my $end = $pos + $len;
    	        	my $strand = ($flag & 0x10 ? '-1' : '1');
    	        	#TODO reverse complement sequence if neg strand?
    	        	$qual = join(' ', map { $_ - $QUAL_ENCODING_OFFSET } unpack("C*", $qual));

                    push(@results, {
                        uniqueID => $qname,
                        "name"     => $qname,
                        "start"    => $start,
                        "end"      => $end,
                        "strand"   => $strand,
                        "score"    => $mapq,
                        "seq"      => $seq,
                        "qual"     => $qual,
                        "Seq length" => $len, #is this used?
                        "cigar"    => $cigar,
                        "md"       => $md
                    });
    	        }
	        }
        }
        else {
        	print STDERR "JBrowse::Experiment::features unknown data type $data_type for experiment $eid\n";
        	next;
        }
    }

    $self->render(json => { "features" => \@results });
}

1;
