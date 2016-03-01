package CoGe::Services::JBrowse::Experiment;
use base 'CGI::Application';

use CoGe::Accessory::histogram;
use CoGe::Accessory::Web;
use CoGe::Core::Experiment;
use CoGe::Core::Storage qw( get_experiment_path );
use CoGeX;
use Data::Dumper;
use JSON::XS;

#TODO: use these from Storage.pm instead of redeclaring them
my $DATA_TYPE_QUANT  = 1; # Quantitative data
my $DATA_TYPE_POLY	 = 2; # Polymorphism data
my $DATA_TYPE_ALIGN  = 3; # Alignments
my $DATA_TYPE_MARKER = 4; # Markers

my $NUM_QUANT_COL   = 6;
my $NUM_VCF_COL     = 9;
my $NUM_GFF_COL     = 7;
my $MAX_EXPERIMENTS = 20;
my $MAX_WINDOW_SIZE = 1000000;#500000;

my $QUAL_ENCODING_OFFSET = 33;

sub setup {
    my $self = shift;
    $self->run_modes(
        'stats_global' => 'stats_global',
        'stats_regionFeatureDensities' => 'stats_regionFeatureDensities',
        'features' => 'features',
        'query_data' => 'query_data',
        'histogram' => 'histogram',
        'snps'  => 'snps',
    );
    $self->mode_param('rm');
}

sub stats_global {
	print STDERR "JBrowse::Experiment::stats_global\n";
    return qq{{
		"scoreMin" : -1,
		"scoreMax" : 1
	}};
}

sub stats_regionFeatureDensities { #FIXME lots of code in common with features()
    my $self     = shift;
    my $eid      = $self->param('eid');
    my $nid      = $self->param('nid');
    my $gid      = $self->param('gid');
    my $chr      = $self->param('chr');
    my $start    = $self->query->param('start');
    my $end      = $self->query->param('end');
    my $bpPerBin = $self->query->param('basesPerBin');

#	print STDERR "JBrowse::Experiment::stats_regionFeatureDensities eid="
#      . ( $eid ? $eid : '' ) . " nid="
#      . ( $nid ? $nid : '' ) . " gid="
#      . ( $gid ? $gid : '' )
#      . " $chr:$start:$end ("
#      . ( $end - $start + 1 )
#      . ") bpPerBin=$bpPerBin\n";

    # Connect to the database
    my ( $db, $user ) = CoGe::Accessory::Web->init;

    # Retrieve experiments
    my @all_experiments;
    if ($eid) {
        my $experiment = $db->resultset('Experiment')->find($eid);
        return unless $experiment;
        push @all_experiments, $experiment;
    }
    elsif ($nid) {
        my $notebook = $db->resultset('List')->find($nid);
        return unless $notebook;
        push @all_experiments, $notebook->experiments;
    }
    elsif ($gid) {
        my $genome = $db->resultset('Genome')->find($gid);
        return unless $genome;
        push @all_experiments, $genome->experiments;
    }

    # Filter experiments based on permissions
    my @experiments;
    foreach my $e (@all_experiments) {
        unless ( $user->has_access_to_experiment($e) ) {
        	print STDERR "JBrowse::Experiment::stats_regionFeatureDensities access denied to experiment $eid for user '", $user->name, "'\n";
        	next;
        }
        push @experiments, $e;
    }
    splice( @experiments, $MAX_EXPERIMENTS, @experiments );

	# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my $results = '';
    my @bins;
    foreach my $exp (@experiments) {
        my $data_type = $exp->data_type;
		if ( !$data_type or $data_type == $DATA_TYPE_QUANT ) {
			next;
		}
        elsif ( $data_type == $DATA_TYPE_POLY || $data_type == $DATA_TYPE_MARKER ) {
            my $pData = CoGe::Core::Experiment::get_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

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
	        my $cmdOut = CoGe::Core::Experiment::get_data(
	            eid   => $eid,
	            data_type  => $exp->data_type,
	            chr   => $chr,
	            start => $start,
	            end   => $end
	        );

	        # Convert SAMTools output into JSON
	        foreach (@$cmdOut) {
	        	chomp;
	        	my ($qname, $flag, $rname, $pos, $mapq, $cigar, undef, undef, undef, $seq, $qual, $tags) = split(/\t/);

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

    return encode_json(
        {
            bins => \@bins,
            stats => { basesPerBin => $bpPerBin, max => $max, mean => $mean }
        }
    );
}

sub histogram {
    my $self  = shift;
    my $eid   = $self->param('eid');

    my $storage_path = get_experiment_path($eid);
    my $hist_file = "$storage_path/value1.hist";
    if (!-e $hist_file) {
		my $result = CoGe::Core::Experiment::query_data(
			eid => $eid,
			col => 'value1',
			chr => $chr
		);
	    my $bins = CoGe::Accessory::histogram::_histogram_bins($result, 20);
	    my $counts = CoGe::Accessory::histogram::_histogram_frequency($result, $bins);
	    open my $fh, ">", $hist_file;
	    print {$fh} encode_json({
	    	first => 0 + $bins->[0][0],
	    	gap => $bins->[0][1] - $bins->[0][0],
	    	counts => $counts
	    });
	    close $fh; 
    }
    open my $fh, $hist_file;
    my $hist = <$fh>;
    close $fh;
    return $hist;
}

sub query_data {
    my $self = shift;
    my $eid = $self->param('eid');
    my $chr = $self->query->param('chr');
    my $type = $self->query->param('type');
    my $gte = $self->query->param('gte');
    my $lte = $self->query->param('lte');

    $chr = undef if $chr eq 'All';
	my $result = CoGe::Core::Experiment::query_data(
		eid => $eid,
		col => 'chr,start,stop,value1',
		chr => $chr,
		type => $type,
		gte => $gte,
		lte => $lte,
	);
	return encode_json($result);
}

sub snp_overlaps_feature {
	my $loc = shift;
	my $features = shift;
	foreach my $feature (@$features) {
		if ($loc >= $feature->[1] && $loc <= $feature->[2]) {
			return 1;
		}
	}
	return 0;
}

sub debug {
	my $data = shift;
	my $new_file = shift;
	my $OUTFILE;
	open $OUTFILE, ($new_file ? ">/tmp/sean" : ">>/tmp/sean");
	print {$OUTFILE} Dumper $data;
	print {$OUTFILE} "\n";
	close $OUTFILE;
}

sub _add_features {
    my ($chr, $type_ids, $dsid, $hits, $dbh) = @_;
    my $query = 'SELECT chromosome,start,stop FROM feature WHERE dataset_id=' . $dsid;
    if ($chr) {
    	$query .= " AND chromosome='" . $chr . "'";
    }
    if ($type_ids) {
    	if (index($type_ids, ',') != -1) {
		    $query .= ' AND feature_type_id IN(' . $type_ids . ')';
    	}
    	else {
    		$query .= ' AND feature_type_id=' . $type_ids;
    	}
    }
    my $sth = $dbh->prepare($query);
    $sth->execute();
    while (my @row = $sth->fetchrow_array) {
    	my $chromosome = $hits->{$row[0]};
    	if (!$chromosome) {
    		$hits->{$row[0]} = [\@row];
    	}
    	else {
	    	push @$chromosome, \@row;
    	}
    }
}

sub snps {
    my $self = shift;
    my $eid = $self->param('eid');
    my $chr = $self->param('chr');

    $chr = undef if $chr eq 'Any';
	my $snps = CoGe::Core::Experiment::query_data(
		eid => $eid,
		col => 'chr,start',
		chr => $chr,
	);

    my ( $db, $user ) = CoGe::Accessory::Web->init;
    my $dbh = $db->storage->dbh;

	my $types = $self->query->param('features');
	my $type_ids;
	if ($types ne 'all') {
		$type_ids = join(',', @{$dbh->selectcol_arrayref('SELECT feature_type_id FROM feature_type WHERE name IN(' . $types . ')')});
	}

    my $features = {};
    my $experiment = $db->resultset('Experiment')->find($eid);
    my $ids = $dbh->selectcol_arrayref('SELECT dataset_id FROM dataset_connector WHERE genome_id=' . $experiment->genome_id);
    foreach my $dsid (@$ids) {
        _add_features $chr, $type_ids, $dsid, $features, $dbh;
    }
	my $hits = [];
    foreach my $snp (@$snps) {
    	my @tokens = split ',', $snp;
    	if (snp_overlaps_feature 0 + $tokens[1], $features->{substr $tokens[0], 1, -1}) {
    		push @$hits, $snp;
    	}
    }
    return encode_json($hits);
}

sub features {
    my $self  = shift;
    my $eid   = $self->param('eid');
    my $nid   = $self->param('nid');
    my $gid   = $self->param('gid');
    my $chr   = $self->param('chr');
    my $start = $self->query->param('start');
    my $end   = $self->query->param('end');
    print STDERR "JBrowse::Experiment::features eid="
      . ( $eid ? $eid : '' ) . " nid="
      . ( $nid ? $nid : '' ) . " gid="
      . ( $gid ? $gid : '' )
      . " $chr:$start:$end ("
      . ( $end - $start + 1 ) . ")\n";
    return unless ( ( $eid or $nid or $gid ) and defined $chr and defined $start and defined $end );

# mdb removed 11/6/15 COGE-678
#    if ( $end - $start + 1 > $MAX_WINDOW_SIZE ) {
#        print STDERR "experiment features maxed\n";
#        return qq{{ "features" : [ ] }};
#    }

    # Connect to the database
    my ( $db, $user ) = CoGe::Accessory::Web->init;

    # Retrieve experiments
    my @all_experiments;
    if ($eid) {
        my $experiment = $db->resultset('Experiment')->find($eid);
        return unless $experiment;
        push @all_experiments, $experiment;
    }
    elsif ($nid) {
        my $notebook = $db->resultset('List')->find($nid);
        return unless $notebook;
        push @all_experiments, $notebook->experiments;
    }
    elsif ($gid) {
        my $genome = $db->resultset('Genome')->find($gid);
        return unless $genome;
        push @all_experiments, $genome->experiments;
    }

    # Filter experiments based on permissions
    my @experiments;
    foreach my $e (@all_experiments) {
        unless ( $user->has_access_to_experiment($e) ) {
        	print STDERR "JBrowse::Experiment::features access denied to experiment $eid for user '", $user->name, "'\n";
        	next;
    	}
        push @experiments, $e;
    }
    splice( @experiments, $MAX_EXPERIMENTS, @experiments );

    if (!@experiments) {
        return qq{{ "features" : [ ] }};
    }

	# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my $results = '';
    foreach my $exp (@experiments) {
        my $eid       = $exp->id;
        my $data_type = $exp->data_type;

        if ( !$data_type || $data_type == $DATA_TYPE_QUANT ) {
            my $pData = CoGe::Core::Experiment::get_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

            # Convert to JSON
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
                $results .= ( $results ? ',' : '') . encode_json(\%result);
            }
        }
        elsif ( $data_type == $DATA_TYPE_POLY ) {
            my $pData = CoGe::Core::Experiment::get_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

            # Convert to JSON
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
                $results .= ( $results ? ',' : '') . encode_json(\%result);
            }
        }
        elsif ( $data_type == $DATA_TYPE_MARKER ) {
            my $pData = CoGe::Core::Experiment::get_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

            # Convert to JSON
            foreach my $d (@$pData) {
                my ($name) = $d->{attr} =~ /ID\=([\w\.\-]+)\;/; # mdb added 4/24/14 for Amanda
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
                $result{'stop'} = $result{'start'} + 1 if ( $result{'stop'} == $result{'start'} ); #FIXME revisit this
                $result{'value1'} = $result{'strand'} * $result{'value1'};
                $results .= ( $results ? ',' : '') . encode_json(\%result);
            }
        }
        elsif ( $data_type == $DATA_TYPE_ALIGN ) {
	        my $cmdOut = CoGe::Core::Experiment::get_data(
	            eid   => $eid,
	            data_type => $exp->data_type,
	            chr   => $chr,
	            start => $start,
	            end   => $end
	        );

	        # Convert SAMTools output into JSON
	        if ($nid) { # return read count if being displayed in a notebook
	           # Bin the read counts by position
                my %bins;
                foreach (@$cmdOut) {
                    chomp;
                    my (undef, undef, undef, $pos, $mapq, undef, undef, undef, undef, $seq) = split(/\t/);
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
                            $results .= ( $results ? ',' : '' )
                                . qq{{ "id": $eid, "start": $start, "end": $stop, "score": $lastCount }};
                        }
                        $start = $pos;
                    }
                    $lastCount = $count;
                    $stop = $pos;
               }
               if (defined $start and $start != $stop) { # do last interval
                    $stop++;
                    $results .= ( $results ? ',' : '' )
                        . qq{{ "id": $eid, "start": $start, "end": $stop, "score": $lastCount }};
                }
	        }
	        else { # else return list reads with qual and seq
    	        foreach (@$cmdOut) {
    	        	chomp;
    	        	my ($qname, $flag, $rname, $pos, $mapq, $cigar, undef, undef, undef, $seq, $qual, $tags) = split(/\t/);

    	        	my $len = length($seq);
    	        	my $start = $pos;
    	        	my $end = $pos + $len;
    	        	my $strand = ($flag & 0x10 ? '-1' : '1');
    	        	#TODO reverse complement sequence if neg strand?
    	        	$qual = join(' ', map { $_ - $QUAL_ENCODING_OFFSET } unpack("C*", $qual));

    	        	$results .= ( $results ? ',' : '' )
                      . qq{{ "uniqueID": "$qname", "name": "$qname", "start": $start, "end": $end, "strand": $strand, "score": $mapq, "seq": "$seq", "qual": "$qual", "Seq length": $len, "CIGAR": "$cigar" }};
    	        }
	        }
        }
        else {
        	print STDERR "JBrowse::Experiment::features unknown data type $data_type for experiment $eid\n";
        	next;
        }
    }

    #print STDERR "{ 'features' : [ $results ] }\n";
    return qq{{ "features" : [ $results ] }};
}

1;
