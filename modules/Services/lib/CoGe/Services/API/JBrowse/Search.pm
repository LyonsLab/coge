package CoGe::Services::API::JBrowse::Search;

use Mojo::Base 'Mojolicious::Controller';

use CoGe::Core::Experiment;
use CoGe::Core::Storage qw( $DATA_TYPE_QUANT $DATA_TYPE_POLY $DATA_TYPE_ALIGN $DATA_TYPE_MARKER get_experiment_path get_upload_path );
use CoGe::Services::Auth;
use CoGeDBI qw( get_dataset_ids feature_type_names_to_id );
use File::Path qw( mkpath );
use File::Spec::Functions qw( catdir catfile );
use File::Temp qw( tempfile tempdir );

my $log10 = log(10);

sub _alignments {
    my $self = shift;
    my $all = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->stash('chr');
    $chr = undef if $chr eq 'Any' || $chr eq 'All';
    my $type_names = $self->param('features');

    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
	my $experiment = $db->resultset('Experiment')->find($eid);
	return unless $self->_can_view($experiment, $user);

    my $data1 = _get_experiment_data($experiment->id, $experiment->data_type, $chr, $all);
    my $data2 = _get_db_data($experiment->genome_id, $type_names, $chr, $db->storage->dbh);
    return _in($data1, $data2);
}

sub _can_view {
	my $self = shift;
	my $experiment = shift;
	my $user = shift;

	if ($experiment->restricted() && (!$user || !($user->is_admin() || $user->has_access_to_experiment($experiment)))) {
		$self->render(json => { error => 'User does not have permission to view this experiment' }, status => 401);
		return 0;
	}
	return 1;
}

sub data { #TODO move this out of this module into Core layer (mdb 8/26/16)
    my $self = shift;
    my $id = int($self->stash('eid'));
    my $chr = $self->stash('chr');
    my $type = $self->param('type');
    my $data_type = $self->param('data_type');
    my $gte = $self->param('gte');
    my $lte = $self->param('lte');
    my $transform = $self->param('transform');
    my $filename = $self->param('filename');
    my $irods_path = $self->param('irods_path');
    my $load_id = $self->param('load_id');
    my $gap_max = $self->param('gap_max');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    # Get experiment
    my $experiment = $db->resultset("Experiment")->find($id);
    unless (defined $experiment) {
        $self->render(json => {
            error => { Error => "Experiment not found" }
        });
        return;
    }
	return unless $self->_can_view($experiment, $user);

    my $exp_data_type = $experiment->data_type;
    $exp_data_type = $DATA_TYPE_MARKER if defined $gap_max;
    my $ext = $exp_data_type == $DATA_TYPE_POLY ? '.vcf' : $exp_data_type == $DATA_TYPE_ALIGN ? '.sam' : $exp_data_type == $DATA_TYPE_MARKER ? '.gff' : '.csv';

    my $fh;
    my $tempfile;
    my $path;
    if ($irods_path) {
        ($fh, $tempfile) = tempfile();
    }
    elsif ($load_id) {
        $path = catdir(get_upload_path($user->name, $load_id), 'upload');
        mkpath($path);
        open $fh, ">", catfile($path, 'search_results' . $ext);
    }

    $filename = 'experiment' unless $filename;
    $filename .= $ext;
    if (!$irods_path && !$load_id) {
        $self->res->headers->content_disposition('attachment; filename=' . $filename . ';');
    }
    $self->_write("##gff-version 3\n", $fh) if $exp_data_type == $DATA_TYPE_MARKER;
    my $comment_char = ($exp_data_type == $DATA_TYPE_ALIGN) ? "\@CO\t" : '# ';
    $self->_write($comment_char . 'experiment: ' . $experiment->name . "\n", $fh);
    $self->_write($comment_char . 'chromosome: ' . $chr . "\n", $fh);

    if ( !$exp_data_type || $exp_data_type == $DATA_TYPE_QUANT ) {
        if ($type) {
            $self->_write('# search: type = ' . $type, $fh);
            $self->_write(", gte = $gte", $fh) if $gte;
            $self->_write(", lte = $lte", $fh) if $lte;
            $self->_write("\n", $fh);
        }
        $self->_write('# transform: ' . $transform . "\n", $fh) if $transform;
        my ($columns, $lines) = _get_quant_lines($id, $chr, $data_type, $type, $gte, $lte);
        $self->_write('# columns: ', $fh);
        for (my $i=0; $i<scalar @{$columns}; $i++) {
            $self->_write(',', $fh) if $i;
            $self->_write($columns->[$i], $fh);
        }
        $self->_write("\n", $fh);

        foreach my $line (@{$lines}) {
            my $tokens = _transform_line($line, $transform);
            $self->_write(join(',', @{$tokens}), $fh);
            $self->_write("\n", $fh);
        }
    }
    elsif ( $exp_data_type == $DATA_TYPE_POLY ) {
        my $type_names = $self->param('features');
        if ($type_names) {
            $self->_write('# search: SNPs in ' . $type_names . "\n", $fh);
        }
        my $snp_type = $self->param('snp_type');
        if ($snp_type) {
            $self->_write('# search: ' . $snp_type . "\n", $fh);
        }
        $self->_write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", $fh);
        my $snps = $self->_snps;
        my $s = (index(@{$snps}[0], ', ') == -1 ? 1 : 2);
        foreach (@{$snps}) {
            my @l = split(',');
            $self->_write(substr($l[0], 1, -1), $fh);
            $self->_write("\t", $fh);
            $self->_write($l[1], $fh);
            $self->_write("\t", $fh);
            $self->_write(substr($l[4], $s, -1), $fh);
            $self->_write("\t", $fh);
            $self->_write(substr($l[5], $s, -1), $fh);
            $self->_write("\t", $fh);
            $self->_write(substr($l[6], $s, -1), $fh);
            $self->_write("\t", $fh);
            $self->_write($l[7], $fh);
            $self->_write("\t.\t", $fh);
            $self->_write(substr($l[8], $s, -1), $fh);
            $self->_write("\n", $fh);
        }
    }
    elsif ( $exp_data_type == $DATA_TYPE_ALIGN) {
        my $type_names = $self->param('features');
        if ($type_names) {
            $self->_write("\@CO\tsearch: Alignments in " . $type_names . "\n", $fh);
        }
        $self->_write("\@HD\tVN:1.5\tSO:coordinate\n", $fh);
        if ($chr eq 'All') {
            my $chromosomes = $experiment->genome->chromosomes_all;
            foreach (@{$chromosomes}) {
                $self->_write("\@SQ\tSN:" . $_->{name} . "\tLN:" . $_->{length} . "\n", $fh);
            }
        }
        else {
            $self->_write("\@SQ\tSN:" . $chr . "\tLN:" . $experiment->genome->get_chromosome_length($chr) . "\n", $fh);
        }

        my $alignments = $self->_alignments(1);
        foreach (@{$alignments}) {
            $self->_write($_, $fh);
            $self->_write("\n", $fh);
        }
    }
    elsif ( $exp_data_type == $DATA_TYPE_MARKER) {
        my $type_names = $self->param('features');
        if ($type_names) {
            $self->_write('# search: Markers in ' . $type_names . "\n", $fh);
        }
        $self->_write("#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n", $fh);
        if (defined $gap_max) {
            my ($columns, $lines) = _get_quant_lines($id, $chr, $data_type, $type, $gte, $lte);
            my $c;
            my @chrs;
            my $markers;
            my $start_f;
            my $start_r;
            my $stop_f;
            my $stop_r;
            foreach my $line (@{$lines}) {
                my $tokens = _transform_line($line, $transform);
                if (!$c) {
                    warn Dumper $tokens;
                    $c = $tokens->[0];
                    $markers = [];
                    push @chrs, [$c, $markers];
                }
                if ($c ne $tokens->[0]) {
                    if ($start_f) {
                        push @$markers, [$start_f, $stop_f, '+'];
                        $start_f = $stop_f = undef;
                    }
                    if ($start_r) {
                        push @$markers, [$start_r, $stop_r, '-'];
                        $start_r = $stop_r = undef;
                    }
                    $c = $tokens->[0];
                    $markers = [];
                    push @chrs, [$c, $markers];
                }
                if ($tokens->[3] eq '1') {
                    if (!$start_f) {
                        $start_f = $tokens->[1];
                        $stop_f = $tokens->[2];
                    }
                    elsif ($tokens->[1] > $stop_f + $gap_max) {
                        push @$markers, [$start_f, $stop_f, '+'];
                        $start_f = $tokens->[1];
                        $stop_f = $tokens->[2];
                    }
                    else {
                        $stop_f = $tokens->[2];
                    }
                }
                else {
                    if (!$start_r) {
                        $start_r = $tokens->[1];
                        $stop_r = $tokens->[2];
                    }
                    elsif ($tokens->[1] > $stop_r + $gap_max) {
                        push @$markers, [$start_r, $stop_r, '-'];
                        $start_r = $tokens->[1];
                        $stop_r = $tokens->[2];
                    }
                    else {
                        $stop_r = $tokens->[2];
                    }
                }
            }
            if ($start_f) {
                push @$markers, [$start_f, $stop_f, '+'];
            }
            if ($start_r) {
                push @$markers, [$start_r, $stop_r, '-'];
            }
            foreach (@chrs) {
                $c = $_->[0];
                my @m = sort { $a->[0] <=> $b->[0] } @{$_->[1]};
                foreach my $l (@m) {
                    my $start = $l->[0];
                    my $stop = $l->[1];
                    $self->_write_marker($c, undef, undef, $start, $stop, 0, $l->[2], undef, undef, $fh);
                }
            }
        }
        else {
            my $markers = $self->_markers;
            my $s = (index(@{$markers}[0], ', ') == -1 ? 1 : 2);
            $chr = undef if $chr eq 'All';
            foreach (@{$markers}) {
                my @l = split(',');
                my $c = substr($l[0], 1, -1);
                next if $chr && $chr ne $c;
                $self->_write_marker($c, undef, substr($l[4], $s, -1), $l[1], $l[2], $l[5], substr($l[3], $s, -1) == 1 ? '+' : '-', undef, substr($l[6], $s, -1), $fh);
            }
        }
    }
    $self->finish();
    if ($fh) {
        close $fh;
        if ($irods_path) {
            irods_iput($tempfile, $irods_path . '/' . $filename);
        }
        elsif ($load_id && $exp_data_type == $DATA_TYPE_ALIGN) {
                my $cmd = $conf->{SAMTOOLS} || 'samtools';
                $cmd .= ' view -bS ' . catfile($path, 'search_results.sam') . ' > ' . catfile($path, 'search_results.bam');
                system($cmd);
        }
    }
}

sub _get_data {
	my ($self, $num, $db) = @_;
	my $type = $self->param('type' . $num);
	if ($type eq 'experiment') {
		my $experiment = $db->resultset('Experiment')->find($self->param('eid' . $num));
		return _get_experiment_data($experiment->id, $experiment->data_type, $self->param('chr'), 0);
	}
	if ($type eq 'features') {
		return _get_db_data($self->param('gid' . $num), $self->param('features' . $num), $self->param('chr'), $db->storage->dbh);
	}
}

# returns a psuedo object with coderefs next() and line()
sub _get_db_data {
	my ($gid, $type_names, $chr, $dbh) = @_;

	my $query = 'SELECT start,stop';
	$query .= ',chromosome' unless $chr;
	$query .= ' FROM feature WHERE dataset_id';
	my $ids = get_dataset_ids($gid, $dbh);
	$query .= (scalar @$ids == 1) ? '=' . $ids->[0] : ' IN(' . join(',', @$ids) . ')';
	$query .= " AND chromosome='" . $chr . "'" if $chr;
	if ($type_names && $type_names ne 'all') {
		my $type_ids = feature_type_names_to_id($type_names, $dbh);
		$query .= ' AND feature_type_id';
		$query .= (index($type_ids, ',') == -1) ? '=' . $type_ids : ' IN(' . $type_ids . ')';
	} else {
        $query .= ' AND feature_type_id!=4'; # ignore chromosomes
    }
	$query .= ' ORDER BY ';
	$query .= 'chromosome,' unless $chr;
	$query .= 'start,stop';
	my $sth = $dbh->prepare($query);
	$sth->execute();
	return {
		sth => $sth,
		line => sub {
				my $self = shift;
				return $self->{row};
			},
		next => sub {
			my $self = shift;
			my $row = $self->{sth}->fetchrow_arrayref;
			return undef unless $row;
			$self->{row} = $row;
			return $row->[0], $row->[1];
		}
	};
}

# returns a psuedo object with coderefs next() and line()
sub _get_experiment_data {
	my ($eid, $data_type, $chr, $all) = @_;

	return {
		all => $all,
		data => CoGe::Core::Experiment::query_data(
				eid => $eid,
				data_type => $data_type,
				chr => $chr,
				all => $all
			),
		data_type => $data_type,
		index => -1,
		line => sub {
				my $self = shift;
				return $self->{data}[$self->{index}];
			},
		next => sub {
				my $self = shift;
				$self->{index}++;
				return undef if $self->{index} == scalar @{$self->{data}};
				my $data_point = $self->{data}[$self->{index}];
				if ($self->{data_type} == $DATA_TYPE_ALIGN && $self->{all}) { # sam format
					my (undef, undef, undef, $pos, undef, undef, undef, undef, undef, $seq, undef, undef) = split(/\t/, $data_point);
					$pos = int($pos);
					return ($pos, $pos + length($seq));
				}
				my $c1 = index($data_point, ',') + 1;
				my $c2 = index($data_point, ',', $c1) + 1;
				my $c3 = index($data_point, ',', $c2);
				return (int(substr($data_point, $c1, $c2 - $c1 - 1)), int(substr($data_point, $c2, $c3 - $c2)));
			}
	};
}

sub _get_quant_lines {
    my ($id, $chr, $data_type, $type, $gte, $lte) = @_;
    my $cols = CoGe::Core::Experiment::get_fastbit_format($id)->{columns};
    my @columns = map { $_->{name} } @{$cols};
    $chr = undef if $chr eq 'All';
    my $lines = CoGe::Core::Experiment::query_data(
        eid => $id,
        data_type => $data_type,
        col => join(',', @columns),
        chr => $chr,
        type => $type,
        gte => $gte,
        lte => $lte,
    );
    return (\@columns, $lines);
}

# pass in two psuedo objects. data points in the first one that overlap those in the second are returned in an array
sub _in {
	my $data1 = shift;
	my $data2 = shift;
	my $hits = [];
	my ($start1, $stop1) = $data1->{next}->($data1);
	my ($start2, $stop2) = $data2->{next}->($data2);
	while ($start1 && $start2) {
		while ($start1 && $stop1 < $start2) {
			($start1, $stop1) = $data1->{next}->($data1);
		}
		last unless $start1;
		while ($start2 && $stop2 < $start1) {
			($start2, $stop2) = $data2->{next}->($data2);
		}
		last unless $start2;
		if ($stop1 >= $start2) {
			push @$hits, $data1->{line}->($data1);
		}
		($start1, $stop1) = $data1->{next}->($data1);
	}
	return $hits;
}

sub _markers {
    my $self = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->stash('chr');
    $chr = undef if $chr eq 'Any';
    my $type_names = $self->param('features');

    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
	my $experiment = $db->resultset('Experiment')->find($eid);
	return unless $self->_can_view($experiment, $user);

    my $data1 = _get_experiment_data($experiment->id, $experiment->data_type, $chr, 0);
    my $data2 = _get_db_data($experiment->genome_id, $type_names, $chr, $db->storage->dbh);
    return _in($data1, $data2);
}

# pass in two psuedo objects. data points in the first one that don't overlap those in the second are returned in an array
sub _not_in {
	my $data1 = shift;
	my $data2 = shift;
	my $hits = [];
	my ($start1, $stop1) = $data1->{next}->($data1);
	my ($start2, $stop2) = $data2->{next}->($data2);
	while ($start1 && $start2) {
		while ($start1 && $stop1 < $start2) {
			push @$hits, $data1->{line}->($data1);
			($start1, $stop1) = $data1->{next}->($data1);
		}
		last unless $start1;
		while ($start2 && $stop2 < $start1) {
			($start2, $stop2) = $data2->{next}->($data2);
		}
		last unless $start2;
		if ($stop1 < $start2) {
			push @$hits, $data1->{line}->($data1);
		}
		($start1, $stop1) = $data1->{next}->($data1);
	}
	return $hits;
}

sub overlaps {
	my $self = shift;
	my $not = $self->param('not');

	my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
	my $data1 = $self->_get_data(1, $db);
	my $data2 = $self->_get_data(2, $db);
	my $hits = $not ? _not_in($data1, $data2) : _in($data1, $data2);
	if ($hits) {
		$self->render(json => $hits);
	}
}

sub query_data {
    my $self = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->param('chr');
    $chr = undef if $chr eq 'All';
    my $type = $self->param('type');
    my $gte = $self->param('gte');
    my $lte = $self->param('lte');

    my ($db, $user) = CoGe::Services::Auth::init($self);
	my $experiment = $db->resultset('Experiment')->find($eid);
	return unless $self->_can_view($experiment, $user);

	my $result = CoGe::Core::Experiment::query_data(
		eid => $eid,
		data_type => $experiment->data_type,
		chr => $chr,
		type => $type,
		gte => $gte,
		lte => $lte,
	);
    $self->render(json => $result);
}

sub snps {
    my $self = shift;
    my $results = $self->_snps;
    if ($results) {
        $self->render(json => $results);
    }
}

sub _snps {
    my $self = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->stash('chr');
    $chr = undef if $chr eq 'Any';

    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
	my $experiment = $db->resultset('Experiment')->find($eid);
	return unless $self->_can_view($experiment, $user);

	my $type_names = $self->param('features');
	if ($type_names) {
        my $data1 = _get_experiment_data($experiment->id, $experiment->data_type, $chr, 0);
        my $data2 = _get_db_data($experiment->genome_id, $type_names, $chr, $db->storage->dbh);
        return _in($data1, $data2);
	}
	elsif ($self->param('snp_type')) {
        my $cmdpath = catfile(CoGe::Accessory::Web::get_defaults()->{BINDIR}, 'snp_search', 'snp_search');
        my $storage_path = get_experiment_path($eid);
        opendir(my $dh, $storage_path);
        my @dir_entries = readdir($dh);
        closedir $dh;
        my @files = grep(/\.processed$/, @dir_entries);
        @files = grep(/\.vcf$/, @dir_entries) unless @files;
        my $cmd = "$cmdpath $storage_path/" . $files[0] . ' ' . $self->stash('chr') . ' "' . $self->param('snp_type') . '"';
        my @cmdOut = qx{$cmd};

        my $cmdStatus = $?;
        if ( $? != 0 ) {
            warn "CoGe::Services::API::JBrowse::Experiment: error $? executing command: $cmd";
        }
        my @lines;
        foreach (@cmdOut) {
            chomp;
            push @lines, $_;
        }
        return \@lines;
	}
    else {
        my $snps;
		$snps = CoGe::Core::Experiment::query_data(
			eid => $eid,
			data_type => $experiment->data_type,
			chr => $chr,
		);
        return $snps;
    }
}

sub _transform_line {
    my $line = shift;
    my $transform = shift;
    my @tokens = split ', ', $line;
    $tokens[0] = substr($tokens[0], 1, -1);
    if ($transform) {
        if ($transform eq 'Inflate') {
            $tokens[4] = 1;
        }
        elsif ($transform eq 'Log2') {
            $tokens[4] = log(1 + $tokens[4]);
        }
        elsif ($transform eq 'Log10') {
            $tokens[4] = log(1 + $tokens[4]) / $log10;
        }
    }
    return \@tokens;
}

sub _write {
    my $self = shift;
    my $s = shift;
    my $fh = shift;
    if ($fh) {
        print $fh $s;
    }
    else {
        $self->write($s);
    }
}

sub _write_marker {
    my ($self, $chr, $source, $method, $start, $stop, $score, $strand, $phase, $group, $fh) = @_;
    $self->_write($chr, $fh);
    $self->_write("\t", $fh);
    $self->_write($source || '.', $fh);
    $self->_write("\t", $fh);
    $self->_write($method || '.', $fh);
    $self->_write("\t", $fh);
    $self->_write($start, $fh);
    $self->_write("\t", $fh);
    $self->_write($stop, $fh);
    $self->_write("\t", $fh);
    $self->_write($score, $fh);
    $self->_write("\t", $fh);
    $self->_write($strand || '.', $fh);
    $self->_write("\t", $fh);
    $self->_write($phase || '.', $fh);
    $self->_write("\t", $fh);
    $self->_write($group || '.', $fh);
    $self->_write("\n", $fh);
}

1;
