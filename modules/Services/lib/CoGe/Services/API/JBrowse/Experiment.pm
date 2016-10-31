package CoGe::Services::API::JBrowse::Experiment;

use Mojo::Base 'Mojolicious::Controller';
use Mojo::JSON;

use CoGeX;
use CoGe::Accessory::histogram;
use CoGe::Accessory::IRODS qw( irods_iput );
use CoGe::Accessory::Utils qw( to_number );
use CoGe::Services::Auth qw( init );
use CoGe::Core::Experiment qw( get_data );
use CoGe::Core::Storage qw( get_experiment_path get_upload_path get_experiment_metadata );
use CoGeDBI qw( get_dataset_ids feature_type_names_to_id );
use Data::Dumper;
use File::Path;
use File::Spec::Functions;
use File::Temp qw( tempfile tempdir );
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

my $log10 = log(10);

sub stats_global { # mdb rewritten 8/26/16 COGE-270
    my $self = shift;
    my $eid  = $self->stash('eid');
    my $nid  = $self->stash('nid');
    my $gid  = $self->stash('gid');
	
	# Authenticate user and connect to the database
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

    # Authenticate user and connect to the database
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

    # Check permissions
    unless (!$experiment->restricted() || ($user && ($user->is_admin() || $user->is_owner_editor(experiment => $id)))) {
        $self->render(json => {
            error => { Auth => "Access denied" }
        }, status => 401);
        return;
    }
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

sub _get_experiments {
	my $db = shift;
	my $user = shift;
	my $eid = shift;
	my $gid = shift;
	my $nid = shift;

    my @all_experiments;
    if ($eid) {
        foreach (split(/,/, $eid)) {
            my $experiment = $db->resultset('Experiment')->find($_);
            push @all_experiments, $experiment if $experiment;
        }
    }
    elsif ($nid) {
        my $notebook = $db->resultset('List')->find($nid);
        push @all_experiments, $notebook->experiments if $notebook;
    }
    elsif ($gid) {
        my $genome = $db->resultset('Genome')->find($gid);
        push @all_experiments, $genome->experiments if $genome;
    }

    # Filter experiments based on permissions
    my @experiments;
    foreach my $e (@all_experiments) {
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

sub query_data {
    my $self = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->param('chr');
    $chr = undef if $chr eq 'All';
    my $type = $self->param('type');
    my $gte = $self->param('gte');
    my $lte = $self->param('lte');

    # Authenticate user and connect to the database
    my ($db, $user) = CoGe::Services::Auth::init($self);

	my $experiments = _get_experiments($db, $user, $eid);
	my $result = [];
	foreach my $experiment (@$experiments) { # this doesn't yet handle multiple experiments (i.e notebooks)
		$result = CoGe::Core::Experiment::query_data(
			eid => $eid,
			data_type => $experiment->data_type,
			chr => $chr,
			type => $type,
			gte => $gte,
			lte => $lte,
		);
	}
    $self->render(json => $result);
}

# returns a psuedo object with coderefs next() and line()
sub get_db_data {
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
sub get_experiment_data {
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
    my $eid = $self->stash('eid');
    my $eid2 = $self->stash('eid2');
    my $chr = $self->stash('chr');
    my $not = $self->param('not');

    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);
    my $experiment = $db->resultset('Experiment')->find($eid);
    my $experiment2 = $db->resultset('Experiment')->find($eid2);
    my $data1 = get_experiment_data($experiment->id, $experiment->data_type, $chr, 0);
    my $data2 = get_experiment_data($experiment2->id, $experiment2->data_type, $chr, 0);
    my $hits = $not ? _not_in($data1, $data2) : _in($data1, $data2);
    if ($hits) {
        $self->render(json => $hits);
    }
}

sub alignments {
    my $self = shift;
    my $results = $self->_alignments;
    if ($results) {
        $self->render(json => $results);
    }
}

sub _alignments {
    my $self = shift;
    my $all = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->stash('chr');
    $chr = undef if $chr eq 'Any' || $chr eq 'All';
    my $type_names = $self->param('features');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    my $experiments = _get_experiments($db, $user, $eid);
    if (scalar @{$experiments} == 0) {
        $self->render(json => { error => 'User does not have permission to view this experiment' });
        return undef;
    }
    
    my $experiment = $experiments->[0]; # this doesn't yet handle multiple experiments (i.e notebooks)
    my $data1 = get_experiment_data($experiment->id, $experiment->data_type, $chr, $all);
    my $data2 = get_db_data($experiment->genome_id, $type_names, $chr, $db->storage->dbh);
    return _in($data1, $data2);
#    return find_overlapping($experiments, $type_names, $chr, $db, $all ? \&get_start_end_sam : \&get_start_end_fastbit, $all);
}

sub markers {
    my $self = shift;
    my $results = $self->_markers;
    if ($results) {
        $self->render(json => $results);
    }
}

sub _markers {
    my $self = shift;
    my $eid = $self->stash('eid');
    my $chr = $self->stash('chr');
    $chr = undef if $chr eq 'Any';
    my $type_names = $self->param('features');

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

    my $experiments = _get_experiments($db, $user, $eid);
    if (scalar @{$experiments} == 0) {
        $self->render(json => { error => 'User does not have permission to view this experiment' });
        return undef;
    }

    my $experiment = $experiments->[0]; # this doesn't yet handle multiple experiments (i.e notebooks)
    my $data1 = get_experiment_data($experiment->id, $experiment->data_type, $chr, 0);
    my $data2 = get_db_data($experiment->genome_id, $type_names, $chr, $db->storage->dbh);
    return _in($data1, $data2);
    # return find_overlapping($experiments, $type_names, $chr, $db, \&get_start_end_fastbit);
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

    # Authenticate user and connect to the database
    my ($db, $user, $conf) = CoGe::Services::Auth::init($self);

	my $experiments = _get_experiments($db, $user, $eid);
    if (scalar @{$experiments} == 0) {
        $self->render(json => { error => 'User does not have permission to view this experiment' });
        return undef;
    }

	my $type_names = $self->param('features');
	if ($type_names) {
        my $experiment = $experiments->[0]; # this doesn't yet handle multiple experiments (i.e notebooks)
        my $data1 = get_experiment_data($experiment->id, $experiment->data_type, $chr, 0);
        my $data2 = get_db_data($experiment->genome_id, $type_names, $chr, $db->storage->dbh);
        return _in($data1, $data2);
        # return find_overlapping($experiments, $type_names, $chr, $db, \&get_start_end_fastbit);
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
        foreach my $experiment (@$experiments) { # this doesn't yet handle multiple experiments (i.e notebooks)
            $snps = CoGe::Core::Experiment::query_data(
                eid => $eid,
                data_type => $experiment->data_type,
                chr => $chr,
            );
        }
        return $snps;
    }
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

    # Authenticate user and connect to the database
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
                        "Seq length" => $len, 
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
