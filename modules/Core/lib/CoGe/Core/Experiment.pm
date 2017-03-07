package CoGe::Core::Experiment;
use strict;

use Data::Dumper;
use Sort::Versions;
use FastBit;
use CoGe::Core::Storage qw( $DATA_TYPE_QUANT $DATA_TYPE_ALIGN $DATA_TYPE_POLY $DATA_TYPE_MARKER $DATA_TYPE_SEQUENCE get_experiment_path );
use CoGe::Accessory::Web qw(get_command_path);
use CoGe::Accessory::IRODS qw($IRODS_METADATA_PREFIX);

our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION, @QUANT_TYPES, @MARKER_TYPES, @POLYMORPHISM_TYPES, @ALIGNMENT_TYPES, @SEQUENCE_TYPES, @SUPPORTED_TYPES );

BEGIN {
    require Exporter;
    $VERSION = 0.1;
    @ISA = qw(Exporter);
    @EXPORT = qw(
        @QUANT_TYPES @MARKER_TYPES @POLYMORPHISM_TYPES @ALIGNMENT_TYPES @SEQUENCE_TYPES @SUPPORTED_TYPES
        delete_experiment undelete_experiment detect_data_type download_data experimentcmp get_data 
        get_fastbit_format get_fastbit_score_column query_data get_irods_metadata
    );
    
    # Setup supported experiment file types
    @QUANT_TYPES        = qw(csv tsv bed wig bw seg);
    @MARKER_TYPES       = qw(gff gtf gff3);
    @POLYMORPHISM_TYPES = qw(vcf gvcf);
    @ALIGNMENT_TYPES    = qw(sam bam);
    @SEQUENCE_TYPES     = qw(fastq fq);
    @SUPPORTED_TYPES = (@QUANT_TYPES, @MARKER_TYPES, @POLYMORPHISM_TYPES, @ALIGNMENT_TYPES, @SEQUENCE_TYPES);
}

sub experimentcmp($$) { # for sorting DBI-X objects or DBI hashes
    my ($a, $b) = @_;

    if ( ref($a) eq 'HASH' && ref($b) eq 'HASH' ) { # DBI
        versioncmp( $b->{version}, $a->{version} )
          || $a->{name} cmp $b->{name}
          || $b->{id} cmp $a->{id};
    }
    else { # DBI-X
        versioncmp( $b->version, $a->version )
          || $a->name cmp $b->name
          || $b->id cmp $a->id;
    }
}

sub delete_experiment {
    my $id = shift;
    my $db = shift;
    my $user = shift;

    my $experiment = $db->resultset('Experiment')->find($id);
    return 'experiment not found' unless $experiment;

    return 'permission denied' unless $user->is_admin or $user->is_owner(experiment => $experiment);

    $experiment->deleted(1);
    $experiment->update;
    return undef;
}

sub undelete_experiment {
    my $id = shift;
    my $db = shift;
    my $user = shift;

    my $experiment = $db->resultset('Experiment')->find($id);
    return 'experiment not found' unless $experiment;

    return 'permission denied' unless $user->is_admin or $user->is_owner(experiment => $experiment);

    $experiment->deleted(0);
    $experiment->update;
    return undef;
}

sub detect_data_type {
    my $filetype = shift;
    my $filepath = shift;
    #print STDOUT "detect_data_type: $filepath\n";

    # Try to determine type based on file extension
    if (!$filetype or $filetype eq 'autodetect') {
        $filepath =~ s/\.gz$//; # remove extension for compressed files -- mdb added 11/22/16
        $filepath =~ s/\.bz2$//; # remove extension for compressed files -- mdb added 3/6/17

        ($filetype) = lc($filepath) =~ /\.([^\.]+)$/;
    }

    $filetype = lc($filetype);

    if ( grep { $_ eq $filetype } @QUANT_TYPES ) {
        print STDOUT "log: Detected a quantitative file ($filetype)\n";
        return ($filetype, $DATA_TYPE_QUANT);
    }
    elsif ( grep { $_ eq $filetype } @ALIGNMENT_TYPES) {
        print STDOUT "log: Detected an alignment file ($filetype)\n";
        return ($filetype, $DATA_TYPE_ALIGN);
    }
    elsif ( grep { $_ eq $filetype } @POLYMORPHISM_TYPES ) {
        print STDOUT "log: Detected a polymorphism file ($filetype)\n";
        return ($filetype, $DATA_TYPE_POLY);
    }
    elsif ( grep { $_ eq $filetype } @MARKER_TYPES ) {
        print STDOUT "log: Detected a marker file ($filetype)\n";
        return ($filetype, $DATA_TYPE_MARKER);
    }
    elsif ( grep { $_ eq $filetype } @SEQUENCE_TYPES ) {
        print STDOUT "log: Detected a sequence file ($filetype)\n";
        return ($filetype, $DATA_TYPE_SEQUENCE);
    }
    else {
        print STDOUT "detect_data_type: unknown file ext '$filetype'\n";
        return ($filetype);
    }
}

sub get_data {
    my %opts = @_;
    my $eid  = $opts{eid};    # required
    my $data_type = $opts{data_type};
    unless ($eid) {
        print STDERR "CoGe::Core::Experiment::get_data: experiment id not specified!\n";
        return;
    }
    my $chr   = $opts{chr};
    my $start = $opts{start};
    my $stop  = $opts{stop};
    $stop = $opts{end} if ( not defined $stop );
    $start = 0 if ($start < 0);
    $stop = 0 if ($stop < 0);

    if (!$data_type ||
        $data_type == $DATA_TYPE_QUANT ||
        $data_type == $DATA_TYPE_POLY ||
        $data_type == $DATA_TYPE_MARKER)
    {
        my $format = get_fastbit_format($eid, $data_type);
        my $columns = join(',', map { $_->{name} } @{$format->{columns}});
        my $lines = FastBit::query("select $columns where 0.0=0.0 and chr='$chr' and start <= $stop and stop >= $start order by start limit 999999999", $eid);
	    my @results = map {_parse_fastbit_line($format, $_, $chr)} @{$lines};
	    return \@results;
    }
    elsif ( $data_type == $DATA_TYPE_ALIGN ) { # FIXME move output parsing from Storage.pm to here
        my $cmdpath = get_command_path('SAMTOOLS');
    	my $storage_path = get_experiment_path($eid);
        my $cmd = "$cmdpath view $storage_path/alignment.bam $chr:$start-$stop 2>&1";
        #print STDERR "$cmd\n";
        my @cmdOut = qx{$cmd};
        #print STDERR @cmdOut;
        my $cmdStatus = $?;
        if ( $? != 0 ) {
            print STDERR "CoGe::Core::Experiment::get_data: error $? executing command: $cmd\n";
            return;
        }
        
        # Return if error message detected (starts with '[')
        map { return if (/^\[/) } @cmdOut; # mdb added 5/6/15 COGE-594
        
        return \@cmdOut;
    }
    else {
        print STDERR "CoGe::Core::Experiment::get_data: unknown data type\n";
        return;
    }
}

sub get_fastbit_format {
    my $eid = shift;
    my $data_type = shift;

    my $storage_path = get_experiment_path($eid);
    my $format_file = $storage_path . '/format.json';

    # Backward compatibility, see issue 352
    # FIXME: remove someday by adding format.json files to old experiments
    if (not -r $format_file) {
        if (!$data_type || $data_type == $DATA_TYPE_QUANT) {
            return {
                columns => [
                    { name => 'chr',    type => 'key' },
                    { name => 'start',  type => 'unsigned long' },
                    { name => 'stop',   type => 'unsigned long' },
                    { name => 'strand', type => 'byte' },
                    { name => 'value1', type => 'double' },
                    { name => 'value2', type => 'double' }
                ]
            };
        }
        elsif ($data_type == $DATA_TYPE_POLY) {
            return {
                columns => [
                    { name => 'chr',   type => 'key' },
                    { name => 'start', type => 'unsigned long' },
                    { name => 'stop',  type => 'unsigned long' },
                    { name => 'type',  type => 'key' },
                    { name => 'id',    type => 'text' },
                    { name => 'ref',   type => 'key' },
                    { name => 'alt',   type => 'key' },
                    { name => 'qual',  type => 'double' },
                    { name => 'info',  type => 'text' }
                ]
            };
        }
        elsif ( $data_type == $DATA_TYPE_MARKER ) {
            return {
                columns => [
                    { name => 'chr',    type => 'key' },
                    { name => 'start',  type => 'unsigned long' },
                    { name => 'stop',   type => 'unsigned long' },
                    { name => 'strand', type => 'key' },
                    { name => 'type',   type => 'key' },
                    { name => 'score',  type => 'double' },
                    { name => 'attr',   type => 'text' }
                ]
            };
        }
        return; # should never happen!
    }

    # Otherwise read format json file
    return CoGe::Accessory::TDS::read($format_file);
}

sub get_fastbit_score_column {
    my $data_type = shift;
	if (!$data_type || $data_type == $DATA_TYPE_QUANT) {
		return 4;
	}
	if ($data_type == $DATA_TYPE_POLY) {
		return 7;
	}
	if ( $data_type == $DATA_TYPE_MARKER ) {
		return 5;
	}
}

sub _parse_fastbit_line {
    my $format = shift;
    my $line = shift;
    my $chr = shift;

    $line =~ s/"//g;
    my @items = split(/,\s*/, $line);
    return if ( $items[0] !~ /^\"?$chr/); # make sure it's a row output line

    my %result;
    foreach (@{$format->{columns}}) {
        my $item = shift @items;
        my $name = $_->{name};
        my $type = $_->{type};
        if ($type =~ /long|byte|double/) { $result{$name} = 0 + $item } # convert to numeric
        else { $result{$name} = '' . $item; }
    }
    return \%result;
}

sub query_data {
    my %opts = @_;
    my $eid  = $opts{eid};    # required
    unless ($eid) {
        warn 'CoGe::Core::Experiment::query_data: experiment id not specified!';
        return;
    }
    my $chr = $opts{chr};
    my $data_type = $opts{data_type};

    # alignment data
    if ($data_type == 3) {
        my $cmdpath = get_command_path('SAMTOOLS');
        my $storage_path = get_experiment_path($eid);
        my $cmd = "$cmdpath view $storage_path/alignment.bam";
        $cmd .= ' ' . $chr if $chr;
        $cmd .= ' 2>&1';
        my @cmdOut = qx{$cmd};
        my $cmdStatus = $?;
        if ( $? != 0 ) {
            warn "CoGe::Core::Experiment::query_data: error $? executing command: $cmd";
            return;
        }
        # Return if error message detected (starts with '[')
        map { return if (/^\[/) } @cmdOut;

        return \@cmdOut if $opts{all};

        my @results;
        foreach (@cmdOut) {
            chomp;
            my ($qname, $flag, $chr, $start, undef, undef, undef, undef, undef, $seq, undef, undef) = split(/\t/);
            my $strand = ($flag & 0x10 ? '-1' : '1');
            push @results, '"' . $chr . '",' . $start . ',' . ($start + length($seq)) . ',' . $strand . ',"' . $qname . '"';
        }
        return \@results;
    }

    my $col = $opts{col};
    my $gte = $opts{gte};
    my $lte = $opts{lte};
    my $order_by = $opts{order_by} || 'chr,start';
    my $type = $opts{type};

	my $where = '0.0=0.0';
	$where .= " and chr='$chr'" if $chr;
	my $value1;
    if ($type eq 'max') {
    	$value1 = FastBit::max($eid, $chr);
    	$where .= ' and value1>' . ($value1 - 0.001);
    }
    elsif ($type eq 'min') {
    	$value1 = FastBit::min($eid, $chr);
    	$where .= ' and value1<' . ($value1 + 0.001);
    }
    elsif ($type eq 'range') {
		$where .= ' and value1>=' . $gte if $gte;
		$where .= ' and value1<=' . $lte if $lte;
    }
    if (!$col) {
        my $format = get_fastbit_format($eid, $data_type);
        $col = join(',', map { $_->{name} } @{$format->{columns}});
    }
    my $results = FastBit::query("select $col where $where order by $order_by limit 999999999", $eid);
    if ($type eq 'max' || $type eq 'min') {
    	my $value_col = get_fastbit_score_column $data_type;
        my @lines = grep { $value1 == (split(',', $_))[$value_col] } @{$results};
        return \@lines;
    }
    return $results;
}

sub get_irods_metadata {
    my $experiment = shift;

    my %md = (
        $IRODS_METADATA_PREFIX.'link'                      => "http://genomevolution.org",
        $IRODS_METADATA_PREFIX.'ExperimentView-link'       => "http://genomevolution.org/coge/ExperimentView.pl?eid=".$experiment->id,
        $IRODS_METADATA_PREFIX.'experiment-id'             => $experiment->id,
        $IRODS_METADATA_PREFIX.'experiment-name'           => $experiment->name,
        $IRODS_METADATA_PREFIX.'experiment-version'        => $experiment->version,
        $IRODS_METADATA_PREFIX.'experiment-data-type'      => $experiment->data_type_desc,
        $IRODS_METADATA_PREFIX.'-genome-id'      => $experiment->genome->id,
        $IRODS_METADATA_PREFIX.'-genome-summary' => $experiment->genome->info(hideRestrictedSymbol => 1)
    );
    $md{$IRODS_METADATA_PREFIX.'experiment-description'} = $experiment->description if $experiment->description;

    # Add sources
    my $i = 1;
    my @sources = $experiment->source;
    foreach my $item (@sources) {
        my $source = $item->name;
        $source .= ": ".$item->description if $item->description;
        my $key = "experiment-source";
        $key .= $i if scalar @sources > 1;
        $md{$IRODS_METADATA_PREFIX.$key} = $source;
        $md{$IRODS_METADATA_PREFIX.$key."-link"} = $item->link if $item->link;
        $i++;
    }

    # Add tags
    $i = 1;
    my @tags = $experiment->tags;
    foreach my $tag (@tags) {
        my $key = $IRODS_METADATA_PREFIX . (scalar(@tags) > 1 ? "experiment-tag-$i" : "experiment-tag");
        $md{$key} = $tag->name;
        $i++;
    }

    # Add custom metadata
    foreach my $a ( $experiment->annotations ) {
        my $group = (
            defined $a->type->group
            ? $a->type->group->name . ',' . $a->type->name
            : $a->type->name
        );

        $md{$IRODS_METADATA_PREFIX.'experiment-'.$group} = $a->info;
    }

    return \%md;
}

1;
