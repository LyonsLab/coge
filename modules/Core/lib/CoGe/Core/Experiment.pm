package CoGe::Core::Experiment;
use strict;

use Data::Dumper;
use Sort::Versions;
use CoGe::Accessory::FastBit;
use CoGe::Core::Storage qw( $DATA_TYPE_QUANT $DATA_TYPE_ALIGN $DATA_TYPE_POLY $DATA_TYPE_MARKER get_experiment_path );

our ( @EXPORT, @EXPORT_OK, @ISA, $VERSION, @QUANT_TYPES, @MARKER_TYPES, 
      @OTHER_TYPES, @SUPPORTED_TYPES );

BEGIN {
    require Exporter;
    $VERSION = 0.1;
    @ISA = qw( Exporter );
    @EXPORT = qw(@QUANT_TYPES @MARKER_TYPES @OTHER_TYPES @SUPPORTED_TYPES);
    @EXPORT_OK = qw( detect_data_type download_data experimentcmp get_data get_fastbit_format get_fastbit_score_column query_data );
    
    # Setup supported experiment file types
    @QUANT_TYPES = qw(csv tsv bed wig);
    @MARKER_TYPES = qw(gff gtf gff3);
    @OTHER_TYPES = qw(bam vcf);
    @SUPPORTED_TYPES = (@QUANT_TYPES, @MARKER_TYPES, @OTHER_TYPES);
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

sub detect_data_type {
    my $filetype = shift;
    my $filepath = shift;
    #print STDOUT "detect_data_type: $filepath\n";

    if (!$filetype or $filetype eq 'autodetect') {
        # Try to determine type based on file extension
        #print STDOUT "log: Detecting file type\n";
        ($filetype) = lc($filepath) =~ /\.([^\.]+)$/;
    }
    
    $filetype = lc($filetype);

    if ( grep { $_ eq $filetype } @QUANT_TYPES ) {
        print STDOUT "log: Detected a quantitative file ($filetype)\n";
        return ($filetype, $DATA_TYPE_QUANT);
    }
    elsif ( $filetype eq 'bam' ) {
        print STDOUT "log: Detected an alignment file ($filetype)\n";
        return ($filetype, $DATA_TYPE_ALIGN);
    }
    elsif ( $filetype eq 'vcf' ) {
        print STDOUT "log: Detected a polymorphism file ($filetype)\n";
        return ($filetype, $DATA_TYPE_POLY);
    }
    elsif ( grep { $_ eq $filetype } @MARKER_TYPES ) {
        print STDOUT "log: Detected a marker file ($filetype)\n";
        return ($filetype, $DATA_TYPE_MARKER);
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
        my $pFormat = get_fastbit_format($eid, $data_type);
        my $columns = join(',', map { $_->{name} } @{$pFormat->{columns}});
        my $lines = CoGe::Accessory::FastBit::query("select $columns where 0.0=0.0 and chr='$chr' and start <= $stop and stop >= $start order by start limit 999999999", $eid);

	    my @results = map {_parse_fastbit_line($pFormat, $_, $chr)} @{$lines};
	    return \@results;
    }
    elsif ( $data_type == $DATA_TYPE_ALIGN ) { # FIXME move output parsing from Storage.pm to here
        my $cmdpath = CoGe::Accessory::Web::get_defaults()->{SAMTOOLS};
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

# FIXME: move to FastBit.pm?
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

# FIXME: move to FastBit.pm?
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

# FIXME: move to FastBit.pm?
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

sub debug {
	my $data = shift;
	my $new_file = shift;
	my $OUTFILE;
	open $OUTFILE, ($new_file ? ">/tmp/sean" : ">>/tmp/sean");
	print {$OUTFILE} Dumper $data;
	print {$OUTFILE} "\n";
	close $OUTFILE;
}

sub query_data {
    my %opts = @_;
    my $eid  = $opts{eid};    # required
    unless ($eid) {
        warn 'CoGe::Core::Experiment::query_data: experiment id not specified!';
        return;
    }
    my $data_type = $opts{data_type};
    my $col = $opts{col};
    my $chr = $opts{chr};
    my $type = $opts{type};
    my $gte = $opts{gte};
    my $lte = $opts{lte};
    my $order_by = $opts{order_by};
    my $limit = $opts{limit};
    if (!$data_type ||
        $data_type == $DATA_TYPE_QUANT ||
        $data_type == $DATA_TYPE_POLY ||
        $data_type == $DATA_TYPE_MARKER)
    {
    	my $where = '0.0=0.0';
    	if ($chr) {
    		$where .= " and chr='$chr'";
    	} 
	    if ($type eq 'max') {
	    	my $max = CoGe::Accessory::FastBit::query('select max(value1) where 0.0=0.0', $eid);
	    	debug $max;
	    	$where .= ' and value1=' . $max->[0];
	    }
	    elsif ($type eq 'min') {
	    	my $min = CoGe::Accessory::FastBit::query('select min(value1) where 0.0=0.0', $eid);
	    	$where .= ' and value1=' . $min->[0];
	    }
	    elsif ($type eq 'range') {
	    	if ($gte) {
	    		$where .= ' and value1>=' . $gte;
	    	}
	    	if ($lte) {
	    		$where .= ' and value1<=' . $lte;
	    	}
	    }
    	my $query = "select $col where $where";
    	if ($order_by) {
    		$query .= " order by $order_by";
    	}
		$query .= ($limit ? " limit $limit" : " limit 999999999");
		debug $query;
        return CoGe::Accessory::FastBit::query($query, $eid);
    }
}

1;
