package CoGe::Services::JBrowse::Experiment;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::Storage qw( get_experiment_data );
use JSON::XS;

#TODO: use these from Storage.pm instead of redeclaring them
my $DATA_TYPE_QUANT = 1;	# Quantitative data
my $DATA_TYPE_POLY	= 2;	# Polymorphism data
my $DATA_TYPE_ALIGN = 3;	# Alignments

my $NUM_QUANT_COL   = 6;
my $NUM_VCF_COL     = 9;
my $MAX_EXPERIMENTS = 20;
#my $MAX_RESULTS = 150000;
my $MAX_WINDOW_SIZE = 500000;

my $QUAL_ENCODING_OFFSET = 33;

sub setup {
    my $self = shift;
    $self->run_modes(
        'stats_global' => 'stats_global',
        'stats_region' => 'stats_region',
        'features'     => 'features',
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

sub stats_region {    #FIXME lots of code in common with features()
    my $self     = shift;
    my $eid      = $self->param('eid');
    my $nid      = $self->param('nid');
    my $gid      = $self->param('gid');
    my $chr      = $self->param('chr');
    my $start    = $self->query->param('start');
    my $end      = $self->query->param('end');
    my $bpPerBin = $self->query->param('bpPerBin');

    print STDERR "JBrowse::Experiment::stats_region eid="
      . ( $eid ? $eid : '' ) . " nid="
      . ( $nid ? $nid : '' ) . " gid="
      . ( $gid ? $gid : '' )
      . " $chr:$start:$end ("
      . ( $end - $start + 1 )
      . ") bpPerBin=$bpPerBin\n";

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;

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
        unless ( !$e->restricted || $user->has_access_to_experiment($e) ) {
        	print STDERR "JBrowse::Experiment::stats_region access denied to experiment $eid\n";
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
        elsif ( $data_type == $DATA_TYPE_POLY ) {
            my $cmdOut = CoGe::Accessory::Storage::get_experiment_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

            # Convert FastBit output into JSON
            foreach (@$cmdOut) {
                chomp;
                if (/^\"/) {    #if (/^\"$chr\"/) { # potential result line
                    s/"//g;
                    my @items = split(/,\s*/);
                    next
                      if ( @items != $NUM_VCF_COL )
                      ; # || $items[0] !~ /^\"?$chr/); # make sure it's a row output line

                    # Bin results
                    my ( undef, $s, $e ) = @items;
                    $e = $s + 1 if ( $e == $s );    #FIXME revisit this
                    my $mid = int( ( $s + $e ) / 2 );
                    my $b   = int( ( $mid - $start ) / $bpPerBin );
                    $b = 0 if ( $b < 0 );
                    $bins[$b]++;
                }
            }
        }
        elsif ( $data_type == $DATA_TYPE_ALIGN ) {
        	next;
        }
        else {
        	print STDERR "JBrowse::Experiment::stats_region unknown data type for experiment $eid\n";
        	next;
        }
    }

    my ( $max, $sum, $count );
    foreach my $x (@bins) {
        $sum += x;
        $max = $x if ( not defined $max or $x > $max );
        $count++;
    }
    my $mean = 0;
    $mean = $sum / $count if ($count);

    return encode_json(
        {
            bins => \@bins,
            stats =>
              [ { basesPerBin => $bpPerBin, max => $max, mean => $mean } ]
        }
    );
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

    if ( $end - $start + 1 > $MAX_WINDOW_SIZE ) {
        print STDERR "experiment features maxed\n";
        return qq{{ "features" : [ ] }};
    }

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;

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
        unless ( !$e->restricted || $user->has_access_to_experiment($e) ) {
        	print STDERR "JBrowse::Experiment::features access denied to experiment $eid\n";
        	next
    	}
        push @experiments, $e;
    }
    splice( @experiments, $MAX_EXPERIMENTS, @experiments );

	# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my $results = '';
    foreach my $exp (@experiments) {
        my $eid       = $exp->id;
        my $data_type = $exp->data_type;

        if ( !$data_type || $data_type == $DATA_TYPE_QUANT ) {
            my $cmdOut = CoGe::Accessory::Storage::get_experiment_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

            # Convert FastBit output into JSON
            foreach (@$cmdOut) {
                chomp;
                if (/^\"/) {    #if (/^\"$chr\"/) { # potential result line
                    s/"//g;
                    my @items = split(/,\s*/);
                    next
                      if ( @items != $NUM_QUANT_COL )
                      ; # || $items[0] !~ /^\"?$chr/); # make sure it's a row output line
                    for ( my $i = 0 ; $i < @items ; $i++ ) {
                        $items[$i] = 1 if $items[$i] !~ /\w/; # what's this for?
                    }

                    #$results .= '[' . join(',', map {"\"$_\""} @items) . ']';

                    my ( $chr, $start, $end, $strand, $value1, $value2 ) =
                      @items;
                    $end = $start + 1 if ( $end == $start ); #FIXME revisit this
                    $strand = -1 if ( $strand == 0 );
                    $value1 = $strand * $value1;
                    $results .= ( $results ? ',' : '' )
                      . qq{{ "id": $eid, "start": $start, "end": $end, "score": $value1, "score2": $value2 }};
                }
            }
        }
        elsif ( $data_type == $DATA_TYPE_POLY ) {
            my $cmdOut = CoGe::Accessory::Storage::get_experiment_data(
                eid   => $eid,
                data_type  => $exp->data_type,
                chr   => $chr,
                start => $start,
                end   => $end
            );

            # Convert FastBit output into JSON
            foreach (@$cmdOut) {
                chomp;
                if (/^\"/) {    #if (/^\"$chr\"/) { # potential result line
                    s/"//g;
                    my @items = split(/,\s*/);
                    next
                      if ( @items != $NUM_VCF_COL )
                      ; # || $items[0] !~ /^\"?$chr/); # make sure it's a row output line

                    my (
                        $chr, $start, $end,  $type, $id,
                        $ref, $alt,   $qual, $info
                    ) = @items;
                    $end = $start + 1 if ( $end == $start ); #FIXME revisit this
                    my $name =
                      ( ( $id && $id ne '.' ) ? "$id " : '' )
                      . "$type $ref > $alt";
                    $type = $type . $ref . 'to' . $alt
                      if ( lc($type) eq 'snp' );
                    $results .= ( $results ? ',' : '' )
                      . qq{{ "id": $eid, "name": "$name", "type": "$type", "start": $start, "end": $end, "ref": "$ref", "alt": "$alt", "score": $qual, "info": "$info" }};
                }
            }
        }
        elsif ( $data_type == $DATA_TYPE_ALIGN ) {
	        my $cmdOut = CoGe::Accessory::Storage::get_experiment_data(
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
	        	my $start = $pos;
	        	my $end = $pos + $len;
	        	my $strand = ($flag & 0x10 ? '-1' : '1');
	        	#TODO reverse complement sequence if neg strand?
	        	my $qual = join(' ', map { $_ - $QUAL_ENCODING_OFFSET } unpack("C*", $qual));
	        	
	        	$results .= ( $results ? ',' : '' )
                  . qq{{ "uniqueID": "$qname", "name": "$qname", "start": $start, "end": $end, "strand": $strand, "score": $mapq, "seq": "$seq", "qual": "$qual", "Seq length": $len, "CIGAR": "$cigar" }};
	        }
        }
        else {
        	print STDERR "JBrowse::Experiment::features unknown data type $data_type for experiment $eid\n";
        	next;
        }
    }

    #	print STDERR "{ 'features' : [ $results ] }\n";
    return qq{{ "features" : [ $results ] }};
}

1;
