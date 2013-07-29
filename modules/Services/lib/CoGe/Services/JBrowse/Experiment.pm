package CoGe::Services::JBrowse::Experiment;
use base 'CGI::Application';

use CoGeX;
use CoGe::Accessory::Web;
use JSON::XS;

my $NUM_QUANT_COL   = 6;
my $NUM_VCF_COL     = 9;
my $MAX_EXPERIMENTS = 20;

#my $MAX_RESULTS = 150000;
my $MAX_WINDOW_SIZE = 500000;

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
    print STDERR "experiment stats_global\n";
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

    print STDERR "experiment stats_region eid="
      . ( $eid ? $eid : '' ) . " nid="
      . ( $nid ? $nid : '' ) . " gid="
      . ( $gid ? $gid : '' )
      . " $chr:$start:$end ("
      . ( $end - $start + 1 )
      . ") bpPerBin=$bpPerBin\n";

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;
    my $cmdpath = $conf->{FASTBIT_QUERY};

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
        next unless ( !$e->restricted || $user->has_access_to_experiment($e) );
        push @experiments, $e;
    }
    splice( @experiments, $MAX_EXPERIMENTS, @experiments );

# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my $results = '';
    my @bins;
    foreach my $exp (@experiments)
    { #TODO need to move this code along with replicate in bin/fastbit_query.pl into CoGe::Web sub-module
        my $storage_path = $exp->storage_path;
        my $data_type    = $exp->data_type;

        if ( !$data_type or $data_type < 2 )
        {    #FIXME hardcoded data_type to "quant"
            next;    # skip this experiment -- is this right?
        }
        elsif ( $data_type == 2 ) {    #FIXME hardcoded data_type to "snp"
             # Call FastBit to do query (see issue 61: query string must contain a "." for fastbit to use consistent output)
            my $cmd =
"$cmdpath -v 1 -d $storage_path -q \"select chr,start,stop,type,id,ref,alt,qual,info where 0.0=0.0 and chr='$chr' and start <= $end and stop >= $start order by start limit 999999999\" 2>&1";

            #print STDERR "$cmd\n";
            my @cmdOut = qx{$cmd};

            #print STDERR @cmdOut;
            my $cmdStatus = $?;
            die "Error executing command $CMDPATH ($cmdStatus)"
              if ( $cmdStatus != 0 );

            # Convert FastBit output into JSON
            foreach (@cmdOut) {
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
    print STDERR "experiment features eid="
      . ( $eid ? $eid : '' ) . " nid="
      . ( $nid ? $nid : '' ) . " gid="
      . ( $gid ? $gid : '' )
      . " $chr:$start:$end ("
      . ( $end - $start + 1 ) . ")\n";
    return unless ( ( $eid or $nid or $gid ) and $chr and $start and $end );

    if ( $end - $start + 1 > $MAX_WINDOW_SIZE ) {
        print STDERR "experiment features maxed\n";
        return qq{{ "features" : [ ] }};
    }

    # Connect to the database
    my ( $db, $user, $conf ) = CoGe::Accessory::Web->init;
    my $cmdpath = $conf->{FASTBIT_QUERY};

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
        next unless ( !$e->restricted || $user->has_access_to_experiment($e) );
        push @experiments, $e;
    }
    splice( @experiments, $MAX_EXPERIMENTS, @experiments );

# Query range for each experiment and build up json response - #TODO could parallelize this for multiple experiments
    my $results = '';
    foreach my $exp (@experiments)
    { #TODO need to move this code along with replicate in bin/fastbit_query.pl into CoGe::Web sub-module
        my $eid          = $exp->id;
        my $storage_path = $exp->storage_path;
        my $data_type    = $exp->data_type;

        if ( !$data_type or $data_type < 2 )
        {    #FIXME hardcoded data_type to "quant"
             # Call FastBit to do query (see issue 61: query string must contain a "." for fastbit to use consistent output)
            my $cmd =
"$cmdpath -v 1 -d $storage_path -q \"select chr,start,stop,strand,value1,value2 where 0.0=0.0 and chr='$chr' and start <= $end and stop >= $start order by start limit 999999999\" 2>&1";

            print STDERR "$cmd\n";
            my @cmdOut = qx{$cmd};

            print STDERR @cmdOut;
            my $cmdStatus = $?;
            die "Error executing command $CMDPATH ($cmdStatus)"
              if ( $cmdStatus != 0 );

            # Convert FastBit output into JSON
            foreach (@cmdOut) {
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
        elsif ( $data_type == 2 ) {    #FIXME hardcoded data_type to "snp"
             # Call FastBit to do query (see issue 61: query string must contain a "." for fastbit to use consistent output)
            my $cmd =
"$CMDPATH -v 1 -d $storage_path -q \"select chr,start,stop,type,id,ref,alt,qual,info where 0.0=0.0 and chr='$chr' and start <= $end and stop >= $start order by start limit 999999999\" 2>&1";

            #print STDERR "$cmd\n";
            my @cmdOut = qx{$cmd};

            #print STDERR @cmdOut;
            my $cmdStatus = $?;
            die "Error executing command $CMDPATH ($cmdStatus)"
              if ( $cmdStatus != 0 );

            # Convert FastBit output into JSON
            foreach (@cmdOut) {
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
                      . qq{{ "id": $eid, "name": "$name", "type": "$type", "start": $start, "end": $end, "score": $qual, "info": "$info" }};
                }
            }
        }
    }

    #	print STDERR "{ 'features' : [ $results ] }\n";
    return qq{{ "features" : [ $results ] }};
}

1;
