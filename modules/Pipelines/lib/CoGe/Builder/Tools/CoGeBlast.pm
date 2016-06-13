package CoGe::Builder::Tools::CoGeBlast;

use Moose;

use CoGe::Accessory::Jex;
use CoGe::Accessory::Web;
use File::Basename;
use JSON::XS;

BEGIN {
    use Exporter 'import';
    our @EXPORT_OK = qw( blast_search create_fasta_file get_blast_db );
}

sub blast_search {
    my $db = shift;
    my $user = shift;
    my $config = shift;
    my %opts = @_;

    my $BLASTDBDIR   = $config->{BLASTDB};
    my $MAX_PROC     = $config->{COGE_BLAST_MAX_PROC};
    my $BLAST_PROGS  = {
        blast_legacy => $config->{BLAST} . " -a $MAX_PROC",
        tblastn      => $config->{TBLASTN} . " -num_threads $MAX_PROC",
        tblastx      => $config->{TBLASTX} . " -num_threads $MAX_PROC",
        blastn       => $config->{BLASTN} . " -num_threads $MAX_PROC -task blastn",
        dcmega       => $config->{BLASTN} . " -num_threads $MAX_PROC -task dc-megablast",
        mega         => $config->{BLASTN} . " -num_threads $MAX_PROC -task megablast",
        lastz        => $config->{LASTZ}
    };
    my $JEX          = CoGe::Accessory::Jex->new( host => $config->{JOBSERVER}, port => $config->{JOBPORT} );
    my $PAGE_TITLE   = "CoGeBlast";
    my $PAGE_NAME    = $PAGE_TITLE . ".pl";
    my $RESULTSLIMIT = 100;
    my $TEMPDIR      = $config->{TEMPDIR} . "CoGeBlast";
    my $TEMPURL      = $config->{TEMPURL} . "CoGeBlast";

    my $color_hsps = $opts{color_hsps};
    my $program    = $opts{program};
    my $expect     = $opts{expect};
    my $job_title  = $opts{job_title};
    my $wordsize   = $opts{wordsize};
    my $type       = $opts{type};

    #$wordsize=11 if $program eq "blastn";
    my $comp         = $opts{comp};
    my $matrix       = $opts{matrix};
    my $gapcost      = $opts{gapcost};
    my $match_score  = $opts{matchscore};
    my $filter_query = $opts{filter_query};
    my $resultslimit = $opts{resultslimit} || $RESULTSLIMIT;
    my $basename     = $opts{basename};
    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    #blastz params
    my $zwordsize      = $opts{zwordsize};
    my $zgap_start     = $opts{zgap_start};
    my $zgap_extension = $opts{zgap_extension};
    my $zchaining      = $opts{zchaining};
    my $zthreshold     = $opts{zthreshold};
    my $zmask          = $opts{zmask};

    my $seq = $opts{seq};
    #this is where the dsgids are stored -- stupid name
    my $blastable = $opts{blastable};

    return encode_json({
        success => JSON::false,
        error => "Please specified a sequence of characters to be blasted."
    }) unless $seq;

    return encode_json({
        success => JSON::false,
        error => "Please select genomes to be blasted."
    }) unless $blastable;

    my @dsg_ids = split( /,/, $blastable );

    my $width = $opts{width};
    my $fid   = $opts{fid};

    my $genomes_url = CoGe::Accessory::Web::get_tiny_link(
        user_id => $user->id,
        page    => "GenomeList",
        url     => $config->{SERVER} . "GenomeList.pl?dsgid=$blastable"
    );

    my $list_link =
        qq{<a href="$genomes_url" target_"blank">}
      . @dsg_ids
      . ' genome'
      . ( @dsg_ids > 1 ? 's' : '' ) . '</a>';

    my $log_msg = 'Blast ' . length($seq) . ' characters against ' . $list_link;

    my $gap;

    if ( $gapcost && $gapcost =~ /^(\d+)\s+(\d+)/ ) {
        $gap = qq[$1,$2];
    }

    my %params = (
        color_hsps   => $color_hsps,
        program      => $program,
        expect       => $expect,
        job_title    => $job_title,
        wordsize     => $wordsize,
        comp         => $comp,
        matrix       => $matrix,
        gapcost      => $gap,
        match_score  => $match_score,
        filter_query => $filter_query,
        resultslimit => $resultslimit,
        basename     => $basename,
        zwordsize    => $zwordsize,
        zgap_start   => $zgap_start,
        zgap_exten   => $zgap_extension,
        zchaining    => $zchaining,
        zthreshold   => $zthreshold,
        zmask        => $zmask,
        type         => $type,

        #Genomes
        dsgid        => $blastable,
    );

    # Optional parameters
    $params{fid} = $fid if $fid;

    my $url = url_for($PAGE_NAME, %params);

    my $link = CoGe::Accessory::Web::get_tiny_link(url => $url);

    my ($tiny_id) = $link =~ /\/(\w+)$/;
    my $workflow = $JEX->create_workflow(
        name    => "cogeblast-$tiny_id",
        id      => 0,
        logfile => $cogeweb->logfile
    );

    CoGe::Accessory::Web::write_log( "process $$", $cogeweb->logfile );

    $width = 400 unless $width =~ /^\d+$/;    #something wrong with how width is calculated in tmpl file

    my $t1 = new Benchmark;
    my ( $fasta_file, $query_seqs_info ) = create_fasta_file($seq, $cogeweb);
    my $opts;
    my $pre_command;
    my $x;
    ( $x, $pre_command ) = CoGe::Accessory::Web::check_taint($pre_command);
    my @results;
    my $count = 1;
    my $t2    = new Benchmark;

    foreach my $dsgid (@dsg_ids) {
        my ( $org, $dbfasta, $dsg ) = get_blast_db($dsgid, $db);
        next unless $dbfasta;
        next unless -s $fasta_file;

        my $name = $dsg->organism->name;
        #my $args = [
        #    ['-i', $dbfasta, 0],
        #    ['-t', qq{"$name"}, 0],
        #    ['-n', $dsgid, 1],
        #];

        #push @$args, ['-p', 'F', 1];

        my $dbpath = File::Spec->catdir(($BLASTDBDIR, $dsgid));
        my $db = File::Spec->catdir(($dbpath, $dsgid));
        #my $outputs = [[$dbpath, 1]]; #["$db.nhr", "$db.nin", "$db.nsq"];

        $workflow->add_job(generate_blastdb_job(
            config  => $config,
            title   => $name,
            out     => $dsgid,
            fasta   => $dbfasta,
            type    => "nucl",
            outdir  => $dbpath,
        ));

        #$workflow->add_job(
        #    cmd     => "mkdir $dsgid && cd $dsgid && $FORMATDB",
        #    script  => undef,
        #    args    => $args,
        #    inputs  => undef,
        #    outputs => $outputs,
        #    description => "Generating blastable database..."
        #);

        my $outfile = $cogeweb->basefile . "-$count.$program";

        my $cmd = $BLAST_PROGS->{$program};
        my $args = [
            [ '', '--adjustment=10', 1 ],
            [ '', $BLAST_PROGS->{$program}, 0 ],
        ];

        if ( $program eq "lastz" ) {
            push @$args, [ '',  $fasta_file, 1 ];
            push @$args, [ '', "W=" . $zwordsize,  1 ] if defined $zwordsize;
            push @$args, [ '', "C=" . $zchaining,  1 ] if defined $zchaining;
            push @$args, [ '', "K=" . $zthreshold, 1 ] if defined $zthreshold;
            push @$args, [ '', "M=" . $zmask,      1 ] if defined $zmask;
            push @$args, [ '', "O=" . $zgap_start, 1 ] if defined $zgap_start;
            push @$args, [ '', "E=" . $zgap_extension, 1 ]
              if defined $zgap_extension;
            push @$args, [ '',  $dbfasta,    0 ];
            push @$args, [ '>', $outfile,    1 ];
        }
        else {
            my ( $nuc_penalty, $nuc_reward, $exist, $extent );
            if ( $gapcost && $gapcost =~ /^(\d+)\s+(\d+)/ ) {
                ( $exist, $extent ) = ( $1, $2 );
            }

            if ($match_score && $match_score =~ /^(\d+)\,(-\d+)/ ) {
                ( $nuc_penalty, $nuc_reward ) = ( $2, $1 );
            }

            push @$args, [ "-comp_based_stats", 1, 1 ] if $program eq "tblastn";
            push @$args, [ '-matrix', $matrix, 1 ] if $program =~ /tblast/i;
            push @$args, [ '-penalty', $nuc_penalty, 1 ]
              unless $program =~ /tblast/i;
            push @$args, [ '-reward', $nuc_reward, 1 ]
              unless $program =~ /tblast/i;
            push @$args, [ '-gapopen', $exist, 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-gapextend', $extent, 1 ]
              unless $program =~ /tblast/i;
            push @$args, [ '-dust', 'no', 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-query',     $fasta_file, 1 ];
            push @$args, [ '-word_size', $wordsize,   1 ];
            push @$args, [ '-evalue',    $expect,     1 ];
            push @$args, [ '-db',        $db,         0 ];
            push @$args, [ '>',          $outfile,    1 ];
        }

        push @results,
          {
            command  => $cmd,
            file     => $outfile,
            organism => $org,
            dsg      => $dsg
          };

        $workflow->add_job({
            cmd     => "/usr/bin/nice",
            script  => undef,
            args    => $args,
            inputs  => [$fasta_file, [$dbpath, 1]],
            outputs => [$outfile],
            description => "Blasting sequence against $name"
        });

        $count++;
    }

    my $response = $JEX->submit_workflow($workflow);

    my $log = CoGe::Accessory::Web::log_history(
        db          => $db,
        user_id     => $user->id,
        page        => $PAGE_TITLE,
        description => $log_msg,
        link        => $link,
        parent_id   => $response->{id},
        parent_type => 7 #FIXME magic number
    ) if $response and $response->{id};

    return encode_json({
        id => $response->{id},
        link => $link,
        logfile => $TEMPURL . "/" . $cogeweb->basefilename . ".log",
        success => $JEX->is_successful($response) ? JSON::true : JSON::false
    })
}

sub build {
    my $self = shift;
    my $genomes = $self->params->{genomes};
    my $notebooks = $self->params->{notebooks};
    my $type = $self->params->{type};
    my $e_value = $self->params->{e_value};
    my $word_size = $self->params->{word_size};
    my $gap_costs = $self->params->{gap_costs};
    my $filter_query = $self->params->{filter_query};
    my $max_results = $self->params->{max_results};
    my $query_seq = $self->params->{query_seq};
    my $matrix = $self->params->{matrix};

    my ( $db, $user, $config ) = CoGe::Accessory::Web->init();

    my @gids;
    @gids = @$genomes if $genomes;
    if ($notebooks) {
        for (@$notebooks) {
            for (@{$_->genomes}) {
                push @gids, $_->id;
            }
            for (@{$_->experiments}) {
                push @gids, $_->genome_id;
            }
        }
    }

    my $resp = blast_search($db, $user, $config,
        workflow     => $self->workflow,
        db           => $self->db,
        config       => $self->config,
        blastable    => join(',', @gids),
        type         => $type,
        expect       => $e_value,
        wordsize     => $word_size,
        gapcost      => $gap_costs->[0] . ' ' . $gap_costs->[1],
        filter_query => $filter_query,
        resultslimit => $max_results,
        seq          => $query_seq,
        matrix       => $matrix
    );
    return 0 if ($resp); # an error occurred
    return 1;
}

sub create_fasta_file {
    my $seq = shift;
    my $cogeweb = shift;
    my %seqs;    #names and lengths
    $seq =~ s/>\s*\n//;
    $seq = ">seq\n" . $seq unless $seq =~ />/;
    if ( $seq =~ />/ ) {
        foreach ( split( /\n>/, $seq ) ) {
            next unless $_;
            my ( $name, $tmp ) = split( /\n/, $_, 2 );
            $name =~ s/^>//;
            next unless $tmp;
            $tmp  =~ s/\n//g;
            $tmp  =~ s/\s//g;
            $name =~ s/\s//g
              ; #need to remove spaces due to how blast breaks query names at spaces or commas
            $seqs{$name} = length($tmp);
        }
    }
    CoGe::Accessory::Web::write_log( "creating user's fasta file",
        $cogeweb->logfile );
    open( NEW, "> " . $cogeweb->basefile . ".fasta" );
    print NEW $seq;
    close NEW;
    return $cogeweb->basefile . ".fasta", \%seqs;
}

sub generate_blastdb_job {
    my %opts = @_;

    # required arguments
    my $config = $opts{config};
    my $title = $opts{title};
    my $fasta = $opts{fasta};
    my $type = $opts{type};
    my $out  = $opts{out};
    my $outdir = $opts{outdir};

    my $logfile = $opts{logfile} || "db.log";
    my $BLASTDB = $config->{MAKEBLASTDB} || "makeblastdb";

    my $args = [
        ["-in", $fasta, 0],
        ["-out", $out, 0],
        ["-dbtype", $type, 0],
        ["-title", qq{"$title"}, 0],
        ["-logfile", $logfile, 0],
    ];

    my $base = basename($outdir);

    return {
        cmd => "mkdir $base && cd $base && $BLASTDB",
        script  => undef,
        args    => $args,
        inputs  => undef,
        outputs => [[$outdir, 1]],
        description => "Generating blastable database..."
    };
}

sub get_blast_db {
    my $dsgid = shift;
    my $db = shift;
    my ($dsg) = $db->resultset('Genome')->search(
        { genome_id => $dsgid },
        {
            join     => [ 'organism', 'genomic_sequence_type' ],
            prefetch => [ 'organism', 'genomic_sequence_type' ],
        }
    );
    unless ($dsg) {
        print STDERR "Problem getting dataset group for dsgid $dsgid\n";
        return;
    }
    my ($ds) = $dsg->datasets;
    my $org_name =
        $dsg->organism->name . " ("
      . $ds->data_source->name . " "
      . $dsg->type->name . " v"
      . $dsg->version . ")";

    #$org_name .= " (".$gst->name.")" if $gst;

    my $file_path      = $dsg->file_path;
    return unless $file_path && -r $file_path;
    return $org_name, $file_path, $dsg;
}

with qw(CoGe::Builder::Buildable);

1;
