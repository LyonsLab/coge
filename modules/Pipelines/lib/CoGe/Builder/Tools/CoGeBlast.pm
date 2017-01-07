package CoGe::Builder::Tools::CoGeBlast;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Accessory::Utils qw(sanitize_name);
use CoGe::Accessory::Web qw(download_url_for get_command_path url_for);
use CoGe::Builder::CommonTasks qw(add_workflow_result);
use CoGe::Core::Storage qw(get_download_path);
use Data::Dumper;
use File::Basename;
use File::Path qw(make_path);
use File::Spec::Functions;
use JSON::XS;

BEGIN {
    use Exporter 'import';
    our @EXPORT_OK = qw( create_fasta_file get_blast_db get_genomes get_tiny_url );
}

sub add_jobs {
    my %opts = @_;
    my $workflow = $opts{workflow};
    my $db = $opts{db};
    my $user = $opts{user};
    my $config = $opts{config};

    my $BLASTDBDIR   = $config->{BLASTDB};
    my $MAX_PROC     = $config->{COGE_BLAST_MAX_PROC} // 8;
    my $BLAST_PROGS  = {
        blast_legacy => get_command_path('BLAST') . " -a $MAX_PROC",
        tblastn      => get_command_path('TBLASTN') . " -num_threads $MAX_PROC",
        tblastx      => get_command_path('TBLASTX') . " -num_threads $MAX_PROC",
        blastn       => get_command_path('BLASTN') . " -num_threads $MAX_PROC -task blastn",
        dcmega       => get_command_path('BLASTN') . " -num_threads $MAX_PROC -task dc-megablast",
        mega         => get_command_path('BLASTN') . " -num_threads $MAX_PROC -task megablast",
        lastz        => get_command_path('LASTZ')
    };

    my $program    = $opts{program};
    my $outfmt    = $opts{outfmt};
    my $expect     = $opts{expect};
    my $wordsize   = $opts{wordsize};
    #$wordsize=11 if $program eq "blastn";
    my $matrix       = $opts{matrix};
    my $gapcost      = $opts{gapcost};
    my $match_score  = $opts{matchscore};
    my $filter_query = $opts{filter_query};
    my $cogeweb      = $opts{cogeweb};
    unless ($cogeweb) {
        $cogeweb = CoGe::Accessory::Web::initialize_basefile(
            basename => $opts{basename},
            tempdir  => $config->{TEMPDIR} . "CoGeBlast"
        );
    }
    $workflow->logfile($cogeweb->logfile);

    #blastz params
    my $zwordsize      = $opts{zwordsize};
    my $zgap_start     = $opts{zgap_start};
    my $zgap_extension = $opts{zgap_extension};
    my $zchaining      = $opts{zchaining};
    my $zthreshold     = $opts{zthreshold};
    my $zmask          = $opts{zmask};

    my $seq = $opts{seq};
    my $blastable = $opts{blastable};

    # add jobs
    my @dsg_ids = split( /,/, $blastable );

    CoGe::Accessory::Web::write_log( "process $$", $cogeweb->logfile );

    my ( $fasta_file, $query_seqs_info ) = create_fasta_file($seq, $cogeweb);
    my $pre_command;
    my $x;
    ( $x, $pre_command ) = CoGe::Accessory::Web::check_taint($pre_command);
    my @results;
    my $count = 1;

    foreach my $dsgid (@dsg_ids) {
        my ( $org, $dbfasta, $dsg ) = get_blast_db($dsgid, $db);
        next unless $dbfasta;
        next unless -s $fasta_file;

        my $name = $dsg->organism->name;
        my $dbpath = File::Spec->catdir($BLASTDBDIR, $dsgid);
        my $BLASTDB = $config->{MAKEBLASTDB} || "makeblastdb";
        my $cmd;
        $cmd = "mkdir $dbpath && " unless -e $dbpath;
        $cmd .= "cd $dbpath && " .
                "$BLASTDB -in $dbfasta -out $dsgid -dbtype nucl -title " . qq{"$name"} . " -logfile db.log && " .
                "touch done";
        $workflow->add_job({
            cmd => $cmd,
            script  => undef,
            args    => [],
            inputs  => undef,
            outputs => ["$dbpath/done"],
            description => "Generating blastable database"
        });

        my $outfile = $cogeweb->basefile . "-$count.$program";

        $cmd = $BLAST_PROGS->{$program};
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
            push @$args, [ '', "E=" . $zgap_extension, 1 ] if defined $zgap_extension;
            push @$args, [ '', '--ambiguous=iupac', 1 ];
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
            push @$args, [ '-penalty', $nuc_penalty, 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-reward', $nuc_reward, 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-gapopen', $exist, 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-gapextend', $extent, 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-dust', 'no', 1 ] unless $program =~ /tblast/i;
            push @$args, [ '-seg', $filter_query ? 'yes' : 'no', 1 ] if $program =~ /tblast/i;
            push @$args, [ '-outfmt', $outfmt, 1] if $outfmt;
            push @$args, [ '-query',     $fasta_file, 1 ];
            push @$args, [ '-word_size', $wordsize,   1 ];
            push @$args, [ '-evalue',    $expect,     1 ];
            push @$args, [ '-db',        File::Spec->catdir($dbpath, $dsgid),   0 ];
            push @$args, [ '>',          $outfile,    1 ];
        }

        $workflow->add_job({
            cmd     => "/usr/bin/nice",
            args    => $args,
            inputs  => [$fasta_file, "$dbpath/done"],
            outputs => [$outfile],
            description => "Blasting sequence against $name"
        });

        if ($opts{link_results}) {
            my $download_path = get_download_path("jobs", $user ? $user->name : 'public', $workflow->id);
            make_path($download_path);
            my $filename = sanitize_name("$org.$program");
            my $outfile_link = catfile($download_path, $filename);
            $workflow->add_job({
                cmd     => "ln -s $outfile \"$outfile_link\"",
                inputs  => [$outfile],
                outputs => [$outfile_link],
                description => "Linking output to downloads"
            });
            $workflow->add_job(add_workflow_result(
                username => $user ? $user->name : 'public',
                wid      => $workflow->id,
                result   => {
                    type => 'url',
                    path => download_url_for(
                        wid => $workflow->id,
                        file => $filename
                    )
                },
                dependency => $outfile_link
            ));
        }

        $count++;
    }
}

sub build {
    my $self = shift;
    my $program = $self->params->{program};
    $program = 'mega' if $program eq 'megablast';
    $program = 'dcmega' if $program eq 'discontinuous_megablast';

    my @gids = get_genomes($self->params->{genomes}, $self->params->{notebooks}, $self->db);

    # set defaults for missing options
    my $expect = $self->params->{e_value};
    if ($program ne 'lastz') {
        $expect = 0.001 unless defined $expect;
    }
    my $wordsize = $self->params->{word_size};
    if ($program eq 'tblastn' || $program eq 'tblastx') {
        $wordsize = 3 unless defined $wordsize;
    }
    elsif ($program eq 'dcmega') {
        $wordsize = 11 unless defined $wordsize;
    }
    else {
        $wordsize = 8 unless defined $wordsize;
    }
    my $gap_costs = $self->params->{gap_costs};
    if ($program eq 'tblastn' || $program eq 'tblastx') {
        $gap_costs = [11, 1] unless $gap_costs;
    }
    elsif ($program ne 'lastz') {
        $gap_costs = [5, 2] unless $gap_costs;
    }
    my $match_score = $self->params->{match_score};
    if ($program eq 'blastn' || $program eq 'mega' || $program eq 'dbmega') {
        $match_score = [1, -2] unless $match_score;
    }
    my $filter_query = $self->params->{filter_query};
    if ($program eq 'tblastn' || $program eq 'tblastx') {
        $filter_query = 1 unless defined $filter_query;
    }

    return 0 if ! scalar @gids;

    my $resp = add_jobs(
        basename     => $self->params->{basename},
        blastable    => join(',', @gids),
        config       => $self->conf,
        db           => $self->db,
        expect       => $expect,
        filter_query => $filter_query,
        gapcost      => $gap_costs->[0] . ' ' . $gap_costs->[1],
        link_results => 1,
        matchscore   => $match_score->[0] . ',' . $match_score->[1],
        matrix       => $self->params->{matrix},
        outfmt       => $self->params->{outfmt},
        program      => $program,
        resultslimit => $self->params->{max_results},
        seq          => $self->params->{query_seq},
        user         => $self->user,
        wordsize     => $wordsize,
        workflow     => $self->workflow,
        zchaining    => $self->params->{zchaining},
        zgap_extension => $self->params->{zgap_extension},
        zgap_start   => $self->params->{zgap_start},
        zmask        => $self->params->{zmask},
        zthreshold   => $self->params->{zthreshold},
        zwordsize    => $self->params->{zwordsize}
    );
    my %opts = %{$self->params};
    $opts{'genomes'} = join(',', @gids);
    $opts{'gapcost'} = $gap_costs->[0] . ',' . $gap_costs->[1];
    $opts{'match_score'} = $match_score->[0] . ',' . $match_score->[1];

    #TODO move into pre_build -- mdb 10/12/16
	$self->{site_url} = get_tiny_url( %opts );
    if ($self->params->{requester}) { # request is from internal web page - external API requests will not have a 'requester' field
        my $page = $self->params->{requester}->{page}; # page name used for logging
        $self->page($page) if $page;
    }

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
    CoGe::Accessory::Web::write_log( "creating user's fasta file", $cogeweb->logfile );
    open( NEW, "> " . $cogeweb->basefile . ".fasta" );
    print NEW $seq;
    close NEW;
    return $cogeweb->basefile . ".fasta", \%seqs;
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

sub get_genomes {
    my $genomes = shift;
    my $notebooks = shift;
    my $db = shift;

    my @gids;
    @gids = @$genomes if $genomes;
    if ($notebooks) {
        for (@$notebooks) {
            for (@{$db->resultset("List")->find($_)->genomes}) {
                push @gids, $_->id;
            }
        }
    }
    return @gids;
}

sub get_name {
	my $self = shift;

    if ($self->request->user) {
        my @gids = @{$self->params->{genomes}};
        my $genomes_url = CoGe::Accessory::Web::get_tiny_link(
            user_id => $self->user->id,
            page    => "GenomeList",
            url     => $self->conf->{SERVER} . "GenomeList.pl?dsgid=" . join(',', @gids)
        );

        my $list_link =
            qq{<a href="$genomes_url" target_"blank">}
        . @gids
        . ' genome'
        . ( @gids > 1 ? 's' : '' ) . '</a>';
        return 'Blast ' . length($self->params->{query_seq}) . ' characters against ' . $list_link;
    }
    return 'CoGeBlast';
}

sub get_tiny_url {
    my %opts = @_;
    my %params = (
        basename     => $opts{basename},
        color_hsps   => $opts{color_hsps},
        comp         => $opts{comp},
        dsgid        => $opts{genomes},
        expect       => $opts{e_value},
        fid          => $opts{fid},
        filter_query => $opts{filter_query},
        gapcost      => $opts{gapcost},
        job_title    => $opts{job_title},
        match_score  => $opts{match_score},
        matrix       => $opts{matrix},
        outfmt       => $opts{outfmt},
        program      => $opts{program},
        resultslimit => $opts{max_results},
        type         => $opts{type},
        wordsize     => $opts{wordsize},
        zchaining    => $opts{zchaining},
        zgap_exten   => $opts{zgap_extension},
        zgap_start   => $opts{zgap_start},
        zmask        => $opts{zmask},
        zthreshold   => $opts{zthreshold},
        zwordsize    => $opts{zwordsize}
    );
    my $url = url_for("CoGeBlast.pl", %params);
    return CoGe::Accessory::Web::get_tiny_link(url => $url);
}

# sub go {
#     my %opts = @_;
#     my $db = $opts{db};
#     my $user = $opts{user};
#     my $config = $opts{config};
#     my $JEX = CoGe::Accessory::Jex->new( host => $config->{JOBSERVER}, port => $config->{JOBPORT} );
#     my $blastable = $opts{blastable};
#     my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
#         basename => $opts{basename},
#         tempdir  => $config->{TEMPDIR} . "CoGeBlast"
#     );

#     my $tiny_url = get_tiny_url(%opts);
#     my $workflow = $JEX->create_workflow(
#         name    => 'cogeblast-' . ($tiny_url =~ /\/(\w+)$/),
#         id      => 0,
#         logfile => $cogeweb->logfile
#     );

#     add_jobs(
#         cogeweb => $cogeweb,
#         workflow => $workflow,
#         %opts
#     );

#     my $response = $JEX->submit_workflow($workflow);

#     my $genomes_url = CoGe::Accessory::Web::get_tiny_link(
#         user_id => $user->id,
#         page    => "GenomeList",
#         url     => $config->{SERVER} . "GenomeList.pl?dsgid=$blastable"
#     );

#     my @dsg_ids = split( /,/, $blastable );
#     my $list_link =
#         qq{<a href="$genomes_url" target_"blank">}
#       . @dsg_ids
#       . ' genome'
#       . ( @dsg_ids > 1 ? 's' : '' ) . '</a>';
#     my $log_msg = 'Blast ' . length($opts{seq}) . ' characters against ' . $list_link;

#     my $log = CoGe::Accessory::Web::log_history(
#         db          => $db,
#         user_id     => $user->id,
#         page        => 'CoGeBlast',
#         description => $log_msg,
#         link        => $tiny_url,
#         parent_id   => $response->{id},
#         parent_type => 7 #FIXME magic number
#     ) if $response and $response->{id};

#     return encode_json({
#         id => $response->{id},
#         link => $tiny_url,
#         logfile => $config->{TEMPURL} . "CoGeBlast/" . $cogeweb->basefilename . ".log",
#         success => $JEX->is_successful($response) ? JSON::true : JSON::false
#     })
# }

1;
