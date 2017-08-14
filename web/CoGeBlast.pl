#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings('redefine');

use CoGeX;
use CoGe::JEX::Jex;
use CoGe::Accessory::Web qw(url_for get_command_path);
use CoGe::Accessory::Utils qw( commify get_link_coords );
use CoGe::Accessory::blast_report;
use CoGe::Accessory::blastz_report;
use CoGe::Builder::Tools::CoGeBlast qw( create_fasta_file get_blast_db get_tiny_url );
use CoGe::Core::Notebook qw(notebookcmp);
use CoGe::Core::Genome qw(fix_chromosome_id genomecmp genomecmp2);
use CoGe::Graphics::GenomeView;
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature::HSP;
use DBI;
use CGI;
use JSON::XS;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use LWP::Simple;
use LWP::Simple::Post qw(post post_xml);
use URI::Escape;
use POSIX;
use File::Temp;
use File::Path;
use File::Spec;
use Spreadsheet::WriteExcel;
use Benchmark qw(:all);
use Parallel::ForkManager;
use Sort::Versions;

use vars qw($P $PAGE_NAME $TEMPDIR $TEMPURL $DATADIR $FASTADIR
  $FORMATDB $FORM $USER $LINK $db $cogeweb $RESULTSLIMIT
  $MAX_PROC $MAX_SEARCH_RESULTS %FUNCTION $PAGE_TITLE $JEX $EMBED);

$PAGE_TITLE = "CoGeBlast";
$PAGE_NAME  = $PAGE_TITLE . ".pl";

$RESULTSLIMIT       = 100;
$MAX_SEARCH_RESULTS = 1000;

$FORM = new CGI;
( $db, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    page_title => $PAGE_TITLE,
    cgi => $FORM
);

$JEX           = CoGe::JEX::Jex->new( host => $P->{JOBSERVER}, port => $P->{JOBPORT} );
$ENV{PATH}     = $P->{COGEDIR};
$ENV{BLASTDB}  = $P->{BLASTDB};
$ENV{BLASTMAT} = $P->{BLASTMATRIX};
$TEMPDIR       = $P->{TEMPDIR} . "CoGeBlast";
$TEMPURL       = $P->{TEMPURL} . "CoGeBlast";
$DATADIR       = $P->{DATADIR};
$FASTADIR      = $P->{FASTADIR};
$FORMATDB      = get_command_path('FORMATDB');

%FUNCTION = (
    source_search            => \&get_data_source_info_for_accn,
    get_types                => \&get_types,
    get_sequence             => \&get_sequence,
    get_url                  => \&get_url,
    check_seq                => \&check_seq,
    set_seq                  => \&set_seq,
    blast_param              => \&blast_param,
    database_param           => \&database_param,
    gen_dsg_menu             => \&gen_dsg_menu,
    generate_feat_info       => \&generate_feat_info,
    get_hsp_info             => \&get_hsp_info,
    generate_overview_image  => \&generate_overview_image,
    overlap_feats_parse      => \&overlap_feats_parse,
    get_nearby_feats         => \&get_nearby_feats,
    export_fasta_file        => \&export_fasta_file,
    export_CodeOn            => \&export_CodeOn,
    export_to_excel          => \&export_to_excel,
    export_top_hits          => \&export_top_hits,
    generate_tab_deliminated => \&generate_tab_deliminated,
    generate_feat_list       => \&generate_feat_list,
    generate_blast           => \&generate_blast,
    export_hsp_info          => \&export_hsp_info,
    export_hsp_query_fasta   => \&export_hsp_query_fasta,
    export_hsp_subject_fasta => \&export_hsp_subject_fasta,
    export_alignment_file    => \&export_alignment_file,
    generate_basefile        => \&generate_basefile,
    get_dsg_for_menu         => \&get_dsg_for_menu,
    get_genome_info          => \&get_genome_info,
    search_lists             => \&search_lists,
    get_list_preview         => \&get_list_preview,
    get_genomes_for_list     => \&get_genomes_for_list,
    get_orgs                 => \&get_orgs,
    read_log                 => \&CoGe::Accessory::Web::read_log,
    save_settings            => \&save_settings,
    get_results              => \&get_results,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $template;
    
    $EMBED = $FORM->param('embed') || 0;
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    
        $template->param( TITLE => 'CoGeBLAST: Perform BLAST Analysis',
                          PAGE_TITLE => 'BLAST',
        				  PAGE_LINK  => $LINK,
        				  SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
        				  HOME       => $P->{SERVER},
                          HELP       => 'CoGeBlast',
                          WIKI_URL   => $P->{WIKI_URL} || '' );
        $template->param( USER => $USER->display_name || '' );
        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
        $template->param( ADMIN_ONLY => $USER->is_admin,
                          CAS_URL    => $P->{CAS_URL} || '',
                          COOKIE_NAME => $P->{COOKIE_NAME} || '' );
    }
    
    $template->param( BODY => gen_body() );
    return $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeBlast.tmpl' );
    my $form   = $FORM;
    my $featid = $form->param('featid');
    my $fid = $form->param('fid');
    $featid .= ",$fid" if $fid;
    $featid =~ s/^,// if $featid;
    $featid =~ s/,$// if $featid;
    my $chr    = $form->param('chr') || 0;
    my $upstream   = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $dsid       = $form->param('dsid') || 0;
    my $dsgid      = $form->param('dsgid') || 0;
    my $gid        = $form->param('gid');
    $dsgid = $gid if $gid;    # alias
    my $lid   = $form->param('lid');
       $lid   = $form->param('nid') unless defined $lid;
    my $rc    = $form->param('rc')    || 0;
    my $seq   = $form->param('seq');
    my $gstid = $form->param('gstid') || 1;
    my $locations = $form->param('locations');

    my $prefs = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $db
    );
    $prefs = {} unless $prefs;

    $template->param(
        MAIN       => 1,
        UPSTREAM   => $upstream,
        DOWNSTREAM => $downstream,
        DSID       => $dsid,
        DSGID      => $dsgid,

        #ORG_LIST		=> get_orgs(html_only => 1),
        RC        => $rc,
        FEATID    => $featid,
        CHR       => $chr,
        GSTID     => $gstid,
        LOCATIONS => $locations
    );

    #$template->param( BETA => 1);

    if ($featid) {
        $template->param( SEQVIEW => 1 );
    }
    elsif ($chr) {
        $template->param( SEQVIEW => 2 );
    }
    elsif ($seq) {
        unless ( $seq =~ />/ ) {
            $columns = 80;
            $seq = wrap( "", "", $seq );
        }
        $template->param(
            SEQVIEW  => 0,
            SEQUENCE => $seq
        );
    }
    elsif ($locations) {
        $template->param(
            SEQVIEW  => 0,
            SEQUENCE => ""
        );
    }
    else {
        $template->param(
            SEQVIEW  => 0,
            SEQUENCE => 'Enter FASTA sequence(s) here'
        );
    }
    $template->param( PAGE_NAME => $PAGE_NAME );
    $template->param( USER_NAME => $USER->user_name );

    #$template->param( REST      => 1 );

    #populate user specified default values
    my $db_list;
    if ( $prefs->{dsgids} ) {
        my $id = get_dsg_for_menu( dsgid => $prefs->{dsgids} );
        $db_list .= qq{add_to_list('$id');} if ($id);
    }

    if ($dsgid) {
        my $id = get_dsg_for_menu( dsgid => $dsgid );
        $db_list .= qq{add_to_list('$id');} if ($id);
    }

    if ($lid) {
        my $list = $db->resultset('List')->find($lid);
        if ($list) {
            my $dsgids = join( ',', map { $_->id } $list->genomes );
            my $id = get_dsg_for_menu( dsgid => $dsgids );
            $db_list .= qq{add_to_list('$id');} if $id;
        }
    }

#set up which columns of results will be displayed by default
#	my $prefs = CoGe::Accessory::Web::load_settings(user=>$USER, page=>$PAGE_NAME);
    unless ( $prefs->{display} && ref( $prefs->{display} ) eq "HASH" ) {
        $prefs->{display} = {
            NumD  => 1,
            EvalD => 1,
            QualD => 1,
            FeatD => 1,
            ChrD  => 1,
            OrgD  => 1,
            QND   => 1,
            PosD  => 1
        };

        #my %orgs;
        #map {$orgs{$_->{HSP_ORG}}++} @hsp;
        #$prefs->{display}{OrgD}=1 if scalar (keys %orgs) >1;
        #my %qseqs;
        #map {$qseqs{$_->{QUERY_SEQ}}++} @hsp;
        #$prefs->{display}{QND}=1 if scalar (keys %qseqs) >1;
        #my %chrs;
        #map {$chrs{$_->{HSP_CHR}}++} @hsp;
        #$prefs->{display}{ChrD}=1 if scalar (keys %qseqs) >1;
    }
    foreach my $key ( keys %{ $prefs->{display} } ) {
        $template->param( $key => "checked" ) if $key;
    }
    $template->param( SAVE_DISPLAY => 1 ) unless $USER->user_name eq "public";

    $template->param( document_ready => $db_list ) if $db_list;
    $template->param( TEMPDIR => $TEMPDIR );
    my $resultslimit = $RESULTSLIMIT;
    $resultslimit = $prefs->{'resultslimit'} if $prefs->{'resultslimit'};
    $template->param( RESULTSLIMIT => $resultslimit );
    $template->param( SAVE_ORG_LIST => 1 ) unless $USER->user_name eq "public";
    return $template->output;
}

sub get_sequence {
    my %opts       = @_;
    my $fids       = $opts{fid};
    my $dsid       = $opts{dsid};
    my $dsgid      = $opts{dsgid};
    my $chr        = $opts{chr};
    my $start      = $opts{start};
    my $stop       = $opts{stop};
    my $blast_type = $opts{blast_type};
    my $upstream   = $opts{upstream};
    my $downstream = $opts{downstream};
    my $locations  = $opts{locations};
    my $gstid      = $opts{gstid};
    my $rc         = $opts{rc};
    my $fasta;
    my $prot = $blast_type =~ /blast_type_p/ ? 1 : 0;

    if ($fids) {
        foreach my $fid ( split( /,/, $fids ) ) {
            my $gstidt;
            if ( $fid =~ /_/ ) {
                ( $fid, $gstidt ) = split( /_/, $fid );
            }
            else {
                $gstidt = $gstid;
            }
            my $feat = $db->resultset('Feature')->find($fid);
            $fasta .=
              ref($feat) =~ /Feature/i
              ? $feat->fasta(
                prot       => $prot,
                rc         => $rc,
                upstream   => $upstream,
                downstream => $downstream,
                gstid      => $gstidt
              )
              : ">Unable to find feature id $fid\n";
        }
    }
    elsif ($dsid) {
        my $ds = $db->resultset('Dataset')->find($dsid);
        $fasta =
          ref($ds) =~ /dataset/i
          ? $ds->fasta(
            start => $start,
            stop  => $stop,
            chr   => $chr,
            prot  => $prot,
            rc    => $rc,
            gstid => $gstid
          )
          : ">Unable to find dataset id $dsid";
    }
    elsif ($dsgid) {
        my $dsg = $db->resultset('Genome')->find($dsgid);
        $fasta =
          ref($dsg) =~ /genome/i
          ? $dsg->fasta(
            start => $start,
            stop  => $stop,
            chr   => $chr,
            prot  => $prot,
            rc    => $rc
          )
          : ">Unable to find genome id $dsgid";
    }
    elsif ($locations) {
        my %genomes;
        foreach ( split( ',', $locations ) ) {
            my ( $gid, $chr, $start, $stop ) = split( ':', $_ );
            next if ( $gid   =~ /\D/ );
            next if ( $start =~ /\D/ );
            next if ( $stop  =~ /\D/ );
            ( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );
            $start -= $upstream;
            $stop += $downstream;
            if ( not defined $genomes{$gid} ) {
                $genomes{$gid} = $db->resultset('Genome')->find($gid);
            }
            my $genome = $genomes{$gid};
            return "Unable to find genome for $gid" unless $genome;
            return "Restricted Access"
              if !$USER->has_access_to_genome($genome);
            $fasta .=
              ref($genome) =~ /genome/i
              ? $genome->fasta(
                start => $start,
                stop  => $stop,
                chr   => $chr,
                prot  => $prot,
                rc    => $rc
              )
              : ">Unable to find genome id $gid";
        }
    }
    return $fasta;
}

sub get_url {
    my %opts    = @_;
    my $program = $opts{program};

    if ( $program eq "blastn" ) {
        return
"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=on&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";
    }
    elsif ( $program eq "blastp" ) {
        return
"http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";
    }
    elsif ( $program eq "blastx" ) {
        return
"http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=blastx&BLAST_PROGRAMS=blastx&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";
    }
    elsif ( $program eq "tblastn" ) {
        return
"http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastn&BLAST_PROGRAMS=tblastn&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";
    }
    elsif ( $program eq "tblastx" ) {
        return
"http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastx&BLAST_PROGRAMS=tblastx&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";
    }
    else {
        return 1;
    }
}

sub blast_param {
    my %opts      = @_;
    my $seq_type  = $opts{blast_type} || "blast_type_n";
    my $translate = $opts{translate};
    my $version   = $opts{version};
    my $pro;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeBlast.tmpl' );

    if ( $seq_type =~ "blast_type_n" ) {
        if ( $version && $version =~ /coge_radio/ ) {
            $template->param( BLAST_NU => 1 );
        }
        else { $template->param( NCBI_BLAST_NU => 1 ); }
    }
    else {
        $pro = 1;
        if ( $version && $version =~ /coge_radio/ ) {
            $template->param( BLAST_PRO => 1 );
        }
        else { $template->param( NCBI_BLAST_PRO => 1 ); }

        unless ($translate) {
            if ( $version =~ /coge_radio/ ) {
                $template->param( BLAST_PRO_COMP => 1 );
            }
            else { $template->param( NCBI_BLAST_PRO_COMP => 1 ); }
        }
    }

    my $html = $template->output;
    return encode_json( { html => $html, version => $version, pro => $pro } );
}

sub database_param {
    my %opts = @_;
    my $program = $opts{program} || "blastn";

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeBlast.tmpl' );
    if ( $program eq "blastn" ) { $template->param( NU_DB => 1 ); }
    elsif ( ( $program eq "blastp" ) || ( $program eq "blastx" ) ) {
        $template->param( PRO_DB => 1 );
    }
    else { $template->param( T_DB => 1 ); }
    my $html = $template->output;
    return $html;
}

sub get_orgs {
    my %opts      = @_;
    my $name_desc = $opts{name_desc};
    my $html_only =
      $opts{html_only};    # optional flag to return html instead of JSON
    my $timestamp = $opts{timestamp};

    #	print STDERR "get_orgs: " . ($name_desc ? $name_desc : '') . "\n";

    my $html;
    my @organisms;
    if ($name_desc) {
        $name_desc = '%' . $name_desc . '%';
        my @organisms = $db->resultset("Organism")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $name_desc ],
                [ 'description', $name_desc ]
            ]
        );

        my @opts;
        foreach
          my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @organisms )
        {
            push @opts,
                "<OPTION value=\""
              . $item->id
              . "\" id=\"o"
              . $item->id . "\">"
              . $item->name
              . "</OPTION>";
        }

#$html .= qq{<FONT class="small" id="org_count">Matching Organisms (} . scalar @opts . qq{)</FONT><BR>};
        if (@opts) {
            if ( @opts <= $MAX_SEARCH_RESULTS ) {

#$html .= qq{<SELECT MULTIPLE id="org_id" SIZE="8" style="min-width:200px;" onchange="gen_dsg_menu(['args__oid','org_id'],['dsgid']);" ondblclick="add_selected_orgs();">\n};
                $html .= join( "\n", @opts );

                #$html .= "\n</SELECT>\n";
                #$html =~ s/OPTION/OPTION SELECTED/;
            }
            else {
                $html .=
"<option id='null_org' style='color:gray;' disabled='disabled'>Too many results to display, please refine your search.</option>";

#$html .= qq{<SELECT MULTIPLE id="org_id" SIZE="8" style="min-width:200px;"><option id='null_org' style='color:gray;'>Too many results to display, please refine your search.</option></SELECT><input type='hidden' id='gstid'>\n};
            }
        }
        else {
            $html .=
"<option id='null_org' style='color:gray;' disabled='disabled'>No results</option>";

#$html .= qq{<SELECT MULTIPLE id="org_id" SIZE="8" style="min-width:200px;"><option id='null_org' style='color:gray;'>No results</option></SELECT><input type='hidden' id='gstid'>\n};
        }
    }
    else {

#		$html .= qq{<FONT class="small" id="org_count">Matching Organisms } . '(' . $db->resultset('Organism')->count() . ')' . qq{</FONT>\n<br>\n};
#		$html .= qq{<SELECT MULTIPLE id="org_id" SIZE="8" style="min-width:200px;"><option id='null_org' style='color:gray;'>Please enter a search term</option></SELECT><input type='hidden' id='gstid'>\n};
        $html .=
"<option id='null_org' style='color:gray;' disabled='disabled'>Please enter a search term</option>";
    }

    return $html if ($html_only);
    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub gen_dsg_menu {
    #my $t1    = new Benchmark;
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $dsgid = $opts{dsgid};
    #print STDERR "gen_dsg_menu: $oid " . (defined $dsgid ? $dsgid : '') . "\n";

    my $favorites = CoGe::Core::Favorites->new(user => $USER);

    my @genomes;
    foreach my $dsg (
        sort { genomecmp2($a, $b, $favorites) }
            $db->resultset('Genome')->search(
                { organism_id => $oid },
                { prefetch    => ['genomic_sequence_type'] }
            )
      )
    {
        next unless $USER->has_access_to_genome($dsg);
        next if $dsg->deleted;

        $dsgid = $dsg->id unless $dsgid;
        my $name;
        $name .= "&#11088; " if ($favorites->is_favorite($dsg));
        $name .= "&#x2705; " if $dsg->certified;
        $name .= "&#x1f512; " if $dsg->restricted; 
        $name .= join( ", ", map { $_->name } $dsg->source ) . ": ";
	    $name .= " (id ". $dsg->id.") ";
        $name .= $dsg->name . ", " if $dsg->name; # : $dsg->datasets->[0]->name;
        $name .= "v"
          . $dsg->version . " "
          . $dsg->type->name . " "
          . commify( $dsg->length ) . "nt";

        push @genomes, [ $dsg->id, $name ];
    }
    my $size = scalar @genomes;
    $size = 5 if $size > 5;

    my $dsg_menu = '';
    if (@genomes) {
#$dsg_menu .= 'Genomes for Organism<br>';
#$dsg_menu .= qq{<select multiple id='dsgid' size='$size' onclick="show_add();" ondblclick="get_dsg_for_menu(['args__dsgid','dsgid'],[add_to_list]);">};
        foreach (@genomes) {
            my ( $numt, $name ) = @$_;
            my $selected = ( $dsgid && $numt == $dsgid ? 'selected' : '' );
            $dsg_menu .= qq{<option value='$numt' $selected>$name</option>};
        }

        #$dsg_menu .= '</select>';
    }

    #my $t2 = new Benchmark;
    #my $time = timestr( timediff( $t2, $t1 ) );
    return $dsg_menu;
}

sub get_dsg_for_menu {
    my %opts   = @_;
    my $dsgids = $opts{dsgid};
    my $orgids = $opts{orgid};
    my %dsgs;

#	print STDERR "get_dsg_for_menu: dsgids=" . ($dsgids ? $dsgids : '') . " orgids=" . ($orgids ? $orgids : '') . "\n";

    if ($orgids) {
        my @orgids = split( /,/, $orgids );
        foreach my $dsg (
            $db->resultset('Genome')->search( { organism_id => [@orgids] } ) )
        {
            next unless $USER->has_access_to_genome($dsg);
            next if $dsg->deleted;

            $dsgs{ $dsg->id } = $dsg;
        }
    }

    if ($dsgids) {
        %dsgs = () if ( $dsgs{$dsgids} );
        foreach my $dsgid ( split( /,/, $dsgids ) ) {
            my $dsg = $db->resultset('Genome')->find($dsgid);
            next unless $USER->has_access_to_genome($dsg);
            next if $dsg->deleted;
            $dsgs{ $dsg->id } = $dsg;
        }
    }

    my $html;
    foreach my $dsg ( values %dsgs ) {
        my ($ds) = $dsg->datasets;
        $html .= ":::" if $html;
        my $org_name = $dsg->organism->name;
        $org_name =~ s/'//g;
        $html .=
            $dsg->id . "::"
          . $org_name . " ("
          . "id " . $dsg->id . " "
          . $ds->data_source->name . " "
          . $dsg->type->name . " v"
          . $dsg->version . ")";
        $html .= ' -- certified' if $dsg->certified;
    }

    return $html;
}

sub generate_basefile {
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    return $cogeweb->basefilename;
}

sub alert {
    my $msg = shift;
    return encode_json({ error => $msg });
}

sub get_results {
    my %opts = @_;

    my $color_hsps = $opts{color_hsps};
    my $program    = $opts{program};
    my $resultslimit = $opts{max_results} || $RESULTSLIMIT;
    my $basename     = $opts{basename};

    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    my $seq = $opts{query_seq};

    my $genomes = $opts{'genomes'};

    my @dsg_ids = split( /,/, $genomes );

    my $width = $opts{width};

    my $genomes_url = CoGe::Accessory::Web::get_tiny_link(
        user_id => $USER->id,
        page    => "GenomeList",
        url     => $P->{SERVER} . "GenomeList.pl?dsgid=$genomes"
    );

#    my $list_link = qq{<a href="$genomes_url" target_"blank">} . @dsg_ids . ' genome' . ( @dsg_ids > 1 ? 's' : '' ) . '</a>';

#    my $log_msg = 'Blast ' . length($seq) . ' characters against ' . $list_link;

    my $tiny_url = get_tiny_url(%opts);
    my $workflow = $JEX->create_workflow(
        name    => 'cogeblast-' . ($tiny_url =~ /\/(\w+)$/),
        id      => 0,
        logfile => $cogeweb->logfile
    );

    CoGe::Accessory::Web::write_log( "process $$", $cogeweb->logfile );

    $width = 400 unless $width =~ /^\d+$/; # something wrong with how width is calculated in tmpl file -- # mdb what does this mean!?

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

        my $outfile = $cogeweb->basefile . "-$count.$program";

        push @results,
          {
            file     => $outfile,
            organism => $org,
            dsg      => $dsg
          };

        $count++;
    }

    my @completed;

    foreach my $item (@results) {
        next unless -r $item->{file};
        my $outfile = $item->{file};
        my $report =
          $outfile =~ /lastz/
          ? new CoGe::Accessory::blastz_report( { file => $outfile } )
          : new CoGe::Accessory::blast_report( { file => $outfile } );

        $item->{report} = $report;
        my $file = $report->file();
        $file =~ s/$TEMPDIR//;
        $file = $TEMPURL . "/" . $file;
        $item->{link} = $file;
        push @completed, $item;
    }

    initialize_sqlite();

    my ( $html, $click_all_links ) = gen_results_page(
        results         => \@completed,
        width           => $width,
        resultslimit    => $resultslimit,
        prog            => $program,
        color_hsps      => $color_hsps,
        query_seqs_info => $query_seqs_info,
        link            => $tiny_url
    );

    return encode_json ({
        html => $html,
        click_all_links => $click_all_links,
        success => JSON::true
    });
}

sub gen_results_page {
    my %opts            = @_;
    my $results         = $opts{results};
#    my $width           = $opts{width};
    my $resultslimit    = $opts{resultslimit};
    my $color_hsps      = $opts{color_hsps};
    my $prog            = $opts{prog};
    my $query_seqs_info = $opts{query_seqs_info};
    my $link            = $opts{link};
    my $click_all_links;
    my $null;
    my %hsp_count;
    my $length;
    my @check;
    my @hsp;

    my $t0 = new Benchmark;
    my ( $chromosome_data, $chromosome_data_large, $genomelist_link ) =
      generate_chromosome_images(
        results      => $results,
        resultslimit => $resultslimit,
        color_hsps   => $color_hsps
      );
    my $t1 = new Benchmark;

    my %query_hit_count;
    foreach my $set (@$results) {
        $hsp_count{ $set->{organism} } = 0
          unless $hsp_count{ $set->{organism} };
        if ( @{ $set->{report}->hsps() } ) {
            foreach my $hsp ( sort { $a->number <=> $b->number }
                @{ $set->{report}->hsps() } )
            {
                my $ttt1 = new Benchmark;
                my $dsg  = $set->{dsg};
                my ($chr) = $hsp->subject_name =~ /\|(\S*)/;
                $chr = $hsp->subject_name unless $chr;
                $chr =~ s/\s+$//;
                my ($ds) = $dsg->datasets( chr => $chr );
                my ($org) = $set->{organism};
                unless ($dsg && defined $chr && $ds) {
                    print STDERR "CoGeBlast::gen_results_page: ERROR, undefined genome/chr/ds for HSP\n"; 
                    next;
                }
                
                $hsp_count{$org}++;
                last if ( $hsp_count{$org} > $resultslimit );
                my $tt1 = new Benchmark;
                populate_sqlite( $hsp, $dsg->id, $org, $chr );
                my $tt2 = new Benchmark;
                my $sqlite_time = timestr( timediff( $tt2, $tt1 ) );

                #print STDERR "sqlite transaction time: $sqlite_time\n";
                my $id = $hsp->number . '_' . $dsg->id;
                $click_all_links .= $id . ',';

                #my $feat_link = qq{<span class="link" onclick="fill_nearby_feats('$id','true')">Click for Closest Feature</span>};
                my $feat_link = qq{<span>Loading...</span>};

                my $qname = $hsp->query_name;
                my $start = $hsp->subject_start;
                my $stop  = $hsp->subject_stop;
                ( $start, $stop ) = ( $stop, $start ) if ( $stop < $start );

                #can we extract a CoGe name for the sequence
                if ( $hsp->query_name =~ /Name: (.*), Type:/ ) {
                    $qname = $1;
                    ($qname) = split( /,/, $qname );
                    $qname =~ s/^\s+//;
                    $qname =~ s/\s+$//;
                    $qname =~ s/,$//;
                }
                $query_hit_count{$qname}{org}{$org}++;
                $query_hit_count{$qname}{orig_name} = $hsp->query_name;
                $query_hit_count{$qname}{color}     = $hsp->{color}
                  if $hsp->{color};
                my $rc = $hsp->strand =~ /-/ ? 1 : 0;
                my $seqview_link =
                  qq{<span class="small link" onclick="window.open('SeqView.pl?dsid=}
                  . $ds->id
                  . qq{;chr=$chr;start=}
                  . $hsp->subject_start
                  . ";stop="
                  . $hsp->subject_stop
                  . qq{;rc=$rc')" target=_new>SeqView</span>};
                push @hsp,
                  {
                    CHECKBOX => $id . "_"
                      . $chr . "_"
                      . $hsp->subject_start . "no",
                    ID        => $id,
                    QUERY_SEQ => $qname,
                    HSP_ORG   => $org,
                    HSP => qq{<span class="link" title="Click for HSP information" onclick="update_hsp_info('table_row$id');">}
                      . $hsp->number
                      . "</span>",
                    HSP_EVAL   => $hsp->pval,
                    HSP_LENGTH => $hsp->length,
                    COVERAGE   => sprintf( "%.1f", $hsp->query_coverage * 100 )
                      . "%",
                    HSP_PID       => $hsp->percent_id . "%",
                    HSP_SCORE     => $hsp->score,
                    HSP_POS_START => ( $hsp->subject_start ),
                    HSP_CHR       => $chr,
                    HSP_LINK      => $feat_link,
                    HSP_QUALITY   => sprintf( "%.1f", $hsp->quality ) . "%",
                    SEQVIEW       => $seqview_link,
                    LOC_VAL       => $dsg->id . ':'
                      . $chr . ':'
                      . $start . ':'
                      . $stop
                  };
            }
        }
        $hsp_count{ $set->{organism} } = scalar @{ $set->{report}->hsps() };
    }
    unless (@hsp) {
        $null = "null";
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeBlast.tmpl' );
    $template->param( RESULT_TABLE => 1 );

    $template->param( NULLIFY => $null ) if $null;
    my $hsp_limit_flag = 0;
    my $hsp_count;
    $hsp_count .=
qq{<table class="small ui-widget-content ui-corner-all"><tr><th>Query Seq<th>}
      . join "<th>", sort keys %hsp_count;
    my $class = "even";
    foreach my $query ( sort keys %query_hit_count ) {
        my $name = $query_hit_count{$query}{orig_name};
        $name =~ s/\s//g;
        my $out_name = $query;
        my $color;
        foreach my $item ( map { sprintf( "%x", $_ ) }
            @{ $query_hit_count{$query}{color} } )
        {
            $item = "0" . $item if length $item == 1;
            $color .= $item;
        }
        $out_name = qq{<span style="color:#} . $color . qq{">$out_name</span>}
          if $color;
        $hsp_count .=
            qq{<tr class="$class"><td>$out_name <span class = species>(}
          . $query_seqs_info->{$name}
          . "nt)</span>";
        foreach my $org ( sort keys %hsp_count ) {
            my $count =
                $query_hit_count{$query}{org}{$org}
              ? $query_hit_count{$query}{org}{$org}
              : 0;
            $hsp_count .= qq{<td align=center>$count};
        }
        $class = $class eq "even" ? "odd" : "even";
    }
    $hsp_count .= qq{<tr class="$class"><th>Total};
    foreach my $org ( sort keys %hsp_count ) {
        my $count = $hsp_count{$org};
        $count = "<span class=alert>$count</span>" if $count > $resultslimit;
        $hsp_count .= qq{<td align=center class=species>$count</td>};

    }
    $hsp_count .= "</table>";
    foreach my $org ( keys %hsp_count ) {
        if ( $hsp_count{$org} > $resultslimit ) {
            $hsp_count .=
"<span class=\"small alert\">Only top $resultslimit HSPs shown for $org.</span><br>";
            $hsp_limit_flag = 1;
        }
    }
    my $total = 0;
    map { $total += $_ } values %hsp_count;
    $hsp_count .= "<span class=\"small\">Total Number of Hits: $total</span>";
    $hsp_count .=
      "<span class=\"small alert\">All results are in the blast report.</span>"
      if $hsp_limit_flag;

    $template->param( HSP_COUNT => $hsp_count );

    if (@hsp) {
        @hsp =
          sort { $a->{HSP_ORG} cmp $b->{HSP_ORG} || $a->{HSP} cmp $b->{HSP} }
          @hsp
          if @hsp;
        $template->param( HSP_TABLE => 1 );
        $template->param( HSPS      => \@hsp );
    }
    my $hsp_results = $template->output;
    $template->param( HSP_TABLE       => 0 );
    $template->param( RESULT_TABLE    => 0 );
    $template->param( HSP_RESULTS     => $hsp_results );
    $template->param( CHROMOSOMES_IF  => 1 );
    $template->param( GENOMELIST_LINK => $genomelist_link );
    $template->param( CHROMOSOME_LOOP => $chromosome_data );
    my $chromosome_element = $template->output;
    $template->param( CHROMOSOMES_IF => 0 );

    if ( scalar(@$chromosome_data_large) > 0 ) {
        $template->param( CHROMOSOMES_LARGE_IF  => 1 );
        $template->param( CHROMOSOME_LOOP_LARGE => $chromosome_data_large );
        my $chr_large_element = $template->output;
        $template->param( CHROMOSOMES_LARGE_IF => 0 );
        $template->param( CHR_LARGE            => $chr_large_element );
    }
    $template->param( CHROMOSOMES   => $chromosome_element );
    $template->param( BLAST_RESULTS => 1 );
    $template->param( DATA_FILES =>
          gen_data_file_summary( prog => $prog, results => $results ) );

    $template->param(LINK => $link);
    my $html = $template->output;

    my $temp = qq{
<span style='padding-left:10px;'>
<span id='link_hidden' class='ui-button ui-corner-all ui-button-icon-right coge-button coge-button-right' onclick="\$('#link_hidden').hide(); \$('#link_shown').fadeIn();"><span class="ui-icon ui-icon-arrowreturnthick-1-w"></span>Link to this</span>
<span id='link_shown' style='display:none;' class='small infobox'>Use this link to return to this page at any time: <span class='link' onclick=window.open('$link');><b>$link</b></span></span>
</span>};

    my $t2      = new Benchmark;

    my $figure_time = timestr( timediff( $t1, $t0 ) );
    my $table_time  = timestr( timediff( $t2, $t1 ) );
    my $render_time = timestr( timediff( $t2, $t0 ) );
    my $benchmark   = qq{
Time to gen tables:              $table_time
Time to gen images:              $figure_time
Time to gen results:             $render_time
};
    CoGe::Accessory::Web::write_log( $benchmark, $cogeweb->logfile );
    return $html, $click_all_links;
}

sub gen_data_file_summary {
    my %opts    = @_;
    my $prog    = $opts{prog};
    my $results = $opts{results};
    my $html    = "<table style='border-top:solid 1px lightgray;'><tr>";
    $html .= qq{<td class="small top">Data Files};
    $html .= "<div class='xsmall'><span class='link' onClick=\"get_all_hsp_data();\">HSP Data</span></DIV>\n";
    $html .= "<div class='xsmall'><span class='link' onClick=\"get_query_fasta();\">Query HSP FASTA File</span></DIV>\n";
    if ( $prog eq "tblastn" ) {
        $html .= "<div class='xsmall'><span class='link' onClick=\"get_subject_fasta();\">Subject HSP Protein FASTA File</span></DIV>\n";
        $html .= "<div class='xsmall'><span class='link' onClick=\"get_subject_fasta('dna');\">Subject HSP DNA FASTA File</span></DIV>\n";
    }
    else {
        $html .= "<div class='xsmall'><span class='link' onClick=\"get_subject_fasta();\">Subject HSP FASTA File</span></DIV>\n";
    }
    $html .= "<div class='xsmall'><span class='link' onClick=\"get_alignment_file();\">Alignment File</span></DIV>\n";
    $html .= qq{<td class='small top'>Analysis Files};
    my $dbname = $TEMPURL . "/" . $cogeweb->basefilename . ".sqlite";
    $html .= "<div class=xsmall><a href=\"$dbname\" target=_new>SQLite DB file</a></div>\n";
    foreach my $item (@$results) {
        my $blast_file = $item->{link};
        my $org        = $item->{organism};
        $html .= qq{<div class='xsmall'><a href="$blast_file" target=_new>Blast file for $org</div>\n};
    }
    $html .= qq{<td class='small top'>Log File};
    my $logfile = $TEMPURL . "/" . $cogeweb->basefilename . ".log";
    $html .= "<div class='xsmall'><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
    $html .= qq{</table>};
}

sub generate_chromosome_images {
    my %opts           = @_;
    my $results        = $opts{results};
    my $hsp_type       = $opts{hsp_type} || "NA";
    my $width          = $opts{width} || 400;
    my $large_width    = $opts{large_width} || 3 * $width;
    my $resultslimit   = $opts{resultslimit};
    my $imagefile_name = $opts{filename} || "null";
    my $color_hsps     = $opts{color_hsps};
    my $height         = ( $width / 16 );
    my $large_height = ( $large_width / 16 ) <= 64 ? ( $large_width / 16 ) : 64;
    my $scale = $opts{scale} || 'linear';
    my %data;
    my $filename;
    my ( @data, @large_data, @no_data );
    my %hsp_count;
    my %query_seqs;

    #colorize!

    if ( $color_hsps eq "query" ) {
        foreach my $set (@$results) {
            map { $query_seqs{ $_->query_name } = 1 }
              @{ $set->{report}->hsps() };
        }
        my $colors = color_pallet( num_seqs => scalar keys %query_seqs );
        my $count = 0;
        foreach my $key ( sort keys %query_seqs ) {
            $query_seqs{$key} = $colors->[$count];
            $count++;
        }
    }
    foreach my $set (@$results) {
        my $org = $set->{organism};
        $org =~ s/\s+$//;
        $data{$org}{file} = $set->{link};
        $filename = $imagefile_name . "_*.png";
        if ( $imagefile_name ne "null" ) {
            if (`ls $filename 2>>/dev/null`) {
                $data{$org}{image} = $imagefile_name;
                $data{$org}{skip}  = 1;
                next;
            }
        }

        #colorize!
        my ( $max_quality, $min_quality, $range_quality );
        if ( $color_hsps eq "quality" ) {
            foreach my $hsp ( @{ $set->{report}->hsps() } ) {
                $max_quality = $hsp->quality unless defined $max_quality;
                $max_quality = $hsp->quality if $hsp->quality > $max_quality;
                $min_quality = $hsp->quality unless defined $min_quality;
                $min_quality = $hsp->quality if $hsp->quality < $min_quality;
            }
            if ( defined $min_quality && defined $max_quality ) {
                $max_quality = sprintf( "%.3f", log($max_quality) )
                  unless $max_quality == 0;
                $min_quality = sprintf( "%.3f", log($min_quality) )
                  unless $min_quality == 0;
                $range_quality = $max_quality - $min_quality;
                $data{$org}{extra} =
qq{<span class=small>Hits colored by Log Quality.  <span style="color:#AA0000">Min: $min_quality</span> <span style="color:#00AA00">Max: $max_quality</span></span>};
            }
        }
        my ( $max_identity, $min_identity, $range_identity );
        if ( $color_hsps eq "identity" ) {
            foreach my $hsp ( @{ $set->{report}->hsps() } ) {
                $max_identity = $hsp->percent_id unless defined $max_identity;
                $max_identity = $hsp->percent_id
                  if $hsp->percent_id > $max_identity;
                $min_identity = $hsp->percent_id unless defined $min_identity;
                $min_identity = $hsp->percent_id
                  if $hsp->percent_id < $min_identity;
            }
            if ( defined $max_identity && $min_identity ) {
                $range_identity = $max_identity - $min_identity;
                $max_identity   = sprintf( "%.1f", $max_identity );
                $min_identity   = sprintf( "%.1f", $min_identity );
                $data{$org}{extra} =
qq{<span class=small>Hits colored by Identity.  <span style="color:#AA0000">Min: $min_identity</span> <span style="color:#00AA00">Max: $max_identity</span></span>};
            }
        }

        $filename = $set->{link};
        if ( @{ $set->{report}->hsps() } ) {
            foreach my $hsp ( @{ $set->{report}->hsps() } ) {

                #first, initialize graphic
                $hsp_count{$org}++;
                next if $hsp_count{$org} > $resultslimit;
                my ($chr) = $hsp->subject_name =~ /\|(\S*)/;
                $chr = $hsp->subject_name unless $chr;
                $chr =~ s/\s+$//;
                $data{$org}{image} = new CoGe::Graphics::GenomeView(
                    {
                        color_band_flag   => 1,
                        image_width       => $width,
                        chromosome_height => $height
                    }
                ) unless $data{$org}{image};
                $data{$org}{large_image} = new CoGe::Graphics::GenomeView(
                    {
                        color_band_flag   => 1,
                        image_width       => $large_width,
                        chromosome_height => $large_height
                    }
                ) unless $data{$org}{large_image};
                my $dsg = $set->{dsg};
                $data{$org}{dsg} = $dsg;

                #add chromosome to graphic
                unless ( $data{$org}{chr}{$chr} ) {
#                    my $last_pos = $dsg->sequence_length($chr);
                    my $last_pos = $dsg->get_chromosome_length($chr);
                    $data{$org}{image}
                      ->add_chromosome( name => "Chr: $chr", end => $last_pos );
                    $data{$org}{chr}{$chr} = 1;
                }
                my $num = $hsp->number . "_" . $dsg->id;
                my $up = $hsp->strand eq "++" ? 1 : 0;
                my $color;
                if ( $color_hsps eq "identity" || $color_hsps eq "quality" ) {
                    my $relative;
                    $relative = sprintf( "%.3f",
                        1 -
                          ( $max_quality - log( $hsp->quality ) ) /
                          ($range_quality) )
                      if $range_quality && $color_hsps eq "quality";
                    $relative = sprintf( "%.3f",
                        1 -
                          ( $max_identity - $hsp->percent_id ) /
                          ($range_identity) )
                      if $range_identity && $color_hsps eq "identity";
                    $relative = 1 unless $relative;
                    my $other = 1 - $relative;
                    my $c1    = 200;
                    $c1 = 200 * ( $relative * 2 ) if $relative < .5;
                    my $c2 = 200;
                    $c2 = 200 * ( $other * 2 ) if $other < .5;
                    $color = [ $c2, $c1, 0 ];
                }
                elsif ( $color_hsps eq "query" ) {
                    $color = $query_seqs{ $hsp->query_name };
                }
                else { $color = [ 0, 200, 0 ]; }
                my $id = $hsp->number . "_" . $set->{dsg}->id;
                $data{$org}{image}->add_feature(
                    name     => $hsp->number,
                    start    => $hsp->sstart,
                    stop     => $hsp->sstop,
                    chr      => "Chr: $chr",
                    imagemap => qq/class="imagemaplink" title="HSP No. /
                      . $hsp->number
                      . qq/" onclick="update_hsp_info('table_row$id');""/,
                    up    => $up,
                    color => $color,
                );
                $hsp->{color} = $color if $color_hsps eq "query";
            }
        }

#let's reverse the order of the image features so that the first one is drawn on top of the latters;
        if ( $data{$org}{image} ) {
            while ( my ( $k, $v ) = each %{ $data{$org}{image}->features } ) {
                $v = [ reverse @$v ];
                $data{$org}{image}->features->{$k} = $v;
            }
        }
    }

    my $count = 1;
    my @dsgids;
    foreach my $org ( sort keys %data ) {
        if ( $data{$org}{image} ) {
            my $image_file;
            my $image_map;
            my $large_image_file;
            my $image_map_large;
            unless ( $data{$org}{skip} ) {
                my $x;
                $large_image_file =
                    $cogeweb->basefile . "_"
                  . $hsp_type
                  . "_$count"
                  . "_large.png";
                ( $x, $large_image_file ) =
                  CoGe::Accessory::Web::check_taint($large_image_file);
                $image_file =
                  $cogeweb->basefile . "_" . $hsp_type . "_$count.png";
                ( $x, $image_file ) =
                  CoGe::Accessory::Web::check_taint($image_file);
                $data{$org}{image}->generate_png( filename => $image_file );
                $image_map =
                  $data{$org}{image}->generate_imagemap(
                    mapname => $cogeweb->basefilename . "_" . $count );
                my $map_file = $cogeweb->basefile . "_$count.$hsp_type.map";
                ( $x, $map_file ) =
                  CoGe::Accessory::Web::check_taint($map_file);
                open( MAP, ">$map_file" );
                print MAP $image_map;
                close MAP;
                $data{$org}{image}->image_width($large_width);
                $data{$org}{image}->chromosome_height($large_height);
                $data{$org}{image}
                  ->generate_png( filename => $large_image_file );
                $image_map_large =
                  $data{$org}{image}
                  ->generate_imagemap( mapname => $cogeweb->basefilename . "_"
                      . $count
                      . "_large" );
                $map_file = $cogeweb->basefile . "_$count.$hsp_type.large.map";
                ( $x, $map_file ) =
                  CoGe::Accessory::Web::check_taint($map_file);
                open( MAP, ">$map_file" );
                print MAP $image_map_large;
                close MAP;
            }
            else {
                my $x;
                $image_file = $data{$org}{image} . "_$count.png";
                $image_map =
                  get_map( $cogeweb->basefile . "_$count.$hsp_type.map" );
                $large_image_file =
                  $data{$org}{image} . "_$count" . "_large.png";
                $image_map_large =
                  get_map( $cogeweb->basefile . "_$count.$hsp_type.large.map" );
                ( $x, $image_file ) =
                  CoGe::Accessory::Web::check_taint($image_file);
                ( $x, $large_image_file ) =
                  CoGe::Accessory::Web::check_taint($large_image_file);
            }

            $image_file       =~ s/$TEMPDIR/$TEMPURL/;
            $large_image_file =~ s/$TEMPDIR/$TEMPURL/;
            my $blast_link = qq{<span class="small link" onclick="window.open('GenomeInfo.pl?gid=}
              . $data{$org}{dsg}->id
              . qq{')">$org</span><br>};
              #. "<a class =\"small\" href=" . $data{$org}{file} . " target=_new>Blast Report</a> "; # mdb removed 8/7/14 - redundant with blast report in Analysis Files section
            $blast_link .= $data{$org}{extra} if $data{$org}{extra};
            $blast_link .= "<br>";
            push @dsgids, $data{$org}{dsg}->id;
            push @large_data,
              {
                DB_NAME_LARGE   => $blast_link,
                CHR_IMAGE_LARGE => "<img src=$large_image_file ismap usemap='#"
                  . $cogeweb->basefilename . "_"
                  . "$count"
                  . "_large' border=0>$image_map_large",
                IMAGE_ID_LARGE => $count,
              };
            push @data,
              {
                DB_NAME => $blast_link,
                CHR_IMAGE => "<img style=\"position:relative; z-index=1;\" src=$image_file ismap usemap='#"
                  . $cogeweb->basefilename . "_"
                  . "$count' border=0>$image_map",
                ENLARGE => "<a href='#' onClick='enlarge_picture_window($count);' class='small large_image'>Enlarge</a>"
              };    #, DIV_STYLE=>'style="position:absolute;"'};

            $count++;
        }
        else {
            push @no_data,
              {
                DB_NAME => "<a href="
                  . $data{$org}{file}
                  . " target=_new>No Hits: $org</a>",
                DIV_STYLE => 'class="small"'
              };
        }
    }
    my $genomelist_link =
      "GenomeList.pl?" . join( ";", map { "dsgid=$_" } @dsgids );

    return [ @data, @no_data ], \@large_data, $genomelist_link;
}

sub get_map {
    my $file = shift;
    my $map;
    open( IN, $file ) || die "$!";
    while (<IN>) {
        $map .= $_;
    }
    close IN;
    return $map;
}

sub generate_fasta {
    my %opts   = @_;
    my $dslist = $opts{dslist};
    my $file   = $opts{file};
    $file = $FASTADIR . "/$file" unless $file =~ /$FASTADIR/;
    CoGe::Accessory::Web::write_log( "creating fasta file.",
        $cogeweb->logfile );
    open( OUT, ">$file" ) || die "Can't open $file for writing: $!";
    foreach my $ds (@$dslist) {
        foreach my $chr ( sort $ds->get_chromosomes ) {
            my $title =
                $ds->organism->name . " (v"
              . $ds->version . ") "
              . "chromosome: $chr"
              . ", CoGe database id: "
              . $ds->id;
            $title =~ s/^>+/>/;
            CoGe::Accessory::Web::write_log( "adding sequence $title",
                $cogeweb->logfile );
            print OUT ">" . $title . "\n";
            print OUT $ds->get_genomic_sequence( chr => $chr ), "\n";
        }
    }
    close OUT;
    CoGe::Accessory::Web::write_log( "Completed fasta creation",
        $cogeweb->logfile );
    return 1 if -r $file;
    CoGe::Accessory::Web::write_log( "Error with fasta file creation",
        $cogeweb->logfile );
    return 0;
}

sub generate_blast_db {
    my %opts    = @_;
    my $fasta   = $opts{fasta};
    my $blastdb = $opts{blastdb};
    my $org     = $opts{org};
    return 1 if -r "$blastdb.nsq";
    my $command = $FORMATDB . " -p F";
    $command .= " -i '$fasta'";
    $command .= " -t '$org'";
    $command .= " -n '$blastdb'";
    CoGe::Accessory::Web::write_log(
        "creating blastdb for $org ($blastdb)\n$command",
        $cogeweb->logfile );
    `$command`;
    return 1 if -r "$blastdb.nsq";
    CoGe::Accessory::Web::write_log(
        "error creating blastdb for $org ($blastdb)",
        $cogeweb->logfile );
    return 0;
}

sub generate_feat_info {
    my %opts   = @_;
    my $featid = $opts{featid};
    $featid =~ s/^table_row//;

    my ( $hsp_num, $dsgid );
    ( $featid, $hsp_num, $dsgid ) = split( /_/, $featid );

    my ($dsg)  = $db->resultset('Genome')->find($dsgid);
    my ($feat) = $db->resultset("Feature")->find($featid);
    unless ( ref($feat) =~ /Feature/i ) {
        return "Unable to retrieve Feature object for id: $featid";
    }

    my $html = $feat->annotation_pretty_print_html( gstid => $dsg->type->id );
    return $html;
}

sub get_hsp_info {
    my %opts     = @_;
    my $hsp_id   = $opts{hsp_id};
    my $filename = $opts{blastfile};
    $filename =~ s/$TEMPDIR//;

    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbfile = $cogeweb->sqlitefile;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "" );

    $hsp_id =~ s/^table_row// if $hsp_id =~ /table_row/;
    $hsp_id =~ s/^\d+_// if $hsp_id =~ tr/_/_/ > 1;

    my (
        $hsp_num, $pval,   $pid,    $psim,      $score,
        $qgap,    $sgap,   $match,  $qmismatch, $smismatch,
        $strand,  $length, $qstart, $qstop,     $sstart,
        $sstop,   $qalign, $salign, $align,     $qname,
        $sname,   $dsgid,  $org,    $chr
    );

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});

    ( $hsp_num, $dsgid ) = split( /_/, $hsp_id );
    my $dsg   = $db->resultset('Genome')->find($dsgid);
    my $gstid = $dsg->type->id;
    $sth->execute($hsp_id) || die "unable to execute";
    while ( my $info = $sth->fetchrow_hashref() ) {
        $pval      = $info->{eval};
        $pid       = $info->{pid};
        $psim      = $info->{psim};
        $score     = $info->{score};
        $qgap      = $info->{qgap};
        $sgap      = $info->{sgap};
        $match     = $info->{match};
        $qmismatch = $info->{qmismatch};
        $smismatch = $info->{smismatch};
        $strand    = $info->{strand};
        $length    = $info->{length};
        $qstart    = $info->{qstart};
        $qstop     = $info->{qstop};
        $sstart    = $info->{sstart};
        $sstop     = $info->{sstop};
        $qalign    = $info->{qalign};
        $salign    = $info->{salign};
        $align     = $info->{align};
        $qname     = $info->{qname};
        $sname     = $info->{sname};
        $chr       = $info->{schr};
        $org       = $info->{org};
    }
    my $qlength = $qstop - $qstart;
    my ($ds) = $dsg->datasets( chr => $chr );
    my $dsid = $ds->id;

    my $query_name = "<pre>" . $qname . "</pre>";
    $query_name = wrap( '', '', $query_name );
    $query_name =~ s/\n/<br>/g;

    my $subject_name = "<pre>" . $sname . "</pre>";
    $subject_name = wrap( '', '', $subject_name );
    $subject_name =~ s/\n/<br>/g;
    my @table1 = (
        {
            HSP_PID_QUERY      => $pid,
            HSP_PSIM_QUERY     => $psim,
            HSP_GAP_QUERY      => $qgap,
            HSP_PID_SUB        => $pid,
            HSP_PSIM_SUB       => $psim,
            HSP_GAP_SUB        => $sgap,
            HSP_MATCH_QUERY    => $match,
            HSP_MISMATCH_QUERY => $qmismatch,
            HSP_MATCH_SUB      => $match,
            HSP_MISMATCH_SUB   => $smismatch,
            HSP_POSITION_QUERY => $qstart . "-" . $qstop,
            HSP_POSITION_SUB   => $sstart . "-" . $sstop,
        }
    );
    my @table2 = (
        {
            HSP_STRAND => $strand,
            HSP_EVAL   => $pval,
            HSP_SCORE  => $score,
            HSP_LENGTH => $length,
            HSP_CHR    => $chr,
        }
    );
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeBlast.tmpl' );
    $template->param( HSP_IF  => 1 );
    $template->param( HSP_NUM => $hsp_num );
    $template->param( HSP_QS  => \@table1 );
    $template->param( HSP_HSP => \@table2 );

    my $align_str = "";
    my $query_seq = $qalign;
    $query_seq = wrap( '', '', $query_seq );
    my @query = split( /\n/, $query_seq );
    $query_seq =~ s/[^atgc]//ig;
    $query_seq = wrap( '', '', $query_seq );
    $query_seq =~ s/\n/<br>/g;
    $query_seq =~ tr/atgc/ATGC/;
    $query_seq = qq{<pre>$query_seq</pre>};

    my $sub_seq = $salign;
    $sub_seq = wrap( '', '', $sub_seq );
    my @sub = split( /\n/, $sub_seq );
    $sub_seq =~ s/[^atgc]//ig;
    $sub_seq = wrap( '', '', $sub_seq );
    $sub_seq =~ s/\n/<br>/g;
    $sub_seq =~ tr/atgc/ATGC/;
    $sub_seq = qq{<pre>$sub_seq</pre>};

    my $alignment = $align;
    $alignment =~ s/ /\./g;
    $alignment = wrap( '', '', $alignment );
    my @align = split( /\n/, $alignment );

    for ( my $i = 0 ; $i < scalar(@sub) ; $i++ ) {
        $align_str .=
          $query[$i] . "<br>" . $align[$i] . "<br>" . $sub[$i] . "<br>";
    }
    $align_str =~ s/<br>$//;
    $align_str =~ s/\./ /g;
    $align_str = "<pre>$align_str</pre>";

    $template->param( QUERY_SEQ =>
qq{<span class="small link" onclick="show_seq('$query_seq','$query_name',1,'seqObj','seqObj','}
          . $qstart . "','"
          . $qstop
          . qq{')">Query Sequence Hit</span>} );
    my $rc = $strand =~ /-/ ? 1 : 0;
    $template->param( SUB_SEQ =>
qq{<span class="small link" onclick="window.open('SeqView.pl?dsid=$dsid;chr=$chr;start=$sstart;stop=$sstop;rc=$rc;gstid=$gstid')" target=_new>Subject Sequence Hit</span>}
    );

    $template->param( ALIGNMENT =>
qq{<span class="link small" onclick="show_seq('$align_str','HSP No. $hsp_num',0,0,0,0)">Alignment</span>}
    );

    my $html = $template->output;
    $template->param( HSP_IF => 0 );

    #get query sequence total length
    my $query = $dbh->selectall_arrayref(
        qq{SELECT * FROM sequence_info WHERE type = "query" AND name = "$qname"}
    );
    my $subject = $dbh->selectall_arrayref(
qq{SELECT * FROM sequence_info WHERE type = "subject" AND name = "$sname"}
    );
    my ( $query_image, $subject_image ) = generate_hit_image(
        hsp_num  => $hsp_num,
        hspdb    => $dbh,
        hsp_name => $hsp_id
    );
    my $query_link =
      qq{<div class=small>Query: $qname</div><img src=$query_image border=0>};
    $query_link =~ s/$TEMPDIR/$TEMPURL/;

	my ($a, $b) = get_link_coords($sstart, $sstop);
    my $subject_link =
qq{<div class=small>Subject: $org, Chromosome: $chr</div>} .
#qq{<a href = 'GenomeView.pl?chr=$chr&ds=$dsid&x=$sstart&z=5;gstid=$gstid' target=_new border=0><img src=$subject_image border=0></a>}; # mdb removed 11/20/13 issue 254
"<a href='GenomeView.pl?gid=$dsgid&loc=$chr:$a..$b' target=_new border='0'><img src='$subject_image' border='0'></a>"; # mdb added 11/20/13 issue 254
    $subject_link =~ s/$TEMPDIR/$TEMPURL/;
    return encode_json(
        {
            html         => $html,
            query_link   => $query_link,
            subject_link => $subject_link
        }
    );
}

sub generate_overview_image {
    my %opts        = @_;
    my $basename    = $opts{basename};
    my $type        = $opts{type};
    my $image_width = $opts{width};
    my @set         = split( /\n/, `ls $TEMPDIR/$basename*.blast` );
    my @reports;
    my $count = 1;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    foreach my $blast (@set) {
        my $report = new CoGe::Accessory::blast_report( { file => $blast } );
        my ($org_name) =
          $report->hsps->[ $count - 1 ]->subject_name =~ /^\s*(.*?)\s*\(/;
        push @reports,
          {
            report   => $report,
            organism => $org_name,
            link     => $TEMPURL . "/" . $basename . "-" . $count . ".blast",
          };
        $count++;
    }
    my $image_filename = $cogeweb->basefile . "_" . $type;
    my ( $chromosome_data, $chromosome_data_large ) =
      generate_chromosome_images(
        results  => \@reports,
        hsp_type => $type,

        #large_width => $image_width,
        filename => $image_filename
      );
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'CoGeBlast.tmpl' );
    $template->param( CHROMOSOMES_IF  => 1 );
    $template->param( CHROMOSOME_LOOP => $chromosome_data );
    my $chromosome_element = $template->output;
    $template->param( CHROMOSOMES_IF => 0 );
    my $chr_large_element;

    if ( scalar(@$chromosome_data_large) > 0 ) {
        $template->param( CHROMOSOMES_LARGE_IF  => 1 );
        $template->param( CHROMOSOME_LOOP_LARGE => $chromosome_data_large );
        $chr_large_element = $template->output;
        $template->param( CHROMOSOMES_LARGE_IF => 0 );
    }
    my $html = $chromosome_element . $chr_large_element;
    return $html;
}

sub generate_hit_image {
    my %opts     = @_;
    my $hsp_name = $opts{hsp_name};
    my $width    = $opts{width} || 400;
    my $dbh      = $opts{dbh};
    my $dbfile   = $cogeweb->sqlitefile;
    $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "" ) unless $dbh;
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});
    $sth->execute($hsp_name) || die "unable to execute";
    my $hsp   = $sth->fetchrow_hashref();
    my $query = $dbh->selectall_arrayref(
            qq{SELECT * FROM sequence_info WHERE type = "query" AND name = "}
          . $hsp->{qname}
          . qq{"} );
    my $subject = $dbh->selectall_arrayref(
            qq{SELECT * FROM sequence_info WHERE type = "subject" AND name = "}
          . $hsp->{sname}
          . qq{"} );
    my ($hsp_num) = $hsp_name =~ /^(\d+)_\d+$/;

    #generate_query_image
    my $cq = new CoGe::Graphics::Chromosome();
    $cq->iw($width);
    $cq->draw_chromosome(1);
    $cq->draw_ruler(1);
    $cq->draw_chr_end(0);
    $cq->minor_tick_labels(0);
    $cq->major_tick_labels(1);
    $cq->draw_hi_qual(0);
    $cq->padding(2);
    $cq->set_region( start => 1, stop => $query->[0][3] );
    $cq->feature_height(10);
    $cq->feature_labels(1);
    my $strand = $hsp->{strand} =~ /-/ ? "-1" : 1;
    my ( $qstart, $qstop ) = ( $hsp->{qstart}, $hsp->{qstop} );
    ( $qstart, $qstop ) = ( $qstop, $qstart ) if $qstart > $qstop;
    my $feat = CoGe::Graphics::Feature::HSP->new(
        {
            start  => $qstart,
            stop   => $qstop,
            strand => $strand,
            label  => $hsp_num,
            type   => "HSP"
        }
    );
    $feat->color( [ 255, 200, 0 ] );
    $cq->add_feature($feat);
    my ($chr) = $hsp->{schr};
    my ( $hspid, $dsgid ) = split( /_/, $hsp_name );
    my $dsg = $db->resultset('Genome')->find($dsgid);
    my ($ds) = $dsg->datasets( chr => $chr );
    my $len = $hsp->{sstop} - $hsp->{sstart} + 1;
    my $start = $hsp->{sstart} - 5000;
    $start = 1 if $start < 1;
    my $stop = $hsp->{sstop} + 5000;
    my $cs   = new CoGe::Graphics::Chromosome();
    $cs->iw($width);
    $cs->draw_chromosome(1);
    $cs->draw_ruler(1);
    $cs->draw_chr_end(0);
    $cs->minor_tick_labels(0);
    $cs->major_tick_labels(1);
    $cs->draw_hi_qual(0);
    $cs->padding(2);
    $cs->set_region( start => $start, stop => $stop );
    $cs->feature_height(10);
    $cs->overlap_adjustment(0);
    $cs->feature_labels(1);
    my ( $sstart, $sstop ) = ( $hsp->{sstart}, $hsp->{sstop} );
    ( $sstart, $sstop ) = ( $sstop, $sstart ) if $sstart > $sstop;
    $feat = CoGe::Graphics::Feature::HSP->new(
        {
            start  => $sstart,
            stop   => $sstop,
            strand => $strand,
            order  => 2,
            label  => $hsp_num,
            type   => "HSP"
        }
    );
    $feat->color( [ 255, 200, 0 ] );
    $cs->add_feature($feat);

    my $graphics = new CoGe::Graphics;
    $graphics->process_features(
        c      => $cs,
        layers => {
            features => { gene => 1, cds => 1, mrna => 1, rna => 1, cns => 1 }
        },
        ds   => $ds,
        chr  => $chr,
        coge => $db
    );
    $cs->overlap_adjustment(1);
    $cq->overlap_adjustment(1);

    #find neighboring hsps
    my $sname     = $hsp->{sname};
    my $qname     = $hsp->{qname};
    my $org       = $hsp->{org};
    my $schr      = $hsp->{schr};
    my $statement = qq{
select * from hsp_data where
((sstart <= $stop AND
sstart >= $start) OR
(sstop <= $stop AND
sstop >= $start)) AND
sname = "$sname" AND
qname = "$qname" AND
org = "$org" AND
schr = "$schr"
};
    $sth = $dbh->prepare($statement);
    $sth->execute;

    while ( my $data = $sth->fetchrow_hashref ) {
        next if $data->{name} eq $hsp_name;
        my ($label) = $data->{name} =~ /^(\d+)_\d+$/;
        $feat = CoGe::Graphics::Feature::HSP->new(
            {
                start  => $data->{sstart},
                stop   => $data->{sstop},
                strand => $data->{strand},
                order  => 2,
                label  => $label,
                type   => "HSP"
            }
        );
        $feat->color( [ 255, 0, 0 ] );
        $cs->add_feature($feat);
        $feat = CoGe::Graphics::Feature::HSP->new(
            {
                start  => $data->{qstart},
                stop   => $data->{qstop},
                strand => $data->{strand},
                order  => 1,
                label  => $label,
                type   => "HSP"
            }
        );
        $feat->color( [ 255, 0, 0 ] );
        $cq->add_feature($feat);
    }

    my $query_file = $cogeweb->basefile . ".q." . $hsp_name . ".png";
    $cq->generate_png( file => $query_file );
    my $sub_file = $cogeweb->basefile . ".s." . $hsp_name . ".png";
    $cs->generate_png( file => $sub_file );
    return $query_file, $sub_file;
}

sub overlap_feats_parse    #Send to GEvo
{
    my %opts      = @_;
    my $accn_list = $opts{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my @list;
    my @no_feats;
    my $url = "GEvo.pl?";
    my ( $chr, $dsgid, $loc );
    my $count = 1;

    foreach my $featid ( split( /,/, $accn_list ) ) {
        if ( $featid =~ /no/ ) {
            ( $dsgid, $chr, $loc ) = $featid =~ /^\d+_(\d+)_(.+?)_(\d+)no$/;
            push @no_feats, { dsgid => $dsgid, chr => $chr, loc => $loc };
        }
        else {
            my ( $hspnum, $dsgid );
            ( $featid, $hspnum, $dsgid ) = $featid =~ m/^(\d+)_(\d+)_(\d+)$/;
            my ($dsg) = $db->resultset('Genome')->find($dsgid);

            if ($dsg->deleted) {
                my $name = $dsg->organism->name;
                my $error_message = "The genome $name was marked as deleted";
                return encode_json({error => $error_message });
            }

            $featid .= "_" . $dsg->type->id if $dsg;
            push @list, $featid;
        }
    }
    my %seen = ();
    @list = grep { !$seen{$_}++ } @list;
    foreach my $featid (@list) {
        $url .= "fid$count=$featid&";
        $count++;
    }
    foreach my $no_feat (@no_feats) {
        $url .=
"dsgid$count=$no_feat->{dsgid}&chr$count=$no_feat->{chr}&x$count=$no_feat->{loc}&";
        $count++;
    }
    $count--;
    $url .= "num_seqs=$count";
    return encode_json( { url => $url, count => $count } );
}

sub initialize_sqlite {
    my $dbfile = $cogeweb->sqlitefile;
    return if -r $dbfile;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "" )
      || die "cant connect to db";
    my $create = qq{
CREATE TABLE hsp_data
(
id INTEGER PRIMARY KEY AUTOINCREMENT,
name varchar(50),
eval varchar(10),
pid varchar(10),
psim varchar(10),
score varchar(10),
qgap integer(10),
sgap integer(10),
match integer(10),
qmismatch integer(10),
smismatch integer(10),
strand varchar(2),
length integer(10),
qstart integer(50),
qstop integer(50),
sstart integer(50),
sstop integer(50),
qalign text,
salign text,
align text,
qname text,
sname text,
hsp_num integer,
org text,
schr text
)
};
    $dbh->do($create);
    my $index = qq{
  ALTER TABLE 'hsp_data' ADD AUTO_INCREMENT 'id';
 };
    $index = qq{
 CREATE INDEX name ON hsp_data (name)
 };
    $dbh->do($index);
    $index = qq{
 CREATE INDEX qname ON hsp_data (qname)
 };
    $dbh->do($index);
    $index = qq{
 CREATE INDEX sname ON hsp_data (sname)
 };
    $dbh->do($index);
    $create = qq{
 CREATE TABLE sequence_info
 (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 name varchar(255),
 type varchar(10),
 length integer(1024)
 )
 };
    $dbh->do($create);
    $index = qq{
  CREATE INDEX seqname ON sequence_info (name)
  };
    $dbh->do($index);
    $index = qq{
  CREATE INDEX type ON sequence_info (type)
  };
    $dbh->do($index);
    system "/bin/chmod +rw $dbfile";
}

sub populate_sqlite {
    my ( $hsp, $dsgid, $org, $chr ) = @_;
    my $pval    = $hsp->pval;
    my $pid     = $hsp->percent_id;
    my $qgap    = $hsp->query_gaps;
    my $sgap    = $hsp->subject_gaps;
    my $score   = $hsp->score;
    my $psim    = $hsp->percent_sim;
    my $length  = $hsp->length;
    my $strand  = $hsp->strand;
    my $match   = $hsp->match;
    my $qalign  = $hsp->query_alignment;
    my $salign  = $hsp->subject_alignment;
    my $align   = $hsp->alignment;
    my $qname   = $hsp->query_name;
    my $sname   = $hsp->subject_name;
    my $name    = $hsp->number . "_" . $dsgid;
    my $hsp_num = $hsp->number;

    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    $dbh->do("begin exclusive transaction") or die $dbh->errstr;
    my $query_length =
      $hsp->query_stop > $hsp->query_start
      ? ( ( $hsp->query_stop ) - ( $hsp->query_start ) + 1 )
      : ( ( $hsp->query_start ) - ( $hsp->query_stop ) + 1 );
    my ( $qstart, $qstop ) =
      $hsp->query_stop > $hsp->query_start
      ? ( $hsp->query_start, $hsp->query_stop )
      : ( $hsp->query_stop, $hsp->query_start );
    my $qmismatch = $query_length - $hsp->match;

    my $subject_length =
      $hsp->subject_stop > $hsp->subject_start
      ? ( ( $hsp->subject_stop ) - ( $hsp->subject_start ) + 1 )
      : ( ( $hsp->subject_start ) - ( $hsp->subject_stop ) + 1 );
    my ( $sstart, $sstop ) = ( $hsp->subject_start, $hsp->subject_stop );
    my $smismatch = $subject_length - $hsp->match;

    my $statement =
qq{INSERT INTO hsp_data (name, eval, pid, psim, score, qgap, sgap,match,qmismatch,smismatch, strand, length, qstart,qstop, sstart, sstop,qalign,salign,align,qname,sname,hsp_num,org,schr) values ('$name', '$pval', '$pid','$psim', '$score', $qgap, $sgap, $match,$qmismatch, $smismatch, '$strand',$length,$qstart, $qstop, $sstart, $sstop,'$qalign','$salign','$align','$qname','$sname',$hsp_num,'$org','$chr')};
    #my $insert_prep = $dbh->prepare('INSERT INTO hsp_data (name, eval, pid, psim, score, qgap, sgap,match,qmismatch,smismatch, strand, length, qstart,qstop, sstart, sstop,qalign,salign,align,qname,sname,hsp_num,org,schr) values (?, ?, ?,?, ?, ?, ?, ?,?, ?, ?,?,?, ?, ?, ?,?,?,?,?,?,?,?,?)');
    print STDERR $statement unless $dbh->do($statement);
    #$insert_prep->execute("$name", "$pval", "$pid","$psim", "$score", $qgap, $sgap, $match,$qmismatch, $smismatch, "$strand",$length,$qstart, $qstop, $sstart, $sstop,"$qalign","$salign","$align","$qname","$sname",$hsp_num,"$org","$chr");

    # Populate sequence_info table
    $statement = "SELECT name FROM sequence_info where name = '$qname'";
    my $val = $dbh->selectall_arrayref($statement);
    unless ( $val->[0][0] ) {
        my $qlength = $hsp->query_length;
        $statement =
qq{INSERT INTO sequence_info (name, type, length) values ('$qname',"query",'$qlength')};
        print STDERR $statement unless $dbh->do($statement);
    }
    $statement = "SELECT name FROM sequence_info where name = '$sname'";
    $val       = $dbh->selectall_arrayref($statement);
    unless ( $val->[0][0] ) {
        my $slength = $hsp->subject_length;
        $statement =
qq{INSERT INTO sequence_info (name, type, length) values ('$sname',"subject",'$slength')};
        print STDERR $statement unless $dbh->do($statement);
    }
    $dbh->do("commit transaction") or die $dbh->errstr;
    $dbh->disconnect();

}
sub get_nearby_feats {
    my %opts     = @_;
    my $hsp_id   = $opts{num};
    my $filename = $opts{basefile};

    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbh = DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    $hsp_id =~ s/^table_row// if $hsp_id =~ /table_row/;
    $hsp_id =~ s/^\d+_// if $hsp_id =~ tr/_/_/ > 1;

    my $name     = "None";
    my $distance = ">250";
    my $sth      = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});
    my ( $sstart, $sstop, $sname, $chr );
    my ( $hsp_num, $dsgid ) = $hsp_id =~ /^(\d+)_(\d+)$/;
    $sth->execute($hsp_id) || die "unable to execute";
    while ( my $info = $sth->fetchrow_hashref() ) {
        $sstart = $info->{sstart};
        $sstop  = $info->{sstop};
        $sname  = $info->{sname};
        $chr    = $info->{schr};
    }
    my $dsg = $db->resultset('Genome')->find($dsgid);
    my @dsids = map { $_->id } $dsg->datasets( chr => $chr );
    my ( $start, $stop ) = ( $sstart, $sstop );
    my @feat;
    my $count = 0;
    my $mid   = ( $stop + $start ) / 2;
    $dbh = $db->storage->dbh;    #DBI->connect( $connstr, $DBUSER, $DBPASS );
    my %fids;
    for (@dsids) {
        my $query = qq{
    select * from (
      (SELECT * FROM
             ((SELECT feature_id,start,stop FROM feature where start<=$mid and dataset_id = $_ and chromosome = '$chr' ORDER BY start DESC LIMIT 10)
        UNION (SELECT feature_id,start,stop FROM feature where start>=$mid and dataset_id = $_ and chromosome = '$chr' ORDER BY start LIMIT 10)) as u)
      UNION
      (SELECT * FROM
             ((SELECT feature_id,start,stop FROM feature where stop<=$mid and dataset_id = $_ and chromosome = '$chr' ORDER BY stop DESC LIMIT 10)
        UNION (SELECT feature_id,start,stop FROM feature where stop>=$mid and dataset_id = $_ and chromosome = '$chr' ORDER BY stop LIMIT 10)) as v)
       ) as w
    order by abs((start + stop)/2 - $mid) LIMIT 10};
        for my $ids ($dbh->selectcol_arrayref($query)) {
            for (@{$ids}) {
                $fids{$_} = 1;
            }
        }
    }

    my $new_checkbox_info;
    my %dist;

    for ( keys %fids ) {
        my $fid = $_;
        my $tmpfeat = $db->resultset('Feature')->find($fid);
        next
          unless $tmpfeat->type->name =~ /gene/i
              || $tmpfeat->type->name =~ /rna/i
              || $tmpfeat->type->name =~ /cds/i;
        my $newmin = abs($tmpfeat->start - $mid) < abs($tmpfeat->stop - $mid) ? abs($tmpfeat->start - $mid) : abs($tmpfeat->stop - $mid);
        $dist{$fid} = { dist => $newmin, feat => $tmpfeat };
    }
    my $feat;
    my $min_dist;
	my @fids = sort keys %dist;
    my @sorted = sort { $dist{$a}{dist} <=> $dist{$b}{dist} } @fids;
    for ( @sorted ) {
        $min_dist = $dist{$_}{dist} unless defined $min_dist;
        $feat     = $dist{$_}{feat} unless $feat;
        last if $feat->type->name eq "CDS";
        my $f = $dist{$_}{feat};
        if ( $f->type->name eq "CDS" ) {
            $feat = $f unless $feat->stop < $f->start || $feat->start > $f->stop;
        }
    }
    if ($feat) {
        if (( $start >= $feat->start && $start <= $feat->stop ) || ( $stop >= $feat->start && $stop <= $feat->stop )) {
            $distance = "overlapping";
        }
        else {
            $distance = abs( ( $stop + $start ) / 2 - ( $feat->stop + $feat->start ) / 2 );
        }
        ($name) = $feat->names;
        $name =
qq{<span class="link" title="Click for Feature Information" onclick=update_info_box('}
          . $feat->id . "_"
          . $hsp_num . "_"
          . $dsgid
          . "')>$name</span>";
        $new_checkbox_info =
            $hsp_id . "_"
          . $chr . "_"
          . $sstart . "no,"
          . $feat->id . "_"
          . $hsp_id;
    }
    else {
        $distance = "No neighboring features found";
    }

    return encode_json(
        {
            name              => $name,
            distance          => $distance,
            hsp_id            => $hsp_id,
            new_checkbox_info => $new_checkbox_info
        }
    );
}

sub export_CodeOn {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "CodeOn.pl?fid=";
    my @list;
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ( $featid, $hspnum, $dsgid ) = $accn =~ m/^(\d+)_(\d+)_(\d+)$/;

        #	my $dsg = $db->resultset('Genome')->find($dsgid);
        #	$featid .= "_".$dsg->type->id if $dsg;
        push @list, $featid;
    }
    my %seen = ();
    @list = grep { !$seen{$_}++ } @list;
    $url .= join( "::", @list );
    $url =~ s/&$//;
    return $url;
}

sub export_fasta_file {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "FastaView.pl?fid=";
    my @list;
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ( $featid, $hspnum, $dsgid ) = $accn =~ m/^(\d+)_(\d+)_(\d+)$/;
        my $dsg = $db->resultset('Genome')->find($dsgid);
        $featid .= "_" . $dsg->type->id if $dsg;
        push @list, $featid;
    }
    my %seen = ();
    @list = grep { !$seen{$_}++ } @list;
    $url .= join ",", @list;
    return $url;
}

sub export_top_hits {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    my $filename  = $opts{filename};
    return unless (defined $accn_list and defined $filename);

    # Setup file path and SQLite DB
    my $web = CoGe::Accessory::Web::initialize_basefile( basename => $filename, tempdir  => $TEMPDIR );
    my $dbh = DBI->connect( "dbi:SQLite:dbname=" . $web->sqlitefile, "", "" );
    my $sth = $dbh->prepare( qq{SELECT * FROM hsp_data WHERE name = ?} );

    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    # Filter top scoring hits
    my %hsp;
    foreach my $accn ( split( /,/, $accn_list ) ) {
	my ($hsp_num, $dsgid);
        if ($accn =~ /no$/) {
		my $accn_with_commas = $accn;
            $accn_with_commas =~ tr/_/,/;
            ( $hsp_num, $dsgid ) = $accn_with_commas =~ /(\d+),(\d+),\w*_?\d+,\d+no$/;
        }
        else {
		if ( $accn =~ tr/_/_/ > 2 ) {
                	my $accn_with_commas = $accn;
                	$accn_with_commas =~ tr/_/,/;
                	( $hsp_num, $dsgid ) =
                  	$accn_with_commas =~
                  		/(\d+),(\d+),(\w*_?\d+),(\d+),(\d+),(\d+.?\d*)/;
            	}
            	else {
                	( undef, $hsp_num, $dsgid ) = $accn =~ /(\d+)_(\d+)_(\d+)/;
            	}
        }
        $sth->execute( $hsp_num . "_" . $dsgid ) || die "unable to execute";
        my ( $qname, $salign, $score );
        while ( my $hsp_info = $sth->fetchrow_hashref() ) { # FIXME: mdb asks "why is this looped?"
            $qname  = $hsp_info->{qname};
            $salign = $hsp_info->{salign};
            $score  = $hsp_info->{score};
        }
        if ( !defined $hsp{$qname} || $score > $hsp{$qname}{score} ) {
            $hsp{$qname}{fasta} = ">$qname\n$salign\n";
            $hsp{$qname}{score} = $score;
        }
    }

    # Generate top hits file
    open( my $fh, "> $TEMPDIR/$filename-tophits.fasta" );
    foreach my $qname (sort keys %hsp) {
        print $fh $hsp{$qname}{fasta};
    }
    close($fh);
    
    return "$TEMPURL/$filename-tophits.fasta";
}

sub export_to_excel {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    my $filename  = $opts{filename};

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( basename => $filename, tempdir  => $TEMPDIR );
    my $dbh = DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});

    my $workbook = Spreadsheet::WriteExcel->new("$TEMPDIR/Excel_$filename.xls");
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $i = 1;
    $worksheet->write( 0, 0, "Query" );
    $worksheet->write( 0, 1, "Organism" );
    $worksheet->write( 0, 2, "Chr" );
    $worksheet->write( 0, 3, "Position" );
    $worksheet->write( 0, 4, "HSP No." );
    $worksheet->write( 0, 5, "Length" );
    $worksheet->write( 0, 6, "E-value" );
    $worksheet->write( 0, 7, "Percent ID" );
    $worksheet->write( 0, 8, "Score" );
    $worksheet->write( 0, 9, "CoGe Feature ID" );
    $worksheet->write( 0, 10,  "Closest Feature" );
    $worksheet->write( 0, 11, "Distance" );

    my (
        $org,   $chr,      $pos,    $hsp_no, $eval,   $pid,
        $score, $distance, $featid, $dsgid,  $length, $qname
    );
    foreach my $accn ( split( /,/, $accn_list ) ) {
        if ( $accn =~ /no/ ) {
            my $accn_with_commas = $accn;
            $accn_with_commas =~ tr/_/,/;
            ( $hsp_no, $dsgid, $chr, $pos ) =
              $accn_with_commas =~ /(\d+),(\d+),(\w*_?\d+),(\d+)no/;
            $sth->execute( $hsp_no . "_" . $dsgid ) || die "unable to execute";
            while ( my $info = $sth->fetchrow_hashref() ) {
                $eval   = $info->{eval};
                $pid    = $info->{pid};
                $score  = $info->{score};
                $pos    = $info->{sstart};
                $org    = $info->{org};
                $chr    = $info->{schr};
                $length = $info->{length};
                $qname  = $info->{qname};
            }
            $worksheet->write( $i, 0, $qname );
            $worksheet->write( $i, 1, $org );
            $worksheet->write( $i, 2, $chr );
            $worksheet->write( $i, 3, $pos );
            $worksheet->write( $i, 4, $hsp_no );
            $worksheet->write( $i, 5, $length );
            $worksheet->write( $i, 6, $eval );
            $worksheet->write( $i, 7, $pid );
            $worksheet->write( $i, 8, $score );
        }
        else {
            # there doesn't appear to code that would pass this in
            # if ( $accn =~ tr/_/_/ > 2 ) {
            #     my $accn_with_commas = $accn;
            #     $accn_with_commas =~ tr/_/,/;
            #     ( $hsp_no, $dsgid, $chr, $pos, $featid, $distance ) =
            #       $accn_with_commas =~
            #       /(\d+),(\d+),(\w*_?\d+),(\d+),(\d+),(\d+.?\d*)/;
            # }
            # else {
                ( $featid, $hsp_no, $dsgid, $distance ) = $accn =~ /(\d+)_(\d+)_(\d+)_(.+)/;
            #     $distance = "overlapping";
            # }
            $sth->execute( $hsp_no . "_" . $dsgid ) || die "unable to execute";
            while ( my $info = $sth->fetchrow_hashref() ) {
                $eval   = $info->{eval};
                $pid    = $info->{pid};
                $score  = $info->{score};
                $pos    = $info->{sstart};
                $org    = $info->{org};
                $chr    = $info->{schr};
                $length = $info->{length};
                $qname  = $info->{qname};
            }
            my ($feat) = $db->resultset("Feature")->find($featid);
            my ($name) = sort $feat->names;
            $worksheet->write( $i, 0, $qname );
            $worksheet->write( $i, 1, $org );
            $worksheet->write( $i, 2, $chr );
            $worksheet->write( $i, 3, $pos );
            $worksheet->write( $i, 4, $hsp_no );
            $worksheet->write( $i, 5, $length );
            $worksheet->write( $i, 6, $eval );
            $worksheet->write( $i, 7, $pid );
            $worksheet->write( $i, 8, $score );
            $worksheet->write( $i, 9, $feat->id );
            $worksheet->write( $i, 10, $P->{SERVER} . "FeatView.pl?accn=$name", $name );
            $worksheet->write( $i, 11, $distance );
        }

        $i++;
    }

    $workbook->close() or die "Error closing file: $!";

    return "$TEMPURL/Excel_$filename.xls";
}

sub generate_tab_deliminated {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    my $filename  = $opts{filename};

    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data WHERE name = ?});

    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    my $str = "Name\tHSP No.\tE-value\tPerc ID\tScore\tOrganism\tCoge Feature ID\n";

    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ( $featid, $hsp_num, $dsgid ) = $accn =~ m/^(\d+)_(\d+)_(\d+)$/;
        my $dsg = $db->resultset("Genome")->find($dsgid);
        my ($feat) = $db->resultset("Feature")->find($featid);

        my ($name) = $feat->names;
        my $org = $dsg->organism->name;
        $sth->execute( $hsp_num . "_" . $dsgid ) || die "unable to execute";
        my ( $pval, $pid, $score );
        while ( my $info = $sth->fetchrow_hashref() ) {
            $pval  = $info->{eval};
            $pid   = $info->{pid};
            $score = $info->{score};
        }
        $str .= "$name\t$hsp_num\t$pval\t$pid\t$score\t$org\t$featid\n";
    }
    $str =~ s/\n$//;

    open( NEW, "> $TEMPDIR/tab_delim$filename.tabbed" );
    print NEW $str;
    close NEW;
    return "$TEMPURL/tab_delim$filename.tabbed";
}

sub generate_feat_list {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    my $filename  = $opts{filename};

    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    my $url = "FeatList.pl?";
    my @list;
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ( $featid, $hspnum, $dsgid ) = $accn =~ m/^(\d+)_(\d+)_(\d+)$/;
        my $dsg = $db->resultset('Genome')->find($dsgid);
        $featid .= "_" . $dsg->type->id if $dsg;
        push @list, $featid;
    }
    my %seen = ();
    @list = grep { !$seen{$_}++ } @list;
    foreach my $featid (@list) {
        $url .= "fid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub generate_blast {
    my %opts      = @_;
    my $accn_list = $opts{accn};
    my $filename  = $opts{filename};

    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    my $url = "CoGeBlast.pl?";
    my @list;
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ( $featid, $hspnum, $dsgid ) = $accn =~ m/^(\d+)_(\d+)_(\d+)$/;
        my $dsg = $db->resultset('Genome')->find($dsgid);
        $featid .= "_" . $dsg->type->id if $dsg;
        push @list, $featid;
    }
    my %seen = ();
    @list = grep { !$seen{$_}++ } @list;
    foreach my $featid (@list) {
        $url .= "featid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub get_genome_info { #FIXME: dup'ed in SynFind.pl
    my %opts  = @_;
    my $dsgid = $opts{dsgid};

    #	print STDERR "get_genome_info: $dsgid\n";
    return " " unless $dsgid;

    my $dsg = $db->resultset("Genome")->find($dsgid);
    return "Unable to create genome object for id: $dsgid" unless $dsg;

    my $html = qq{<table class='small'>}
      ;    # = qq{<div style="overflow:auto; max-height:78px">};
    $html .= qq{<tr valign='top'><td style='white-space:nowrap'>Name:<td><span class='link' onclick=window.open('GenomeInfo.pl?gid=$dsgid')>}
      . $dsg->organism->name
      . "</span>";
    $html .= "<tr valign='top'><td style='white-space:nowrap'>Description:<td>" . join(
        "; ",
        map {
qq{<span class=link onclick="\$('#org_desc').val('$_').focus();">$_</span>}
          } split( /;\s+/, $dsg->organism->description )
      ) if $dsg->organism->description;

    my $chr_num = $dsg->chromosome_count;
    my $total_length = $dsg->length;
    my ($ds) = $dsg->datasets;
    my $link = $ds->data_source->link;
    $link = "http://" . $link unless $link =~ /^http/;
    $html .=
        "<tr valign='top'><td>Source:"
      . '<td>' . ($link ? "<a class='link' href='$link' target=_new>" : '')
      . $ds->data_source->name . "</a>"
      . qq{<tr valign='top'><td style='white-space:nowrap'>Chromosomes:<td>$chr_num}
      . qq{<tr valign='top'><td style='white-space:nowrap'>Total length:<td>}
      . commify($total_length) . " bp"
      . qq{<tr valign='top'><td style='white-space:nowrap'>Sequence Type:<td>}
      . $dsg->genomic_sequence_type->name
      . qq{<input type='hidden' id='gstid' value=}
      . $dsg->genomic_sequence_type->id;

    #    foreach my $gs (@gs)
    #      {
    #	my $chr = $gs->chromosome;
    #	my $length = $gs->sequence_length;
    #	$length = commify($length);
    #	$html .= qq{$chr:  $length bp<br>};
    #     }
    #    $html .= "</table>";
    return $html;
}

sub export_hsp_info {
    my %opts = @_;

    #	my $accn = $opts{accn};
    my $filename = $opts{filename};

    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data ORDER BY org,hsp_num});
    $sth->execute();

    my $str =
qq{Org\tChr\tPosition\tStrand\tHSP No.\tPercent ID\tAlignment length\tE-value\tScore\tMatch\tQuery Mismatch\tQuery Gap\tQuery Length\tSubject Mismatch\tSubject Gap\tSubject Length\tQuery HSP Sequence\tSubject HSP Sequence\n};

    my (
        $org,            $chr,              $pos,
        $hsp_num,        $pid,              $align_length,
        $eval,           $score,            $match,
        $strand,         $query_mismatch,   $query_gap,
        $query_length,   $subject_mismatch, $subject_gap,
        $subject_length, $query_seq,        $subject_seq,
        $align_seq,      $sstart,           $sstop,
        $sname,          $align
    );

    while ( my $row = $sth->fetchrow_hashref() ) {
        ($hsp_num) = $row->{name} =~ /(\d+)_\d+/;
        $org              = $row->{org};
        $eval             = $row->{eval};
        $pid              = $row->{pid};
        $score            = $row->{score};
        $query_gap        = $row->{qgap};
        $subject_gap      = $row->{sgap};
        $match            = $row->{match};
        $query_mismatch   = $row->{qmismatch};
        $subject_mismatch = $row->{smismatch};
        $strand           = $row->{strand};
        $align_length     = $row->{length};
        $sstart           = $row->{sstart};
        $sstop            = $row->{sstop};
        $query_seq        = $row->{qalign};
        $subject_seq      = $row->{salign};
        $align            = $row->{align};
        $sname            = $row->{sname};

        ($chr) = $row->{schr};
        $pos       = $sstart . " - " . $sstop;
        $align_seq = $query_seq . "<br>" . $align . "<br>" . $subject_seq;

        $query_seq =~ s/-//g;
        $query_length = length $query_seq;

        $subject_seq =~ s/-//g;
        $subject_length = length $subject_seq;

        $str .=
"$org\t$chr\t$pos\t$strand\t$hsp_num\t$pid\t$align_length\t$eval\t$score\t$match\t$query_mismatch\t$query_gap\t$query_length\t$subject_mismatch\t$subject_gap\t$subject_length\t$query_seq\t$subject_seq\n";
    }

    open( NEW, "> $TEMPDIR/tab_delim_$filename.txt" );
    print NEW $str;
    close NEW;

    return "$TEMPURL/tab_delim_$filename.txt";
}

sub export_hsp_query_fasta {
    my %opts     = @_;
    my $filename = $opts{filename};

    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data});
    $sth->execute();

    my ( $fasta, $query_seq, $name, $qstart, $qstop, $hsp_num );
    while ( my $row = $sth->fetchrow_hashref() ) {
        $query_seq = $row->{qalign};
        $name      = $row->{qname};
        $name =~ s/^\s+//;
        $name =~ s/\s+$//;
        $qstart  = $row->{qstart};
        $qstop   = $row->{qstop};
        $hsp_num = $row->{hsp_num};
        $query_seq =~ s/-//g;
        $fasta .=
            ">HSP"
          . $hsp_num . "_"
          . $name
          . ", Subsequence: "
          . $qstart . "-"
          . $qstop . "\n"
          . $query_seq . "\n";
    }
    open( NEW, "> $TEMPDIR/query_fasta_$filename.txt" );
    print NEW $fasta;
    close NEW;
    return "$TEMPURL/query_fasta_$filename.txt";
}

sub export_hsp_subject_fasta {
    my %opts     = @_;
    my $filename = $opts{filename};
    my $dna      = $opts{dna};        #1 if shift;

    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );

    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data});
    $sth->execute();

    my (
        $fasta, $subject_seq, $sstart, $sstop, $hsp_num,
        $chr,   $org,         $dsgid,  $strand
    );
    while ( my $row = $sth->fetchrow_hashref() ) {
        $org         = $row->{org};
        $chr         = $row->{schr};
        $strand      = $row->{strand} =~ /\+/ ? 1 : -1;
        $sstart      = $row->{sstart};
        $sstop       = $row->{sstop};
        $hsp_num     = $row->{hsp_num};
        $subject_seq = $row->{salign};
        $subject_seq =~ s/-//g;
        $fasta .=
            ">HSP "
          . $hsp_num . " "
          . $org
          . ", Chromsome: "
          . $chr
          . ", Location: "
          . $sstart . "-"
          . $sstop
          . "($strand)" . "\n"
          . $subject_seq . "\n";
    }
    open( NEW, "> $TEMPDIR/subject_fasta_$filename.txt" );
    print NEW $fasta;
    close NEW;
    return "$TEMPURL/subject_fasta_$filename.txt";
}

sub export_alignment_file {
    my %opts     = @_;
    my $filename = $opts{filename};

    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $filename,
        tempdir  => $TEMPDIR
    );
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    my $sth = $dbh->prepare(qq{SELECT * FROM hsp_data ORDER BY org,hsp_num});
    $sth->execute();

    my (
        $sname, $qname, $qstop, $qstart,  $sstart, $sstop, $align,
        $qseq,  $sseq,  $str,   $hsp_num, $chr,    $org
    );
    while ( my $row = $sth->fetchrow_hashref() ) {
        $hsp_num = $row->{hsp_num};
        $sseq    = $row->{salign};
        $sname   = $row->{sname};
        $sstart  = $row->{sstart};
        $sstop   = $row->{sstop};
        $qseq    = $row->{qalign};
        $qname   = $row->{qname};
        $qname =~ s/^\s+//;
        $qname =~ s/\s+$//;
        $sname =~ s/^\s+//;
        $sname =~ s/\s+$//;
        $qstart = $row->{qstart};
        $qstop  = $row->{qstop};
        $align  = $row->{align};
        $chr    = $row->{schr};
        $org    = $row->{org};
        $str .=
            "HSP: "
          . $hsp_num
          . "\n>Query: "
          . $qname
          . ", Subsequence: "
          . $qstart . "-"
          . $qstop
          . "\n>Subject: "
          . $org
          . ", Chromosome"
          . $chr
          . ", Location: "
          . $sstart . "-"
          . $sstop . "\n"
          . $qseq . "\n"
          . $align . "\n"
          . $sseq . "\n\n";
    }
    $str =~ s/\n+$//;
    open( NEW, "> $TEMPDIR/alignment_file_$filename.txt" );
    print NEW $str;
    close NEW;
    return "$TEMPURL/alignment_file_$filename.txt";
}

sub save_settings {
    my %opts = @_;

    my $prefs = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $db
    );
    delete $prefs->{display} unless ref( $prefs->{display} ) eq "HASH";
    foreach my $key ( keys %opts ) {
        if ( $key eq "display" ) {
            delete $prefs->{$key};
            my %settings = (
                1  => 'QND',
                2  => 'OrgD',
                3  => 'ChrD',
                4  => 'PosD',
                5  => 'NumD',
                6  => 'LengthD',
                7  => 'QSCD',
                8  => 'EvalD',
                9  => 'PIDD',
                10 => 'ScoreD',
                11 => 'QualD',
                12 => 'FeatD',
                13 => 'DistD',
                14 => 'SeqViewD',
            );
            foreach my $index ( split( /,/, $opts{$key} ) ) {
                $prefs->{$key}{ $settings{$index} } = 1;
            }
        }
        else {
            $prefs->{$key} = $opts{$key} if $opts{$key};
        }
        delete $prefs->{$key} unless $opts{$key};
    }
    my $item = CoGe::Accessory::Web::save_settings(
        opts => $prefs,
        user => $USER,
        page => $PAGE_NAME,
        coge => $db
    );
}

sub color_pallet {
    my %opts     = @_;
    my $start    = $opts{start} || [ 20, 200, 20 ];
    my $offset   = $opts{offset} || 75;
    my $num_seqs = $opts{num_seqs};

    my @colors;
    my %set = (
        HSP_NUM => 1,
        RED     => $start->[0],
        GREEN   => $start->[1],
        BLUE    => $start->[2]
    );

    my $temp = [@$start];
    for ( my $i = 1 ; $i <= $num_seqs ; $i++ ) {
        foreach (@$temp) {
            if ( $_ < 20 ) {
                my $color = $i > 18 ? int( 255 * rand() ) / 2 + 125 : 255;
                $_ = $color;
            }
        }
        my @color;
        @color = @$temp;
        push @colors, [ $color[0], $color[1], $color[2], ];
        unless ( $i % 3 ) {
            $temp = [ map { int( $_ / 2 ) } @color ];
        }
        else {
            $temp = [ $temp->[2], $temp->[0], $temp->[1] ];
        }
        $temp = [ map { $_ - 25 } @$temp ] unless ( $i % 6 ) || $i < 3;

    }

    #print STDERR map {join ("\t", @$_)."\n"} @colors;
    return wantarray ? @colors : \@colors;
}

sub search_lists {   # FIXME this coded is dup'ed in User.pl and NotebookView.pl and SynFind.pl
    my %opts = @_;
    return if ( $USER->user_name eq 'public' );
    my $search_term = $opts{search_term};
    my $timestamp   = $opts{timestamp};
    #print STDERR "$search_term $timestamp\n";

    my @notebooks;
    my $num_results = 0;
    my $group_str = join( ',', map { $_->id } $USER->groups );

    # Try to get all items if blank search term - mdb removed 3/6/14, this is too slow
#    if ( !$search_term ) {
#        my $sql = "locked=0"; # AND restricted=0 OR user_group_id IN ( $group_str ))"; # FIXME
#        $num_results = $db->resultset("List")->count_literal($sql);
#        if ( $num_results < $MAX_SEARCH_RESULTS ) {
#            foreach my $notebook ( $db->resultset("List")->search_literal($sql) )
#            {
#                next unless $USER->has_access_to_notebook($notebook);
#                push @notebooks, $notebook;
#            }
#        }
#    }

    # Perform search
    if ($search_term) {
        # Get public lists and user's private lists
        $search_term = '%' . $search_term . '%';
        foreach my $notebook ($db->resultset("List")->search_literal("locked=0 AND (name LIKE '$search_term' OR description LIKE '$search_term')"))
        {
            next unless $USER->has_access_to_notebook($notebook);
            push @notebooks, $notebook;
        }
        $num_results = @notebooks;
    }

    # Limit number of results display
    if ( $num_results > $MAX_SEARCH_RESULTS ) {
        return encode_json({
            timestamp => $timestamp,
            html => "<option>$num_results matches, please refine your search.</option>"
        });
    }

    # Build select items out of results
    my $html;
    foreach my $n ( sort notebookcmp @notebooks ) {
        my $item_spec = 1 . ':' . $n->id; #FIXME magic number for item_type
        $html .= "<option value='$item_spec'>" . $n->info . "</option><br>\n";
    }
    if (!$html) {
        if (!$search_term) {
            $html = "<option disabled='disabled'></option>";
        }
        else {
            $html = "<option disabled='disabled'>No matches</option>";
        }
    }

    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub get_list_preview {
    my %opts      = @_;
    my $item_spec = $opts{item_spec};

    my ( undef, $lid ) = split( ':', $item_spec );

    my $list = $db->resultset('List')->find($lid);
    my $html = '';
    if ($list) {
        $html .=
          "List <b>'" . $list->name . "'</b> (id" . $list->id . ') contains ';
        $html .=
            @{ $list->genomes }
          . ' genomes <ul>'
          . join( '', map { '<li>' . $_->info . '</li>' } $list->genomes )
          . '</ul>';
    }

    return $html;
}

sub get_genomes_for_list {
    my %opts      = @_;
    my $item_spec = $opts{item_spec};

    my ( undef, $lid ) = split( ':', $item_spec );

    my $list    = $db->resultset('List')->find($lid);
    my $genomes = '';
    if ($list) {
        my $dsgids = join( ',', map { $_->id } $list->genomes );
        $genomes = get_dsg_for_menu( dsgid => $dsgids );
    }

    return $genomes;
}
