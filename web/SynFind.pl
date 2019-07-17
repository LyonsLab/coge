#! /usr/bin/perl -w
use v5.10;
use strict;
use CoGeX;
use CoGeX::Result::Genome qw(ERROR LOADING);
use DBIxProfiler;
use CoGe::Accessory::Utils qw( commify html_escape);
use CoGe::JEX::Jex;
use CoGe::Accessory::Web;
use CoGe::Core::Notebook qw(notebookcmp);
use CoGe::Core::Genome qw(genomecmp2);
use CoGe::Core::Favorites;
use CGI;
use JSON::XS;
use HTML::Template;
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use POSIX;
use File::Path;
use File::Basename qw(basename fileparse);
use File::Spec::Functions;
use Benchmark qw(:all);
use Parallel::ForkManager;
use DBI;
use Sort::Versions;

no warnings 'redefine';

#example URL: http://genomevolution.org/coge/SynFind.pl?fid=34519245;qdsgid=3;dsgid=4241,6872,7084,7094,7111

use vars qw($config $DBNAME $DBHOST $DBPORT $DBUSER $DBPASS
    $PAGE_TITLE $PAGE_NAME $DIR $LINK $TEMPDIR $TEMPURL $DATADIR $FASTADIR
    $BLASTDBDIR $DIAGSDIR $BEDDIR $LASTZ $LAST $LASTDB $CONVERT_BLAST $BLAST2BED
    $BLAST2RAW $SYNTENY_SCORE $DATASETGROUP2BED $PYTHON26 $FORM $USER $DATE
    $coge $cogeweb $RESULTSLIMIT $MAX_SEARCH_RESULTS $MAX_PROC $SERVER $COOKIE_NAME
    $JEX $GEN_FASTA $MAX_FEATURES_IN_LIST $SCRIPTDIR);

$PAGE_TITLE = "SynFind";
$PAGE_NAME  = $PAGE_TITLE . ".pl";
$RESULTSLIMIT         = 100;
$MAX_SEARCH_RESULTS   = 1000;
$MAX_FEATURES_IN_LIST = 100;

$FORM = new CGI;
( $coge, $USER, $config, $LINK ) = CoGe::Accessory::Web->init(
    page_title => $PAGE_TITLE,
    cgi => $FORM
);

$JEX = CoGe::JEX::Jex->new( host => $config->{JOBSERVER}, port => $config->{JOBPORT} );

$ENV{PATH}     = $config->{COGEDIR};
$TEMPDIR       = $config->{TEMPDIR} . $PAGE_TITLE;
$TEMPURL       = $config->{TEMPURL} . $PAGE_TITLE;
$SERVER        = $config->{SERVER};
$ENV{BLASTDB}  = $config->{BLASTDB};
$ENV{BLASTMAT} = $config->{BLASTMATRIX};
$DIR           = $config->{COGEDIR};
$DATADIR       = $config->{DATADIR};
$DIAGSDIR      = $config->{DIAGSDIR};

$FASTADIR   = $config->{FASTADIR};
$BLASTDBDIR = $config->{BLASTDB};
$BEDDIR     = $config->{BEDDIR};
mkpath( $BEDDIR, 0, 0777 ) unless -d $BEDDIR;

$MAX_PROC = $config->{MAX_PROC};
$LASTZ = 'nice ' . get_command_path('PYTHON') . ' ' . $config->{MULTI_LASTZ} . " -A $MAX_PROC --path=" . get_command_path('LASTZ');
#$LAST  = 'nice ' . $config->{PYTHON} . ' ' . $config->{MULTI_LAST} . " -a $MAX_PROC --path=" . $config->{LAST_PATH}; $ mdb removed /3/21/16 for LAST v731 -- wrapper no longer needed
$LAST = $config->{LASTAL} // 'lastal'; $LAST .= " -u 0 -P $MAX_PROC -i3G -f BlastTab"; # mdb added 3/21/16 for LAST v731
$LASTDB = $config->{LASTDB2} // 'lastdb'; # mdb added 3/21/16 for LAST v731

$SCRIPTDIR        = catdir($config->{SCRIPTDIR}, 'synmap');
$GEN_FASTA        = 'nice ' . catfile($SCRIPTDIR, 'generate_fasta.pl');
$CONVERT_BLAST    = 'nice ' . $config->{CONVERT_BLAST};
$BLAST2BED        = 'nice ' . catfile($SCRIPTDIR, 'blast2bed.pl');
$BLAST2RAW        = 'nice ' . $config->{BLAST2RAW};
$SYNTENY_SCORE    = 'nice ' . $config->{SYNTENY_SCORE};
$PYTHON26         = $config->{PYTHON};
$DATASETGROUP2BED = 'nice ' . $config->{DATASETGROUP2BED};
$COOKIE_NAME      = $config->{COOKIE_NAME};

if ( $FORM->param('get_master') ) {
    get_master_syn_sets();
    exit;
}
if ( $FORM->param('ug') ) { # generating unique gene lists
    get_unique_genes();
    exit;
}

my %FUNCTION = (
    get_orgs                => \&get_orgs,
    gen_dsg_menu            => \&gen_dsg_menu,
    get_genome_info         => \&get_genome_info,
    get_dsg_for_menu        => \&get_dsg_for_menu,
#    generate_basefile       => \&generate_basefile,
    search_lists            => \&search_lists,
    get_list_preview        => \&get_list_preview,
    get_genomes_for_list    => \&get_genomes_for_list,
    get_types               => \&get_types,
    cogefeatsearch          => \&cogefeatsearch,
    get_anno                => \&get_anno,
    get_orgs_feat           => \&get_orgs_feat,
    source_search           => \&source_search,
    generate_feat_info      => \&generate_feat_info,
    save_settings           => \&save_settings,
    go_synfind              => \&go_synfind,
    get_results             => \&get_results,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my ($body) = gen_body();

    my $template = HTML::Template->new( filename => $config->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'SynFind: Syntenic Compiler',
                      PAGE_TITLE => 'SynFind',
                      PAGE_LINK  => $LINK,
                      SUPPORT_EMAIL => $config->{SUPPORT_EMAIL},
                      HOME       => $config->{SERVER},
                      HELP       => 'SynFind',
                      WIKI_URL   => $config->{WIKI_URL} || '',
		              USER       => $USER->display_name || '' );

    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin,
                      CAS_URL    => $config->{CAS_URL} || '',
                      COOKIE_NAME => $config->{COOKIE_NAME} || '' );

    my $prebox = HTML::Template->new( filename => $config->{TMPLDIR} . 'SynFind.tmpl' );
    $prebox->param( RESULTS_DIV => 1 );
    $template->param( PREBOX => $prebox->output );

    return $template->output;
}

sub gen_body {
    my $template = HTML::Template->new( filename => $config->{TMPLDIR} . 'SynFind.tmpl' );
    #$template->param( BETA => 1 );

    #comparison algorithm
    my $algo;
    $algo = $FORM->param('algo') if $FORM->param('algo');
    $algo = "last" unless $algo;
    if ( $algo eq "last" ) {
        $template->param( 'LAST' => "selected" );
    }
    elsif ( $algo eq "lastz" ) {
        $template->param( 'LASTZ' => "selected" );
    }

    #synteny_score parameters
    #synteny_score window size
    my $ws;
    $ws = $FORM->param('ws') if $FORM->param('ws');
    $ws = 40 unless defined $ws;
    $template->param( WS => $ws );

    #synteny_score cutoff
    my $co;
    $co = $FORM->param('co') if $FORM->param('co');
    $co = 4 unless defined $co;
    $template->param( CO => $co );

    #synteny_score scoring function
    my $sf;
    $sf = $FORM->param('sf') if $FORM->param('sf');
    $sf = 1 unless defined $sf;
    $template->param( SF_COLLINEAR => "selected" ) if $sf == 1;
    $template->param( SF_DENSITY   => "selected" ) if $sf == 2;

    #syntenic depth reporting maximum
    my $sd;
    $sd = $FORM->param('sd') if $FORM->param('sd');
    $template->param( 'SD' => $sd ) if $sd;
    $template->param( JAVASCRIPT => 1,
                      FRONT_PAGE => 1 );

    #set up org search
    $template->param( ORG_LIST => get_orgs() );

    #set up feat search
    my $accn;
    $accn = $FORM->param('accn');
    $template->param( ACCN => $accn );

    #$template->param(FEAT_TYPE=> get_feature_types());
    $template->param( ORG_LIST_FEAT => get_orgs_feat() );
    my $doc_ready;
    $doc_ready .= qq{search_chain(1);\n} if $accn;

    # Load saved dsgids
    my $prefs = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    $prefs = {} unless $prefs;

    my $gid = $FORM->param('dsgid') || $FORM->param('gid');
    if ($gid) {
        $template->param(TARGETS => $gid);
        foreach my $item ($gid) {
            foreach my $dsgid ( split( /,/, $item ) ) {
                my $id = get_dsg_for_menu( dsgid => $dsgid );
                $doc_ready .= qq{add_to_list('$id');};
            }
        }
    }
    elsif ( $prefs->{dsgids} ) {
        foreach my $dsgid ( split( /,/, $prefs->{dsgids} ) ) {
            my $id = get_dsg_for_menu( dsgid => $dsgid );
            $doc_ready .= qq{add_to_list('$id');};
        }
    }

    if ( my $fid = $FORM->param('fid') ) {
        my ( $fid, $seq_type_id ) = split( /_/, $fid );
        $seq_type_id = 1 unless $seq_type_id;
        my $qdsgid = $FORM->param('qdsgid');
        unless ($qdsgid) {
            my $feat = $coge->resultset('Feature')->find($fid);
            foreach my $dsg ( $feat->genomes ) {
                if ( $dsg->type->id eq $seq_type_id ) # find the unmasked version of the data
                {
                    $qdsgid = $dsg->id;
                    last;
                }
            }
            $qdsgid = $feat->genomes->[0]->id unless $qdsgid;
        }
        $doc_ready .= qq{get_anno_chain('', $qdsgid, $fid);}; #qq{get_anno(['args__fid','args__$fid','args__dsgid','args__$qdsgid'],[show_anno]);};
        $template->param( FEAT_DSGID => qq{<select MULTIPLE type='hidden' id='feat_dsgid' name='feat_dsgid'><option value='$qdsgid' selected></option></select>} )
          if $qdsgid;

        $template->param( FID => $fid);
    }

    $template->param( DOCUMENT_READY => $doc_ready ) if $doc_ready;
    $template->param( SAVE_ORG_LIST => 1 ) unless $USER->user_name eq "public";

    $template->param( RUN => 1) if $FORM->param('run');

    return $template->output;
}

#sub generate_basefile {
#    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
#    return $cogeweb->basefilename;
#}

sub get_orgs { #FIXME: dup'ed in CoGeBlast.pl
    my %opts      = @_;
    my $name_desc = $opts{name_desc};
    my $html_only = $opts{html_only}; # optional flag to return html instead of JSON
    my $timestamp = $opts{timestamp};
#    print STDERR "get_orgs: " . ($name_desc ? $name_desc : '') . "\n";

    my $html;
    my @organisms;
    if ($name_desc) {
        $name_desc = '%' . $name_desc . '%';
        my @organisms = $coge->resultset("Organism")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $name_desc ],
                [ 'description', $name_desc ]
            ]
        );

        my @opts;
        foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @organisms ) {
            push @opts,
                "<OPTION value=\""
              . $item->id
              . "\" id=\"o"
              . $item->id . "\">"
              . $item->name
              . "</OPTION>";
        }

        if (@opts) {
            if ( @opts <= $MAX_SEARCH_RESULTS ) {
                $html .= join( "\n", @opts );
            }
            else {
                $html .= "<option id='null_org' style='color:gray;' disabled='disabled'>Too many results to display, please refine your search.</option>";
            }
        }
        else {
            $html .= "<option id='null_org' style='color:gray;' disabled='disabled'>No results</option>";
        }
    }
    else {
        $html .= "<option id='null_org' style='color:gray;' disabled='disabled'>Please enter a search term</option>";
    }

    return $html if ($html_only);
    return encode_json( { timestamp => $timestamp, html => $html } );
}

sub gen_dsg_menu {
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $dsgid = $opts{dsgid};

    my $favorites = CoGe::Core::Favorites->new(user => $USER);

    my @genomes;
    foreach my $dsg (
        sort { genomecmp2($a, $b, $favorites) } $coge->resultset('Genome')->search(
            { organism_id => $oid },
            { prefetch    => ['genomic_sequence_type'] }
        )
      )
    {
        next unless $USER->has_access_to_genome($dsg);
        next if $dsg->deleted; #skip deleted genomes
        
        my $name;
        $name .= "&#11088; " if ($favorites->is_favorite($dsg));
        $name .= "&#x2705; " if $dsg->certified;
        $name .= "&#x1f512; " if $dsg->restricted;
        
        my $has_cds = has_cds( $dsg->id );
        $name .= " NO CDS ANNOTATIONS - CAN'T BE USED: " unless $has_cds;
        $name .= " Genome is still being loaded - CAN'T BE USED: " if $dsg->status && $dsg->status == LOADING;
        $name .= " There was an error while loading this genome - CAN'T BE USED: " if $dsg->status && $dsg->status == ERROR;
        
        $dsgid = $dsg->id unless $dsgid;
	    $name .= " (id ". $dsg->id.") ";
        $name .= $dsg->name . ", " if $dsg->name; # : $dsg->datasets->[0]->name;
        $name .= "v"
          . $dsg->version . " "
          . $dsg->type->name . " "
          . commify( $dsg->length ) . "nt";
        $name .= " (Source: ".join( ", ", map { $_->name } $dsg->source ).")";
        push @genomes, [ $dsg->id, $name ];
    }
    my $size = scalar @genomes;
    $size = 5 if $size > 5;

    my $dsg_menu = '';
    if (@genomes) {
        foreach (@genomes) {
            my ( $numt, $name ) = @$_;
            my $selected = ( $dsgid && $numt == $dsgid ? 'selected' : '' );
            $dsg_menu .= qq{<option value='$numt' $selected>$name</option>};
        }
    }

    return $dsg_menu;
}

sub get_dsg_for_menu { #FIXME: dup'ed in CoGeBlast.pl
    my %opts   = @_;
    my $dsgids = $opts{dsgid};
    my $orgids = $opts{orgid};
    my %dsgs;
#   print STDERR "get_dsg_for_menu: dsgids=" . ($dsgids ? $dsgids : '') . " orgids=" . ($orgids ? $orgids : '') . "\n";

    if ($orgids) {
        my @orgids = split( /,/, $orgids );
        foreach my $dsg ( $coge->resultset('Genome')->search( { organism_id => [@orgids] } ) )
        {
            next unless $USER->has_access_to_genome($dsg);
	        next if $dsg->deleted;
            $dsgs{ $dsg->id } = $dsg;
        }
    }

    if ($dsgids) {
        %dsgs = () if ( $dsgs{$dsgids} );
        foreach my $dsgid ( split( /,/, $dsgids ) ) {
            my $dsg = $coge->resultset('Genome')->find($dsgid);
            next unless $USER->has_access_to_genome($dsg);
            next if $dsg->deleted;
            $dsgs{ $dsg->id } = $dsg;
        }
    }

    my $html;
    foreach my $dsg ( values %dsgs ) {
        next unless has_cds( $dsg->id ); # skip if it has no CDS annotations
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
    }

    return html_escape($html);
}

sub get_genome_info { #FIXME: dup'ed in CoGeBlast.pl
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    #print STDERR "get_genome_info: $dsgid\n";
    return " " unless $dsgid;

    my $dsg = $coge->resultset("Genome")->find($dsgid);
    return "Unable to create genome object for id: $dsgid" unless $dsg;

    my $html = qq{<table class='small'>}; # = qq{<div style="overflow:auto; max-height:78px">};
    $html .= qq{<tr valign='top'><td style='white-space:nowrap'>Name:<td><span class='link' onclick=window.open('OrganismView.pl?dsgid=$dsgid')>}
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
    $link = "http://" . $link unless ($link && $link =~ /^http/);
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
#   my $chr = $gs->chromosome;
#   my $length = $gs->sequence_length;
#   $length = commify($length);
#   $html .= qq{$chr:  $length bp<br>};
#     }
#    $html .= "</table>";

    return $html;
}

sub get_types {
    my %opts = @_;
    #my ($dsid, $gstid) = split(/_/, $opts{dsid});
    my $dsgid = $opts{dsgid};
    my $accn  = $opts{accn};
    my $ftid  = $opts{ftid};
    #print STDERR "get_types: dsgid=$dsgid accn=$accn ftid=$ftid\n";

    my $html;
    my $blank = qq{<font class='small note'>(0)</font><br><SELECT id="type_name" size="10" style="min-width:60px"></SELECT>};#<input type="hidden" id="type_name">};
    my %seen;
    my $search = {
        'feature_names.name'           => $accn,
        'dataset_connectors.genome_id' => $dsgid,
    };
    $search->{feature_type_id} = 3;

    my @opts =
      sort map { "<OPTION>$_</OPTION>" }
      grep     { !$seen{$_}++ }
      map      { $_->type->name }
          $coge->resultset('Feature')->search(
            $search,
            ,
            { join => [ 'feature_names', { 'dataset' => 'dataset_connectors' } ] }
          );

    if (@opts) {
        $html .= "<font class='small note'>(" . scalar @opts . ")</font>\n<BR>\n"
            . qq{<SELECT id="type_name" size="10" style="min-width:60px" MULTIPLE }
            . qq{onChange="get_anno_chain('', $dsgid)" >\n}
            . join( "\n", @opts )
            . "\n</SELECT>\n";
        $html =~ s/OPTION/OPTION SELECTED/;
    }
    else {
        return encode_json({ html => $blank });
    }

    #return $blank unless $html =~ /OPTION/;
    return encode_json({ html => $html, dsgid => $dsgid });
}

sub search_lists {   # FIXME this coded is dup'ed in User.pl and NotebookView.pl and CoGeBlast.pl
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
#        $num_results = $coge->resultset("List")->count_literal($sql);
#        if ( $num_results < $MAX_SEARCH_RESULTS ) {
#            foreach my $notebook ( $coge->resultset("List")->search_literal($sql) )
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
        foreach my $notebook ($coge->resultset("List")->search_literal("locked=0 AND (name LIKE '$search_term' OR description LIKE '$search_term')"))
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

    my $list = $coge->resultset('List')->find($lid);
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

    my $list    = $coge->resultset('List')->find($lid);
    my $genomes = '';
    if ($list) {
        my $dsgids = join( ',', map { $_->id } $list->genomes );
        $genomes = get_dsg_for_menu( dsgid => $dsgids );
    }

    return $genomes;
}

sub cogefeatsearch {
    my %opts = @_;
    my $accn = $opts{accn};
    $accn =~ s/^\s+//;
    $accn =~ s/\s+$//;
    my $anno = $opts{anno};
    $anno =~ s/^\s+//;
    $anno =~ s/\s+$//;
    my $fid      = $opts{fid};
    my $type     = $opts{type};
    my $org_id   = $opts{org_id};
    my $org_name_desc = $opts{org_name_desc};
    my $blank = qq{<input type="hidden" id="accn_select"><select MULTIPLE type='hidden' id='feat_dsgid'></select>};
    my $weak_query = "Query needs to be better defined.";
    #print STDERR "cogefeatsearch: accn=$accn anno=$anno fid=$fid\n";

    if ( !$accn && !$anno && !$fid ) {
        return $weak_query . $blank unless $org_id && $type;
    }
    my $html;
    my %seen;
    my $search = {};
    $search->{feature_type_id} = 3;    #$type if $type;
    $search->{'organism.name'} = { like => "%" . $org_name_desc . "%" } if $org_name_desc;
    $search->{'organism.description'} = { like => "%" . $org_name_desc . "%" } if $org_name_desc;
    my $join = {
        join => [
            { 'feature' => {
                    'dataset' =>
                      { 'dataset_connectors' => { 'genome' => 'organism' } }
                }
            }
        ]
    };
    push @{ $join->{join} }, 'feature_annotation' if $anno; #=>['annotation',]};

    #trying to get fulltext to work (and be fast!)
    my @names;
    if ($accn) {
        push @names,
          $coge->resultset('FeatureName')->search( $search, $join )->search_literal( 'MATCH(me.name) AGAINST (?)', $accn );
    }
    if ($anno) {
        push @names,
          $coge->resultset('FeatureName')->search( $search, $join )->search_literal( 'MATCH(annotation) AGAINST (?)', $anno );
    }
    my @opts;
    foreach my $name (@names) {
        my $item = $name->name;
        next if $seen{ uc($item) };
        $seen{ uc($item) }++;
        push @opts, "<OPTION>$item</OPTION>" if $name->feature->dataset;
    }
    if ( @opts > 10000 ) {
        return $blank . "Search results over 10000, please refine your search.\n";
    }
    $html .= "<font class='small note'>(" . scalar @opts . ")</font>\n<BR>\n"
        . qq{<SELECT id="accn_select" SIZE="10" MULTIPLE onChange="source_search_chain();">\n}
        . join( "\n", @opts )
        . "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank . "No results found\n" unless $html =~ /OPTION/;
    return $html;
}

sub get_anno {
    my %opts  = @_;
    my $accn  = $opts{accn};
    my $fid   = $opts{fid};
    my $type  = $opts{type};
    my $dsgid = $opts{dsgid};
    return unless $accn || $fid;
    print STDERR "get_anno: " . ($accn ? "accn=$accn " : '') . ($fid ? "fid=$fid " : '') . ($type ? "type=$type " : '') . "dsgid=$dsgid\n";

    my $genome = $coge->resultset('Genome')->find($dsgid);
    return unless ( $USER->has_access_to_genome($genome) );

    my @feats;
    if ($accn) {
        foreach my $feat (
            $coge->resultset('Feature')->search(
                {   'feature_names.name'           => $accn,
                    "dataset_connectors.genome_id" => $dsgid
                },
                {   join => [ 'feature_names', { 'dataset' => 'dataset_connectors' } ] }))
        {
            push @feats, $feat if ( $feat->type->name eq $type );
        }
    }
    elsif ($fid) {
        my $feat = $coge->resultset('Feature')->find($fid);
        push @feats, $feat if $feat;
    }
    return unless @feats;

    my $anno;
    $anno .= "<div style='border-bottom:1px solid gray;'><strong>Annotations</strong> <font class='small note'>(" . scalar @feats . ")</font></div><br>\n"
        if scalar @feats;
    my $i   = 0;
    my $dsg = $coge->resultset('Genome')->find($dsgid);
    unless ($dsg) {
        foreach my $tdsg ( $feats[0]->genomes ) {
            $dsg = $tdsg if $tdsg->type->id eq "1"; # default to unmasked sequence if possible
        }
        ($dsg) = $feats[0]->genomes unless $dsg;
    }
    my $gstid = $dsg->type->id;
    foreach my $feat (@feats) {
        my ($genome) = $feat->genomes;
        return encode_json({ anno => 'restricted', fid => 0 }) unless ($USER->has_access_to_genome($genome));
        #next if $feat->dataset->restricted && !$USER->has_access_to_dataset($feat->dataset);
        $i++;

        #my $featid = $feat->id;
        #my $chr = $feat->chr;
        #my $rc = 0;
        #my $pro = 0;
        #my $ds = $feat->dataset->id;
        #my $x = $feat->start;
        #my $z = 4;
        #$anno .= qq{<span class="ui-button ui-corner-all" onClick="window.open('FastaView.pl?featid=$featid&gstid=$gstid');">Get Sequence</span>};
        #$anno .= qq{<span class="ui-button ui-corner-all" onClick="window.open('CoGeBlast.pl?featid=$featid;gstid=$gstid');">CoGeBlast</span>};
        #$anno .= qq{<span class="ui-button ui-corner-all" onClick="window.open('GenomeView.pl?chr=$chr&ds=$ds&x=$x&z=$z;gstid=$gstid');">Genome Browser</span>};
        $anno .= '<span class="small">Organism: </span>';
        foreach my $item ( split( /;/, $feat->organism->description ) ) {
            $anno .= qq{<span class='small note link' onclick="\$('#org_desc').val('$item').focus()">$item;</span>};
        }
        $anno .= '<br><br>';
        $anno .= join("\n<BR><HR><BR>\n", $feat->annotation_pretty_print_html( gstid => $gstid ));
    }
    $anno = "<font class=\"annotation\">No annotations for this entry</font>" unless $anno;
    return encode_json({ anno => $anno, fid => $feats[0]->id });
}

sub get_orgs_feat {
    my %opts    = @_;
    my $search  = $opts{search};
    my $id_only = $opts{id_only};

    #print STDERR "get_orgs_feat: search=$search\n";

    my @organisms;
    my $count;

    if ($search) {
        $search = '%' . $search . '%';
        @organisms = $coge->resultset("Organism")->search(
            \[
                'name LIKE ? OR description LIKE ?',
                [ 'name',        $search ],
                [ 'description', $search ]
            ]
        );
        $count = @organisms;
    }
    else {
        $count = $coge->resultset("Organism")->count();
    }
    return map { $_->id } @organisms if $id_only;

    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @organisms ) {
        push @opts,
            "<OPTION value=\"" . $item->id . "\" id=\"o" . $item->id . "\">" . $item->name . "</OPTION>";
    }
    my $html;
    $html .= qq{<div class="note" id="org_count">Matches: } . $count . qq{</div>\n};
    if ( $search && !@opts ) {
        $html .= qq{<input type = hidden name="org_id_feat" id="org_id_feat"><br>} . "No results";
        return $html;
    }
    unshift( @opts, "<OPTION value=\"all\" id=\"all\">All Listed Organisms</OPTION>" );
    my $size = scalar @opts;
    $size = 8 if $size > 8;
    $html .= qq{<SELECT id="org_id_feat" SIZE="$size" MULTIPLE >\n}
         . join( "\n", @opts )
         . "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;

    return $html;
}

sub source_search {
    my %opts     = @_;
    my $accn     = $opts{accn};
    my $org_id   = $opts{org_id};
    my $org_name_desc = $opts{org_name_desc};
    $org_id = "all" unless $org_id;
    #print STDERR "source_search: accn=$accn org_id=$org_id org_name=$org_name org_desc=$org_desc\n";

    my %org_ids;
    if ( $org_id eq "all" ) {
        %org_ids = map { $_ => 1 }
            get_orgs( id_only => 1, search => $org_name_desc )
                if $org_name_desc;
    }
    else {
        $org_ids{$org_id} = 1;
    }
    my $blank = qq{<input type="hidden" id="dsid">------------};
    return $blank unless $accn;
    my @feats = $coge->resultset('Feature')->search(
        { 'feature_type_id' => 3, 'feature_names.name' => $accn },
        {
            join       => 'feature_names',
            'prefetch' => {
                'dataset' => [
                    'data_source',
                    {
                        'dataset_connectors' => {
                            'genome' => [ 'organism', 'genomic_sequence_type' ]
                        }
                    }
                ]
            }
        }
    );
    my %sources;
    foreach my $feat (@feats) {

#	my ($genome) = $feat->genomes;
#	next unless $USER->has_access_to_genome($genome);
#	next if $feat->dataset->restricted && !$USER->has_access_to_dataset($feat->dataset);
        foreach my $dsg ( $feat->dataset->genomes ) {
            next if $dsg->deleted;
            next unless $USER->has_access_to_genome($dsg);
            my $org = $dsg->organism->name if $dsg->organism;
            if ( keys %org_ids ) { next unless $org_ids{ $dsg->organism->id }; }
            my $name  = $dsg->name;
            my $ver   = $dsg->version;
            my $desc  = $dsg->description;
            my $sname = $feat->dataset->data_source->name
              if $feat->dataset->data_source;
            my $source_id = $feat->dataset->data_source_id
              if $feat->dataset->data_source;
            my $dsgid   = $dsg->id;
            my $gstname = $dsg->sequence_type->name;
            my $title   = "$org: ";
            $title .= $name if $name;
            $title .= "($sname, v$ver, $gstname, dsgid$dsgid)";
            $sources{$title} = {
                id       => $dsgid,
                v        => $ver,
                gstid    => $dsg->sequence_type->id,
                sourceid => $source_id
            };
        }
    }
    my $html;
    $html .= qq{<SELECT MULTIPLE name="feat_dsgid" id="feat_dsgid" size="10" onChange="get_types_chain();">};
    my $count = 0;
    foreach my $title (
        sort {
                 $sources{$b}{v} <=> $sources{$a}{v}
              || $sources{$a}{gstid} <=> $sources{$b}{gstid}
              || $sources{$a}{sourceid} <=> $sources{$b}{sourceid}
              || $a cmp $b
        } keys %sources
      )
    {
        my $id = $sources{$title}{id};

        #	my $gstid = $sources{$title}{gstid};
        #	my $val = $id."_".$gstid;
        $html .= qq{  <option value="$id">$title\n};
        $html =~ s/option/option selected/ unless $count;
        $count++;
    }
    $html .= qq{</SELECT>\n};
    return ("<font class='small note'>(" . $count . ")</font>\n<BR>\n" . $html );
}

sub get_feature_types {
    my %args  = @_;
    my $orgid = $args{orgid};
    my $html  = qq{<select  name="type" id="type" />
<OPTION VALUE=0>All</OPTION>
};
    if ($orgid) {
        my @tmp = $coge->resultset('FeatureType')->search(
            { organism_id => $orgid, },
            {
                distinct => ['me.name'],
                join     => { 'features' => 'dataset' }
            }
        );
        $html .= join(
            "\n",
            map { "<OPTION value=\"" . $_->id . "\">" . $_->name . "</OPTION>" }
              sort { uc( $a->name ) cmp uc( $b->name ) } @tmp
        );
    }
    else {
        $html .= join(
            "\n",
            map { "<OPTION value=\"" . $_->id . "\">" . $_->name . "</OPTION>" }
              sort { uc( $a->name ) cmp uc( $b->name ) }
              $coge->resultset('FeatureType')->all
        );
    }
    return $html;
}

sub go_synfind {
    my %opts             = @_;
    my $dsgids           = $opts{dsgids};
    my $fid              = $opts{fid};
    my $source_dsgid     = $opts{qdsgid};
    my $algo             = $opts{algo};
    #my $basename         = $opts{basename}; # mdb removed 3/6/14
    my $window_size      = $opts{window_size};
    my $cutoff           = $opts{cutoff};
    my $scoring_function = $opts{scoring_function};
    my $depth            = $opts{depth}; # max syntenic depth to report

    $window_size = 40  unless defined $window_size;
    $cutoff      = 0.1 unless defined $cutoff;
    $cutoff      = "0" . $cutoff if $cutoff =~ /^\./; # add the prepending 0 for filename consistence
    $scoring_function = 1 unless defined $scoring_function;

    #need to set this link before the scoring function changes to a different name type
    my $synfind_link = $SERVER . "SynFind.pl?fid=$fid;ws=$window_size;co=$cutoff;sf=$scoring_function;algo=$algo";
    my $unique_gene_link = $synfind_link;
    $synfind_link .= ";dsgid=$dsgids;qdsgid=$source_dsgid";
    $synfind_link .= ";sd=$depth" if $depth;
    $synfind_link .= ";run=1";

    my @ids = split( /,/, $dsgids );

    my $genomes_url = CoGe::Accessory::Web::get_tiny_link(
        user_id => $USER->id,
        page    => "GenomeList",
        url     => $config->{SERVER} . "/GenomeList.pl?dsgid=$dsgids"
    );

    my $list_link = qq{<a href="$genomes_url" target="_blank">} . @ids . ' genome' . ( @ids > 1 ? 's' : '' ) . '</a>';

    my $feat = $coge->resultset('Feature')->find($fid);
    my $feat_url = "<a href='FeatView.pl?fid=$fid' target='_blank'>" . $feat->name . "</a>";

    my ( $source_name, $titleA ) = gen_org_name( dsgid => $source_dsgid, write_log => 0 );

    #my $name = "<a href='OrganismView.pl?dsgid=$source_dsgid' target='_blank'>$source_name</a>";
    my $log_msg = 'Searched ' . $list_link . '  for feature ' . $feat_url;

    my $tiny_synfind_link = CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $synfind_link
    );

    my ($tiny_id) = $tiny_synfind_link =~ /\/(\w+)$/;
    my $workflow_name = "synfind-$tiny_id";

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );

    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Creating Workflow", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Link to Regenerate Analysis", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "$tiny_synfind_link", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "",           $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Created Workflow: $workflow_name", $cogeweb->logfile );

    #convert numerical codes for different scoring functions to appropriate types
    if ( $scoring_function eq '2' ) {
        $scoring_function = "density";
    }
    else {
        $scoring_function = "collinear";
    }

    #need to blast source_dsg against each dsgids
    my @blast_results;
    my $html;

    ###########################################################################
    # Setup workflow
    ###########################################################################
    my $workflow = $JEX->create_workflow(
        id      => 0,
        name    => $workflow_name,
        logfile => $cogeweb->logfile
    );

    my ( $query_info, @target_info );
    foreach my $dsgid ( $source_dsgid, split( /,/, $dsgids ) ) {
        #check for taint
        ($dsgid) = $dsgid =~ /(\d+)/;

        my $has_cds   = has_cds($dsgid);
        my $feat_type = $has_cds ? "CDS" : "genomic";
        my $fasta     = $FASTADIR . "/$dsgid-$feat_type.fasta";

        if ( $dsgid && $has_cds ) {
            my ( $org_name, $title ) = gen_org_name(
                dsgid     => $dsgid,
                feat_type => $feat_type,
                write_log => 1
            );

            push @target_info,
              {
                dsgid     => $dsgid,
                feat_type => $feat_type,
                fasta     => $fasta,
                org_name  => $org_name,
              };
        }
        else {
            my $dsg = $coge->resultset('Genome')->find($dsgid);

            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "WARNING:", $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( $dsg->organism->name . " does not have CDS sequences.  Can't process in SynFind.\n", $cogeweb->logfile );
            next;
        }

        my $fasta_args = [
            [ "--config",       $config->{_CONFIG_PATH}, 0 ],
            [ "--genome_id",    $dsgid,     1 ],
            [ "--feature_type", $feat_type, 1 ],
            [ "--fasta",        $fasta,     1 ]
        ];

        $workflow->add_job({
            cmd         => $GEN_FASTA,
            description => "Generating fasta file...",
            script      => undef,
            args        => $fasta_args,
            inputs      => undef,
            outputs     => [$fasta]
        });

        my $bed_args = [
            [ " -cf ",  $config->{_CONFIG_PATH}, 0 ],
            [ '-dsgid', $dsgid,                    1 ],
            [ '>',      $BEDDIR . $dsgid . ".bed", 1 ]
        ];

        $workflow->add_job({
            cmd         => $DATASETGROUP2BED,
            description => "Creating bed files...",
            script      => undef,
            args        => $bed_args,
            inputs      => undef,
            outputs     => [ $BEDDIR . $dsgid . ".bed" ]
        });
    }

    #query is the first item on this list.
    $query_info = shift @target_info;

    foreach my $target (@target_info) {
        my ( $org1, $org2 ) = ( $query_info->{org_name}, $target->{org_name} );
        my ( $dsgid1, $dsgid2 ) = ( $query_info->{dsgid}, $target->{dsgid} );
        my ( $feat_type1, $feat_type2 ) = ( $query_info->{feat_type}, $target->{feat_type} );
        my ( $fasta1, $fasta2 ) = ( $query_info->{fasta}, $target->{fasta} );

        #my ($db1, $db2) = ($query_info->{blastdb},$target->{blastdb});
        ( $org1, $org2, $dsgid1, $dsgid2, $feat_type1, $feat_type2, $fasta1, $fasta2 )
            = ($org2, $org1, $dsgid2, $dsgid1, $feat_type2, $feat_type1, $fasta2, $fasta1
        ) if ( $dsgid2 lt $dsgid1 );

        my $basedir = $DIAGSDIR . "/" . $dsgid1 . "/" . $dsgid2;
        mkpath( $basedir, 0, 0777 ) unless -d $basedir;
        my $basename = $dsgid1 . "_" . $dsgid2 . "." . $feat_type1 . "-" . $feat_type2;
        my $blastfile = $basedir . "/" . $basename . ".$algo";
        my $bedfile1  = $BEDDIR . $dsgid1 . ".bed";
        my $bedfile2  = $BEDDIR . $dsgid2 . ".bed";
        $target->{synteny_score_db}    = $basedir . "/" . $basename . "_" . $window_size . "_" . $cutoff . "_" . $scoring_function . ".$algo" . ".v053" . ".db"; #v053 is the version of synteny_score
        $target->{basedir}             = $basedir;
        $target->{basename}            = $basename;
        $target->{blastfile}           = $blastfile;
        $target->{converted_blastfile} = $blastfile . ".short_names";
        $target->{filtered_blastfile}  = $target->{converted_blastfile} . ".filtered";
        $target->{bedfile1}    = $bedfile1;
        $target->{bedfile2}    = $bedfile2;
        $target->{dsgid1}      = $dsgid1;
        $target->{dsgid2}      = $dsgid2;
        $target->{query_fasta} = $fasta1
          ; #need to determine the correct query/target order so output is compatibable with SynMap

        #$target->{target_db} = $db2;#need to determine the correct query/target order so output is compatibable with SynMap
        $target->{target_fasta} = $fasta2;
    }

    #blast and conquer
    foreach my $target (@target_info) {
        #######################################################################
        # Blast
        #######################################################################
        if ( $algo =~ /lastz/i ) {
            $workflow->add_job({
                cmd         => $LASTZ,
                description => "Running blast ($algo) algorithm...",
                script      => undef,
                args        => [
                    [ "-i", $target->{query_fasta},  1 ],
                    [ "-d", $target->{target_fasta}, 1 ],
                    [ "-o", $target->{blastfile},    1 ]
                ],
                inputs      => [ $target->{query_fasta}, $target->{target_fasta} ],
                outputs     => [ $target->{blastfile} ]
            });
        }
        else { # lastal
            my $dbpath = catdir($config->{LASTDB}, $target->{dsgid2});
            mkpath($dbpath, 0, 0777);

# mdb removed 3/21/16 for LAST v731
#            $blast_args = [
#                [ "--dbpath", $dbpath, 0 ],
#                [ "",   $target->{target_fasta}, 1 ],
#                [ "",   $target->{query_fasta},  1 ],
#                [ "-o", $target->{blastfile},    1 ]
#            ];
#            $blast_cmd = $LAST;

            
            # mdb added 3/21/16 for LAST v731
            #$dbpath = catfile($dbpath, basename($target->{query_fasta})); #Eric is changing all query to target 4/10/16
            $dbpath = catfile($dbpath, basename($target->{target_fasta}));
            my @dbfiles = map { $dbpath . '.' . $_ } ('bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis'); 
            $workflow->add_job({
                cmd         => "$LASTDB",
                description => "Creating database for $algo...",
                script      => undef,
                args        => [
                    ['', $dbpath, 0],
                    ['', $target->{target_fasta}, 1]
                ],
                inputs      => [ $target->{target_fasta} ],
                outputs     => \@dbfiles
            });
            #Eric is changing all target to query 4/10/16
            $workflow->add_job({
                cmd         => "$LAST $dbpath $target->{query_fasta} > $target->{blastfile} ; touch $target->{blastfile}.done ;",
                description => "Running blast ($algo) algorithm...",
                script      => undef,
                args        => [],
                inputs      => [ $target->{query_fasta}, @dbfiles],
                outputs     => [ $target->{blastfile}, "$target->{blastfile}.done" ]
            });
        }

        #######################################################################
        # Convert Blast
        #######################################################################
        my $convert_args = [
            [ '<', $target->{blastfile},           1 ],
            [ '>', $target->{converted_blastfile}, 1 ]
        ];
        my $convert_inputs = [ $target->{blastfile}, "$target->{blastfile}.done" ];
        my $convert_outputs = [ $target->{converted_blastfile}, ];

        $workflow->add_job({
            cmd         => $CONVERT_BLAST,
            description => "Converting blast file to short names...",
            script      => undef,
            args        => $convert_args,
            inputs      => $convert_inputs,
            outputs     => $convert_outputs
        });

        #######################################################################
        # Blast 2 Raw
        #######################################################################
        my $raw_args = [
            [ "",              $target->{converted_blastfile}, 1 ],
            [ "--qbed",        $target->{bedfile1},            1 ],
            [ "--sbed",        $target->{bedfile2},            1 ],
            [ "--tandem_Nmax", 10,                             1 ],
            [ ">",             $target->{filtered_blastfile},  1 ],
        ];

        my $raw_inputs = [
            $target->{bedfile1}, $target->{bedfile2},
            $target->{converted_blastfile}
        ];

        my $raw_outputs = [
            $target->{filtered_blastfile},
            #$target->{filtered_blastfile} . ".q.localdups",
            #$target->{filtered_blastfile} . ".s.localdups",
        ];

        $workflow->add_job({
            cmd         => $BLAST2RAW,
            description => "Finding and removing local duplications...",
            script      => undef,
            args        => $raw_args,
            inputs      => $raw_inputs,
            outputs     => $raw_outputs
        });

        #######################################################################
        # Synteny Score
        #######################################################################
        $cutoff = sprintf( "%.2f", $cutoff / $window_size ) if ( $cutoff >= 1 );
        my $synteny_score_args = [
            [ '',          $target->{filtered_blastfile}, 1 ],
            [ '--qbed',    $target->{bedfile1},           1 ],
            [ '--sbed',    $target->{bedfile2},           1 ],
            [ '--window',  $window_size,                  1 ],
            [ '--cutoff',  $cutoff,                       1 ],
            [ '--scoring', $scoring_function,             1 ],
            [ '--qnote',   $target->{dsgid1},             1 ],
            [ '--snote',   $target->{dsgid2},             1 ],
	        #added new commands for synteny_score v0.5.3 Eric Lyons 4/29/2015
            [ '--qbedlift',$target->{bedfile1},           1 ],
            [ '--sbedlift',$target->{bedfile2},           1 ],
            [ '--lift',    $target->{converted_blastfile},1 ],
	        ####
            [ '--sqlite',  $target->{synteny_score_db},   1 ],
        ];

        my $synteny_score_inputs = [
            $target->{filtered_blastfile},
            $target->{bedfile1}, $target->{bedfile2},
	    $target->{converted_blastfile}
        ];

        my $synteny_score_outputs = [ $target->{synteny_score_db}, ];

        $workflow->add_job({
            cmd         => $SYNTENY_SCORE,
            description => "Running Synteny Score...",
            script      => undef,
            args        => $synteny_score_args,
            inputs      => $synteny_score_inputs,
            outputs     => $synteny_score_outputs,
        });
    }

    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Running Workflow", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    my $response = $JEX->submit_workflow($workflow);
    my $success = JSON::true;
    $success = JSON::false if lc($response->{status}) eq "error";

    my $log = CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        description => $log_msg,
        page        => $PAGE_TITLE,
        link        => $tiny_synfind_link,
        parent_id   => $response->{id},
        parent_type => 7 #FIXME magic numbe
    ) if $response and $response->{id};

    my $log_url = $cogeweb->logfile;
    $log_url =~ s/$TEMPDIR/$TEMPURL/;

    return encode_json({
        success => $success,
        logfile => $log_url,
        link    => $tiny_synfind_link,
        request => 'jex/status/' . $response->{id}
    });
}

sub get_results {
    my %opts             = @_;
    my $dsgids           = $opts{dsgids};
    my $fid              = $opts{fid};
    my $source_dsgid     = $opts{qdsgid};
    my $algo             = $opts{algo};
    #my $basename         = $opts{basename}; # mdb removed 3/6/14
    my $window_size      = $opts{window_size};
    my $cutoff           = $opts{cutoff};
    my $scoring_function = $opts{scoring_function};
    my $depth            = $opts{depth}; # max syntenic depth to report
    my $logfile          = $opts{logfile};
    $logfile =~ s/$TEMPURL/$TEMPDIR/ if $logfile;

    $window_size = 40  unless defined $window_size;
    $cutoff      = 0.1 unless defined $cutoff;
    $cutoff      = "0" . $cutoff if $cutoff =~ /^\./; # add the prepending 0 for filename consistence
    $scoring_function = 1 unless defined $scoring_function;

    #need to set this link before the scoring function changes to a different name type
    my $synfind_link = $SERVER . "SynFind.pl?fid=$fid;ws=$window_size;co=$cutoff;sf=$scoring_function;algo=$algo";
    my $unique_gene_link = $synfind_link;
    $synfind_link .= ";dsgid=$dsgids;qdsgid=$source_dsgid";
    $synfind_link .= ";sd=$depth" if $depth;

    my @ids = split( /,/, $dsgids );

    my $genomes_url = CoGe::Accessory::Web::get_tiny_link(
        user_id => $USER->id,
        page    => "GenomeList",
        url     => $config->{SERVER} . "/GenomeList.pl?dsgid=$dsgids"
    );

    my $list_link = qq{<a href="$genomes_url" target="_blank">} . @ids . ' genome' . ( @ids > 1 ? 's' : '' ) . '</a>';

    my $feat = $coge->resultset('Feature')->find($fid);
    my $feat_url = "<a href='FeatView.pl?fid=$fid' target='_blank'>" . $feat->name . "</a>";

    my ( $source_name, $titleA ) = gen_org_name( dsgid => $source_dsgid, write_log => 0 );

    #my $name = "<a href='OrganismView.pl?dsgid=$source_dsgid' target='_blank'>$source_name</a>";
    my $log_msg = 'Searched ' . $list_link . '  for feature ' . $feat_url;

    my $tiny_synfind_link = CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $synfind_link
    );

    #convert numerical codes for different scoring functions to appropriate types
    if ( $scoring_function eq '2' ) {
        $scoring_function = "density";
    }
    else {
        $scoring_function = "collinear";
    }

    my ($name, $path, $ext) = fileparse($logfile, qr{\.log}) if $logfile;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR, basename => $name);

    #need to blast source_dsg against each dsgids
    my @blast_results;
    my $html;

    my ( $query_info, @target_info, @invalid );
    foreach my $dsgid ( $source_dsgid, split( /,/, $dsgids ) ) {
        #check for taint
        ($dsgid) = $dsgid =~ /(\d+)/;

        my $has_cds   = has_cds($dsgid);
        my $feat_type = $has_cds ? "CDS" : "genomic";
        my $fasta     = $FASTADIR . "/$dsgid-$feat_type.fasta";
        my ( $org_name, $title ) = gen_org_name(
            dsgid     => $dsgid,
            feat_type => $feat_type,
        );

        if ( $dsgid && $has_cds ) {
            return encode_json({
                success => JSON::false,
                message => "A genome is missing a fasta file.",
            }) unless -r $fasta;

            push @target_info,
              {
                dsgid     => $dsgid,
                feat_type => $feat_type,
                fasta     => $fasta,
                org_name  => $org_name,
              };
        }
        else {
            push @invalid,
              {
                dsgid     => $dsgid,
                feat_type => $feat_type,
                org_name  => $org_name,
                fasta     => undef
              };
        }
    }

    #query is the first item on this list.
    $query_info = shift @target_info;

       
    foreach my $target (@target_info) {
        my ( $org1, $org2 ) = ( $query_info->{org_name}, $target->{org_name} );
        my ( $dsgid1, $dsgid2 ) = ( $query_info->{dsgid}, $target->{dsgid} );
        my ( $feat_type1, $feat_type2 ) = ( $query_info->{feat_type}, $target->{feat_type} );
        my ( $fasta1, $fasta2 ) = ( $query_info->{fasta}, $target->{fasta} );

        #my ($db1, $db2) = ($query_info->{blastdb},$target->{blastdb});
        ( $org1, $org2, $dsgid1, $dsgid2, $feat_type1, $feat_type2, $fasta1, $fasta2 )
            = ($org2, $org1, $dsgid2, $dsgid1, $feat_type2, $feat_type1, $fasta2, $fasta1
        ) if ( $dsgid2 lt $dsgid1 );

        my $basedir = $DIAGSDIR . "/" . $dsgid1 . "/" . $dsgid2;
        my $basename = $dsgid1 . "_" . $dsgid2 . "." . $feat_type1 . "-" . $feat_type2;
        my $blastfile = $basedir . "/" . $basename . ".$algo";
        my $bedfile1  = $BEDDIR . $dsgid1 . ".bed";
        my $bedfile2  = $BEDDIR . $dsgid2 . ".bed";
        $target->{synteny_score_db}    = $basedir . "/" . $basename . "_" . $window_size . "_" . $cutoff . "_" . $scoring_function . ".$algo" . ".v053" . ".db"; #v053 is the version of synteny_score
        #$target->{synteny_score_db}    = $basedir . "/" . $basename . "_" . $window_size . "_" . $cutoff . "_" . $scoring_function . ".$algo" . ".db";

        return encode_json({
            success => JSON::false,
            message => "A genome is missing a database."
        }) unless -r $target->{synteny_score_db};

        $target->{basedir}             = $basedir;
        $target->{basename}            = $basename;
        $target->{blastfile}           = $blastfile;
        $target->{converted_blastfile} = $blastfile . ".short_names";
        $target->{filtered_blastfile}  = $target->{converted_blastfile} . ".filtered";
        $target->{bedfile1}    = $bedfile1;
        $target->{bedfile2}    = $bedfile2;
        $target->{dsgid1}      = $dsgid1;
        $target->{dsgid2}      = $dsgid2;
        $target->{query_fasta} = $fasta1
          ; #need to determine the correct query/target order so output is compatibable with SynMap

        #$target->{target_db} = $db2;#need to determine the correct query/target order so output is compatibable with SynMap
        $target->{target_fasta} = $fasta2;
    }

    my ( $gevo_link, $matches ) = gen_gevo_link(
        fid         => $fid,
        dbs         => [ map { $_->{synteny_score_db} } @target_info ],
        window_size => $window_size,
        depth       => $depth
    );

    my $tiny_gevo_link = CoGe::Accessory::Web::get_tiny_link( url => $gevo_link );
    CoGe::Accessory::Web::write_log( "#TINY GEVO LINK: $tiny_gevo_link", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Finished!", $cogeweb->logfile );

    # make table of results
    # table to look them up later;
    my %dsgids = map { $_->{dsgid}, => 1 } @target_info;
    $dsgids{ $query_info->{dsgid} } = 1;

    my $synmap_algo_num;    #from SynMap's lookup table on algorithms

    if ( $algo eq "lastz" ) {
        $synmap_algo_num = 4;    #lastz
    }
    elsif ( $algo eq "last" ) {
        $synmap_algo_num = 6;    #last
    }
    my $synmap_link = "SynMap.pl?autogo=1;b=$synmap_algo_num;dsgid1=$source_dsgid;dsgid2=";
    my @res;
    my %open_all_synmap;
    my @homologs;
    my $count = 0;
    my $last_tdsgid;

    my $table_header .= qq{<table id='syntelog_table' class="ui-widget-content ui-corner-all">}
        . qq{<THEAD><tr>}
        . qq{<th>Organism}
        . qq{<th>Genome}
        . qq{<th>Type}
        . qq{<th>Name}
        . qq{<th>Chr}
        . qq{<th>Synteny Score}
        . qq{<th>SynMap}
        . qq{</tr></THEAD><TBODY>};

    my $table_content;

    foreach my $item ( [ $fid, "query" ], @$matches ) {
        my ( $tfid, $match_type, $synteny_score ) = @$item;
        my $feat = $coge->resultset('Feature')->find($tfid);
        my $dsg;

        foreach my $dsgt ( $feat->genomes ) {
            if ( $dsgids{ $dsgt->id } ) {
                #delete $dsgids{$dsgt->id};
                $dsg = $dsgt;
                last;
            }
        }

        next unless $dsg;

        $last_tdsgid = $dsg->id unless $last_tdsgid;
        $count = 0 unless $dsg->id eq $last_tdsgid;
        $last_tdsgid = $dsg->id;
        next if $depth && $count > $depth;
        $table_content .= qq{<tr><td>};
        my $name;
        if ( $match_type eq "S" ) {
            $match_type = "syntelog";
            ($name) = $feat->names;
            push @homologs, $feat;
        }
        elsif ( $match_type eq "query" ) {
            ($name) = $feat->names;
            push @homologs, $feat;
        }
        else {
            $match_type = "proxy for region";
            $name       = "pos " . $feat->start;
        }
        $synteny_score = "0" unless defined $synteny_score;
        my $dsg_name;
        $dsg_name = $dsg->name . ": " if $dsg->name;
        $dsg_name .= " (v" . $dsg->version . " " . $dsg->type->name . ")";
        my $synmap_open_link =
            qq{window.open("}
          . $synmap_link
          . $dsg->id
          . ";fid1=$fid;fid2=$tfid" . qq{");};
        $open_all_synmap{$synmap_open_link} = 1;
        my $oid   = $feat->organism->id;
        my $dsgid = $dsg->id;
        $table_content .= join( "<td>", "<span class='link' onclick=window.open('OrganismView.pl?oid=$oid')>"
              . $feat->organism->name
              . "</span>", "<span class='link' onclick=window.open('OrganismView.pl?dsgid=$dsgid')>"
              . $dsg_name
              . "</span>",
            $match_type,
            qq{<span class='link' onclick='generate_feat_info("$tfid} . "_"
              . $dsg->id
              . qq{")'>}
              . $name
              . "</span>",
            $feat->chromosome,
            $synteny_score,
            qq{<span class="link" onclick='$synmap_open_link'>Dotplot</span>},
        );
        $count++;
    }

    if ($table_content) {
        $html .= $table_header . $table_content . "</tbody></table>" if $table_content;
    } else {
        #XXX Move this error checking up
        return encode_json({
            success => JSON::false,
            message => "The feature selected does not exist in any of the genomes selected."
        });
    }

    my $featlist_link = gen_featlist_link( fids => [ $fid, map { $_->[0] } @$matches ] );

    $html .= '<br><strong>Links</strong><br>'
        . qq{<span class="small">Regenerate this analysis: </span><a href='$tiny_synfind_link' class='small link' target=_new_synfind>$tiny_synfind_link</a><br>}
        . qq{<a class="small link" href='$tiny_gevo_link' target=_blank)">Compare and visualize region in GEvo</a><br>};

    my $open_all_synmap = join( "\n", keys %open_all_synmap );
    $html .= qq{<a onclick='$open_all_synmap' class='small link'>Generate all dotplots</a><br>};
    my $feat_list_link = qq{FeatList.pl?fid=} . join( ",", map { $_->id } @homologs );
    $html .= qq{<a onclick='window.open("$feat_list_link")' class='small link'>Generate Feature List</a><br>};
    my $fasta_link = qq{FastaView.pl?fid=} . join( ",", map { $_->id } @homologs );
    $html .= qq{<a onclick='window.open("$fasta_link")' class='small link'>Generate Fasta Sequences</a><br>};

    my $master_list_link = $synfind_link . ";get_master=1";

    #$html .= "<a onclick=window.open('$master_list_link') class='ui-button ui-corner-all' target=_new_synfind>Generate master gene set table</a>";
    $html .= qq{<span onclick="get_master('$master_list_link')" class='small link' target=_new_synfind>Generate master gene set table</span><br>};
    $master_list_link .= ";limit=1";

    #$html .= "<a onclick=window.open('$master_list_link') class='ui-button ui-corner-all' target=_new_synfind>Generate master gene set table (top one syntenlog per organism)</a>";
    $html .= qq{<span onclick="get_master('$master_list_link')" class='small link' target=_new_synfind>Generate master gene set table (top one syntenlog per organism)</span>};
    #$html .= qq{<span class="small">Pad Sequence in GEvo <input type="text" size=11 id=pad size=11 value=0></span>};
    my $log_file = $cogeweb->logfile;
    $log_file =~ s/$TEMPDIR/$TEMPURL/;

    $html .= qq{<br><br><strong>Syntenic Depth</strong> }
        . qq{<a href="/wiki/index.php/SynFind#Syntenic_Depth" target=_new class="ui-icon ui-icon-help" style='border:1px solid gray;border-radius:50%;'></a>}
        . qq{ <span class="small note">Query genome depth coverage</span><br>}
        . get_master_histograms(
                target_dbs       => \@target_info,
                query_dsgid      => $source_dsgid,
                unique_gene_link => $unique_gene_link );
    $html .= "<br><strong>Downloads</strong><br>"
        . qq{<a href="$log_file" target=_new class="small">Log File</a><br>};

    foreach my $item (@target_info) {
        # mdb added 10/8/13
        my $blastfile_link = $item->{blastfile};
        $blastfile_link =~ s/$config->{COGEDIR}//;

        $html .= '<span class="small">' . $item->{org_name} . ':</span> ' . qq{<a href="$blastfile_link" class="small" target=_new>Raw Blast</a>, };

        # mdb added 10/8/13
        $blastfile_link = $item->{filtered_blastfile};
        $blastfile_link =~ s/$config->{COGEDIR}//;

        $html .= qq{<a href="$blastfile_link" class="small" target=_new>Filtered Blast, </a>};
	my $db = $item->{synteny_score_db};
        $db =~ s/$config->{COGEDIR}//;
	$html .= qq{<a href="$db" class="small" target=_new> Synteny_Score SQLite Database</a><br>};
	
    }

    return encode_json({
        success => JSON::true,
        html => $html
    });
}

sub gen_fasta {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};

    my $feat_type = $opts{feat_type};
    my $write_log = $opts{write_log} || 0;
    my ( $org_name, $title );
    ( $org_name, $title ) = gen_org_name(
        dsgid     => $dsgid,
        feat_type => $feat_type,
        write_log => $write_log
    );
    my $file = $FASTADIR . "/$dsgid-$feat_type.fasta";
    my $res;
    CoGe::Accessory::Web::write_log( "#FASTA#", $cogeweb->logfile )
      if $write_log;

    while ( -e "$file.running" ) {
        print STDERR "detected $file.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $file ) {
        CoGe::Accessory::Web::write_log(
            "fasta file for *" . $org_name . "* ($file) exists",
            $cogeweb->logfile )
          if $write_log;
        $res = 1;
    }
    else {
        system "/usr/bin/touch $file.running"; #track that a blast anlaysis is running for this
        $res =
          generate_fasta( dsgid => $dsgid, file => $file, type => $feat_type )
          unless -r $file;
        system "/bin/rm $file.running" if -r "$file.running"; #remove track file
    }
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile ) if $write_log;
    return $file, $org_name, $title if $res;
    return 0;
}

sub gen_org_name {
    my %opts      = @_;
    my $dsgid     = $opts{dsgid};
    my $feat_type = $opts{feat_type} || 1;
    my $write_log = $opts{write_log} || 0;

    my ($dsg) = $coge->resultset('Genome')->search(
        { genome_id => $dsgid },
        { join => 'organism', prefetch => 'organism' }
    );
    return unless $dsg;

    my $org_name = $dsg->organism->name;
    my $title = $org_name . " (v" . $dsg->version . ", dsgid" . $dsgid . ")" . $feat_type;
    $title =~ s/(`|')//g;

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile )
      if $write_log;
    CoGe::Accessory::Web::write_log( "ORGANISM: " . $title, $cogeweb->logfile )
      if $write_log;

    return ( $org_name, $title );
}

sub generate_fasta {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    my $file  = $opts{file};
    my $type  = $opts{type};

    my ($dsg) =
      $coge->resultset('Genome')->search( { "me.genome_id" => $dsgid } );
#      $coge->resultset('Genome')->search( { "me.genome_id" => $dsgid },
#        { join => 'genomic_sequences', prefetch => 'genomic_sequences' } );
    $file = $FASTADIR . "/$file" unless $file =~ /$FASTADIR/;
    CoGe::Accessory::Web::write_log( "creating fasta file", $cogeweb->logfile );
    open( OUT, ">$file" ) || die "Can't open $file for writing: $!";
    if ( $type eq "CDS" ) {
        my $count = 1;
        foreach my $feat (
            sort {
                     $a->chromosome cmp $b->chromosome
                  || $a->start <=> $b->start
            } $coge->resultset('Feature')->search(
                {
                    feature_type_id => [ 3, 5, 8 ],
                    genome_id       => $dsgid
                },
                {
                    join => [ { dataset => 'dataset_connectors' } ],
                    prefetch => ['feature_names']
                }
            )
          )
        {
            my ($chr) = $feat->chromosome;    #=~/(\d+)/;
            my $name;
            foreach my $n ( $feat->names ) {
                $name = $n;
                last unless $name =~ /\s/;
            }
            $name =~ s/\s+/_/g;
            my $title = join( "||",
                $chr, $feat->start, $feat->stop, $name, $feat->strand,
                $feat->type->name, $feat->id, $count );
            my $seq = $feat->genomic_sequence( dsgid => $dsg );
            next unless $seq;

            #skip sequences that are only 'x' | 'n';
            next unless $seq =~ /[^x|n]/i;
            print OUT ">" . $title . "\n";
            print OUT $seq, "\n";
            $count++;
        }
    }
    else {
        foreach my $chr ( sort $dsg->get_chromosomes ) {
#		    my $title = join ("||",$chr, 1, $ds->last_chromosome_position($chr), "Chr_$chr",1, "genomic", "N/A");
            my $seq = $dsg->get_genomic_sequence( chr => $chr );
            next unless $seq;
            print OUT ">" . $chr . "\n";
            print OUT $seq, "\n";
        }
    }
    close OUT;
    return 1 if -r $file;
    CoGe::Accessory::Web::write_log( "Error with fasta file creation", $cogeweb->logfile );
    return 0;
}

sub get_feature_count {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    my $count = $coge->resultset('Feature')->count(
        {
            feature_type_id => [ 3, 5, 8 ],
            genome_id       => $dsgid
        },
        { join => [ { dataset => 'dataset_connectors' } ], }
    );
    return $count;
}

sub run_blast {
    my %opts    = @_;
    my $fasta1  = $opts{fasta1};
    my $fasta2  = $opts{fasta2};
    my $outfile = $opts{outfile};
    my $prog    = $opts{prog};
    $prog = "last" unless $prog;
    CoGe::Accessory::Web::write_log( "#RUN BLAST: $prog#", $cogeweb->logfile );
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $outfile || -r "$outfile.gz" ) {
        unless ( -s $outfile || -s "$outfile.gz" ) {
            CoGe::Accessory::Web::write_log(
                "WARNING: Blast output file ($outfile) contains no data!",
                $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
            return 0;
        }
        CoGe::Accessory::Web::write_log( "blastfile $outfile already exists",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        return 1;
    }
    my $pre_command;
    if ( $prog eq "lastz" ) {
        $pre_command .= "$LASTZ -i $fasta1 -d $fasta2 -o $outfile";
    }
    else {
        $pre_command .= "$LAST $fasta2 $fasta1 -o $outfile";
    }
    my $x;
    system "/usr/bin/touch $outfile.running"; #track that a blast anlaysis is running for this
    ( $x, $pre_command ) = CoGe::Accessory::Web::check_taint($pre_command);
    CoGe::Accessory::Web::write_log( "running $pre_command", $cogeweb->logfile );
    `$pre_command`;
    system "/bin/rm $outfile.running" if -r "$outfile.running"; #remove track file
    unless ( -s $outfile ) {
        CoGe::Accessory::Web::write_log("WARNING: Problem running $pre_command command.  Blast output file contains no data!", $cogeweb->logfile);
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        return 0;
    }
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    return 1 if -r $outfile;
}

sub make_bed {
    my %opts    = @_;
    my $dsgid   = $opts{dsgid};
    my $outfile = $opts{outfile};
    CoGe::Accessory::Web::write_log( "#BED FILES#", $cogeweb->logfile );
    my $cmd = $DATASETGROUP2BED . " -dsgid $dsgid > $outfile";
    if ( -r $outfile && -s $outfile ) {
        CoGe::Accessory::Web::write_log( "bed file $outfile already exists", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        return $outfile;
    }
    CoGe::Accessory::Web::write_log( "Creating bedfiles: $cmd", $cogeweb->logfile );
    `$cmd`;
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    return $outfile;
}

sub blast2bed {
    my %opts     = @_;
    my $infile   = $opts{infile};
    my $outfile1 = $opts{outfile1};
    my $outfile2 = $opts{outfile2};
    CoGe::Accessory::Web::write_log( "#BLAST 2 BED#", $cogeweb->logfile );
    if (
        ( -r $outfile1 && -s $outfile1 && -r $outfile2 && -s $outfile2 )
        || (   -r "$outfile1.gz"
            && -s "$outfile1.gz"
            && -r "$outfile2.gz"
            && -s "$outfile2.gz" )
      )
    {
        CoGe::Accessory::Web::write_log(
            ".bed files $outfile1 and $outfile2 already exist.",
            $cogeweb->logfile );
        return;
    }
    CoGe::Accessory::Web::gunzip("$infile");
    my $cmd =
        $PYTHON26 . " "
      . $BLAST2BED
      . " -infile $infile -outfile1 $outfile1 -outfile2 $outfile2";
    CoGe::Accessory::Web::write_log( "Creating bed files: $cmd", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    `$cmd`;
}

sub run_convert_blast {
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};
    CoGe::Accessory::Web::write_log( "#CONVERT BLAST#", $cogeweb->logfile );
    my $cmd = $CONVERT_BLAST . " < $infile > $outfile";
    if (   ( -r $outfile && -s $outfile )
        || ( -r "$outfile.gz" || -s "$outfile.gz" ) )
    {
        CoGe::Accessory::Web::write_log(
            "converted blast file with short names exists: $outfile",
            $cogeweb->logfile );
        return $outfile;
    }
    CoGe::Accessory::Web::gunzip("$infile");
    CoGe::Accessory::Web::write_log(
        "convering blast file to short names: $cmd",
        $cogeweb->logfile );
    `$cmd`;
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    return $outfile;
}

sub run_blast2raw {
    my %opts      = @_;
    my $blastfile = $opts{blastfile};
    my $bedfile1  = $opts{bedfile1};
    my $bedfile2  = $opts{bedfile2};
    my $outfile   = $opts{outfile};
    CoGe::Accessory::Web::write_log( "#BLAST 2 RAW#", $cogeweb->logfile );
    if (   ( -r $outfile && -s $outfile )
        || ( -r "$outfile.gz" && -s "$outfile.gz" ) )
    {
        CoGe::Accessory::Web::write_log(
"Filtered blast file found where tandem dups have been removed: $outfile",
            $cogeweb->logfile
        );
        return $outfile;
    }
    CoGe::Accessory::Web::gunzip("$blastfile");
    CoGe::Accessory::Web::gunzip("$bedfile1");
    CoGe::Accessory::Web::gunzip("$bedfile2");
    unless ( -r $blastfile ) {
        warn "can't read $blastfile\n";
        return;
    }
    my $tandem_distance = $opts{tandem_distance};
    $tandem_distance = 10 unless defined $tandem_distance;
    my $cmd =
        $PYTHON26 . " "
      . $BLAST2RAW
      . " $blastfile --qbed $bedfile1 --sbed $bedfile2 --tandem_Nmax $tandem_distance > $outfile";
    CoGe::Accessory::Web::write_log(
        "BLAST2RAW (finding and removing local duplications): running $cmd",
        $cogeweb->logfile );
    `$cmd`;
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    return $outfile;
}

sub run_synteny_score {
    my %opts        = @_;
    my $blastfile   = $opts{blastfile};
    my $bedfile1    = $opts{bedfile1};
    my $bedfile2    = $opts{bedfile2};
    my $window_size = $opts{window_size};
    my $cutoff      = $opts{cutoff};
    if ( $cutoff >= 1 ) {
        $cutoff = sprintf( "%.2f", $cutoff / $window_size );
    }
    my $outfile          = $opts{outfile};
    my $scoring_function = $opts{scoring_function};
    my $dsgid1           = $opts{dsgid1};
    my $dsgid2           = $opts{dsgid2};
    CoGe::Accessory::Web::write_log( "#SYNTENY SCORE#", $cogeweb->logfile );
    while ( -e "$outfile.running" ) {
        print STDERR "detected $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if (   ( -r $outfile && -s $outfile )
        || ( -r "$outfile.gz" && -s "$outfile.gz" ) )
    {
        CoGe::Accessory::Web::write_log(
            "synteny_score database ($outfile) exists",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        return $outfile;
    }
    else {
        system "/usr/bin/touch $outfile.running"
          ;    #track that a blast anlaysis is running for this
    }
    CoGe::Accessory::Web::gunzip("$blastfile");    #turned on debugging
    CoGe::Accessory::Web::gunzip("$bedfile1");
    CoGe::Accessory::Web::gunzip("$bedfile2");
    unless ( -r $blastfile ) {
        warn "can't read $blastfile\n";
        return;
    }
    my $cmd = $SYNTENY_SCORE
      . " $blastfile --qbed $bedfile1 --sbed $bedfile2 --window $window_size --cutoff $cutoff --scoring $scoring_function --qnote $dsgid1 --snote $dsgid2 --sqlite $outfile";

    CoGe::Accessory::Web::write_log( "Synteny Score:  running $cmd",
        $cogeweb->logfile );
    print STDERR $cmd, "\n";
    system("$PYTHON26 $cmd");
    system "/bin/rm $outfile.running"
      if -r "$outfile.running";    #remove track file
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    return $outfile;
}

sub gen_gevo_link {
    my %opts = @_;
    my $fid  = $opts{fid};
    my $dbs  = $opts{dbs};
    my $depth =
      $opts{depth}; #syntenic depth is the maximum number of syntenic regions per organism requested
    my $window_size =
      $opts{window_size}; #this is the the number of genes searched around the feature, half up and half down

    return "no feature id specified" unless $fid;

    #determine distance of $window_size/2 genes up and down from query feature
    my ( $up, $down ) =
      get_neighboring_region( fid => $fid, window_size => $window_size );
    my $link = "$SERVER/GEvo.pl?fid1=$fid;dr1up=$up;dr1down=$down";
    my @matched_fids;
    my $count = 2;
    my %seen_fids;
    foreach my $db (@$dbs) {
        next unless -s $db;
        my $dbh   = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );
        my $query = "SELECT * FROM  synteny where query = $fid";
        my $sth   = $dbh->prepare($query);
        next unless $sth;
        $sth->execute();
        my $depth_count = 0;
        while ( my $data = $sth->fetchrow_arrayref ) {
            last if $depth && $depth_count >= $depth;
            next if $seen_fids{ $data->[1] };
            push @matched_fids, [ $data->[1], $data->[2], $data->[3] ]
              ;    #fid, match_type, synteny_score
            $link .= ";fid$count" . "=" . $data->[1];
            $seen_fids{ $data->[1] } = 1;
            $link .= ";ref$count" . "=0";
            $link .= ";dr$count" . "up=" . $data->[4];
            $link .= ";dr$count" . "down=" . $data->[4];

            if ( $data->[5] =~ /-/ ) {
                $link .= ";rev$count" . "=1";
            }
            $count++;
            $depth_count++;
        }
    }
    $count--;
    $link .= ";num_seqs=$count;autogo=1";
    CoGe::Accessory::Web::write_log( "#GEVO LINK: $link", $cogeweb->logfile );
    return $link, \@matched_fids;
}

sub gen_featlist_link {
    my %opts = @_;
    my $fids = $opts{fids};
    my $link = "$SERVER/FeatList.pl?fid=" . join( ";fid=", @$fids );
    return $link;
}

sub generate_feat_info {
    my %opts = @_;
    my $featid = $opts{featid};
    ( $featid, my $dsgid ) = split( /_/, $featid );
    my ($dsg)  = $coge->resultset('Genome')->find($dsgid);
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless ( $dsg && $feat && ref($feat) =~ /Feature/i ) {
        return "Unable to retrieve Feature object for id: $featid";
    }
    my $html = $feat->annotation_pretty_print_html( gstid => $dsg->type->id );
    return $html;
}

sub has_cds {
    my $dsgid   = shift;
    my $has_cds = 0;

#add check to make sure that if type is "CDS" that CDS annotations exist!  Otherwise return error!
    foreach my $ft (
        $coge->resultset('FeatureType')->search(
            {
                genome_id            => $dsgid,
                'me.feature_type_id' => 3
            },
            {
                join => { features => { dataset => 'dataset_connectors' } },
                rows => 1,
            }
        )
      )
    {
        $has_cds = 1;
    }
    return $has_cds;
}

sub get_master_syn_sets {
    my $window_size      = $FORM->param('ws');
    my $cutoff           = $FORM->param('co');
    my $scoring_function = $FORM->param('sf');
    my $pad              = $FORM->param('pad');
    my $algo             = $FORM->param('algo');
    if ( $scoring_function == 2 ) {
        $scoring_function = "density";
    }
    else {
        $scoring_function = "collinear";
    }

    my $qdsgid = $FORM->param('qdsgid');
    my $limit  = $FORM->param('limit')
      ;    #limit the number of syntelogs returned per organism searched
    my $qdsg = $coge->resultset('Genome')->find($qdsgid);
    my %qdsids = map { $_->id, 1 } $qdsg->datasets;
    my @dsgs;
    foreach my $item ( $FORM->param('dsgid') ) {
        foreach my $dsgid ( split( /,/, $item ) ) {
            push @dsgs, $coge->resultset('Genome')->find($dsgid);
        }
    }

    my $header = "Content-disposition: attachement; filename=";  #test.gff\n\n";
    $header .= join( "_",
        ( map { $_->id } ( $qdsg, @dsgs ) ),
        $window_size, $cutoff, $scoring_function, $algo );
    $header .= ".txt\n\n";
    print $header;

    my %data;
    my %lookup;
    foreach my $dsg (@dsgs) {
        my $org1 = $qdsg->organism->name;
        my $org2 = $dsg->organism->name;

        my $dsgid1 = $qdsg->id;
        my $dsgid2 = $dsg->id;
        ( $org1, $org2, $dsgid1, $dsgid2 ) = ( $org2, $org1, $dsgid2, $dsgid1 )
          if ( $dsgid2 lt $dsgid1 );
        my $basedir  = $DIAGSDIR . "/" . $dsgid1 . "/" . $dsgid2;
        my $basename = $dsgid1 . "_" . $dsgid2 . "." . "CDS-CDS";
        my $db =
            $basedir . "/"
          . $basename . "_"
          . $window_size . "_"
          . $cutoff . "_"
          . $scoring_function
          . ".$algo"
	  . ".v053"  . ".db";

        unless (-r $db) { # mdb added 3/4/14, issue 324
            print STDERR qq{SQLite database not found: $db\n};
            next;
        }

        my $dbh   = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );
        unless ($dbh) {
            print STDERR qq{Problem connecting to $db\n};
            next;
        }

        my $query = "SELECT * FROM synteny";
        my $sth   = $dbh->prepare($query);
        unless ($sth) {
            print STDERR qq{Problem querying $db\n};
            next;
        }

        $sth->execute();
        while ( my @data = $sth->fetchrow_array ) {
            next unless $data[6] == $qdsgid;
            my $id = $data[0];
            unless ( $data{$id} ) {
                my ($feat) = $coge->resultset('Feature')->find($id);
                my $chr    = $feat->chromosome;
                my $start  = $feat->start;
                $lookup{$chr}{$start}{$id} = 1;
            }
            my $sdsgid = $data[7];

            push @{ $data{$id}{$sdsgid} },
              \@data;    #data{query_feature_id}{subject_genome_id}
        }
    }

    print "#",
      join(
        "\t", "COUNTS",
        (
            map { "ORG: " . $_->organism->name." (".$_->id.")" } (
                $qdsg, sort { $a->organism->name cmp $b->organism->name } @dsgs
            )
        ),
        (
            map { "CHR: " . $_->organism->name } (
                $qdsg, sort { $a->organism->name cmp $b->organism->name } @dsgs
            )
        )
      ),
      "\tGEvo link", "\n";

    my $total;
    foreach my $sort_chr ( sort keys %lookup ) {
        foreach my $sort_start ( sort keys %{ $lookup{$sort_chr} } ) {
            foreach my $id ( keys %{ $lookup{$sort_chr}{$sort_start} } ) {
                my $link = $SERVER . "/GEvo.pl?";

#@data contains:
# element 1:  array_ref of hard coded stuff to fake data from query genome
# elements 2 -- N+1: data from synteny database for each of the subject genomes of which there is N
                my @data = (
                    [ [ 0, $id, "S", 100000 ] ],
                    map    { $data{$id}{ $_->id } }
                      sort { $a->organism->name cmp $b->organism->name } @dsgs
                );
                my @names;
                my @chr;
                my $count = 1;
                my $max;
                my @syntelog_count;
              SET:

                foreach my $set (
                    @data) #iterate through each genome -- first is query genome
                {

                    unless ($set) {
                        push @names,          "-";
                        push @chr,            "-";
                        push @syntelog_count, 0;
                        next;
                    }
                    my $name;
                    my $chr;
                    my $limit_count = 0;
                    foreach my $data ( sort { $b->[3] <=> $a->[3] }
                        @$set )    #sort by synteny score
                    {
                        if ($limit) {
                            last if $limit_count >= $limit;
                        }
                        my $fid = $data->[1];
                        unless ($fid) {
                            $name .= "-,";
                            $chr  .= "-,";
                            next;
                        }
                        my ($feat) = $coge->resultset('Feature')->find($fid);
                        $chr .= $feat->chromosome . ",";
                        if ( $count == 1 ) {

#		my ($up, $down) = get_neighboring_region(fid=>$fid, window_size=>$window_size);
#		$link .= ";dr$count"."up=".$up;
#		$link .= ";dr$count"."down=".$down
                        }
                        else {
                            $link .= ";ref$count=0";
                            $link .= ";dr$count" . "up=" . $data->[4];
                            $link .= ";dr$count" . "down=" . $data->[4];
                            $max = $data->[4] unless $max;
                            $max = $data->[4] if $max < $data->[4];
                            if ( $data->[5] =~ /-/ ) {
                                $link .= ";rev$count" . "=1";
                            }
                        }
                        if ( $data->[2] eq "S" ) {
                            my $rs = $coge->resultset('FeatureName');
                            $rs->result_class(
                                'DBIx::Class::ResultClass::HashRefInflator');

#		    my ($name_hash) = sort {$b->{primary_name} <=> $a->{primary_name} || $a->{name} cmp $b->{name}} $rs->search({feature_id=>$fid});
                            my ($name_hash) = sort {
                                     $b->primary_name <=> $a->primary_name
                                  || $a->name cmp $b->name
                              } $coge->resultset('FeatureName')
                              ->search( { feature_id => $fid } );
                            $name .=  $name_hash?$name_hash->name . ",":$fid.",";
                            $link .= ";fid$count=$fid";
                        }
                        else {
                            $name .= "proxy" . ",";
                            my $rs = $coge->resultset('Feature');
                            $rs->result_class(
                                'DBIx::Class::ResultClass::HashRefInflator');
                            my ($feat_hash) = $rs->find($fid);
                            $link .=
                                ";x$count="
                              . $feat_hash->{start}
                              . ";chr$count="
                              . $feat_hash->{chromosome}
                              . ";dsgid$count="
                              . $data->[7];
                        }
                        $limit_count++;
                        $count++;
                    }
                    push @syntelog_count, $limit_count;
                    $name =~ s/,$//;    #trim trailing ','
                    push @names, $name;
                    $chr =~ s/,$//;     #trim trailing ','
                    push @chr, $chr;
                }

#this is a short cut for specifying the dr up/down of query sequence.  Much faster than looking it up in the database
                $link .= ";dr1up=" . $max;
                $link .= ";dr1down=" . $max;
                $count--;
                $link .= ";num_seqs=$count;autogo=1";
                $link .= ";apply_all=$pad" if $pad;
                print join( "\t",
                    join( ",", @syntelog_count ),
                    @names, @chr, $link ),
                  "\n";
                $total++;
            }
        }
    }
}

sub get_unique_genes {
    my $window_size      = $FORM->param('ws');
    my $cutoff           = $FORM->param('co');
    my $scoring_function = $FORM->param('sf');
    my $pad              = $FORM->param('pad');
    my $algo             = $FORM->param('algo');
    my $fid              = $FORM->param('fid');
    my $qdsgid           = $FORM->param('qdsgid');
    my $sdsgid           = $FORM->param('sdsgid');
    my $synfind_link     = $SERVER . "SynFind.pl?fid=$fid;qdsgid=$qdsgid;ws=$window_size;co=$cutoff;sf=$scoring_function;algo=$algo;dsgid=$sdsgid";

    if ( $scoring_function == 2 ) {
        $scoring_function = "density";
    }
    else {
        $scoring_function = "collinear";
    }

    my ($qdsg) = $coge->resultset('Genome')->find($qdsgid);
    my ($sdsg) = $coge->resultset('Genome')->find($sdsgid);
    print "Error" unless $qdsg && $sdsg;
    my $org1 = $qdsg->organism->name;
    my $org2 = $sdsg->organism->name;

    my $dsgid1 = $qdsg->id;
    my $dsgid2 = $sdsg->id;
    ( $org1, $org2, $dsgid1, $dsgid2 ) = ( $org2, $org1, $dsgid2, $dsgid1 )
        if ( $dsgid2 lt $dsgid1 );
    my $basedir  = $DIAGSDIR . "/" . $dsgid1 . "/" . $dsgid2;
    my $basename = $dsgid1 . "_" . $dsgid2 . "." . "CDS-CDS";
    my $db = $basedir . "/" . $basename . "_" . $window_size . "_" . $cutoff . "_" . $scoring_function . ".$algo" . ".v053" . ".db";
    my $dbh   = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );
    my $query = "SELECT * FROM synteny";
    my $sth   = $dbh->prepare($query);
    unless ($sth) {
        print STDERR qq{Problem connecting to $db\n};
        next;
    }
    $sth->execute();

    my %qdata;
    my %sdata;
    while ( my @data = $sth->fetchrow_array ) {
        next unless $data[6] == $qdsgid;
        my $idq = $data[0];
        my $ids = $data[1];
        $qdata{$idq}{$ids} = 1;
        $sdata{$ids}{$idq} = 1;
    }
    $dbh->disconnect();

    my @qunique;
#    my @sunique;

    my $rs = $coge->resultset('Feature')->search(
        {   feature_type_id => 3,
            genome_id       => $qdsgid
        },
        { join => [ { dataset => 'dataset_connectors' } ] }
    );
    foreach my $item ( $rs->get_column('feature_id')->all ) {
        push @qunique, $item unless $qdata{$item};
    }

#    $rs = $coge->resultset('Feature')->search(
#        {   feature_type_id => 3,
#            genome_id       => $sdsgid
#        },
#        { join => [ { dataset => 'dataset_connectors' } ], }
#    );
#    foreach my $item ( $rs->get_column('feature_id')->all ) {
#        push @sunique, $item unless $sdata{$item};
#    }

    # Create a list of these items and add them to the user's lists
    # mdb added 3/6/14 - limit max number of items
    if (@qunique <= $MAX_FEATURES_IN_LIST) {
# mdb removed 12/14/16 COGE-800
#        my $list_type = $coge->resultset('ListType')->find_or_create(
#            {   name        => 'SynFind Unique Genes',
#                description => 'Auto-generated by SynFind'
#            }
#        );

        my $list = $coge->resultset('List')->create({
            name         => "Unique genes in " . $qdsg->organism->name,
            description  => "Compared to " . $sdsg->organism->name,
            #list_type_id => $list_type->id, # mdb removed 12/14/16 COGE-800
            creator_id   => $USER->id,
            restricted   => 1
        });
        return unless $list;

        my $conn = $coge->resultset('UserConnector')->create({
            parent_id   => $USER->id,
            parent_type => 5,           # FIXME hardcoded to "user"
            child_id    => $list->id,
            child_type  => 1,           # FIXME hardcoded to "list"
            role_id     => 2,           # FIXME hardcoded to "owner"
        });
        return unless $conn;

        my $annotation_type = $coge->resultset('AnnotationType')->find_or_create( { name => "CoGe Note" } );
        $list->add_to_list_annotations({
            annotation_type_id => $annotation_type->id,
            annotation         => "Auto created by SynFind",
            locked => 1
        });
        $list->add_to_list_annotations({
            annotation_type_id => $annotation_type->id,
            annotation         => "Regenerate Analysis in SynFind",
            link => CoGe::Accessory::Web::get_tiny_link( url => $synfind_link ),
            locked => 1
        });
        $list->add_to_list_annotations({
            annotation_type_id => $annotation_type->id,
            annotation         => "Query Genome: " . $qdsg->info,
            link => $SERVER . "OrganismView.pl?dsgid=" . $qdsg->id,
            locked => 1
        });
        $list->add_to_list_annotations({
            annotation_type_id => $annotation_type->id,
            annotation         => "Target Genome: " . $sdsg->info,
            link => $SERVER . "OrganismView.pl?dsgid=" . $sdsg->id,
            locked => 1
        });

        foreach my $fid (@qunique) {
            $list->add_to_list_connectors({
                parent_id => $list->id, child_id => $fid, child_type => 4 # FIXME: hardcoded type
            });
        }

        my $listview = $SERVER . "NotebookView.pl?lid=" . $list->id;
        print $FORM->header, qq{<meta HTTP-EQUIV="REFRESH" content="0; url=$listview">};
    }
}

sub get_master_histograms {
    my %opts             = @_;
    my $target_dbs       = $opts{target_dbs};
    my $query_dsgid      = $opts{query_dsgid};
    my $unique_gene_link = $opts{unique_gene_link};

    my $total_feature_count = get_feature_count( dsgid => $query_dsgid );

    my $html;
    my @all;
    foreach my $item (@$target_dbs) {
        next unless -s $item->{synteny_score_db};
        my $db    = $item->{synteny_score_db};
        my $dbh   = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );
        my $query = "SELECT * FROM synteny";
        my $sth   = $dbh->prepare($query);
        unless ($sth) {
            print STDERR qq{Problem connecting to $db\n};
            next;
        }
        $sth->execute();
        my %data;
        while ( my @data = $sth->fetchrow_array ) { # data base contains reciprocal entries for each item (a=>b; b=>a)
            next unless $data[6] == $query_dsgid;
            my $id     = $data[0];    #query feature ID
            my $sdsgid = $data[7];    #target genome ID
            $data{$id}++;             #data{query_feature_id}{subject_genome_id}
        }

        my $genome = $coge->resultset('Genome')->find( $item->{dsgid2} );
        next unless $genome;

        $html .= qq{<span class="small link" onclick='window.open("OrganismView.pl?dsgid=} . $item->{dsgid2} . "\")'>" . $genome->info . "</span>";
        $html .= '<div>';
        if ( $query_dsgid eq $item->{dsgid2} ) {
            $html .= "<span class='small'>(self-self comparison: self-self syntenic regions ignored)</span><br>";
        }
        my $link = $unique_gene_link . ";ug=1;sdsgid=" . $item->{dsgid2} . ";qdsgid=$query_dsgid";
        $html .= qq{<span class="small link" onclick="window.open('$link')">Query genome unique genes</span><br>};
        $link = $unique_gene_link . ";ug=1;qdsgid=" . $item->{dsgid2} . ";sdsgid=$query_dsgid";
        $html .= qq{<span class="small link" onclick="window.open('$link')">Target genome unique genes</span><br>};
        my $depths;
        $html .= depth_table(
            depth_data => \%data,
            total      => $total_feature_count,
            depths     => \$depths
        );
        push @all, [ $depths, $item->{dsgid2} ];
    }
    $html .= master_depth_table( data => \@all );
    $html .= '</div>';
    return $html;
}

sub master_depth_table {
    my %opts = @_;
    my $data = $opts{data};
    my $html = "<table class='small' style='border-top:1px solid lightgray;border-bottom:1px solid lightgray;'>";
    my $max  = 0;
    foreach
      my $item ( sort { scalar @{ $b->[0] } <=> scalar @{ $a->[0] } } @$data )
    {
        unless ($max) {
            $max = scalar @{ $item->[0] };
            $html .= "<tr>" . "<th style='text-align:left'>Organism";
            foreach ( my $i = 0 ; $i < $max ; $i++ ) {
                $html .= "<th>Depth: $i";
            }
            $html .= "</tr>\n";
        }
        $html .= "<tr>";
        my $genome = $coge->resultset('Genome')->find( $item->[1] );
        $html .= "<td>" . $genome->info;
        foreach ( my $i = 0 ; $i < $max ; $i++ ) {
            my $val = $item->[0][$i] ? $item->[0][$i] : 0;
            $html .= "<td>$val";
        }
        $html .= "</tr>\n";
    }
    $html .= "</table>";
    return $html;
}

sub depth_table {
    my %opts   = @_;
    my $data   = $opts{depth_data};
    my $total  = $opts{total};        #total number of features;
    my $depths = $opts{depths};
    my %depths;
    map { $depths{$_}++ } values %$data;
    my $total_w_depth = 0;
    map { $total_w_depth += $_ } values %depths;
    $depths{0} = ( $total - $total_w_depth );
    my $html = "<table class='small' style='border-top:1px solid lightgray;border-bottom:1px solid lightgray;'>";
    my @depths;

    foreach my $depth ( sort { $a <=> $b } keys %depths ) {
        $html .= "<tr>"
         . "<td>Depth $depth"
         . "<td>" . $depths{$depth}
         . "<td>of"
         . "<td>$total"
         . "<td>" . sprintf( "%.2f", $depths{$depth} / $total * 100 ) . "%"
         . "</tr>";
        push @depths, $depths{$depth};
    }
    $html .= "</table><br>";
    $$depths = \@depths;
    return $html;
}

sub save_settings {
    my %opts   = @_;

    my $prefs  = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );

    foreach my $key ( keys %opts ) {
        if ($opts{$key}) {
            $prefs->{$key} = $opts{$key};
        }
        else {
            delete $prefs->{$key};
        }
    }

    my $item = CoGe::Accessory::Web::save_settings(
        opts => $prefs,
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
}

sub get_neighboring_region {
    my %opts        = @_;
    my $fid         = $opts{fid};
    my $window_size = $opts{window_size};
    return "no feature id specified" unless $fid;

    #determine distance of $window_size/2 genes up and down from query feature
    my $feat = $coge->resultset('Feature')->find($fid);

    my $count = 0;
    my %seen;
    my $last_item;
  item: foreach my $item (
        $coge->resultset('Feature')->search(
            {
                dataset_id      => $feat->dataset_id,
                feature_type_id => 3,
                chromosome      => $feat->chromosome,
                start           => { '>', $feat->stop }
            },
            {

                #							  join =>'feature_names',
                #							  prefetch =>'feature_names',
                order_by => 'start ASC',
                limit    => $window_size
                , #search for more than we need as some are alternative spliced transcripts.
            }
        )
      )
    {
        foreach my $name ( $item->names ) {
            next item if $seen{$name};
            $seen{$name} = 1;
        }
        $last_item = $item;
        $count++;
        last if $count >= $window_size / 2;
    }
    my $down = $last_item->stop - $feat->stop if $last_item && $feat;
    $down = 10000 unless defined $down;
    $down  = 10000 if $down < 0;
    $count = 0;
    %seen  = ();
  item: foreach my $item (
        $coge->resultset('Feature')->search(
            {
                dataset_id      => $feat->dataset_id,
                feature_type_id => 3,
                chromosome      => $feat->chromosome,
                start           => { '<', $feat->start }
            },
            {

                #							  join =>'feature_names',
                prefetch => 'feature_names',
                order_by => 'start DESC',
                limit    => $window_size
                , #search for more than we need as some are alternative spliced transcripts.
            }
        )
      )
    {
        foreach my $name ( $item->names ) {
            next item if $seen{$name};
            $seen{$name} = 1;
        }
        $last_item = $item;
        $count++;
        last if $count >= $window_size / 2;
    }
    my $up = $feat->start - $last_item->start if $feat && $last_item;
    $up = 10000 unless defined $up;
    $up = 10000 if $up < 0;
    return ( $up, $down );
}
