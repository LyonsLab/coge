#! /usr/bin/perl -w
use strict;
use CoGeX;
use DBIxProfiler;

use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use CGI;
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use POSIX;
use File::Path;
use Benchmark qw(:all);
use Parallel::ForkManager;
use DBI;
use Sort::Versions;

no warnings 'redefine';

#example URL: http://toxic.berkeley.edu/CoGe/SynFind.pl?fid=34519245;qdsgid=3;dsgid=4241,6872,7084,7094,7111

use vars
  qw($P $PAGE_TITLE $PAGE_NAME $DIR $URL $TEMPDIR $TEMPURL $DATADIR $FASTADIR $BLASTDBDIR
  $DIAGSDIR $BEDDIR $LASTZ $LAST $CONVERT_BLAST $BLAST2BED $BLAST2RAW
  $SYNTENY_SCORE $DATASETGROUP2BED $PYTHON26 $FORM $USER $coge
  $cogeweb $RESULTSLIMIT $MAX_PROC $SERVER);

$PAGE_TITLE = 'SynFind_old';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$ENV{PATH}     = $P->{COGEDIR};
$TEMPDIR       = $P->{TEMPDIR} . "SynFind";
$TEMPURL       = $P->{TEMPURL} . "SynFind";
$SERVER        = $P->{SERVER};
$ENV{BLASTDB}  = $P->{BLASTDB};
$ENV{BLASTMAT} = $P->{BLASTMATRIX};
$DIR           = $P->{COGEDIR};
$URL           = $P->{URL};
$DATADIR       = $P->{DATADIR};
$DIAGSDIR      = $P->{DIAGSDIR};
$FASTADIR      = $P->{FASTADIR};
$BLASTDBDIR    = $P->{BLASTDB};
$BEDDIR        = $P->{BEDDIR};
mkpath( $BEDDIR, 0, 0777 ) unless -d $BEDDIR;

$MAX_PROC = $P->{MAX_PROC};
$LASTZ =
    $P->{PYTHON} . " "
  . $P->{MULTI_LASTZ}
  . " -A $MAX_PROC --path="
  . $P->{LASTZ};
$LAST =
    $P->{PYTHON} . " "
  . $P->{MULTI_LAST}
  . " -a $MAX_PROC --path="
  . $P->{LAST_PATH}
  . " --dbpath="
  . $P->{LASTDB};

$CONVERT_BLAST = $P->{CONVERT_BLAST};
$BLAST2BED     = $P->{BlAST2BED};
$BLAST2RAW     = $P->{BLAST2RAW};
$SYNTENY_SCORE = $P->{SYNTENY_SCORE};
$PYTHON26      = $P->{PYTHON};
$DATASETGROUP2BED =
  $P->{DATASETGROUP2BED} . " -cf " . $ENV{COGE_HOME} . 'coge.conf';

my %ajax = CoGe::Accessory::Web::ajax_func();

if ( $FORM->param('get_master') ) {
    get_master_syn_sets();
}
if ( $FORM->param('ug') )    #generating unique gene lists
{
    get_unique_genes();
}
my $pj = new CGI::Ajax(
    %ajax,
    gen_html                => \&gen_html,
    get_orgs                => \&get_orgs,
    gen_dsg_menu            => \&gen_dsg_menu,
    get_genome_info         => \&get_genome_info,
    get_dsg_for_search_menu => \&get_dsg_for_search_menu,
    generate_basefile       => \&generate_basefile,
    save_orglist_synfind    => \&save_orglist_synfind,
    get_types               => \&get_types,
    cogefeatsearch          => \&cogefeatsearch,
    get_anno                => \&get_anno,
    get_orgs_feat           => \&get_orgs_feat,
    source_search           => \&get_data_source_info_for_accn,
    generate_feat_info      => \&generate_feat_info,
    go                      => \&go_synfind,
);
$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );

#print $FORM->header, gen_html();

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'Syntenic Compiler' );
    $template->param( PAGE_TITLE => 'SynFind' );
    $template->param( HELP       => '/wiki/index.php?title=SynFind' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER     => $name );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    my $prebox =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SynFind.tmpl' );
    $prebox->param( RESULTS_DIV => 1 );
    $template->param( PREBOX => $prebox->output );
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SynFind.tmpl' );
    my $form = $FORM;

    #comparison algorithm
    my $algo = $FORM->param('algo') if $FORM->param('algo');
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

    $template->param( JAVASCRIPT => 1 );
    $template->param( FRONT_PAGE => 1 );

    #set up org search
    $template->param( ORG_LIST => get_orgs() );

    #set up feat search
    my $accn;
    $accn = $form->param('accn');
    $template->param( ACCN => $accn );

    #    $template->param(FEAT_TYPE=> get_feature_types());
    $template->param( ORG_LIST_FEAT => get_orgs_feat( type => "none" ) );
    my $doc_ready;
    $doc_ready .= qq{search_chain(1);\n} if $accn;

#This was changed because of an error it was causing in some browsers on Windows where the "()" were causing JS errors.  This bug was reported by johannes.hofberger and this fixed his problem (Eric L.)
    my $prefs
      ; # = CoGe::Accessory::Web::load_settings(user=>$USER, page=>$PAGE_NAME, coge=>$coge);
    ####
    if ( $form->param('dsgid') ) {
        foreach my $item ( $form->param('dsgid') ) {
            foreach my $dsgid ( split( /,/, $item ) ) {
                my $id = get_dsg_for_search_menu( dsgid => $dsgid );
                $doc_ready .= qq{
  add_to_list('$id');
};
            }
        }
    }
    elsif ( $prefs->{dsgids} ) {
        foreach my $dsgid ( split( /,/, $prefs->{dsgids} ) ) {
            my $id = get_dsg_for_search_menu( dsgid => $dsgid );
            $doc_ready .= qq{
  add_to_list('$id');
};
        }
    }
    if ( my $fid = $form->param('fid') ) {
        my ( $fid, $seq_type_id ) = split( /_/, $fid );
        $seq_type_id = 1 unless $seq_type_id;
        my $qdsgid = $form->param('qdsgid');
        unless ($qdsgid) {
            my $feat = $coge->resultset('Feature')->find($fid);
            foreach my $dsg ( $feat->genomes ) {
                if ( $dsg->type->id eq
                    $seq_type_id )    #find the unmasked version of the data
                {
                    $qdsgid = $dsg->id;
                    last;
                }
            }
            $qdsgid = $feat->genomes->[0]->id unless $qdsgid;
        }
        $doc_ready .= qq{
 get_anno(['args__fid','args__$fid','args__dsgid','args__$qdsgid'],[show_anno]);
};
        $template->param(
            FEAT_DSGID => qq{<input type=hidden id=feat_dsgid value=$qdsgid>} )
          if $qdsgid;
    }
    $template->param( DOCUMENT_READY => $doc_ready ) if $doc_ready;
    $template->param( SAVE_ORG_LIST => 1 ) unless $USER->user_name eq "public";
    $template->param( TEMPDIR => $TEMPDIR );
    return $template->output;
}

sub generate_basefile {
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    return $cogeweb->basefilename;
}

sub get_orgs {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    $name = "" if $name && $name =~ /^Search$/i;
    $desc = "" if $desc && $desc =~ /^Search$/i;
    my $html;
    my @db;
    if ($name) {
        @db =
          $coge->resultset('Organism')
          ->search( { name => { like => "%" . $name . "%" } } );
    }
    elsif ($desc) {
        @db =
          $coge->resultset("Organism")
          ->search( { description => { like => "%" . $desc . "%" } } );
    }
    else {
        $html .=
            qq{<FONT CLASS ="small" id="org_count">Organism count: }
          . $coge->resultset('Organism')->count()
          . qq{</FONT>\n<BR>\n};
        $html .=
qq{<SELECT id="org_id" SIZE="8" MULTIPLE"><option id=null_org>Please search</option></SELECT><input type=hidden id=gstid>\n};
        return $html;
    }
    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
        push @opts,
            "<OPTION value=\""
          . $item->id
          . "\" id=\"o"
          . $item->id . "\">"
          . $item->name
          . "</OPTION>";
    }

    $html .=
        qq{<FONT CLASS ="small" id="org_count">Organism count: }
      . scalar @opts
      . qq{</FONT>\n<BR>\n};
    unless (@opts) {
        $html .= qq{<input type = hidden name="org_id" id="org_id"><br>};
        $html .= "No results";
        return $html;
    }
    $html .=
qq{<SELECT id="org_id" SIZE="8" MULTIPLE onclick="show_add()" onchange="gen_dsg_menu(['args__oid','org_id'],['org_seq_types']);" ondblclick="add_selected_orgs();">\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
}

sub gen_dsg_menu {
    my $t1    = new Benchmark;
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $dsgid = $opts{dsgid};
    my @dsg_menu;
    foreach my $dsg (
        sort {
                 versioncmp( $b->version, $a->version )
              || $a->type->id <=> $b->type->id
              || $a->name cmp $b->name
              || $b->id cmp $a->id
        } $coge->resultset('Genome')->search(
            { organism_id => $oid },
            { prefetch    => ['genomic_sequence_type'] }
        )
      )
    {
        next if $dsg->deleted;
        next if $dsg->restricted && !$USER->has_access_to_genome($dsg);
        $dsgid = $dsg->id unless $dsgid;
        my $name    = $dsg->info;
        my $has_cds = has_cds( $dsg->id );
        $name .= " NO CDS ANNOTATIONS.  CAN'T BE USED." unless $has_cds;
        push @dsg_menu, [ $dsg->id, $name ];
    }
    my $size = scalar @dsg_menu;
    $size = 5 if $size > 5;

    my $dsg_menu = qq{
   <select id=dsgid multiple size=$size onclick="show_add();" ondblclick="get_dsg_for_search_menu(['args__dsgid','dsgid'],[add_to_list]);">>
};
    foreach (@dsg_menu) {
        my ( $numt, $name ) = @$_;
        my $selected = " selected" if $dsgid && $numt == $dsgid;
        $selected = " " unless $selected;
        $dsg_menu .= qq{
   <OPTION VALUE=$numt $selected>$name</option>
};
    }
    $dsg_menu .= "</select>";
    my $t2 = new Benchmark;
    my $time = timestr( timediff( $t2, $t1 ) );
    return ($dsg_menu);

}

sub get_genome_info {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return " " unless $dsgid;
    my $dsg = $coge->resultset("Genome")->find($dsgid);
    return "Unable to create genome object for id: $dsgid" unless $dsg;
    my $html;    # = qq{<div style="overflow:auto; max-height:78px">};
    my $total_length;
    my @gs = sort { $a->chromosome cmp $b->chromosome } $dsg->genomic_sequences;
    my $chr_num = scalar @gs;

    #    $html .=qq{<table class=small ><tr class=small valign=top>};
    my $count = 0;
    foreach my $gs (@gs) {
        my $chr    = $gs->chromosome;
        my $length = $gs->sequence_length;
        $total_length += $length;
        $length = commify($length);
        $html .= qq{$chr:  $length bp<br>};
        $count++;
    }
    my ($ds) = $dsg->datasets;
    my $link = $ds->data_source->link;
    $link = "http://" . $link unless $link =~ /^http/;
    $html =
        qq{<span class='link' onclick="window.open('OrganismView.pl?dsgid=}
      . $dsgid
      . qq{')">View in OrganismView</span><br>}
      . "Source:  <a class = 'link' href="
      . $link
      . " target=_new>"
      . $ds->data_source->name
      . "</a><br>"
      . qq{Chromosome count: $chr_num<br>}
      . qq{<div style="float:left;">Total length: }
      . commify($total_length) . " bp" . "<br>"
      . qq{Sequence Type: }
      . $dsg->genomic_sequence_type->name
      . qq{<input type=hidden id=gstid value=}
      . $dsg->genomic_sequence_type->id . qq{>}
      . qq{<br>}
      . qq{-----------<br>Chr:   (length)<br>}
      . $html;
    return $html;
}

sub get_dsg_for_search_menu {
    my %opts   = @_;
    my $dsgids = $opts{dsgid};
    my $orgids = $opts{orgid};
    my %dsgs;
    if ($orgids) {
        my @orgids = split( /,/, $orgids );
        foreach my $dsg (
            $coge->resultset('Genome')->search( { organism_id => [@orgids] } ) )
        {
            next if $dsg->deleted;
            $dsgs{ $dsg->id } = $dsg;
        }
    }
    if ($dsgids) {
        %dsgs = () if ( $dsgs{$dsgids} );
        foreach my $dsgid ( split( /,/, $dsgids ) ) {
            my $dsg = $coge->resultset('Genome')->find($dsgid);
            next unless $dsg;
            $dsgs{ $dsg->id } = $dsg;
        }
    }
    my $opts;
    foreach my $dsg ( values %dsgs ) {
        my ($ds) = $dsg->datasets;
        $opts .= ":::" if $opts;
        next unless has_cds( $dsg->id );    #skip if it has no CDS annotations
        my $item =
            $dsg->id . "::"
          . $dsg->organism->name . " ("
          . $ds->data_source->name . " "
          . $dsg->type->name . " v"
          . $dsg->version . ")";
        $opts .= $item;
    }
    return $opts;
}

sub get_types {
    my %opts = @_;

    #    my ($dsid, $gstid) = split(/_/, $opts{dsid});
    my $dsgid = $opts{dsgid};
    my $accn  = $opts{accn};
    my $ftid  = $opts{ftid};
    my $html;
    my $blank = qq{<input type="hidden" id="Type_name">-------};
    my %seen;
    my $search = {
        'feature_names.name'           => $accn,
        'dataset_connectors.genome_id' => $dsgid,
    };
    $search->{feature_type_id} = 3;

    my @opts =
      sort map { "<OPTION>$_</OPTION>" }
      grep     { !$seen{$_}++ }
      map      { $_->type->name } $coge->resultset('Feature')->search(
        $search,
        ,
        { join => [ 'feature_names', { 'dataset' => 'dataset_connectors' } ] }
      );

    if (@opts) {
        $html .=
          "<font class=small>Type count: " . scalar @opts . "</font>\n<BR>\n";

        $html .=
qq{<SELECT id="Type_name" SIZE="10" MULTIPLE onChange="get_anno(['args__accn','accn_select','args__type','Type_name', 'args__dsgid','dsgid'],[show_anno])" >\n};
        $html .= join( "\n", @opts );
        $html .= "\n</SELECT>\n";
        $html =~ s/OPTION/OPTION SELECTED/;
    }
    else { return $blank; }

    #    return $blank unless $html =~ /OPTION/;
    return ( $html, $dsgid );
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
    my $org_name = $opts{org_name};
    $org_name = undef if $org_name =~ /^search$/i;
    my $org_desc = $opts{org_desc};
    $org_desc = undef if $org_desc =~ /^search$/i;
    my $blank =
qq{<input type="hidden" id="accn_select"><input type="hidden" id="feat_dsgid">};
    my $weak_query = "Query needs to be better defined.";

    if ( !$accn && !$anno && !$fid ) {
        return $weak_query . $blank unless $org_id && $type;
    }
    my $html;
    my %seen;
    my $search = {};
    $search->{feature_type_id} = 3;    #$type if $type;
    $search->{'organism.name'} = { like => "%" . $org_name . "%" } if $org_name;
    $search->{'organism.description'} = { like => "%" . $org_desc . "%" }
      if $org_desc;
    my $join = {
        join => [
            {
                'feature' => {
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
          $coge->resultset('FeatureName')->search( $search, $join )
          ->search_literal( 'MATCH(me.name) AGAINST (?)', $accn );
    }
    if ($anno) {
        push @names,
          $coge->resultset('FeatureName')->search( $search, $join )
          ->search_literal( 'MATCH(annotation) AGAINST (?)', $anno );
    }
    my @opts;
    foreach my $name (@names) {
        my $item = $name->name;
        next if $seen{ uc($item) };
        $seen{ uc($item) }++;
        push @opts, "<OPTION>$item</OPTION>" if $name->feature->dataset;
    }
    if ( @opts > 10000 ) {
        return $blank
          . "Search results over 10000, please refine your search.\n";
    }
    $html .=
      "<font class=small>Name count: " . scalar @opts . "</font>\n<BR>\n";
    $html .=
qq{<SELECT id="accn_select" SIZE="10" MULTIPLE onChange="source_search_chain(); " >\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
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

    my @feats;
    my $genome = $coge->resultset('Genome')->find($dsgid);
    next if $genome->restricted && !$USER->has_access_to_genome($genome);

    if ($accn) {
        foreach my $feat (
            $coge->resultset('Feature')->search(
                {
                    'feature_names.name'           => $accn,
                    "dataset_connectors.genome_id" => $dsgid
                },
                {
                    join => [
                        'feature_names', { 'dataset' => 'dataset_connectors' }
                    ]
                }
            )
          )
        {
            push @feats, $feat if ( $feat->type->name eq $type );
        }
    }
    else {
        push @feats, $coge->resultset('Feature')->find($fid);
    }
    return unless @feats;
    my $anno;
    $anno .=
        "<font class=small>Annotation count: "
      . scalar @feats
      . "</font>\n<BR><hr>\n"
      if scalar @feats;
    my $i   = 0;
    my $dsg = $coge->resultset('Genome')->find($dsgid);
    unless ($dsg) {
        foreach my $tdsg ( $feats[0]->genomes ) {
            $dsg = $tdsg
              if $tdsg->type->id eq
                  "1";    #default to unmasked sequence if possible
        }
        ($dsg) = $feats[0]->genomes unless $dsg;
    }
    my $gstid = $dsg->type->id;
    foreach my $feat (@feats) {
        my ($genome) = $feat->genomes;
        return ( "restricted", 0 )
          if $genome->restricted && !$USER->has_access_to_genome($genome);

#	next if $feat->dataset->restricted && !$USER->has_access_to_dataset($feat->dataset);
        $i++;

#	my $featid = $feat->id;
#	my $chr = $feat->chr;
#	my $rc = 0;
#	my $pro = 0;
#	my $ds = $feat->dataset->id;
#	my $x = $feat->start;
#	my $z = 4;
#	$anno .= qq{<span class="ui-button ui-corner-all" onClick="window.open('FastaView.pl?featid=$featid&gstid=$gstid');">Get Sequence</span>};
#	$anno .= qq{<span class="ui-button ui-corner-all" onClick="window.open('CoGeBlast.pl?featid=$featid;gstid=$gstid');">CoGeBlast</span>};
#	$anno .= qq{<span class="ui-button ui-corner-all" onClick="window.open('GenomeView.pl?chr=$chr&ds=$ds&x=$x&z=$z;gstid=$gstid');">Genome Browser</span>};
        foreach my $item ( split( /;/, $feat->organism->description ) ) {
            $anno .=
qq{<span class=link onclick="\$('#org_desc').val('$item').focus()">$item</span>;};
        }
        $anno .= join "\n<BR><HR><BR>\n",
          $feat->annotation_pretty_print_html( gstid => $gstid );
    }
    $anno = "<font class=\"annotation\">No annotations for this entry</font>"
      unless $anno;
    return ( $anno, $feats[0]->id );
}

sub get_orgs_feat {
    my %opts    = @_;
    my $type    = $opts{type};
    my $search  = $opts{search};
    my $id_only = $opts{id_only};
    my @db;
    my $count;
    if ( $type && $type eq "name" ) {
        @db =
          $coge->resultset("Organism")
          ->search( { name => { like => "%" . $search . "%" } } );
        $count = scalar @db;
    }
    elsif ( $type && $type eq "desc" ) {
        @db =
          $coge->resultset("Organism")
          ->search( { description => { like => "%" . $search . "%" } } );
        $count = scalar @db;
    }
    else {
        $count = $coge->resultset("Organism")->count();

        #	@db = $coge->resultset("Organism")->all;
    }
    return map { $_->id } @db if $id_only;

#my @db = $name ? $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}})
#  : $coge->resultset('Organism')->all();
    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
        push @opts,
            "<OPTION value=\""
          . $item->id
          . "\" id=\"o"
          . $item->id . "\">"
          . $item->name
          . "</OPTION>";
    }
    my $html;
    $html .=
        qq{<FONT CLASS ="small" id="org_count">Organism count: }
      . $count
      . qq{</FONT>\n<BR>\n};
    if ( $search && !@opts ) {
        $html .=
          qq{<input type = hidden name="org_id_feat" id="org_id_feat"><br>};
        $html .= "No results";
        return $html;
    }
    unshift( @opts,
        "<OPTION value=\"all\" id=\"all\">All Listed Organisms</OPTION>" );
    my $size = scalar @opts;
    $size = 8 if $size > 8;
    $html .= qq{<SELECT id="org_id_feat" SIZE="$size" MULTIPLE >\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
}

sub get_data_source_info_for_accn {
    my %opts     = @_;
    my $accn     = $opts{accn};
    my $org_id   = $opts{org_id};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    $org_id = "all" unless $org_id;
    my %org_ids;
    if ( $org_id eq "all" ) {
        my ( $type, $search );
        ( $type, $search ) = ( "name", $org_name )
          if $org_name && $org_name ne "Search";
        ( $type, $search ) = ( "desc", $org_desc )
          if $org_desc && $org_desc ne "Search";
        %org_ids =
          map { $_ => 1 }
          get_orgs( id_only => 1, type => $type, search => $search )
          if $type && $search;
    }
    else {
        $org_ids{$org_id} = 1;
    }
    my $blank = qq{<input type="hidden" id="dsid">------------};
    return $blank unless $accn;
    my @feats = $coge->resultset('Feature')->search(
        { 'feature_names.name' => $accn },
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
#	next if $genome->restricted && !$USER->has_access_to_genome($genome);
#	next if $feat->dataset->restricted && !$USER->has_access_to_dataset($feat->dataset);
        foreach my $dsg ( $feat->dataset->genomes ) {
            next if $dsg->deleted;
            next if $dsg->restricted && !$USER->has_access_to_genome($dsg);
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
    $html .= qq{
<SELECT name = "feat_dsgid" id="feat_dsgid" MULTIPLE SIZE="10" onChange="get_types_chain();">
};
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
        $html .= qq{  <option value="$id" >$title\n};
        $html =~ s/option/option selected/ unless $count;
        $count++;
    }
    $html .= qq{</SELECT>\n};
    return ("<font class=small>Dataset count: "
          . $count
          . "</font>\n<BR>\n"
          . $html );
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
    my $basename         = $opts{basename};
    my $window_size      = $opts{window_size};
    my $cutoff           = $opts{cutoff};
    my $scoring_function = $opts{scoring_function};
    my $depth            = $opts{depth};           #max syntenic depth to report

    $window_size = 40  unless defined $window_size;
    $cutoff      = 0.1 unless defined $cutoff;
    $cutoff      = "0" . $cutoff
      if $cutoff =~ /^\./;    #add the prepending 0 for filename consistence
    $scoring_function = 1 unless defined $scoring_function;

#need to set this link before the scoring function changes to a different name type
    my $synfind_link = $SERVER
      . "SynFind.pl?fid=$fid;ws=$window_size;co=$cutoff;sf=$scoring_function;algo=$algo";
    my $unique_gene_link = $synfind_link;
    $synfind_link .= ";dsgid=$dsgids;qdsgid=$source_dsgid";
    $synfind_link .= ";sd=$depth" if $depth;

    my @ids = split( /,/, $dsgids );

    my $list_link =
        "<a href='GenomeList.pl?dsgid=$dsgids' target='_blank'>"
      . @ids
      . ' genome'
      . ( @ids > 1 ? 's' : '' ) . '</a>';

    my $feat = $coge->resultset('Feature')->find($fid);
    my $feat_url = "<a href='FeatView.pl?fid=$fid' target='_blank'>" . $feat->name . "</a>";

    my ( $source_name, $titleA ) = gen_org_name(
        dsgid     => $source_dsgid,
        write_log => 0
    );

    #my $name = "<a href='OrganismView.pl?dsgid=$source_dsgid' target='_blank'>$source_name</a>";

    my $log_msg = 'Searched ' . $list_link . '  for feature ' . $feat_url;

    my $tiny_synfind_link = CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $synfind_link
    );

    my $log = CoGe::Accessory::Web::log_history(
        db      => $coge,
        user_id => $USER->id,
        description => $log_msg,
        page    => $PAGE_TITLE,
        link => $tiny_synfind_link,
    );

   #convert numerical codes for different scoring functions to appropriate types
    if ( $scoring_function == 2 ) {
        $scoring_function = "density";
    }
    else {
        $scoring_function = "collinear";
    }

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );

    #need to blast source_dsg against each dsgids
    my @blast_results;
    my $html;

    #check for taint
    ($source_dsgid) = $source_dsgid =~ /(\d+)/;
    my $source_type = has_cds($source_dsgid);
    $source_type = $source_type ? "CDS" : "genomic";
    my @dsgids;
    foreach my $dsgid ( split( /,/, $dsgids ) ) {
        ($dsgid) = $dsgid =~ /(\d+)/;
        my $has_cds = has_cds($dsgid);
        my $feat_type = $has_cds ? "CDS" : "genomic";
        push @dsgids, [ $dsgid, $feat_type ] if $dsgid && $has_cds;
        if ( !$has_cds ) {
            my $dsg = $coge->resultset('Genome')->find($dsgid);
            CoGe::Accessory::Web::write_log( "#WARNING:#", $cogeweb->logfile );
            CoGe::Accessory::Web::write_log(
                $dsg->organism->name
                  . " does not have CDS sequences.  Can't process in SynFind.\n",
                $cogeweb->logfile
            );
        }
    }

    my $pm = new Parallel::ForkManager($MAX_PROC);

    #Generate fasta files and blastdbs
    my @to_process;
    push @to_process, [ $source_dsgid, $source_type ] if $source_dsgid;
    push @to_process, @dsgids;
    foreach my $item (@to_process) {
        $pm->start and next;
        my ( $dsgid, $feat_type ) = @$item;

        my ( $fasta, $org_name ) =
          gen_fasta( dsgid => $dsgid, feat_type => $feat_type, write_log => 1 );

#	my $blastdb = gen_blastdb(dbname=>"$dsgid-$feat_type-new",fasta=>$fasta,org_name=>$org_name);
        make_bed( dsgid => $dsgid, outfile => $BEDDIR . $dsgid . ".bed" );
        $pm->finish;
    }
    $pm->wait_all_children();

    #Generate fasta files and blastdbs
    my @target_info;    #store all the stuff about a genome
    foreach my $item (@to_process) {
        my ( $dsgid, $feat_type ) = @$item;
        my ( $fasta, $org_name ) =
          gen_fasta( dsgid => $dsgid, feat_type => $feat_type, write_log => 0 );

#	my $blastdb = gen_blastdb(dbname=>"$dsgid-$feat_type-new",fasta=>$fasta,org_name=>$org_name, write_log=>0);
        push @target_info, {
            dsgid     => $dsgid,
            feat_type => $feat_type,
            fasta     => $fasta,
            org_name  => $org_name,

            #			    blastdb=>$blastdb,
        };
    }
    my $query_info = shift @target_info
      if $source_dsgid;    #query is the first item on this list.
     #need to create blastfile name.  Must be alphabetized on query and target names.
    foreach my $target (@target_info) {
        my ( $org1, $org2 ) = ( $query_info->{org_name}, $target->{org_name} );
        my ( $dsgid1, $dsgid2 ) = ( $query_info->{dsgid}, $target->{dsgid} );
        my ( $feat_type1, $feat_type2 ) =
          ( $query_info->{feat_type}, $target->{feat_type} );
        my ( $fasta1, $fasta2 ) = ( $query_info->{fasta}, $target->{fasta} );

        #	my ($db1, $db2) = ($query_info->{blastdb},$target->{blastdb});

        (
            $org1, $org2, $dsgid1, $dsgid2, $feat_type1, $feat_type2, $fasta1,
            $fasta2
          )
          = (
            $org2, $org1, $dsgid2, $dsgid1, $feat_type2, $feat_type1, $fasta2,
            $fasta1
          ) if ( $org2 lt $org1 );

        #prep names for file system
        foreach my $tmp ( $org1, $org2 ) {
            $tmp =~ s/\///g;
            $tmp =~ s/\s+/_/g;
            $tmp =~ s/\(//g;
            $tmp =~ s/\)//g;
            $tmp =~ s/://g;
            $tmp =~ s/;//g;
            $tmp =~ s/#/_/g;
            $tmp =~ s/'//g;
            $tmp =~ s/"//g;
        }
        my $basedir = $DIAGSDIR . "/" . $org1 . "/" . $org2;
        mkpath( $basedir, 0, 0777 ) unless -d $basedir;
        my $basename =
          $dsgid1 . "_" . $dsgid2 . "." . $feat_type1 . "-" . $feat_type2;
        my $blastfile = $basedir . "/" . $basename . ".$algo";
        my $bedfile1  = $BEDDIR . $dsgid1 . ".bed";
        my $bedfile2  = $BEDDIR . $dsgid2 . ".bed";
        $target->{synteny_score_db} =
            $basedir . "/"
          . $basename . "_"
          . $window_size . "_"
          . $cutoff . "_"
          . $scoring_function
          . ".$algo" . ".db";
        $target->{basedir}             = $basedir;
        $target->{basename}            = $basename;
        $target->{blastfile}           = $blastfile;
        $target->{converted_blastfile} = $blastfile . ".short_names";
        $target->{filtered_blastfile} =
          $target->{converted_blastfile} . ".filtered";
        $target->{bedfile1}    = $bedfile1;
        $target->{bedfile2}    = $bedfile2;
        $target->{dsgid1}      = $dsgid1;
        $target->{dsgid2}      = $dsgid2;
        $target->{query_fasta} = $fasta1
          ; #need to determine the correct query/target order so output is compatibable with SynMap

#	$target->{target_db} = $db2;#need to determine the correct query/target order so output is compatibable with SynMap
        $target->{target_fasta} = $fasta2;
    }

    #blast and conquer
    foreach my $target (@target_info) {
        $pm->start and next;
        my $blastfile = $target->{blastfile};
        my $success   = run_blast(
            fasta1  => $target->{query_fasta},
            fasta2  => $target->{target_fasta},
            outfile => $blastfile,
            prog    => $algo
        );    #for last/lastz with sequence databases
        CoGe::Accessory::Web::write_log( "failed blast run for " . $blastfile,
            $cogeweb->logfile )
          unless $success;

        $blastfile = run_convert_blast(
            infile  => $blastfile,
            outfile => $target->{converted_blastfile}
        );

#	blast2bed(infile=>$blastfile, outfile1=>$target->{bedfile1}, outfile2=>$target->{bedfile2});
        run_blast2raw(
            blastfile => $blastfile,
            bedfile1  => $target->{bedfile1},
            bedfile2  => $target->{bedfile2},
            outfile   => $target->{filtered_blastfile}
        );
        unless ( -r $target->{filtered_blastfile}
            || -r $target->{filtered_blastfile} . ".gz" )
        {
            return "error creating " . $target->{filtered_blastfile}, "\n";
        }
        run_synteny_score(
            blastfile        => $target->{filtered_blastfile},
            bedfile1         => $target->{bedfile1},
            bedfile2         => $target->{bedfile2},
            outfile          => $target->{synteny_score_db},
            window_size      => $window_size,
            cutoff           => $cutoff,
            scoring_function => $scoring_function,
            dsgid1           => $target->{dsgid1},
            dsgid2           => $target->{dsgid2}
        );
        $pm->finish;
    }
    $pm->wait_all_children;
    my ( $gevo_link, $matches ) = gen_gevo_link(
        fid         => $fid,
        dbs         => [ map { $_->{synteny_score_db} } @target_info ],
        window_size => $window_size,
        depth       => $depth
    );

    my $tiny_gevo_link =
      CoGe::Accessory::Web::get_tiny_link( url => $gevo_link );
    CoGe::Accessory::Web::write_log( "#TINY GEVO LINK: $tiny_gevo_link",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Finished!", $cogeweb->logfile );

    foreach my $item (@target_info) {
        foreach my $file (
            $item->{blastfile},
            $item->{converted_blastfile},
            $item->{filtered_blastfile},
          )
        {
            $pm->start and next;

            #	    $file = CoGe::Accessory::Web::gzip($file);
            $pm->finish;
        }
    }
    $pm->wait_all_children();

    #make table of results
    my %dsgids = map { $_->[0], => 1 } @dsgids;    #table to look them up later;
    $dsgids{$source_dsgid} = 1;

    $html .=
      qq{<table id=syntelog_table class="ui-widget-content ui-corner-all">};
    $html .= qq{<THEAD><tr>};
    $html .= qq{<th>Organism};
    $html .= qq{<th>Genome};
    $html .= qq{<th>Type};
    $html .= qq{<th>Name};
    $html .= qq{<th>Chr};
    $html .= qq{<th>Synteny Score};
    $html .= qq{<th>SynMap};
    $html .= qq{</tr></THEAD><TBODY>};
    my $synmap_algo_num;    #from SynMap's lookup table on algorithms

    if ( $algo eq "lastz" ) {
        $synmap_algo_num = 4;    #lastz
    }
    elsif ( $algo eq "last" ) {
        $synmap_algo_num = 6;    #last
    }
    my $synmap_link =
      "SynMap.pl?autogo=1;b=$synmap_algo_num;dsgid1=$source_dsgid;dsgid2=";
    my @res;
    my %open_all_synmap;
    my @homologs;
    my $count = 0;
    my $last_tdsgid;
    foreach my $item ( [ $fid, "query" ], @$matches ) {
        my ( $tfid, $match_type, $synteny_score ) = @$item;
        my $feat = $coge->resultset('Feature')->find($tfid);
        my $dsg;
        foreach my $dsgt ( $feat->genomes ) {
            if ( $dsgids{ $dsgt->id } ) {

                #		delete $dsgids{$dsgt->id};
                $dsg = $dsgt;
                last;
            }
        }
	#next unless $dsg;
        $last_tdsgid = $dsg->id unless $last_tdsgid;
        $count = 0 unless $dsg->id eq $last_tdsgid;
        $last_tdsgid = $dsg->id;
        next if $depth && $count > $depth;
        $html .= qq{<tr><td>};
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
        $html .= join( "<td>",
"<span class='link' onclick=window.open('OrganismView.pl?oid=$oid')>"
              . $feat->organism->name
              . "</span>",
"<span class='link' onclick=window.open('OrganismView.pl?dsgid=$dsgid')>"
              . $dsg_name
              . "</span>",
            $match_type,
            qq{<span class='link' onclick='update_info_box("$tfid} . "_"
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
    $html .= "</tbody></table>";
    $html .=
"<br><a style='font-size: 1em' href='$tiny_gevo_link' onclick=window.open('$tiny_gevo_link') class='ui-button ui-corner-all' target=_new_gevo>Compare and visualize region in GEvo: $tiny_gevo_link</a><br><br>";

    my $featlist_link =
      gen_featlist_link( fids => [ $fid, map { $_->[0] } @$matches ] );
    $html .=
"<a href='$tiny_synfind_link' class='ui-button ui-corner-all' target=_new_synfind>Regenerate this analysis: $tiny_synfind_link</a>";
    my $open_all_synmap = join( "\n", keys %open_all_synmap );
    $html .=
qq{<a onclick='$open_all_synmap' class='ui-button ui-corner-all'>Generate all dotplots</a><br><br>};
    my $feat_list_link =
      qq{FeatList.pl?fid=} . join( ",", map { $_->id } @homologs );
    $html .=
qq{<a onclick='window.open("$feat_list_link")' class='ui-button ui-corner-all'>Generate Feature List</a>};
    my $fasta_link =
      qq{FastaView.pl?fid=} . join( ",", map { $_->id } @homologs );
    $html .=
qq{<a onclick='window.open("$fasta_link")' class='ui-button ui-corner-all'>Generate Fasta Sequences</a><br><br>};

    my $master_list_link = $synfind_link . ";get_master=1";

#    $html .= "<a onclick=window.open('$master_list_link') class='ui-button ui-corner-all' target=_new_synfind>Generate master gene set table</a>";
    $html .=
"<span onclick=get_master('$master_list_link') class='link ui-button ui-corner-all' target=_new_synfind>Generate master gene set table</span>";
    $master_list_link .= ";limit=1";

#    $html .= "<a onclick=window.open('$master_list_link') class='ui-button ui-corner-all' target=_new_synfind>Generate master gene set table (top one syntenlog per organism)</a>";
    $html .=
"<a onclick=get_master('$master_list_link') class='ui-button ui-corner-all' target=_new_synfind>Generate master gene set table (top one syntenlog per organism)</a>";
    $html .=
qq{<span class="small">Pad Sequence in GEvo <input type="text" size=11 id=pad size=11 value=0></span>};
    CoGe::Accessory::Web::write_log( "#SYNFIND LINK $synfind_link",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#TINY SYNFIND LINK $tiny_synfind_link",
        $cogeweb->logfile );
    my $log_file = $cogeweb->logfile;
    $log_file =~ s/$TEMPDIR/$TEMPURL/;

    $html .=
qq{<br><br><a href="/wiki/index.php/SynFind#Syntenic_Depth" target=_new>Syntenic Depth:</a> <span class=small>(Query genome depth coverage)</span><br>}
      . get_master_histograms(
        target_dbs       => \@target_info,
        query_dsgid      => $source_dsgid,
        unique_gene_link => $unique_gene_link
      );
    $html .= "<br>Downloads:";
    $html .= qq{<table class="small ui-widget-content ui-corner-all">};
    $html .=
      qq{<tr><td><a href=$log_file target=_new class="small">Log File</a>};

    foreach my $item (@target_info) {
    	# mdb added 10/8/13
    	my $blastfile_link = $item->{blastfile};
    	$blastfile_link =~ s/$P->{COGEDIR}//;

        $html .= qq{<tr>};
        $html .= qq{<td>} . $item->{org_name};
        $html .=
            qq{<td><a href="$blastfile_link" target=_new>Raw Blast File</a>};

        # mdb added 10/8/13
    	$blastfile_link = $item->{filtered_blastfile};
    	$blastfile_link =~ s/$P->{COGEDIR}//;

        $html .=
            qq{<td><a href="$blastfile_link" target=_new>Filtered Blast File</a>};
    }
    $html .= qq{</table>};
    return $html;
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
    my $file = $FASTADIR . "/old_synmap/$dsgid-$feat_type.fasta";
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
        system "/usr/bin/touch $file.running"
          ;    #track that a blast anlaysis is running for this
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
    my ($dsg) = $coge->resultset('Genome')->search( { genome_id => $dsgid },
        { join => 'organism', prefetch => 'organism' } );
    my $org_name = $dsg->organism->name;
    my $title =
      $org_name . " (v" . $dsg->version . ", dsgid" . $dsgid . ")" . $feat_type;
    $title =~ s/(`|')//g;
    CoGe::Accessory::Web::write_log( "#INITIALIZING#", $cogeweb->logfile )
      if $write_log;
    CoGe::Accessory::Web::write_log( "ORGANISM: " . $title, $cogeweb->logfile )
      if $write_log;
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile ) if $write_log;
    return ( $org_name, $title );
}

sub generate_fasta {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    my $file  = $opts{file};
    my $type  = $opts{type};

    my ($dsg) =
      $coge->resultset('Genome')->search( { "me.genome_id" => $dsgid },
        { join => 'genomic_sequences', prefetch => 'genomic_sequences' } );
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

#		my $title = join ("||",$chr, 1, $ds->last_chromosome_position($chr), "Chr_$chr",1, "genomic", "N/A");
            my $seq = $dsg->get_genomic_sequence( chr => $chr );
            next unless $seq;
            print OUT ">" . $chr . "\n";
            print OUT $seq, "\n";
        }
    }
    close OUT;
    return 1 if -r $file;
    CoGe::Accessory::Web::write_log( "Error with fasta file creation",
        $cogeweb->logfile );
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
    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    ( $x, $pre_command ) = CoGe::Accessory::Web::check_taint($pre_command);
    CoGe::Accessory::Web::write_log( "running $pre_command",
        $cogeweb->logfile );
    `$pre_command`;
    system "/bin/rm $outfile.running"
      if -r "$outfile.running";    #remove track file
    unless ( -s $outfile ) {
        CoGe::Accessory::Web::write_log(
"WARNING: Problem running $pre_command command.  Blast output file contains no data!",
            $cogeweb->logfile
        );
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
        CoGe::Accessory::Web::write_log( "bed file $outfile already exists",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        return $outfile;
    }
    CoGe::Accessory::Web::write_log( "Creating bedfiles: $cmd",
        $cogeweb->logfile );
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
    CoGe::Accessory::Web::write_log( "Creating bed files: $cmd",
        $cogeweb->logfile );
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
    print STDERR "In sub run_convert_blast\n";
    CoGe::Accessory::Web::gunzip($infile);
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
    print STDERR "IN SUB run_blast2raw\n";
    CoGe::Accessory::Web::gunzip($blastfile);
    CoGe::Accessory::Web::gunzip($bedfile1);
    CoGe::Accessory::Web::gunzip($bedfile2);
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
    print STDERR "Path: " . $ENV{COGE_HOME}, "\n";
    CoGe::Accessory::Web::gunzip($blastfile);    #turned on debugging
    CoGe::Accessory::Web::gunzip($bedfile1);
    CoGe::Accessory::Web::gunzip($bedfile2);
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
        my $dbh   = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );
        my $query = "SELECT * FROM  synteny where query = $fid";
        my $sth   = $dbh->prepare($query);
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
    my $featid = shift;
    my ($dsgid);
    ( $featid, $dsgid ) = split( /_/, $featid );
    my ($dsg)  = $coge->resultset('Genome')->find($dsgid);
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless ( ref($feat) =~ /Feature/i ) {
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
    my $form             = $FORM;
    my $window_size      = $form->param('ws');
    my $cutoff           = $form->param('co');
    my $scoring_function = $form->param('sf');
    my $pad              = $form->param('pad');
    my $algo             = $form->param('algo');
    if ( $scoring_function == 2 ) {
        $scoring_function = "density";
    }
    else {
        $scoring_function = "collinear";
    }

    my $qdsgid = $form->param('qdsgid');
    my $limit  = $form->param('limit')
      ;    #limit the number of syntelogs returned per organism searched
    my $qdsg = $coge->resultset('Genome')->find($qdsgid);
    my %qdsids = map { $_->id, 1 } $qdsg->datasets;
    my @dsgs;
    foreach my $item ( $form->param('dsgid') ) {
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
        foreach my $tmp ( $org1, $org2 ) {
            $tmp =~ s/\///g;
            $tmp =~ s/\s+/_/g;
            $tmp =~ s/\(//g;
            $tmp =~ s/\)//g;
            $tmp =~ s/://g;
            $tmp =~ s/;//g;
            $tmp =~ s/#/_/g;
            $tmp =~ s/'//g;
            $tmp =~ s/"//g;
        }
        my $dsgid1 = $qdsg->id;
        my $dsgid2 = $dsg->id;
        ( $org1, $org2, $dsgid1, $dsgid2 ) = ( $org2, $org1, $dsgid2, $dsgid1 )
          if ( $org2 lt $org1 );
        my $basedir  = $DIAGSDIR . "/" . $org1 . "/" . $org2;
        my $basename = $dsgid1 . "_" . $dsgid2 . "." . "CDS-CDS";
        my $db =
            $basedir . "/"
          . $basename . "_"
          . $window_size . "_"
          . $cutoff . "_"
          . $scoring_function
          . ".$algo" . ".db";
        my $dbh   = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );
        my $query = "SELECT * FROM synteny";
        my $sth   = $dbh->prepare($query);

        unless ($sth) {
            print STDERR qq{Problem connecting to $db\n};
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
            map { "ORG: " . $_->organism->name } (
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
                            $name .= $name_hash->name . ",";
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
    exit;
}

sub get_unique_genes {
    my $form             = $FORM;
    my $window_size      = $form->param('ws');
    my $cutoff           = $form->param('co');
    my $scoring_function = $form->param('sf');
    my $pad              = $form->param('pad');
    my $algo             = $form->param('algo');
    my $fid              = $form->param('fid');
    my $qdsgid           = $form->param('qdsgid');
    my $sdsgid           = $form->param('sdsgid');
    my $synfind_link     = $SERVER
      . "SynFind.pl?fid=$fid;qdsgid=$qdsgid;ws=$window_size;co=$cutoff;sf=$scoring_function;algo=$algo;dsgid=$sdsgid";

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
    foreach my $tmp ( $org1, $org2 ) {
        $tmp =~ s/\///g;
        $tmp =~ s/\s+/_/g;
        $tmp =~ s/\(//g;
        $tmp =~ s/\)//g;
        $tmp =~ s/://g;
        $tmp =~ s/;//g;
        $tmp =~ s/#/_/g;
        $tmp =~ s/'//g;
        $tmp =~ s/"//g;
    }
    my $dsgid1 = $qdsg->id;
    my $dsgid2 = $sdsg->id;
    ( $org1, $org2, $dsgid1, $dsgid2 ) = ( $org2, $org1, $dsgid2, $dsgid1 )
      if ( $org2 lt $org1 );
    my $basedir  = $DIAGSDIR . "/" . $org1 . "/" . $org2;
    my $basename = $dsgid1 . "_" . $dsgid2 . "." . "CDS-CDS";
    my $db =
        $basedir . "/"
      . $basename . "_"
      . $window_size . "_"
      . $cutoff . "_"
      . $scoring_function
      . ".$algo" . ".db";
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
    my @sunique;

    my $rs = $coge->resultset('Feature')->search(
        {
            feature_type_id => 3,
            genome_id       => $qdsgid
        },
        { join => [ { dataset => 'dataset_connectors' } ], }
    );
    foreach my $item ( $rs->get_column('feature_id')->all ) {
        push @qunique, $item unless $qdata{$item};
    }

    $rs = $coge->resultset('Feature')->search(
        {
            feature_type_id => 3,
            genome_id       => $sdsgid
        },
        { join => [ { dataset => 'dataset_connectors' } ], }
    );
    foreach my $item ( $rs->get_column('feature_id')->all ) {
        push @sunique, $item unless $sdata{$item};
    }

    my $list_type = $coge->resultset('ListType')->find_or_create(
        {
            name        => 'SynFind Unique Genes',
            description => 'Auto-generated by SynFind'
        }
    );

    #let's create a list of these items and add them to the user's lists
    my $list = $coge->resultset('List')->create(
        {
            name         => "Unique genes in " . $qdsg->organism->name,
            description  => "Compared to " . $sdsg->organism->name,
            list_type_id => $list_type->id,
            restricted   => 1
        }
    );
    return unless $list;

    my $conn = $coge->resultset('UserConnector')->create(
        {
            parent_id   => $USER->id,
            parent_type => 5,           # FIXME hardcoded to "user"
            child_id    => $list->id,
            child_type  => 1,           # FIXME hardcoded to "list"
            role_id     => 2,           # FIXME hardcoded to "owner"
        }
    );
    return unless $conn;

    my $annotation_type =
      $coge->resultset('AnnotationType')
      ->find_or_create( { name => "CoGe Note" } );
    $list->add_to_list_annotations(
        {
            annotation_type_id => $annotation_type->id,
            annotation         => "Auto created by SynFind",
        }
    );
    $list->add_to_list_annotations(
        {
            annotation_type_id => $annotation_type->id,
            annotation         => "Regenerate Analysis in SynFind",
            link => CoGe::Accessory::Web::get_tiny_link( url => $synfind_link ),
        }
    );
    $list->add_to_list_annotations(
        {
            annotation_type_id => $annotation_type->id,
            annotation         => "Query Genome: " . $qdsg->info,
            link => $SERVER . "OrganismView.pl?dsgid=" . $qdsg->id,
        }
    );
    $list->add_to_list_annotations(
        {
            annotation_type_id => $annotation_type->id,
            annotation         => "Target Genome: " . $sdsg->info,
            link => $SERVER . "OrganismView.pl?dsgid=" . $sdsg->id,
        }
    );

    #add feature_ids
    #FIXME:  hardcoded value for item_type:  id == 4 for features
    my $count = 0;
    foreach my $fid (@qunique) {
        $count++;
        $list->add_to_list_connectors_as_parent(
            { child_id => $fid, child_type => 4 } );
    }
    print $form->header;
    my $listview = $SERVER . "NotebookView.pl?lid=" . $list->id;
    print qq{<meta HTTP-EQUIV="REFRESH" content="0; url=$listview">};
    exit;
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
        while ( my @data = $sth->fetchrow_array
          )    #data base contains reciprocal entries for each item (a=>b; b=>a)
        {
            next unless $data[6] == $query_dsgid;
            my $id     = $data[0];    #query feature ID
            my $sdsgid = $data[7];    #target genome ID
            $data{$id}++;             #data{query_feature_id}{subject_genome_id}
        }
        my $genome = $coge->resultset('Genome')->find( $item->{dsgid2} );

        $html .=
            qq{<span class=link onclick='window.open("OrganismView.pl?dsgid=}
          . $item->{dsgid2} . "\")'>"
          . $genome->info
          . "</span>";
        if ( $query_dsgid eq $item->{dsgid2} ) {
            $html .=
"<br><span class=small>(self-self comparison: self-self syntenic regions ignored)</span>";
        }
        my $link =
            $unique_gene_link
          . ";ug=1;sdsgid="
          . $item->{dsgid2}
          . ";qdsgid=$query_dsgid";
        $html .=
qq{<br><span class="small link" onclick=window.open('$link')>Query genome unique genes</span>};
        $link =
            $unique_gene_link
          . ";ug=1;qdsgid="
          . $item->{dsgid2}
          . ";sdsgid=$query_dsgid";
        $html .=
qq{<br><span class="small link" onclick=window.open('$link')>Target genome unique genes</span>};
        my $depths;
        $html .= depth_table(
            depth_data => \%data,
            total      => $total_feature_count,
            depths     => \$depths
        );
        push @all, [ $depths, $item->{dsgid2} ];
    }
    $html .= master_depth_table( data => \@all );
    return $html;
}

sub master_depth_table {
    my %opts = @_;
    my $data = $opts{data};
    my $html = "<table class='small  ui-widget-content ui-corner-all'>";
    my $max  = 0;
    foreach
      my $item ( sort { scalar @{ $b->[0] } <=> scalar @{ $a->[0] } } @$data )
    {
        unless ($max) {
            $max = scalar @{ $item->[0] };
            $html .= "<tr>" . "<th>Org";
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
    my $html = "<table class='small  ui-widget-content ui-corner-all'>";
    my @depths;

    foreach my $depth ( sort { $a <=> $b } keys %depths ) {
        $html .= "<tr>";
        $html .= "<td>Depth $depth";
        $html .= "<td>" . $depths{$depth};
        $html .= "<td>of";
        $html .= "<td>$total";
        $html .=
          "<td>" . sprintf( "%.2f", $depths{$depth} / $total * 100 ) . "%";
        $html .= "</tr>";
        push @depths, $depths{$depth};
    }
    $html .= "</table>";
    $$depths = \@depths;
    return $html;
}

sub save_orglist_synfind {
    my %opts   = @_;
    my $dsgids = $opts{dsgids};
    my $prefs  = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    $prefs->{dsgids} = $dsgids;
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
