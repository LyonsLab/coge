#! /usr/bin/perl -w
use strict;
use CoGeX;
use CoGe::Accessory::Web;
use CoGe::Accessory::genetic_code;
use CoGe::Accessory::Utils qw( get_link_coords );
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use POSIX;
use Storable qw(dclone);
no warnings 'redefine';

use vars qw($P $USER $FORM $ACCN $FID $coge $PAGE_NAME $PAGE_TITLE $LINK $EMBED);

$PAGE_TITLE = 'FeatView';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$| = 1;    # turn off buffering

$FORM = new CGI;
$ACCN = $FORM->param('accn');
$FID  = $FORM->param('fid');

( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my %FUNCTION = (
    db_lookup          => \&db_lookup,
    source_search      => \&get_data_source_info_for_accn,
    get_types          => \&get_types,
    cogesearch         => \&cogesearch,
    cogesearch_featids => \&cogesearch_featids,
    get_anno           => \&get_anno,
    show_express       => \&show_express,
    gen_data           => \&gen_data,
    get_feature_types  => \&get_feature_types,
    codon_table        => \&codon_table,
    protein_table      => \&protein_table,
    gc_content         => \&gc_content,
    update_featlist    => \&update_featlist,
    parse_for_FeatList => \&parse_for_FeatList,
    get_orgs           => \&get_orgs,
    codon_aa_alignment => \&codon_aa_alignment,
);
my $pj = new CGI::Ajax(%FUNCTION);
$pj->JSDEBUG(0);
$pj->DEBUG(0);

#$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );

#print $FORM->header, "HERE!\n",gen_html();

sub gen_data {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
}

sub get_types {
    my %opts = @_;
    my ( $dsid, $gstid ) = split /_/, $opts{dsid};
    my $accn = $opts{accn};
    my $ftid = $opts{ftid};
    
    my $html;
    my %seen;
    my $search = {
        'feature_names.name' => $accn,
        dataset_id           => $dsid,
    };
    $search->{feature_type_id} = $ftid if $ftid;
    
    my @opts =
      sort map { "<OPTION>$_</OPTION>" }
      grep     { !$seen{$_}++ }
      map      { $_->type->name }
      $coge->resultset('Feature')->search( $search,, { join => 'feature_names' } );

    if (@opts) {
        $html .= "<font class='small info'>[" . scalar @opts . "]</font>\n<BR>\n";
        $html .= qq{<SELECT id="Type_name" SIZE="10" MULTIPLE onChange="get_anno(['args__accn','accn_select','args__type','Type_name', 'args__dsid','dsid', 'args__gstid','args__$gstid'],[show_anno])" >\n};
        $html .= join( "\n", @opts );
        $html .= "\n</SELECT>\n";
        $html =~ s/OPTION/OPTION SELECTED/; # set the first option to be selected
    }
    else { 
        return qq{<input type="hidden" id="Type_name">-------};
    }

    return ( $html, $gstid );
}

sub update_featlist {
    my %opts = @_;
    my $accn = $opts{accn};
    return unless $accn;
    my $type   = $opts{type};
    my $featid = $opts{fid};
    my $gstid  = $opts{gstid};
    
    my $gst    = $coge->resultset('GenomicSequenceType')->find($gstid);
    $accn .= " ($type";
    $accn .= ", " . $gst->name if $gst;
    $accn .= ")";
    return $accn, $featid, $gstid;
}

sub parse_for_FeatList {
    my $featlist = shift;
    my $url = "FeatList.pl?";
    foreach my $featid ( split(/,/, $featlist) ) {
        next unless $featid;
        $url .= "fid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub cogesearch {
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
    my $blank      = qq{<input type="hidden" id="accn_select">};
    my $weak_query = "Query needs to be better defined.";

    if ( !$accn && !$anno && !$fid ) {
        return $weak_query . $blank unless $org_id && $type;
    }
    my $html;
    my %seen;
    my $search = {};
    $search->{feature_type_id} = $type if $type;
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
        my $search1 = dclone($search);
        $search1->{"me.name"} = $accn;
        push @names, $coge->resultset('FeatureName')->search( $search1, $join );
        unless (@names) {
            push @names,
              $coge->resultset('FeatureName')->search( $search, $join )
              ->search_literal( 'MATCH(me.name) AGAINST (?)', $accn );
        }
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
        push @opts, "<OPTION>$item</OPTION>";
    }
    if ( @opts > 10000 ) {
        return $blank
          . "Search results over 10000, please refine your search.\n";
    }
    $html .=
      "<font class='small info'>[" . scalar @opts . "]</font>\n<BR>\n";
    $html .= qq{<SELECT id="accn_select" SIZE="10" MULTIPLE onChange="source_search_chain();">\n}; # 'Matches' select
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/<OPTION/<OPTION SELECTED/;
    return $blank . "No results found\n" unless $html =~ /OPTION/;
    return $html;
}

sub cogesearch_featids {
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
    my $blank      = qq{<input type="hidden" id="accn_select">};
    my $weak_query = "Query needs to be better defined.";

    if ( !$accn && !$anno && !$fid ) {
        return;
    }
    my %seen;
    my $search = {};
    $search->{feature_type_id} = $type if $type;
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
    push @{ $join->{join} }, 'annotation' if $anno;    #=>['annotation',]};
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
        my $key =
          $name->feature_id . "_" . $name->feature->dataset->sequence_type->id;
        $seen{$key} = $name->name . " (" . $name->feature->type->name;
        $seen{$key} .= ", " . $name->feature->dataset->sequence_type->name;
        $seen{$key} .= ")";
    }
    return join "||", map { $_ . "::" . $seen{$_} } keys %seen;
}

sub cogesearch_featids_old {
    my %opts = @_;
    my $accn = $opts{accn};
    $accn =~ s/^\s+//;
    $accn =~ s/\s+$//;
    my $anno = $opts{anno};
    $anno =~ s/^\s+//;
    $anno =~ s/\s+$//;
    my $type     = $opts{type};
    my $org_id   = $opts{org_id};
    my $org_name = $opts{org_name};
    $org_name = undef if $org_name =~ /^search$/i;
    my $org_desc = $opts{org_desc};
    $org_desc = undef if $org_desc =~ /^search$/i;
    my @org_ids;
    $org_id = "all" unless $org_id;

    if ( $org_id eq "all" ) {
        my ( $otype, $search ) = ( "name", $org_name )
          if $org_name && $org_name ne "Search";
        ( $otype, $search ) = ( "desc", $org_desc )
          if $org_desc && $org_desc ne "Search";
        @org_ids = get_orgs( id_only => 1, type => $otype, search => $search );
    }
    else {
        push @org_ids, $org_id;
    }
    my $feat_accn_wild = $opts{feat_name_wild};
    my $feat_anno_wild = $opts{feat_anno_wild};
    my $blank          = qq{<input type="hidden" id="accn_select">};
    my $weak_query     = "Query needs to be better defined.";
    if ( !$accn && !$anno ) {
        return $weak_query . $blank unless $org_id && $type;
    }
    my $html;
    my %seen;
    my @opts;
    $accn = "%" . $accn
      if $accn && ( $feat_accn_wild eq "both" || $feat_accn_wild eq "left" );
    $accn = $accn . "%"
      if $accn && ( $feat_accn_wild eq "both" || $feat_accn_wild eq "right" );
    $anno = "%" . $anno
      if $anno && ( $feat_anno_wild eq "both" || $feat_anno_wild eq "left" );
    $anno = $anno . "%"
      if $anno && ( $feat_anno_wild eq "both" || $feat_anno_wild eq "right" );
    my $search = {};

    if ( $accn =~ /%/ ) {
        $search->{'me.name'} = { like => $accn } if $accn;
    }
    else {
        $search->{'me.name'} = $accn if $accn;
    }
    if ( $anno =~ /%/ ) {
        $search->{'annotation'} = { like => $anno } if $anno;
    }
    else {
        $search->{'annotation'} = $anno if $anno;
    }
    $search->{feature_type_id} = $type if $type;
    $search->{organism_id}{-in} = [@org_ids] if @org_ids;
    my $join =
      { 'feature' => { 'dataset' => { 'dataset_connectors' => 'genome' } } };
    $join->{'feature'} = [ 'dataset', 'annotations' ] if $anno;
    foreach my $name (
        $coge->resultset('FeatureName')->search(
            $search,
            {
                join     => $join,
                prefetch => {
                    'feature' => {
                        'dataset' => {
                            'dataset_connectors' =>
                              { 'genome' => 'genomic_sequence_type' }
                        }
                    }
                }
            }
        )
      )
    {
        my $key =
          $name->feature_id . "_" . $name->feature->dataset->sequence_type->id;
        $seen{$key} = $name->name . " (" . $name->feature->type->name;
        $seen{$key} .= ", " . $name->feature->dataset->sequence_type->name;
        $seen{$key} .= ")";
    }
    return join "||", map { $_ . "::" . $seen{$_} } keys %seen;
}

sub get_anno {
    my %opts            = @_;
    my $accn            = $opts{accn};
    my $featlist_choice = $opts{featlist_choice};

    my $fid        = $opts{fid};
    my $type       = $opts{type};
    my $dataset_id = $opts{dsid};
    my $gstid      = $opts{gstid};
    ( $fid, $gstid ) = split /_/, $featlist_choice if $featlist_choice;
    return unless $accn || $fid;

    my @feats;
    if ($accn) {
        foreach my $feat (
            $coge->resultset('Feature')->search(
                {
                    'feature_names.name' => $accn,
                    dataset_id           => $dataset_id
                },
                { join => 'feature_names' }
            ))
        {
            push @feats, $feat if ( $feat->type->name eq $type );
        }
    }
    else {
        push @feats, $coge->resultset('Feature')->find($fid);
    }
        
    my $anno;
    $anno = qq{<h4>Annotation count: <span class="note">(} . scalar @feats . qq{)</span></h4><hr>}
      if (scalar(@feats) && !$EMBED );

    my $i = 0;
    foreach my $feat (@feats) {
        next unless $feat;
        my ($dsg) = $feat->dataset->genomes;
        return "Restricted Access" unless $USER->has_access_to_genome($dsg);
        return "Deleted" if $dsg->deleted;
        #next if ($feat->dataset->restricted && !$USER->has_access_to_dataset($feat->dataset));
        
        $i++;
        my $featid = $feat->id;
        my $chr    = $feat->chr;
        my $rc     = 0;
        my $pro    = 0;
        my $ds     = $feat->dataset->id;
        my $x      = $feat->start;
        my $z      = 4;
        my $gid    = $dsg->id;
        $anno .= qq{<div class="coge-buttonset"><span class="coge-button" onClick="window.open('FastaView.pl?fid=$featid&gstid=$gstid');">Get Sequence</span>};
        $anno .= qq{<span class="coge-button" onClick="window.open('CoGeBlast.pl?featid=$fid;gstid=$gstid');">CoGeBlast</span>};
        my ($a, $b) = get_link_coords($feat->start, $feat->stop);
        $anno .= #qq{<span class="coge-button" onClick="window.open('GenomeView.pl?chr=$chr&ds=$ds&x=$x&z=$z;gstid=$gstid');">Genome Browser</span>}; # mdb removed 11/20/13 issue 254
            qq{<span class="coge-button" onClick="window.open('GenomeView.pl?gid=$gid&loc=$chr:$a..$b');">Genome Browser</span>}; # mdb added 11/20/13 issue 254
        $anno .= qq{<span class="coge-button" onClick="window.open('SynFind.pl?fid=$featid');">SynFind</span>};

        #$anno .= qq{<DIV id="exp$i"><input type="button" value = "Click for expression tree" onClick="gen_data(['args__Generating expression view image'],['exp$i']);show_express(['args__}.$accn.qq{','args__}.'1'.qq{','args__}.$i.qq{'],['exp$i']);"></DIV>};
        $anno .= qq{<span class="coge-button" onClick="update_featlist(['args__accn', 'args__$accn','args__type', 'args__$type','args__fid', 'args__$featid', 'args__gstid','args__$gstid'],[add_to_featlist]);\$('#feat_list').dialog('option', 'width', 'auto').dialog('open');">Add to list</span></div>}
          if $accn;

        eval {
            $anno .= join "\n<hr>\n",
            '<div class="padded">' . $feat->annotation_pretty_print_html( gstid => $gstid ) . '</div>';
        };

        if ( $feat->type->name eq "CDS" ) {
            $anno .= qq{<div class="coge-buttonset">};
            $anno .= qq{<span class="coge-button" onClick="codon_table(['args__featid','args__$featid', 'args__gstid','args__$gstid'],['codon_table']); \$('#codon_table').dialog('option', 'width', 600).dialog('open');">Codon Usage</span>};
            $anno .= qq{<span class="coge-button" onClick="protein_table(['args__featid','args__$featid', 'args__gstid','args__$gstid'],['aa_table']);\$('#aa_table').dialog('open');">Amino Acid Usage</span>};
            $anno .= qq{<span class="coge-button" onClick="codon_aa_alignment(['args__featid','args__$featid', 'args__gstid','args__$gstid'],['codon_aa_alignment']); \$('#codon_aa_alignment').dialog('option', 'width', 650).dialog('open');">Codon/AA alignment</span></div>};
        }
    }
    $anno = "<h4 class=\"annotation\">No annotations for this entry</h4>" unless $anno;
    
    return ($anno);
}

sub show_express {
    my %opts = @_;
    my $accn = shift;
    my $log  = shift;
    my $div  = shift;
    $log = 1 unless defined $log;
    $accn =~ s/\..*//;
    my $link = qq{<img src="expressiontree.pl?locus=$accn&label=1&rw=80&rh=8&name=1&legend=1&mean=1&log_trans=$log">\n};
    $log = $log ? 0 : 1;
    my $type = $log ? "log transformed" : "normal";

    #    print STDERR $link;
    $link .= qq{<br><input type="button" value = "Click for $type expression tree" onClick="gen_image([],['exp$div']);show_express(['args__}
      . $accn
      . qq{','args__}
      . $log
      . qq{','args__}
      . $div
      . qq{'],['exp$div']);">};
    return $link;
}

sub gen_html {
    my $template;

    $EMBED = $FORM->param('embed') || 0;
    if ($EMBED) {
        $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'embedded_page.tmpl' );
    }
    else {
        $template = HTML::Template->new( filename => $P->{TMPLDIR}.'generic_page.tmpl' );

        $template->param(
            PAGE_TITLE    => 'FeatView',
            TITLE         => 'FeatView: Search Features Across Organisms',
            PAGE_LINK     => $LINK,
            SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
            HOME          => $P->{SERVER},
            HELP          => 'FeatView',
            WIKI_URL      => $P->{WIKI_URL} || '',
            USER          => $USER->display_name || ''
        );

        $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
        $template->param(
            ADMIN_ONLY  => $USER->is_admin,
            CAS_URL     => $P->{CAS_URL} || '',
            COOKIE_NAME => $P->{COOKIE_NAME} || ''
        );
    }

    $template->param( BODY => gen_body() );

    return $template->output;
}

sub gen_body {
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'FeatView.tmpl' );
    $template->param( EMBED     => $EMBED );
    $template->param( ACCN      => $ACCN );
    $template->param( FEAT_TYPE => get_feature_types() );
    $template->param( ORG_LIST  => get_orgs( type => "none" ) );
    
    if ($FID) {
        my $anno = get_anno( fid => $FID );
        $template->param( FID_ANNO => $anno );
    }

    return $template->output;
}

sub get_feature_types {
    my %args  = @_;
    my $orgid = $args{orgid};
    my $html  = qq{<select name="type" id="type" /><option value=0>All</option>};
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

sub get_data_source_info_for_accn {
    my %opts     = @_;
    my $accn     = $opts{accn};
    return qq{<input type="hidden" id="dsid">------------} unless $accn;
    my $org_id   = $opts{org_id};
    my $org_name = $opts{org_name};
    my $org_desc = $opts{org_desc};
    $org_id = "all" unless $org_id;
    
    my %org_ids;
    if ( $org_id eq "all" ) {
        my ( $type, $search );
        ( $type, $search ) = ( "name", $org_name ) if $org_name && $org_name ne "Search";
        ( $type, $search ) = ( "desc", $org_desc ) if $org_desc && $org_desc ne "Search";
        %org_ids =
          map { $_ => 1 }
          get_orgs( id_only => 1, type => $type, search => $search )
          if $type && $search;
    }
    else {
        $org_ids{$org_id} = 1;
    }
    
    my @feats = $coge->resultset('Feature')->search(
        { 'feature_names.name' => $accn },
        {  join       => 'feature_names',
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
        my $val = $feat->dataset;
        next unless $val;

        my ($dsg) = $feat->dataset->genomes;
        next unless ($dsg && !$dsg->deleted);

        unless ( $USER->has_access_to_genome($dsg) ) {
            $sources{"Restricted Access"} = {
                v     => 0,
                gstid => 0,
                id    => 0,
            };
            next;
        }
        unless ($val) {
            warn "error with feature: " . $feat->id . " for name $accn\n";
            next;
        }
        
        my $org = $val->organism->name if $val->organism;
        if ( keys %org_ids ) { 
            next unless $org_ids{ $val->organism->id }; 
        }
        my $name    = $val->name;
        my $ver     = $val->version;
        my $desc    = $val->description;
        my $sname   = $val->data_source->name if $val->data_source;
        my $ds_name = $val->name;
        my $dsid    = $val->id;
        my @gstypes = $val->sequence_type;

        foreach my $type (@gstypes) {
            my $gstname = $type->name;
            my $title   = "$org (gid ".$dsg->id. " v".$dsg->version."): $ds_name (dsid: $dsid, $sname, v$ver, $gstname)";
            #my $title = "$org: $ds_name (v$ver, $type)";
            $sources{$title}{id}    = $val->id;
            $sources{$title}{v}     = $ver;
            $sources{$title}{gstid} = $type->id;
        }
    }
    
    my $html = qq{<SELECT name="dsid" id="dsid" MULTIPLE SIZE="10" onChange="get_types_chain();">}; # 'Genomes' select
    my $count = 0;
    foreach my $title (
        sort { $sources{$b}{v} <=> $sources{$a}{v} || $a cmp $b }
        keys %sources
      )
    {
        my $id    = $sources{$title}{id};
        my $gstid = $sources{$title}{gstid};
        my $val   = $id . "_" . $gstid;
        $html .= qq{  <option value="$val" >$title\n};
        $html =~ s/option/option selected/ unless $count;
        $count++;
    }
    $html .= qq{</SELECT>\n};
    
    return ("<font class='small info'>[" . $count . "]</font>\n<BR>\n" . $html );
}

sub get_orgs {
    my %opts    = @_;
    my $type    = $opts{type};
    my $search  = $opts{search};
    my $id_only = $opts{id_only};
    
    my @db;
    my $count;
    if ( $type && $type eq "name" ) {
        @db = $coge->resultset("Organism")->search( { name => { like => "%".$search."%" } } );
        $count = scalar @db;
    }
    elsif ( $type && $type eq "desc" ) {
        @db = $coge->resultset("Organism")->search( { description => { like => "%".$search."%" } } );
        $count = scalar @db;
    }
    else {
        $count = $coge->resultset("Organism")->count();
        #@db = $coge->resultset("Organism")->all;
    }
    return map { $_->id } @db if $id_only;

    #my @db = $name ? $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}})
    #  : $coge->resultset('Organism')->all();
    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
        push @opts,
            "<OPTION value=\"" . $item->id . "\" id=\"o" . $item->id . "\">" . $item->name . "</OPTION>";
    }
    
    my $html =
        qq{<FONT CLASS ="small" id="org_count">Organism count: }
      . $count
      . qq{</FONT>\n<BR>\n};
    if ( $search && !@opts ) {
        $html .= qq{<input type = hidden name="org_id" id="org_id"><br>};
        $html .= "No results";
        return $html;
    }
    unshift( @opts, "<OPTION value=\"all\" id=\"all\">All Listed Organisms</OPTION>" );
    my $size = scalar @opts;
    $size = 8 if $size > 8;
    $html .= qq{<SELECT id="org_id" SIZE="$size" MULTIPLE>\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
}

sub gc_content {
    my %opts   = @_;
    my $featid = $opts{featid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ( $gc, $at ) = $feat->gc_content;
    my $html = "GC:" . ( 100 * $gc ) . "%" . ", AT:" . ( 100 * $at ) . "%";
    return $html;
}

sub codon_table {
    my %opts   = @_;
    my $featid = $opts{featid};
    my $gstid  = $opts{gstid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ( $codon, $code_type ) =
      $feat->codon_frequency( counts => 1, gstid => $gstid );
    my %aa;
    my ($code) = $feat->genetic_code;

    foreach my $tri ( keys %$code ) {
        $aa{ $code->{$tri} } += $codon->{$tri};
    }
    my $html = "Codon Usage: $code_type";
    $html .= CoGe::Accessory::genetic_code->html_code_table(
        data   => $codon,
        code   => $code,
        counts => 1
    );

    #    $html .= "Predicted amino acid usage for $code_type genetic code:";
    #    $html .= CoGe::Accessory::genetic_code->html_aa(data=>\%aa, counts=>1);
    return $html;
}

sub protein_table {
    my %opts   = @_;
    my $featid = $opts{featid};
    my $gstid  = $opts{gstid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $aa     = $feat->aa_frequency( counts => 1, gstid => $gstid );
    my $html   = "Amino Acid Usage";
    $html .= CoGe::Accessory::genetic_code->html_aa( data => $aa, counts => 1 );
    return $html;
}

sub codon_aa_alignment {
    my %opts   = @_;
    my $featid = $opts{featid};
    my $gstid  = $opts{gstid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $seq =
      join( " ", $feat->genomic_sequence( gstid => $gstid ) =~ /(...)/g );
    my $aa =
      join( "   ", split //, $feat->protein_sequence( gstid => $gstid ) );
    my @seq = $seq =~ /(.{1,80})/g;
    my @aa  = $aa  =~ /(.{1,80})/g;

    my $aln = "<pre>";
    for ( my $i = 0 ; $i < @seq ; $i++ ) {
        $aln .= $aa[$i] . "\n" . $seq[$i] . "\n\n";
    }
    $aln .= "</pre>";
    return $aln;
}
