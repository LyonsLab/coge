#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CGI::Carp 'fatalsToBrowser';
use HTML::Template;
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use DBI;
use Data::Dumper;
no warnings 'redefine';

delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
use vars qw($P $PAGE_TITLE $PAGE_NAME $DEBUG $USER $FORM $coge $LINK);

$PAGE_TITLE = 'GenomeView_old';
$PAGE_NAME  = "$PAGE_TITLE.pl";
$DEBUG      = 0;
$FORM       = new CGI;

( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my $pj = new CGI::Ajax(
    gen_html        => \&gen_html,
    grab_sequence   => \&grab_sequence,
    save_options    => \&save_options,
    get_genome_info => \&get_genome_info,
);

#$pj->js_encode_function('escape');
print $pj->build_html( $FORM, \&gen_html );

#print $FORM->header; gen_html();

sub gen_html {
    my $html;    #=  "Content-Type: text/html\n\n";
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );

    #$template->param(TITLE=>'Genome Viewer');
    $template->param( PAGE_TITLE => 'Genome Viewer',
    				  PAGE_LINK  => $LINK,
    				  HOME       => $P->{SERVER},
                      HELP       => 'GenomeView',
                      WIKI_URL   => $P->{WIKI_URL} || '' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    #$template->param(BOX_WIDTH=>"100%");
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    my ( $body, $org_name ) = gen_body();
    $template->param( BODY     => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $form = shift || $FORM;
    my ( $chr, $dsid, $z, $loc, $gstid, $fid, $dsgid, $show_legend, $prefs );
    $chr  = $form->param('chr')  if defined $form->param('chr');
    $dsid = $form->param('ds')   if $form->param('ds');
    $dsid = $form->param('dsid') if $form->param('dsid');
    $z    = $form->param('z')    if defined( $form->param('z') );
    $loc  = $form->param('x')    if $form->param('x');
    $loc =~ s/,|\.//g if $loc;    #remove commas and points if present
    $gstid       = $form->param('gstid') if $form->param('gstid');
    $fid         = $form->param('fid')   if $form->param('fid');
    $dsgid       = $form->param('dsgid') if $form->param('dsgid');
    $dsgid       = $form->param('gid')   if $form->param('gid');
    $show_legend = $form->param('sl')    if $form->param('sl');
    $prefs       = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    my ( @ds, $dsg, $gst );

    if ($fid) {
        my $feat = $coge->resultset('Feature')->find($fid);
        $chr = $feat->chromosome;
        push @ds, $feat->dataset;
        $loc = $feat->start;
    }
    if ($dsgid) {
        $dsg = $coge->resultset('Genome')->find($dsgid);
        return "unable to find genome for $dsgid" unless $dsg;
        $gst = $dsg->type;
        ($chr) = $dsg->chromosomes unless ($chr);
        @ds = $dsg->datasets( chr => $chr );
    }
    if ($dsid) {
        my $dataset = $coge->resultset('Dataset')->find($dsid);
        push @ds, $dataset;
        ($chr) = $dataset->chromosomes unless ($chr);
    }
    unless ($dsg) {
        foreach my $dsgt (
            sort {
                $a->genomic_sequence_type_id <=> $b->genomic_sequence_type_id
            } $ds[0]->genomes
          )
        {
            last if $dsgid;
            if ( $gstid && $dsgt->genomic_sequence_type_id == $gstid ) {
                $dsg = $dsgt;
                last;
            }
            $dsg = $dsgt;
        }
    }
    if ( !$USER->has_access_to_genome($dsg) ) {
        return "Permission denied";
    }

    $gst   = $dsg->type;
    $dsgid = $dsg->id;
    $gstid = $gst->id;
    my $ver        = $dsg->version;
    my $org        = $ds[0]->organism->name;
    my $chr_length = $ds[0]->last_chromosome_position($chr);

    my @feat_types;
    my $in_clause = "IN (" . ( join ",", map { $_->id } @ds ) . ")";

#my $query = qq{select distinct(feature_type_id) from feature where dataset_id = $dsid};
    my $query =
qq{select distinct(feature_type_id) from feature where dataset_id $in_clause};
    my $dbh = $coge->storage->dbh;  #DBI->connect( $connstr, $DBUSER, $DBPASS );
    my $sth = $dbh->prepare($query);
    $sth->execute;
    while ( my $row = $sth->fetchrow_arrayref ) {
        my $ftid = $row->[0];
        push @feat_types, $coge->resultset('FeatureType')->find($ftid);
    }
    $loc = 1 unless $loc;
    $z   = 5 unless defined $z;
    $z = 0 if $z < 0;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );

#$template->param( JBROWSE => $P->{SERVER}."/js/jbrowse/index.html?gid=$dsgid&loc=1%3A0..105000");

    #set tiler program
    #$template->param(TILE_SERVER=>$P->{SERVER}."tiler.pl?");
    my @servers = ( $P->{SERVER} . 'tiler.pl?' );
    my $tiler_servers = join( ",", map { "'" . $_ . "'" } @servers );

    #'http://b.coge.iplantcollaborative.org/coge/tiler.pl?',
    #'http://c.coge.iplantcollaborative.org/coge/tiler.pl?',
    #'http://d.coge.iplantcollaborative.org/coge/tiler.pl?',
    #'http://e.coge.iplantcollaborative.org/coge/tiler.pl?',
    #'http://f.coge.iplantcollaborative.org/coge/tiler.pl?',
    #};
    $template->param( TILE_SERVER => $tiler_servers );

    #set layers
    my $count = 0;
    my $base1;
    my $base2;
    my @colors = (
        [ 200, 50,  200 ],
        [ 0,   0,   200 ],
        [ 0,   200, 0 ],
        [ 200, 200, 50 ],
        [ 200, 50,  0 ],
    );
    my $color = 0;
    foreach my $ft (@feat_types) {

        if ( $ft->name eq "gene_space" ) {
            $template->param( GENE_SPACE_LAYER => 1 );
        }
        elsif ( $ft->name =~ /duplicate/i ) {
            $template->param( LOCAL_DUP_LAYER => 1 );
        }
        elsif ( $ft->name =~ /domain/i ) {
            $template->param( FUNC_DOMAIN_LAYER => 1 );
        }
        elsif ( $ft->name =~ /repeat/i ) {
            $template->param( REPEATS_LAYER => 1 );
        }
        elsif ( $ft->name =~ /transposable/i ) {
            $template->param( TE_LAYER => 1 );
        }
        elsif ( $ft->name =~ /quantitation/i ) {

            #	    $template->param(QUANT_LAYER=>1);
            my $ftid       = $ft->id;
            my $layer_name = "layer$count";
            my $name       = $ft->name;
            my ( $r, $g, $b ) = @{ $colors[$color] };
            $color++;
            $color = 0 if $color == @colors;
            $base1 .= qq{
$layer_name = new OpenLayers.Layer.Genomic( "<span id='gene_space'>$name</span>" , server
                 ,{'ds':ds,'dsg':dsg,  'gstid':gstid,'chr':chr,layers:"quant",'ftid': '$ftid', 'r':'$r', 'g':'$g', 'b':'$b'}
                 ,{isBaseLayer:false}
             );
};

            $base2 .= qq{
map.addLayer($layer_name);\n
};
        }
        $count++;
    }
    $count = 0;
    foreach my $exp ( $dsg->experiments
      ) #this genome has experiments, let's create layers for those experiments.  This will need to be greatly expanded for the final EPIC-CoGe project so that users can do fancy stuff
    {
        next
          ; #skipping loading these right now due to slowness when loading many of then
        next if $exp->deleted;
        my $expid      = $exp->id;
        my $name       = $exp->name;
        my $layer_name = "exp_layer$count";
        my ( $r, $g, $b ) = @{ $colors[$color] };
        $color++;
        $color = 0 if $color == @colors;
        $base1 .= qq{
$layer_name = new OpenLayers.Layer.Genomic( "<span id='gene_space'>$name</span>" , server
                 ,{'ds':ds,'dsg':dsg,  'gstid':gstid,'chr':chr,layers:"quant",'expid': '$expid', 'r':'$r', 'g':'$g', 'b':'$b'}
                 ,{isBaseLayer:false}
             );
};

        $base2 .= qq{
map.addLayer($layer_name);
$layer_name.setVisibility(0);
};
        $count++;
    }
    $template->param( DYNAMIC_LAYERS     => $base1 );
    $template->param( DYNAMIC_LAYERS_ADD => $base2 );
    my ($gevo_group) =
      $coge->resultset('AnnotationTypeGroup')
      ->search( { name => "gevo link" } );
    if ($gevo_group) {
        my ($anno) = $coge->resultset('FeatureAnnotation')->count(
            {
                'feature.dataset_id' => [ map { $_->id } @ds ],
                'annotation_type.annotation_type_group_id' => $gevo_group->id
            },
            { join => [ 'feature', 'annotation_type' ], limit => 1 }
        );
        $template->param( GEVO_LINK_LAYER => 1 ) if ($anno);
    }

    my ($tandem_group) =
      $coge->resultset('AnnotationTypeGroup')
      ->search( { name => "Tandem duplicates" } );
    if ($tandem_group) {
        my ($anno) = $coge->resultset('FeatureAnnotation')->count(
            {
                'feature.dataset_id' => [ map { $_->id } @ds ],
                'annotation_type.annotation_type_group_id' => $tandem_group->id
            },
            { join => [ 'feature', 'annotation_type' ], limit => 1 }
        );
        $template->param( LOCAL_DUP_LAYER => 1 ) if $anno;
    }
    $template->param( CHR           => $chr );
    $template->param( VER           => $ver );
    $template->param( ORG           => $org );
    $template->param( DS            => $dsid ) unless $dsgid;
    $template->param( DSG           => $dsgid );
    $template->param( LOC           => $loc );
    $template->param( ZOOM          => $z );
    $template->param( GSTID         => $gstid );
    $template->param( CHR_LENGTH    => $chr_length );
    $template->param( SAVE_SETTINGS => 1 ) unless $USER->user_name eq "public";
    $template->param( FLAT          => "checked" )
      if $prefs->{'flat'} && $prefs->{'flat'} eq "true";
    $template->param( EXPAND => "checked" )
      if $prefs->{'expand'} && $prefs->{'expand'} eq "true";
    $template->param( POPUPANNO => "checked" )
      if $prefs->{'popupanno'} && $prefs->{'popupanno'} eq "true";
    $template->param( SHOW_LEGEND => 1 ) if $show_legend;

    #stuff for people coming from maizegdb
    if ( $ENV{HTTP_REFERER} && $ENV{HTTP_REFERER} =~ /maizegdb/ ) {
        $template->param( MAIZEGDB => 1 );
        $template->param( EXPAND   => "checked" );
    }
    my %default_true = ( gc => 'true', genes => 'true' );
    foreach my $item (
        qw(gc gaga gbox genes wobblegc wobble50gc localdup funcdomain prot repeats other gevo_link)
      )
    {
        my $show = $prefs->{$item};
        $show = $default_true{$item} unless $show;
        $show = 'false' unless $show;
        $template->param( uc($item) => $show );
    }
    my $html = $template->output;
    my $org_name =
"<span class=link onclick=window.open('OrganismView.pl?dsgid=$dsgid')>$org (v$ver),";
    $org_name .= " " . $dsg->name         if $dsg->name;
    $org_name .= ": " . $dsg->description if $dsg->description;
    $org_name .=
        " Chromosome: $chr "
      . $gst->name
      . " (dsgid$dsgid dsid"
      . join( ",", map { $_->id } @ds )
      . ")</span>";
    return $html, $org_name;
}

sub grab_sequence {
    my $new_value   = shift;
    my $first_value = shift;
    my @vals        = ( $new_value, $first_value );
    @vals = sort { $a <=> $b } (@vals);
    return $vals[0], $vals[1];
}

sub save_options {
    my %opts = @_;
    my $opts = Dumper \%opts;
    my $item = CoGe::Accessory::Web::save_settings(
        opts => $opts,
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
}

sub get_genome_info {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return " " unless $dsgid;
    my $dsg = $coge->resultset("Genome")->find($dsgid);
    return "Unable to create genome object for id: $dsgid" unless $dsg;
    my $html = "<span class=species>Datasets:</span><br><table class='small'>";
    $html .= "<tr align='left'><th>ID<th>Name<th>Description";
    foreach my $ds ( $dsg->datasets ) {
        $html .=
            "<tr><td>"
          . $ds->id
          . "<td><a href=OrganismView.pl?dsid="
          . $ds->id
          . " target=_new>"
          . $ds->name
          . "</a><td>"
          . $ds->description;
    }
    $html .= "</table>";
    return $html;
}
