#! /usr/bin/perl -w

#TODO -- FIX FTID SO THAT NO _ AT END OF ID NUM., EX. ID='12345_'

use strict;
use CGI;
use HTML::Template;
use Spreadsheet::WriteExcel;
use CoGeX;
use CoGeX::Result::Feature;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );

use vars qw($P $PAGE_TITLE $PAGE_NAME $LINK
  $TEMPDIR $TEMPURL $USER $DATE $BASEFILE $coge $cogeweb $FORM %FUNCTION);

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$PAGE_TITLE = 'FeatList';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$ENV{PATH} = $P->{COGEDIR};
$TEMPDIR   = $P->{TEMPDIR} . 'FeatList/' . $USER->id;
$TEMPURL   = $P->{TEMPURL} . 'FeatList/' . $USER->id;

%FUNCTION = (
    send_to_gevo           => \&send_to_gevo,
    send_to_blast          => \&send_to_blast,
    send_to_fasta          => \&send_to_fasta,
    send_to_xls            => \&send_to_xls,
    codon_table            => \&codon_table,
    protein_table          => \&protein_table,
    gc_content             => \&gc_content,
    gen_data               => \&gen_data,
    send_to_featmap        => \&send_to_featmap,
    send_to_msa            => \&send_to_msa,
    send_to_featlist       => \&send_to_featlist,
    send_to_SynFind        => \&send_to_SynFind,
    get_anno               => \&get_anno,
    get_gc                 => \&get_gc,
    get_wobble_gc          => \&get_wobble_gc,
    save_FeatList_settings => \&save_FeatList_settings,
    add_to_user_history    => \&add_to_user_history,
    send_to_CodeOn         => \&send_to_CodeOn,
    send_to_list           => \&send_to_list,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );

    $template->param( TITLE => 'FeatList',
                      PAGE_TITLE => 'FeatList',
                      PAGE_LINK  => $LINK,
                      SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
                      HOME       => $P->{SERVER},
                      HELP       => 'FeatList',
                      WIKI_URL   => $P->{WIKI_URL} || '',
                      USER     => $USER->display_name || '' );
    $template->param( LOGON    => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin,
                      CAS_URL    => $P->{CAS_URL} || '' );
    return $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'FeatList.tmpl' );

    my $link = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
    $link = CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $link
    );
    $template->param( LINK => $link );

    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $sort_by_type     = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $prefs            = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    $prefs = {} unless $prefs;
    $prefs->{display} = {
        FeatNameD => 1,
        TypeD     => 1,
        ChrD      => 1,
        StartD    => 1,
        StopD     => 1,
        StrandD   => 1,
        LengthD   => 1,
        OrgD      => 1,
        AnnoD     => 1,
      }
      unless $prefs->{display};

    foreach my $key ( keys %{ $prefs->{display} } ) {
        $template->param( $key => "checked" );
    }

    $template->param( SAVE_DISPLAY => 1 ) unless $USER->user_name eq "public";

#Don't show button to save list data unless valid user
#$template->param( SAVE_DATA => 1 ) unless $USER->user_name eq "public"; # mdb removed 8/30/12

    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;    #: $opts{feature_list};
          #feat ids may be in the format of <fid>_<gstid>
    foreach my $item ( $form->param('fid') ) {
        foreach my $item2 ( split /(::)|(,)/, $item ) {
            push @$feat_list, $item2 if $item2 =~ /^\d+_?\d*$/;
        }
    }
    foreach my $item ( $form->param('featid') ) {
        foreach my $item2 ( split /(::)|(,)/, $item ) {
            push @$feat_list, $item2 if $item2 =~ /^\d+_?\d*$/;
        }
    }
    my $dsid  = $form->param('dsid')  if $form->param('dsid');
    my $dsgid = $form->param('dsgid') if $form->param('dsgid');
    my $gid = $form->param('gid') if !$dsgid && $form->param('gid');
    $dsgid = $gid if $gid;
    my $chr   = $form->param('chr')   if $form->param('chr');
    my $ftid  = $form->param('ftid')  if $form->param('ftid');
    my $start = $form->param('start') if defined $form->param('start');
    my $stop  = $form->param('stop')  if $form->param('stop');
    my $gstid = $form->param('gstid')
      if $form->param('gstid');    # genomic_sequence_type_id
    my $lid = $form->param('lid') if $form->param('lid');    # list_id
    $template->param( 'GSTID' => $gstid ) if $gstid;
    push @$feat_list,
      @{
        get_fids(
            lid   => $lid,
            dsid  => $dsid,
            dsgid => $dsgid,
            ftid  => $ftid,
            chr   => $chr,
            start => $start,
            stop  => $stop
        )
      }
      if ( $dsid || $dsgid || $lid );
    my ( $table, $feat_types, $count ) = generate_table(
        feature_list => $feat_list,
        ftid         => $ftid,
        gstid        => $gstid
    );
    $template->param( 'CDS_COUNT'             => $feat_types->{CDS} );
    $template->param( 'FEAT_COUNT'            => $count );
    $template->param( 'SHOW_ALL_CODON_TABLES' => $feat_types->{CDS} )
      if $feat_types->{CDS};

    my $type = qq{<SELECT ID="feature_type">};
    $type .= join( "\n",
        map { "<OPTION value=$_>" . $_ . "</option>" } sort keys %$feat_types )
      . "\n";
    $type .= "</select>";
    $template->param( 'FEAT_TYPES' => $type );
    if ($table) {
        $template->param( INFO => $table );
        return $template->output;
    }
    else {
        return "No feature ids were specified.";
    }
}

sub add_to_user_history {
    my %opts = @_;

    #my $url_params = $ENV{'QUERY_STRING'};
    my $url = $opts{url};
    if ( $opts{archive} ) {
        $USER->add_to_works(
            {
                'name'        => $opts{work_name},
                'archive'     => $opts{archive},
                'page'        => $PAGE_NAME,
                'parameter'   => $url,
                'description' => $opts{description},
                'note'        => $opts{note},
            }
        );
    }
    else {
        $USER->add_to_works(
            {
                'name'        => 'FeatList-' . $DATE,
                'archive'     => 0,
                'page'        => $PAGE_NAME,
                'parameter'   => $url,
                'description' => 'Feature List created on ' . $DATE,
            }
        );
    }
    return ();
}

sub get_fids {
    my %opts  = @_;
    my $lid   = $opts{lid};
    my $dsid  = $opts{dsid};
    my $dsgid = $opts{dsgid};
    my $ftid  = $opts{ftid};
    my $chr   = $opts{chr};
    my $start = $opts{start};
    my $stop  = $opts{stop};
    my $search;
    my @dsids;
    push @dsids, $dsid if $dsid;
    my @ids;

    if ($lid) {
        my $list = $coge->resultset('List')->find($lid);
        if ($list) {
            foreach ( $list->features ) {
                push @ids, $_->id;
            }
        }
        else {
            warn "unable to create list object for id $lid\n";
        }
        return \@ids;
    }

    if ($dsgid) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        if ($dsg) {
            foreach my $ds ( $dsg->datasets ) {
                push @dsids, $ds->id;
            }
        }
        else {
            warn "unable to create dsg object for id $dsgid\n";
        }
    }
    if (@dsids) {
        $search->{-or} = [ dataset_id => [@dsids] ];
    }
    $search->{feature_type_id} = $ftid if $ftid;
    $search->{chromosome}      = $chr  if $chr;

    if ( defined $start ) {
        @ids = map { $_->id } $coge->get_features_in_region(
            dataset => $dsid,
            chr     => $chr,
            start   => $start,
            stop    => $stop
        );
        foreach my $item (@dsids) {
            push @ids, map { $_->id } $coge->get_features_in_region(
                dataset => $item,
                chr     => $chr,
                start   => $start,
                stop    => $stop
            );
        }
    }
    else {
        @ids = map { $_->id } $coge->resultset('Feature')->search($search);
    }

    return \@ids;
}

sub generate_table {
    my %opts      = @_;
    my $feat_list = $opts{feature_list};
    my $ftid      = $opts{ftid};
    my $gstid     = $opts{gstid};
    $gstid = 1 unless defined $gstid;
    return unless @$feat_list;
    my @table;
    my %feat_types;
    my $count = 1;
    my %feats;

    foreach my $item (@$feat_list) {
        my ( $fid, $gstidt );
        if ( $item =~ /_/ ) {
            ( $fid, $gstidt ) = split /_/, $item;
        }
        else {
            $fid = $item;
            $gstidt = $gstid if $gstid;
        }
        my ($feat) = $coge->resultset("Feature")->find(
            { 'me.feature_id' => $fid },
            {
                join => [
                    'feature_names',
                    'locations',
                    'feature_type',
                    {
                        'dataset' =>
                          { 'dataset_connectors' => { genome => 'organism' } }
                    }
                ],
                prefetch => [
                    'feature_names',
                    'locations',
                    'feature_type',
                    {
                        'dataset' =>
                          { 'dataset_connectors' => { genome => 'organism' } }
                    }
                ],
            }
        );
        next unless $feat;
        $feats{$item} = {
            fid   => $fid,
            feat  => $feat,
            gstid => $gstidt,
        };
    }
    foreach my $item (
        sort {
            $feats{$a}{feat}->organism->name
              cmp $feats{$b}{feat}->organism->name
              || $feats{$a}{feat}->type->name cmp $feats{$b}{feat}->type->name
              || $feats{$a}{feat}->chromosome cmp $feats{$b}{feat}->chromosome
              || $feats{$a}{feat}->start <=> $feats{$b}{feat}->start
        } keys %feats
      )
    {
        my $feat = $feats{$item}{feat};
        unless ($feat) {

        #	  warn "feature id $featid failed to return a valid feature object\n";
            next;
        }
        if ($ftid) {
            next unless $feat->type->id eq $ftid;
        }
        $feat_types{ $feat->type->name }++;
        my $featid = $feat->id;
        $item = $featid . "_" . $feats{$item}{gstid} unless $item =~ /_/;
        my ($name) = $feat->names;
        my $hpp =
qq{<div id="anno_$count" class="link" onclick="get_anno('$item', $count);">Get Annotation</div>};
        my $other;
        my $cds_count = $feat_types{CDS};
        $other .=
qq{<div class='link' id='codon_usage$cds_count'><DIV onclick="\$('#codon_usage$cds_count').removeClass('link'); codon_table('$item',$cds_count);">}
          . "Click for codon usage"
          . "</DIV></DIV><input type='hidden' id='CDS$cds_count' value='$item'>"
          if $feat->type->name eq "CDS";
        my $gc = qq{Get GC};
        my $at = qq{Get GC};
        my ( $wat, $wgc );

        if ( $feat->type->name eq "CDS" ) {
            $wgc = qq{Get GC};
            $wat = qq{Get GC};
        }
        push @table, {
            COUNT  => $count,
            FEATID => $item,
            NAME   => $name,
            TYPE   => $feat->type->name,
            CHR    => $feat->chr,
            STRAND => $feat->strand,
            START  => commify( $feat->start ),
            STOP   => commify( $feat->stop ),
            ORG    => $feat->organism->name . " (v" . $feat->version . ")",
            HPP    => $hpp,
            LENGTH => commify( $feat->length() ),
            OTHER  => $other,
            AT     => $at,
            GC     => $gc,
            WAT    => $wat,
            WGC    => $wgc,

        };
        $count++;
    }
    $count--;
    return \@table, \%feat_types, $count;
}

sub gen_data {
    my $message = shift;
    return qq{<font class="loading">$message. . .</font>};
}

sub read_file {
    my $file = "$TEMPDIR/$BASEFILE.featlist";
    my @featlist;
    unless ( -r $file ) {
        warn "unable to read file $file for feature ids\n";
        return \@featlist;
    }
    open( IN, $file ) || die "can't open $file for reading: $!";
    while (<IN>) {
        chomp;
        push @featlist, $_;
    }
    close IN;
    return \@featlist;
}

sub send_to_gevo    #Send to GEvo
{

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    #my $url = $URL . "GEvo.pl?";
    my $url   = "GEvo.pl?";
    my $count = 1;

    #print STDERR $url,"\n";
    foreach my $featid ( split( /,/, $accn_list ) ) {
        $url .= "fid$count=$featid&";
        $count++;
    }
    $count--;
    return ( "alert", $count ) if $count > 20;
    $url .= "num_seqs=$count";
    return $url;
}

sub send_to_SynFind {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my ($fid) = split( /,/, $accn_list );

    #my $url = $URL . "SynFind.pl?fid=$fid";
    my $url = "SynFind.pl?fid=$fid";
    return $url;
}

sub send_to_featmap {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    #my $url = $URL . "FeatMap.pl?";
    my $url = "FeatMap.pl?";
    foreach my $featid ( split( /,/, $accn_list ) ) {
        $url .= "fid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub send_to_featlist {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    #my $url = $URL . "FeatList.pl?";
    my $url = "FeatList.pl?";
    foreach my $featid ( split( /,/, $accn_list ) ) {
        $url .= "fid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub send_to_msa {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    #my $url = $URL . "CoGeAlign.pl?";
    my $url = "CoGeAlign.pl?";
    foreach my $featid ( split( /,/, $accn_list ) ) {
        $url .= "fid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub send_to_blast    #send to cogeblast
{

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    #my $url = $URL . "CoGeBlast.pl?fid=$accn_list";
    my $url = "CoGeBlast.pl?fid=$accn_list";
    return $url;
}

sub send_to_list    #send to list
{
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    # Create the new list
    my $list = $coge->resultset('List')->create(
        {
            name         => 'featlist',
            description  => 'Created by FeatList',
            #list_type_id => 4, # FIXME hardcoded type! # mdb removed 12/14/16 COGE-800
            creator_id   => $USER->id,
            restricted   => 1
        }
    );
    return unless $list;

    # Make this user the owner of the new list
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

    # Add each feature to the new list
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ($fid) = split( '_', $accn );
        my $conn =
          $coge->resultset('ListConnector')
          ->create(
            { parent_id => $list->id, child_id => $fid, child_type => 4 } )
          ;    #FIXME hardcoded type!
        return unless $conn;
    }

    # Record in the log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => $PAGE_NAME,
        description => 'created feature list ' . $list->info_html,
        parent_id   => $list->id,
        parent_type => 1 #FIXME magic number
    );

    my $url = "NotebookView.pl?lid=" . $list->id;
    return $url;
}

sub send_to_fasta {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "FastaView.pl?";
    foreach my $featid ( split( /,/, $accn_list ) ) {
        $url .= "fid=$featid&";
    }
    $url =~ s/&$//;
    return $url;
}

sub send_to_CodeOn {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    my $url = "CodeOn.pl?fid=";
    my @list;
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ( $featid, $hspnum, $dsgid ) = $accn =~ m/^(\d+)_(\d+)?_?(\d+)?$/;
        push @list, $featid;
    }
    my %seen = ();
    @list = grep { !$seen{$_}++ } @list;
    $url .= join( "::", @list );
    $url =~ s/&$//;
    return $url;
}

sub send_to_xls {

    #my $accn_list = shift;
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $basename = $cogeweb->basefile;
    my ($filename) = $basename =~ /_(.*)$/;
    $filename = 'Excel_FeatList_' . $filename . '.xls';

    my $workbook = Spreadsheet::WriteExcel->new("$TEMPDIR/$filename");
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i         = 1;

    $worksheet->write( 0, 0,  "Feature Name" );
    $worksheet->write( 0, 1,  "Type" );
    $worksheet->write( 0, 2,  "Location" );
    $worksheet->write( 0, 3,  "Strand" );
    $worksheet->write( 0, 4,  "Chromosome" );
    $worksheet->write( 0, 5,  "Length" );
    $worksheet->write( 0, 6,  "Percent GC" );
    $worksheet->write( 0, 7,  "Percent AT" );
    $worksheet->write( 0, 8,  "Percent Wobble GC" );
    $worksheet->write( 0, 9,  "Percent Wobble AT" );
    $worksheet->write( 0, 10, "Organism (version)" );
    $worksheet->write( 0, 11, "More information" );
    $worksheet->write( 0, 12, "Link" );
    $worksheet->write( 0, 13, "Sequence DNA" );
    $worksheet->write( 0, 14, "Sequence Protein" );

    foreach my $item ( split /,/, $accn_list ) {
        my ( $featid, $gstid ) = split /_/, $item;
        my ($feat) = $coge->resultset("Feature")->find($featid);

        next unless $feat;
        my ($name) = sort $feat->names;
        my $app = $feat->annotation_pretty_print();
        $app =~ s/(<\/?span(\s*class=\"\w+\")?\s*>)?//ig;

        #my ($anno) = $app =~ /annotation:<td>(.+)?/i;
        #($anno) = split (/<BR/, $anno);
        #$anno =~ s/;/;\n/g;
        my ( $at, $gc ) = $feat->gc_content;
        $at *= 100;
        $gc *= 100;
        my ( $wat, $wgc ) = $feat->wobble_content;
        $wat *= 100;
        $wgc *= 100;
        $worksheet->write( $i, 0, $P->{SERVER} . "FeatView.pl?accn=$name",
            $name );
        $worksheet->write( $i, 1, $feat->type->name );
        $worksheet->write( $i, 2, $feat->start . "-" . $feat->stop );
        $worksheet->write( $i, 3, $feat->strand );
        $worksheet->write( $i, 4, $feat->chr );
        $worksheet->write( $i, 5, $feat->length );
        $worksheet->write( $i, 6, $gc );
        $worksheet->write( $i, 7, $at );
        $worksheet->write( $i, 8, $wgc );
        $worksheet->write( $i, 9, $wat );
        $worksheet->write( $i, 10,
            $feat->organism->name . "(v " . $feat->version . ")" );
        $worksheet->write( $i, 11, $app );

        if ( my ($geneid) = $app =~ /geneid.*?(\d+)/i ) {
            my $link =
"http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=full_report&list_uids="
              . $geneid;
            $worksheet->write( $i, 12, $link );
        }
        my $seq = $feat->genomic_sequence( gstid => $gstid );
        $worksheet->write( $i, 13, $seq );
        if ( $feat->type->name eq "CDS" ) {
            $worksheet->write( $i, 14, $feat->protein_sequence() );
        }
        $i++;
    }
    $workbook->close() or die "Error closing file: $!";
    return "$TEMPURL/$filename";
}

sub gc_content {
    my %args   = @_;
    my $featid = $args{featid};
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ( $gc, $at ) = $feat->gc_content;
    my $html = "GC:" . ( 100 * $gc ) . "%" . ", AT:" . ( 100 * $at ) . "%";
    return $html;
}

sub codon_table {
    my %args   = @_;
    my $featid = $args{fid};
    my $gstid  = $args{gstid};
    ( $featid, $gstid ) = split /_/, $featid if $featid =~ /_/;
    return unless $featid;
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my ( $codon, $code_type ) =
      $feat->codon_frequency( counts => 1, gstid => $gstid );
    my %aa;
    my ($code) = $feat->genetic_code;
    my $count = 0;

    foreach my $tri ( keys %$code ) {
        $aa{ $code->{$tri} } += $codon->{$tri};
        $count += $codon->{$tri};
    }
    my $html;
    $html .= "<table><tr valign=top><td>";
    $html .= "Codon Usage: $code_type<br>";
    my ( $at, $gc ) = $feat->gc_content( gstid => $gstid );
    $at *= 100;
    $gc *= 100;
    my ( $wat, $wgc ) = $feat->wobble_content( gstid => $gstid );
    $wat *= 100;
    $wgc *= 100;
    $html .=
      "Codon Count: $count" . ", GC: $at% $gc%" . ", Wobble GC: $wat% $wgc%";
    $html .= "<td>";
    $html .= "Predicted amino acid usage";
    $html .= "<tr valign=top><td>";
    $html .= CoGe::Accessory::genetic_code->html_code_table(
        data   => $codon,
        code   => $code,
        counts => 1
    );

    #$html .= "</div>";
    #$html .= "Predicted amino acid usage for $code_type genetic code:";
    $html .= "<td>";
    $html .= CoGe::Accessory::genetic_code->html_aa(
        data   => \%aa,
        counts => 1,
        split  => 1
    );
    $html .= "</table>";
    return $html;
}

sub protein_table {
    my %args   = @_;
    my $featid = $args{featid};
    my ($feat) = $coge->resultset('Feature')->find($featid);
    my $aa     = $feat->aa_frequency( counts => 1 );
    my $html   = "Amino Acid Usage";
    $html .= CoGe::Accessory::genetic_code->html_aa( data => $aa, counts => 1 );
    return $html;
}

sub get_anno {
    my %opts  = @_;
    my $fid   = $opts{fid};
    my $gstid = $opts{gstid};
    ( $fid, $gstid ) = split /_/, $fid if ( $fid =~ /_/ );
    return unless $fid;
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return "No feature for id $fid" unless $feat;
    return $feat->annotation_pretty_print_html( gstid => $gstid );
}

sub get_gc {
    my %opts  = @_;
    my $fid   = $opts{fid};
    my $gstid = $opts{gstid};
    return unless $fid;
    ( $fid, $gstid ) = split /_/, $fid if ( $fid =~ /_/ );
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return "No feature for id $fid" unless $feat;
    my ( $gc, $at, $n ) = $feat->gc_content( gstid => $gstid );
    $at *= 100;
    $gc *= 100;
    return encode_json( [ $gc, $at ] );
}

sub get_wobble_gc {
    my %opts  = @_;
    my $fid   = $opts{fid};
    my $gstid = $opts{gstid};
    ( $fid, $gstid ) = split /_/, $fid if ( $fid =~ /_/ );
    return unless $fid;
    my ($feat) = $coge->resultset('Feature')->find($fid);
    return unless $feat->type->name eq "CDS";
    return "No feature for id $fid" unless $feat;
    my ( $wgc, $wat ) = $feat->wobble_content( gstid => $gstid );
    $wat *= 100;
    $wgc *= 100;
    return encode_json( [ $wgc, $wat ] );
}

sub save_FeatList_settings {
    my %opts    = @_;
    my $display = $opts{display};
    my %save;
    if ($display) {
        my %settings = (
            1  => 'FeatNameD',
            2  => 'TypeD',
            3  => 'ChrD',
            4  => 'StartD',
            5  => 'StopD',
            6  => 'StrandD',
            7  => 'LengthD',
            8  => 'GCD',
            9  => 'ATD',
            10 => 'WGCD',
            11 => 'WATD',
            12 => 'OrgD',
            13 => 'AnnoD',
            14 => 'OtherD',
        );
        foreach my $index ( split( /,/, $display ) ) {
            $save{display}{ $settings{$index} } = 1;
        }
    }
    CoGe::Accessory::Web::save_settings(
        opts => \%save,
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
}
