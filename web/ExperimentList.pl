#! /usr/bin/perl -w

use strict;
use CGI;
use CoGeX;
use DBI;

use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use Digest::MD5 qw(md5_base64);
use DBIxProfiler;
use File::Path;
no warnings 'redefine';

use vars qw( $P $PAGE_TITLE $PAGE_NAME $TEMPDIR $USER $DATE $BASEFILE $LINK
  $coge $cogeweb $FORM $URL $TEMPURL %FUNCTION );

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

$PAGE_TITLE = 'ExperimentList';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$TEMPDIR   = $P->{TEMPDIR} . "$PAGE_TITLE/";
$TEMPURL   = $P->{TEMPURL} . "$PAGE_TITLE/";
$ENV{PATH} = $P->{COGEDIR};
$URL       = $P->{URL};

%FUNCTION = (
    gen_html               => \&gen_html,
    get_feature_counts     => \&get_feature_counts,
    send_to_xls            => \&send_to_xls,
    send_to_csv            => \&send_to_csv,
    send_to_ExperimentList => \&send_to_ExperimentList,
    send_to_list           => \&send_to_list,
    gen_data               => \&gen_data,
    save_FeatList_settings => \&save_FeatList_settings,
    add_to_user_history    => \&add_to_user_history,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => 'ExperimentList',
				      TITLE	     => 'ExperimentList',
    				  PAGE_LINK  => $LINK,
    				  SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
    				  HOME       => $P->{SERVER},
                      HELP       => 'ExperimentList',
                      WIKI_URL   => $P->{WIKI_URL} || '' );
    $template->param( USER       => $USER->display_name || '' );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE       => $DATE );
    $template->param( BODY       => $body );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );
    $html .= $template->output;
}

sub gen_body {
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'ExperimentList.tmpl' );

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
        NameDescD   => 1,
        TypeD       => 1,
        AnnotationD => 1,
        SourceD     => 1,
        VerD        => 1
      }
      unless $prefs->{display};

    foreach my $key ( keys %{ $prefs->{display} } ) {
        $template->param( $key => "checked" );
    }

    $template->param( SAVE_DISPLAY => 1 ) unless $USER->user_name eq "public";

    #Don't show button to save list data unless valid user
    $template->param( SAVE_DATA => 1 ) unless $USER->user_name eq "public";
    my %expids;

    #	$dsgids = read_file() if $BASEFILE;    #: $opts{feature_list};
    foreach ( $form->param('dsgid') ) {
        foreach my $dsgid ( split /(::)|(,)/ ) {

            #next if $dsgid =~ /^\d+_?\d*$/;
            my $dsg = $coge->resultset('Genome')->find($dsgid);
            next unless ($dsg);
            foreach ( $dsg->experiments ) {
                $expids{ $_->id }++ unless $_->deleted;
            }
        }
    }

    #	$expids = read_file() if $BASEFILE;    #: $opts{feature_list};
    foreach ( $form->param('eid') ) {
        foreach my $expid ( split /(::)|(,)/ ) {
            $expids{$expid}++;    # if $expid =~ /^\d+_?\d*$/;
        }
    }

    foreach ( $form->param('lid') ) {
        foreach my $lid ( split /(::)|(,)/ ) {
            my $list = $coge->resultset('List')->find($lid);
            next unless ($list);
            foreach ( $list->experiments ) {
                $expids{ $_->id }++;
            }
        }
    }

    my ( $table, $count ) =
      generate_table( expids => [ sort { $a <=> $b } keys %expids ] );

    #	$template->param( 'EXPERIMENT_COUNT' => $count );

    if ($table) {
        $template->param( INFO => $table );
        return $template->output;
    }
    else {
        return "No genomes/experiments were specified.";
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
    return $opts{work_name};
}

sub generate_table {
    my %opts   = @_;
    my $expids = $opts{expids};
    return unless @$expids;
    my @table;
    my $count = 1;

    foreach my $expid (@$expids) {
        my $exp = $coge->resultset('Experiment')->find($expid);
        next unless $exp;

        # Build sub-table of types
        my %types;
        foreach ( $exp->types ) {
            my $key = $_->name;
            $types{$key}++;
        }
        my $type_tbl = join( '<br>', sort keys %types );

        # Build sub-table of annotations
        %types = ();
        foreach ( $exp->annotations ) {
            my $type      = $_->annotation_type;
            my $group     = $type->group;
            my $groupname = ( $group ? $group->name : '' );
            my $typename  = $type->name;
            push @{ $types{$groupname}{$typename} }, $_->annotation;
        }
        my $annot_tbl;
        foreach my $groupname ( sort keys %types ) {
            my $first1 = 1;
            foreach my $typename ( sort keys %{ $types{$groupname} } ) {
                my $first2 = 1;
                foreach my $annot ( sort @{ $types{$groupname}{$typename} } ) {
                    $annot_tbl .= ( $first1 ? "$groupname: " : '' )
                      if ($groupname);
                    $annot_tbl .= ( $first2 ? "$typename: " : '' );
                    $annot_tbl .= $annot . '<br>';
                    $first1 = $first2 = 0;
                }
            }
        }

        # Build source link
        my $src  = $exp->source;
        my $link = $src->link;
        my $source;
        $source .= "<span class='link' onclick=window.open('$link');>"
          if ($link);
        $source .= $src->name . ( $src->desc ? ': ' . $src->desc : '' );
        $source .= "</span>" if ($link);

        push @table,
          {
            COUNT      => $count,
            EXPID      => $expid,
            NAME       => $exp->name . " (id" . $exp->id . ")",
            DESC       => $exp->desc,
            TYPE       => $type_tbl,
            ANNOTATION => $annot_tbl,
            SOURCE     => $source,
            VER        => $exp->version
          };
        $count++;
    }
    $count--;

    return \@table, $count;
}

sub get_feature_counts {
    my %opts  = @_;
    my $dsgid = $opts{dsgid};
    return "No specified dataset group id" unless $dsgid;
    my $dsg   = $coge->resultset('Genome')->find($dsgid);
    my $query = qq{
SELECT count(distinct(feature_id)), ft.name, ft.feature_type_id
  FROM feature
  JOIN feature_type ft using (feature_type_id)
  JOIN dataset_connector dc using (dataset_id)
 WHERE genome_id = $dsgid
  GROUP BY ft.name

};
    my $dbh = $coge->storage->dbh;  #DBI->connect( $connstr, $DBUSER, $DBPASS );
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $feats = {};

    while ( my $row = $sth->fetchrow_arrayref ) {
        my $name = $row->[1];
        $name =~ s/\s+/&nbsp/g;
        $feats->{$name} = {
            count => $row->[0],
            id    => $row->[2],
            name  => $name,
        };
    }
    my $feat_string .= qq{<table class="small ui-widget-content ui-corner-all">
<thead></thead><tbody>};
    $feat_string .= "<tr valign=top>" . join(
        "\n<tr valign=top>",
        map {
                "<td valign=top><div id=$_>"
              . $feats->{$_}{name}
              . "</div>"
              . "<td valign=top align=right>"
              . $feats->{$_}{count}
          } sort { $a cmp $b } keys %$feats
    );
    $feat_string .= "</tbody></table>";
    $feat_string .= "None" unless keys %$feats;
    return $feat_string;
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

sub send_to_ExperimentList    #send to ExperimentList
{
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    #my $url = $URL . "ExperimentList.pl?eid=$accn_list";
    my $url = "ExperimentList.pl?eid=$accn_list";
    return $url;
}

sub send_to_list              #send to list
{
    my %opts      = @_;
    my $accn_list = $opts{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;

    # Create the new list
    my $list = $coge->resultset('List')->create(
        {
            name         => 'experimentlist',
            description  => 'Created by ExperimentList',
            list_type_id => 2, # FIXME hardcoded type!
            creator_id => $USER->id,
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

    # Add each experiment to the new list
    foreach my $accn ( split( /,/, $accn_list ) ) {
        next if $accn =~ /no$/;
        my ($eid) = split( '_', $accn );
        my $conn =
          $coge->resultset('ListConnector')
          ->create(
            { parent_id => $list->id, child_id => $eid, child_type => 3 } )
          ;    #FIXME hardcoded type!
        return unless $conn;
    }

    # Record in the log
    CoGe::Accessory::Web::log_history(
        db          => $coge,
        user_id     => $USER->id,
        page        => $PAGE_NAME,
        description => 'created experiment list ' . $list->info_html,
        parent_id   => $list->id,
        parent_type => 1 #FIXME magic number
    );

    my $url = "NotebookView.pl?lid=" . $list->id;
    return $url;
}

sub send_to_xls {    # FIXME mdb
    my %args      = @_;
    my $accn_list = $args{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $basename = $cogeweb->basefilename;
    my $file     = "$TEMPDIR/Excel_$basename.xls";
    my $workbook = Spreadsheet::WriteExcel->new($file);
    $workbook->set_tempdir("$TEMPDIR");
    my $worksheet = $workbook->add_worksheet();
    my $i         = 1;

    $worksheet->write( 0, 0,  "Name" );
    $worksheet->write( 0, 1,  "Description" );
    $worksheet->write( 0, 2,  "Source" );
    $worksheet->write( 0, 3,  "Provenance" );
    $worksheet->write( 0, 4,  "Sequence Type" );
    $worksheet->write( 0, 5,  "Chr Count" );
    $worksheet->write( 0, 6,  "Length (bp)" );
    $worksheet->write( 0, 7,  "Percent GC" );
    $worksheet->write( 0, 8,  "Percent AT" );
    $worksheet->write( 0, 9,  "Percent N|X" );
    $worksheet->write( 0, 10, "OrganismView Link" );

    foreach my $dsgid ( split( /,/, $accn_list ) ) {
        my ($dsg) = $coge->resultset("Genome")->find($dsgid);
        next unless $dsg;

        my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
        my $desc =
          $dsg->description ? $dsg->description : $dsg->organism->description;
        my ($ds_source) = $dsg->source;
        my $source = $ds_source->name;
        my $provenance = join( " ", map { $_->name } $dsg->datasets );
        my $length     = $dsg->length;
        my $chr_count  = $dsg->chromosome_count;
        my $type       = $dsg->type->name;
        my ( $gc, $at, $n ) = $dsg->percent_gc();
        $at *= 100;
        $gc *= 100;

        #my ($wgc, $wat) = $dsg->wobble_content;
        #$wat*=100;
        #$wgc*=100;

        $worksheet->write( $i, 0, $name );
        $worksheet->write( $i, 1, $desc );
        $worksheet->write( $i, 2, $source );
        $worksheet->write( $i, 3, $provenance );
        $worksheet->write( $i, 4, $type );
        $worksheet->write( $i, 5, $chr_count );
        $worksheet->write( $i, 6, $length );
        $worksheet->write( $i, 7, $gc . '%' );
        $worksheet->write( $i, 8, $at . '%' );
        $worksheet->write( $i, 9, $n . '%' );
        $worksheet->write( $i, 10,
            $P->{SERVER} . 'OrganismView.pl?dsgid=' . $dsgid );

        $i++;
    }
    $workbook->close() or die "Error closing file: $!";
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
}

sub send_to_csv {    # FIXME mdb
    my %args      = @_;
    my $accn_list = $args{accn_list};
    $accn_list =~ s/^,//;
    $accn_list =~ s/,$//;
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    my $basename = $cogeweb->basefilename;
    my $file     = "$TEMPDIR/$basename.csv";
    open( OUT, ">$file" );
    print OUT join( "\t",
        "CoGe Genome ID",
        "Name",
        "Description",
        "Source",
        "Provenance",
        "Sequence Type",
        "Chr Count",
        "Length (bp)",
        "Percent GC",
        "Percent AT",
        "Percent N|X",
        "OrganismView Link" ),
      "\n";

    foreach my $dsgid ( split( /,/, $accn_list ) ) {
        next unless $dsgid;
        my ($dsg) = $coge->resultset("Genome")->find($dsgid);
        next unless $dsg;
        my $name = $dsg->name ? $dsg->name : $dsg->organism->name;
        my $desc =
          $dsg->description ? $dsg->description : $dsg->organism->description;
        my ($ds_source) = $dsg->source;
        my $source = $ds_source->name;
        my $provenance = join( "||", map { $_->name } $dsg->datasets );
        my $chr_count  = $dsg->chromosome_count;
        my $length     = $dsg->length;
        my $type       = $dsg->type->name;
        my ( $gc, $at, $n ) = $dsg->percent_gc();
        $at *= 100;
        $gc *= 100;
        print OUT join( "\t",
            $dsgid,      $name,
            $desc,       $source,
            $provenance, $type,
            $chr_count,  $length,
            $gc,         $at,
            $n,          $P->{SERVER} . 'OrganismView.pl?dsgid=' . $dsgid ),
          "\n";
    }
    close OUT;
    $file =~ s/$TEMPDIR/$TEMPURL/;
    return $file;
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
