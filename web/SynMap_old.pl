#! /usr/bin/perl -w
use v5.10;

use strict;
use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Algos::KsCalc;
use CoGe::Accessory::Utils qw( commify );
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use DBIxProfiler;
use Data::Dumper;
use HTML::Template;
use LWP::Simple;
use Parallel::ForkManager;
use GD;
use File::Path;
use File::Spec;
use Benchmark;
use DBI;
use POSIX;
use Sort::Versions;

no warnings 'redefine';

umask(0);
our (
    $P,             $DBNAME,         $DBHOST,      $DBPORT,
    $DBUSER,        $DBPASS,         $connstr,     $DATE,
    $DEBUG,         $DIR,            $URL,         $SERVER,
    $USER,          $FORM,           $coge,        $cogeweb,
    $PAGE_NAME,     $FORMATDB,       $BLAST,       $TBLASTX,
    $BLASTN,        $BLASTP,         $LASTZ,       $LAST,
    $DATADIR,       $FASTADIR,       $BLASTDBDIR,  $DIAGSDIR,
    $MAX_PROC,      $DAG_TOOL,       $PYTHON,      $PYTHON26,
    $TANDEM_FINDER, $RUN_DAGCHAINER, $EVAL_ADJUST, $FIND_NEARBY,
    $DOTPLOT,       $SVG_DOTPLOT,    $NWALIGN,     $QUOTA_ALIGN,
    $CLUSTER_UTILS, $BLAST2RAW,      $BASE_URL,    $BLAST2BED,
    $SYNTENY_SCORE, $TEMPDIR,        $TEMPURL,     $ALGO_LOOKUP,
    $GZIP,          $GUNZIP,         $COOKIE_NAME, %FUNCTIONS,
    $CONFIG
);

$P = CoGe::Accessory::Web::get_defaults();
$ENV{PATH} = join ":",
  (
    $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin"
  );
$ENV{BLASTDB}    = $P->{BLASTDB};
$ENV{BLASTMAT}   = $P->{BLASTMATRIX};
$ENV{PYTHONPATH} = File::Spec->catdir( $P->{BINDIR}, "dagchainer_bp" );

$DEBUG    = 0;
$BASE_URL = $P->{SERVER};
$DIR      = $P->{COGEDIR};
$URL      = $P->{URL};
$SERVER   = $P->{SERVER};
$TEMPDIR  = $P->{TEMPDIR} . "SynMap";
$TEMPURL  = $P->{TEMPURL} . "SynMap";
$FORMATDB = $P->{FORMATDB};
$MAX_PROC = $P->{MAX_PROC};
$BLAST    = $P->{BLAST} . " -a " . $MAX_PROC . " -K 80 -m 8 -e 0.0001";
my $blast_options = " -num_threads $MAX_PROC -evalue 0.0001 -outfmt 6";
$TBLASTX = $P->{TBLASTX} . $blast_options;
$BLASTN  = $P->{BLASTN} . $blast_options;
$BLASTP  = $P->{BLASTP} . $blast_options;
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
$GZIP   = $P->{GZIP};
$GUNZIP = $P->{GUNZIP};

#in the web form, each sequence search algorithm has a unique number.  This table identifies those and adds appropriate options
$ALGO_LOOKUP = {
    0 => {
        algo => $BLASTN . " -task megablast",    #megablast
        opt             => "MEGA_SELECT",  #select option for html template file
        filename        => "megablast",
        displayname     => "MegaBlast",
        html_select_val => 0,
        formatdb        => 1,
    },
    1 => {
        algo     => $BLASTN . " -task dc-megablast",   #discontinuous megablast,
        opt      => "DCMEGA_SELECT",
        filename => "dcmegablast",
        displayname     => "Discontinuous MegaBlast",
        html_select_val => 1,
        formatdb        => 1,
    },
    2 => {
        algo            => $BLASTN . " -task blastn",    #blastn
        opt             => "BLASTN_SELECT",
        filename        => "blastn",
        displayname     => "BlastN",
        html_select_val => 2,
        formatdb        => 1,
    },
    3 => {
        algo            => $TBLASTX,                     #tblastx
        opt             => "TBLASTX_SELECT",
        filename        => "tblastx",
        displayname     => "TBlastX",
        html_select_val => 3,
        formatdb        => 1,
    },
    4 => {
        algo            => $LASTZ,                       #lastz
        opt             => "LASTZ_SELECT",
        filename        => "lastz",
        displayname     => "(B)lastZ",
        html_select_val => 4,
    },
    5 => {
        algo            => $BLASTP . " -task blastp",    #blastn
        opt             => "BLASTP_SELECT",
        filename        => "blastp",
        displayname     => "BlastP",
        html_select_val => 5,
        formatdb        => 1,
    },
    6 => {
        algo            => $LAST,                        #last
        opt             => "LAST_SELECT",
        filename        => "last",
        displayname     => "Last",
        html_select_val => 6,
    },
};

$DATADIR  = $P->{DATADIR};
$DIAGSDIR = $P->{DIAGSDIR};
$FASTADIR = $P->{FASTADIR};

mkpath( $TEMPDIR,     1, 0777 );
mkpath( $FASTADIR,    1, 0777 );
mkpath( $DIAGSDIR,    1, 0777 );    # mdb added 7/9/12
mkpath( $P->{LASTDB}, 1, 0777 );    # mdb added 7/9/12
$BLASTDBDIR = $P->{BLASTDB};

$PYTHON        = $P->{PYTHON};         #this was for python2.5
$PYTHON26      = $P->{PYTHON};
$DAG_TOOL      = $P->{DAG_TOOL};
$BLAST2BED     = $P->{BLAST2BED};
$TANDEM_FINDER = $P->{TANDEM_FINDER}
  . " -d 5 -s -r"
  ; #-d option is the distance (in genes) between dups -- not sure if the -s and -r options are needed -- they create dups files based on the input file name

#$RUN_DAGHAINER = $DIR."/bin/dagchainer/DAGCHAINER/run_DAG_chainer.pl -E 0.05 -s";
$RUN_DAGCHAINER = $P->{DAGCHAINER};
$EVAL_ADJUST    = $P->{EVALUE_ADJUST};

$FIND_NEARBY = $P->{FIND_NEARBY}
  . " -d 20"
  ; #the parameter here is for nucleotide distances -- will need to make dynamic when gene order is selected -- 5 perhaps?

#programs to run Haibao Tang's quota_align program for merging diagonals and mapping coverage
$QUOTA_ALIGN   = $P->{QUOTA_ALIGN};     #the program
$CLUSTER_UTILS = $P->{CLUSTER_UTILS};   #convert dag output to quota_align input
$BLAST2RAW     = $P->{BLAST2RAW};       #find local duplicates
$SYNTENY_SCORE = $P->{SYNTENY_SCORE};
$CONFIG = $ENV{COGE_HOME} . 'coge.conf';
$DOTPLOT     = $P->{DOTPLOT} . " -cf " . $CONFIG;
$SVG_DOTPLOT = $P->{SVG_DOTPLOT};

#$CONVERT_TO_GENE_ORDER = $DIR."/bin/SynMap/convert_to_gene_order.pl";
#$NWALIGN = $DIR."/bin/nwalign-0.3.0/bin/nwalign";
$NWALIGN = $P->{NWALIGN};
$|       = 1;                           # turn off buffering
$DATE    = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);
$FORM = new CGI;

$PAGE_NAME = "SynMap_old.pl";

my %ajax = CoGe::Accessory::Web::ajax_func();

#$ajax{read_log}=\&read_log_test;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect( $connstr, $DBUSER, $DBPASS );

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    cookie_name => $COOKIE_NAME,
    ticket      => $cas_ticket,
    coge        => $coge,
    this_url    => $FORM->url()
) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge
) unless $USER;

#print $pj->build_html( $FORM, \&gen_html );

#print "Content-Type: text/html\n\n";print gen_html($FORM);

%FUNCTIONS = (
    go                     => \&go,
    get_orgs               => \&get_orgs,
    get_genome_info        => \&get_genome_info,
    get_previous_analyses  => \&get_previous_analyses,
    get_pair_info          => \&get_pair_info,
    check_address_validity => \&check_address_validity,
    generate_basefile      => \&generate_basefile,
    get_dotplot            => \&get_dotplot,
    gen_dsg_menu           => \&gen_dsg_menu,
    get_dsg_gc             => \&get_dsg_gc,
    get_query_link         => \&get_query_link,
    read_log               => \&CoGe::Accessory::Web::read_log,
    %ajax,
);

my $pj = new CGI::Ajax(%FUNCTIONS);
if ( $FORM->param('jquery_ajax') ) {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};

    #print STDERR Dumper \%args;
    if ( $fname and defined $FUNCTIONS{$fname} ) {
        if ( $args{args} ) {
            my @args_list = split( /,/, $args{args} );
            print $FORM->header, $FUNCTIONS{$fname}->(@args_list);
        }
        else {
            print $FORM->header, $FUNCTIONS{$fname}->(%args);
        }
    }
}
else {
    $pj->js_encode_function('escape');
    print $pj->build_html( $FORM, \&gen_html );

    #	print $FORM->header; print gen_html();
}

sub read_log_test {
    my %args    = @_;
    my $logfile = $args{logfile};
    my $prog    = $args{prog};
    return unless $logfile;
    $logfile .= ".log" unless $logfile =~ /log$/;
    $logfile = $TEMPDIR . "/$logfile" unless $logfile =~ /^$TEMPDIR/;
    return unless -r $logfile;
    my $str;
    open( IN, $logfile );

    while (<IN>) {
        $str .= $_;
    }
    close IN;
    return $str;
}

sub gen_html {
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => 'SynMap' );
    $template->param( TITLE      => 'Whole Genome Synteny' );
    $template->param( HEAD       => qq{} );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( DATE => $DATE );

    #$template->param(ADJUST_BOX=>1);
    $template->param( LOGO_PNG => "SynMap-logo.png" );
    $template->param( BODY     => $body );
    $template->param( HELP     => "/wiki/index.php?title=SynMap" );
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $form = shift || $FORM;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SynMap_old.tmpl' );

    $template->param( MAIN => 1 );

    #$template->param( EMAIL       => $USER->email )  if $USER->email;

    my $master_width = $FORM->param('w') || 0;
    $template->param( MWIDTH => $master_width );

    #set search algorithm on web-page
    if ( defined( $FORM->param('b') ) ) {
        $template->param(
            $ALGO_LOOKUP->{ $FORM->param('b') }{opt} => "selected" );
    }
    else {
        $template->param( $ALGO_LOOKUP->{6}{opt} => "selected" );
    }
    my ( $D, $A, $Dm, $gm, $dt, $dupdist, $cscore );
    $D       = $FORM->param('D');
    $A       = $FORM->param('A');
    $Dm      = $FORM->param('Dm');
    $gm      = $FORM->param('gm');
    $gm      = 40 unless defined $gm;
    $dt      = $FORM->param('dt');
    $cscore  = $FORM->param('csco');
    $dupdist = $FORM->param('tdd');

#	$cvalue = $FORM->param('c');       #different c value than the one for cytology.  But if you get that, you probably shouldn't be reading this code

    my $display_dagchainer_settings;
    if ( $D && $A && $dt ) {
        my $type;
        if ( $dt =~ /gene/i ) {
            $type = " genes";
            $template->param( 'DAG_GENE_SELECT' => 'checked' );
        }
        else {
            $type = " bp";
            $template->param( 'DAG_DISTANCE_SELECT' => 'checked' );
        }
        $display_dagchainer_settings =
          qq{display_dagchainer_settings([$D,$A, '$gm', $Dm],'$type');};
    }
    else {
        $template->param( 'DAG_GENE_SELECT' => 'checked' );
        $display_dagchainer_settings = qq{display_dagchainer_settings();};
    }

    #	$cvalue = 4 unless defined $cvalue;
    #	$template->param( 'CVALUE'                      => $cvalue );
    $dupdist = 10 unless defined $dupdist;
    $template->param( 'DUPDIST' => $dupdist );
    $template->param( 'CSCORE' => $cscore ) if defined $cscore;
    $template->param(
        'DISPLAY_DAGCHAINER_SETTINGS' => $display_dagchainer_settings );
    $template->param( 'MIN_CHR_SIZE' => $FORM->param('mcs') )
      if $FORM->param('mcs');

    #will the program automatically run?
    my $autogo = $FORM->param('autogo');
    $autogo = 0 unless defined $autogo;
    $template->param( AUTOGO => $autogo );

    #populate organism menus
    for ( my $i = 1 ; $i <= 2 ; $i++ ) {
        my $dsgid = $form->param( 'dsgid' . $i ) || 0;
        my $feattype_param = $FORM->param( 'ft' . $i )
          if $FORM->param( 'ft' . $i );
        my $name = $FORM->param( 'name' . $i ) if $FORM->param( 'name' . $i );
        my $org_menu = gen_org_menu(
            dsgid          => $dsgid,
            num            => $i,
            feattype_param => $feattype_param,
            name           => $name
        );
        $template->param( "ORG_MENU" . $i => $org_menu );
    }

    #set ks for coloring syntenic dots
    if ( $FORM->param('ks') ) {
        if ( $FORM->param('ks') eq 1 ) {
            $template->param( KS1 => "selected" );
        }
        elsif ( $FORM->param('ks') eq 2 ) {
            $template->param( KS2 => "selected" );
        }
        elsif ( $FORM->param('ks') eq 3 ) {
            $template->param( KS3 => "selected" );
        }
    }
    else {
        $template->param( KS0 => "selected" );
    }

    #set color_scheme
    my $cs = 1;
    $cs = $FORM->param('cs') if defined $FORM->param('cs');
    $template->param( "CS" . $cs => "selected" );

    #set codeml min and max
    my $codeml_min;
    $codeml_min = $FORM->param('cmin') if defined $FORM->param('cmin');
    my $codeml_max;
    $codeml_max = $FORM->param('cmax') if defined $FORM->param('cmax');
    $template->param( 'CODEML_MIN' => $codeml_min ) if defined $codeml_min;
    $template->param( 'CODEML_MAX' => $codeml_max ) if defined $codeml_max;
    my $logks;
    $logks = $FORM->param('logks') if defined $FORM->param('logks');
    $logks = 1 unless defined $logks;    #turn on by default if not specified
    $template->param( 'LOGKS' => "checked" ) if defined $logks && $logks;

    #show non syntenic dots:  on by default
    my $snsd = 0;
    $snsd = $FORM->param('snsd') if ( defined $FORM->param('snsd') );
    $template->param( 'SHOW_NON_SYN_DOTS' => 'checked' ) if $snsd;

    #are the axes flipped?
    my $flip = 0;
    $flip = $FORM->param('flip') if ( defined $FORM->param('flip') );
    $template->param( 'FLIP' => 'checked' ) if $flip;

    #are the chromosomes labeled?
    my $clabel = 1;
    $clabel = $FORM->param('cl') if ( defined $FORM->param('cl') );
    $template->param( 'CHR_LABEL' => 'checked' ) if $clabel;

    #are the chromosomes labeled?
    my $skip_rand = 1;
    $skip_rand = $FORM->param('sr') if ( defined $FORM->param('sr') );
    $template->param( 'SKIP_RAND' => 'checked' ) if $skip_rand;

    #what is the sort for chromosome display?
    my $chr_sort_order = "S";
    $chr_sort_order = $FORM->param('cso') if ( defined $FORM->param('cso') );
    if ( $chr_sort_order =~ /N/i ) {
        $template->param( 'CHR_SORT_NAME' => 'selected' );
    }
    elsif ( $chr_sort_order =~ /S/i ) {
        $template->param( 'CHR_SORT_SIZE' => 'selected' );
    }

    #set axis metric for dotplot
    if ( $FORM->param('ct') ) {
        if ( $FORM->param('ct') eq "inv" ) {
            $template->param( 'COLOR_TYPE_INV' => 'selected' );
        }
        elsif ( $FORM->param('ct') eq "diag" ) {
            $template->param( 'COLOR_TYPE_DIAG' => 'selected' );
        }
    }
    else {
        $template->param( 'COLOR_TYPE_NONE' => 'selected' );
    }
    if ( $FORM->param('am') && $FORM->param('am') =~ /g/i ) {
        $template->param( 'AXIS_METRIC_GENE' => 'selected' );
    }
    else {
        $template->param( 'AXIS_METRIC_NT' => 'selected' );
    }

    #axis relationship:  will dotplot be forced into a square?
    if ( $FORM->param('ar') && $FORM->param('ar') =~ /s/i ) {
        $template->param( 'AXIS_RELATIONSHIP_S' => 'selected' );
    }
    else {
        $template->param( 'AXIS_RELATIONSHIP_R' => 'selected' );
    }

    #merge diags algorithm
    if ( $FORM->param('ma') ) {
        $template->param( QUOTA_MERGE_SELECT => 'selected' )
          if $FORM->param('ma') eq "1";
        $template->param( DAG_MERGE_SELECT => 'selected' )
          if $FORM->param('ma') eq "2";
    }
    if ( $FORM->param('da') ) {
        if ( $FORM->param('da') eq "1" ) {
            $template->param( QUOTA_ALIGN_SELECT => 'selected' );
        }
    }
    my $depth_org_1_ratio = 1;
    $depth_org_1_ratio = $FORM->param('do1') if $FORM->param('do1');
    $template->param( DEPTH_ORG_1_RATIO => $depth_org_1_ratio );
    my $depth_org_2_ratio = 1;
    $depth_org_2_ratio = $FORM->param('do2') if $FORM->param('do2');
    $template->param( DEPTH_ORG_2_RATIO => $depth_org_2_ratio );
    my $depth_overlap = 40;
    $depth_overlap = $FORM->param('do') if $FORM->param('do');
    $template->param( DEPTH_OVERLAP => $depth_overlap );

    $template->param( 'BOX_DIAGS' => "checked" ) if $FORM->param('bd');
    my $spa = $FORM->param('sp') if $FORM->param('sp');
    $template->param( 'SYNTENIC_PATH' => "checked" ) if $spa;
    $template->param( 'SHOW_NON_SYN' => "checked" ) if $spa && $spa =~ /2/;
    $template->param( 'SPA_FEW_SELECT'  => "selected" ) if $spa && $spa > 0;
    $template->param( 'SPA_MORE_SELECT' => "selected" ) if $spa && $spa < 0;
    my $file = $form->param('file');
    if ($file) {
        my $results = read_file($file);
        $template->param( RESULTS => $results );
    }

#place to store fids that are passed into SynMap to highlight that pair in the dotplot (if present)
    my $fid1 = 0;
    $fid1 = $FORM->param('fid1') if $FORM->param('fid1');

    $template->param( 'FID1' => $fid1 );
    my $fid2 = 0;
    $fid2 = $FORM->param('fid2') if $FORM->param('fid2');
    $template->param( 'FID2'      => $fid2 );
    $template->param( 'PAGE_NAME' => $PAGE_NAME );
    $template->param( 'TEMPDIR'   => $TEMPDIR );
    return $template->output;
}

sub gen_org_menu {
    my %opts           = @_;
    my $oid            = $opts{oid};
    my $num            = $opts{num};
    my $name           = $opts{name};
    my $desc           = $opts{desc};
    my $dsgid          = $opts{dsgid};
    my $feattype_param = $opts{feattype_param};
    $feattype_param = 1 unless $feattype_param;
    my $org;

    if ($dsgid) {
        $org = $coge->resultset('Genome')->find($dsgid)->organism;
        $oid = $org->id;
    }
    $name = "Search" unless $name;
    $desc = "Search" unless $desc;
    my $menu_template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'SynMap_old.tmpl' );
    $menu_template->param( ORG_MENU   => 1 );
    $menu_template->param( NUM        => $num );
    $menu_template->param( 'ORG_NAME' => $name );
    $menu_template->param( 'ORG_DESC' => $desc );
    $menu_template->param(
        'ORG_LIST' => get_orgs( name => $name, i => $num, oid => $oid ) );
    my ($dsg_menu) = gen_dsg_menu( oid => $oid, dsgid => $dsgid, num => $num );
    $menu_template->param( DSG_MENU => $dsg_menu );

    if ($dsgid) {
        my ( $dsg_info, $feattype_menu, $message ) = get_genome_info(
            dsgid    => $dsgid,
            org_num  => $num,
            feattype => $feattype_param
        );
        $menu_template->param( DSG_INFO       => $dsg_info );
        $menu_template->param( FEATTYPE_MENU  => $feattype_menu );
        $menu_template->param( GENOME_MESSAGE => $message );
    }
    return $menu_template->output;
}

sub gen_dsg_menu {
    my $t1    = new Benchmark;
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $num   = $opts{num};
    my $dsgid = $opts{dsgid};
    my @dsg_menu;
    my $message;
    my $org_name;

    #	print STDERR"here!\n";

    foreach my $dsg (
        $coge->resultset('Genome')->search(
            { organism_id => $oid },
            { prefetch    => ['genomic_sequence_type'] }
        )
      )
    {
        my $name;
        my $has_cds = 0;
        if ( $dsg->restricted && !$USER->has_access_to_genome($dsg) ) {
            next unless $dsgid && $dsg->id == $dsgid;
            $name = "Restricted";
        }
        else {
            $name .= $dsg->name . ": " if $dsg->name;
            $name .=
              $dsg->type->name . " (v" . $dsg->version . ",id" . $dsg->id . ")";
            $org_name = $dsg->organism->name unless $org_name;
            foreach my $ft (
                $coge->resultset('FeatureType')->search(
                    {
                        genome_id            => $dsg->id,
                        'me.feature_type_id' => 3
                    },
                    {
                        join =>
                          { features => { dataset => 'dataset_connectors' } },
                        rows => 1,
                    }
                )
              )
            {
                $has_cds = 1;
            }
        }
        push @dsg_menu, [ $dsg->id, $name, $dsg, $has_cds ];

    }

    my $dsg_menu = qq{
   <select id=dsgid$num onChange="\$('#dsg_info$num').html('<div class=dna_small class=loading class=small>loading. . .</div>'); get_genome_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],[handle_dsg_info])">
};
    foreach (
        sort {
                 versioncmp( $b->[2]->version, $a->[2]->version )
              || $a->[2]->type->id <=> $b->[2]->type->id
              || $b->[3] cmp $a->[3]
        } @dsg_menu
      )
    {
        my ( $numt, $name ) = @$_;
        my $selected = " selected" if $dsgid && $numt == $dsgid;
        $selected = " " unless $selected;
        $numt = 0 if $name eq "Restricted";
        $dsg_menu .= qq{
   <OPTION VALUE=$numt $selected>$name</option>
};
    }
    $dsg_menu .= "</select>";
    my $t2 = new Benchmark;
    my $time = timestr( timediff( $t2, $t1 ) );

    #    print STDERR qq{
    #-----------------
    #sub gen_dsg_menu runtime:  $time
    #-----------------
    #};
    return ( $dsg_menu, $message );

}

sub read_file {
    my $file = shift;

    my $html;
    open( IN, $TEMPDIR . $file ) || die "can't open $file for reading: $!";
    while (<IN>) {
        $html .= $_;
    }
    close IN;
    return $html;
}

sub get_orgs {
    my %opts = @_;
    my $name = $opts{name};
    my $desc = $opts{desc};
    my $oid  = $opts{oid};
    my $i    = $opts{i};
    my @db;

    #get rid of trailing white-space
    $name =~ s/^\s+//g if $name;
    $name =~ s/\s+$//g if $name;
    $desc =~ s/^\s+//g if $desc;
    $desc =~ s/\s+$//g if $desc;

    $name = ""
      if $name && $name =~ /Search/;    #need to clear to get full org count
    $desc = ""
      if $desc && $desc =~ /Search/;    #need to clear to get full org count
    my $org_count;
    if ($oid) {
        my $org = $coge->resultset("Organism")->find($oid);
        $name = $org->name if $org;
        push @db, $org if $name;
    }
    elsif ($name) {
        @db =
          $coge->resultset("Organism")
          ->search( { name => { like => "%" . $name . "%" } } );
    }
    elsif ($desc) {
        @db =
          $coge->resultset("Organism")
          ->search( { description => { like => "%" . $desc . "%" } } );
    }
    else {
        $org_count = $coge->resultset("Organism")->count;
    }

    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @db ) {
        my $option = "<OPTION value=\"" . $item->id . "\"";
        $option .= " selected" if $oid && $oid == $item->id;
        $option .= ">" . $item->name . " (id" . $item->id . ")</OPTION>";
        push @opts, $option;

    }
    $org_count = scalar @opts unless $org_count;
    my $html;
    $html .=
        qq{<FONT CLASS ="small">Organism count: }
      . $org_count
      . qq{</FONT>\n<BR>\n};
    unless ( @opts && ( $name || $desc ) ) {
        $html .= qq{<input type = hidden name="org_id$i" id="org_id$i">};
        return $html;
    }

    $html .=
qq{<SELECT id="org_id$i" SIZE="5" MULTIPLE onChange="get_genome_info_chain($i)" >\n};
    $html .= join( "\n", @opts );
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/ unless $oid;
    return $html;
}

sub get_genome_info {
    my $t1       = new Benchmark;
    my %opts     = @_;
    my $dsgid    = $opts{dsgid};
    my $org_num  = $opts{org_num};
    my $feattype = $opts{feattype};
    $feattype = 1 unless defined $feattype;
    return " ", " ", " " unless $dsgid;
    my $html_dsg_info;

#my ($dsg) = $coge->resultset("Genome")->find({genome_id=>$dsgid},{join=>['organism','genomic_sequences'],prefetch=>['organism','genomic_sequences']});
    my ($dsg) = $coge->resultset("Genome")->find( { genome_id => $dsgid },
        { join => 'organism', prefetch => 'organism' } );
    return " ", " ", " " unless $dsg;
    my $org     = $dsg->organism;
    my $orgname = $org->name;
    $orgname =
        "<a href=\"OrganismView.pl?oid="
      . $org->id
      . "\" target=_new>$orgname</a>";
    my $org_desc;
    if ( $org->description ) {
        $org_desc = join(
            "; ",
            map {
                    qq{<span class=link onclick="\$('#org_desc}
                  . qq{$org_num').val('$_').focus();search_bar('org_desc$org_num'); timing('org_desc$org_num')">$_</span>}
              } split /\s*;\s*/,
            $org->description
        );
    }
    $html_dsg_info .=
      qq{<div><span>Organism: </span><span class="small">$orgname</span></div>};
    $html_dsg_info .=
qq{<div><span>Description:</span><span class="small">$org_desc</span></div>};
    $html_dsg_info .=
qq{<div><span class="link" onclick=window.open('OrganismView.pl?dsgid=$dsgid')>Genome Information: </span><br>};
    my $i = 0;

    my (
        $percent_gc, $percent_at, $percent_n, $percent_x, $chr_length,
        $chr_count,  $plasmid,    $contig,    $scaffold
    ) = get_dsg_info($dsg);
    my ($ds) = $dsg->datasets;
    my $link = $ds->data_source->link;
    $link = $BASE_URL unless $link;
    $link = "http://" . $link unless $link && $link =~ /^http/;
    $html_dsg_info .= qq{<table class=small>};
    $html_dsg_info .= "<tr><td>Name: <td>" . $dsg->name if $dsg->name;
    $html_dsg_info .= "<tr><td>Description: <td>" . $dsg->description
      if $dsg->description;
    $html_dsg_info .=
        "<tr><td>Source:  <td><a href="
      . $link
      . " target=_new>"
      . $ds->data_source->name . "</a>";

    #$html_dsg_info .= $dsg->chr_info(summary=>1);
    $html_dsg_info .= "<tr><td>Dataset: <td>" . $ds->name;
    $html_dsg_info .= ": " . $ds->description if $ds->description;
    $html_dsg_info .= "<tr><td>Chromosome count: <td>" . commify($chr_count);
    if ( $percent_gc > 0 ) {
        $html_dsg_info .=
"<tr><td>DNA content: <td>GC: $percent_gc%, AT: $percent_at%, N: $percent_n%, X: $percent_x%";
    }
    else {
        $html_dsg_info .=
qq{<tr><td>DNA content: <td id=gc_content$org_num class='link' onclick="get_gc($dsgid, 'gc_content$org_num')">Click to retrieve};
    }
    $html_dsg_info .= "<tr><td>Total length: <td>" . commify($chr_length);
    $html_dsg_info .= "<tr><td>Contains plasmid" if $plasmid;
    $html_dsg_info .= "<tr><td>Contains contigs" if $contig;
    $html_dsg_info .= "<tr><td>Contains scaffolds" if $scaffold;
    $html_dsg_info .= "</table>";
    if ( $dsg->restricted && !$USER->has_access_to_genome($dsg) ) {
        $html_dsg_info = "Restricted";
    }
    my $t2 = new Benchmark;
    my $time = timestr( timediff( $t2, $t1 ) );

    #    print STDERR qq{
    #-----------------
    #sub get_genome_info runtime:  $time
    #-----------------
    #};

    my $message;

    #create feature type menu
    my $has_cds;
    foreach my $ft (
        $coge->resultset('FeatureType')->search(
            {
                genome_id            => $dsg->id,
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
    my ( $cds_selected, $genomic_selected ) = ( " ", " " );
    $cds_selected     = "selected" if $feattype eq 1 || $feattype eq "CDS";
    $genomic_selected = "selected" if $feattype eq 2 || $feattype eq "genomic";

    my $feattype_menu = qq{
  <select id="feat_type$org_num" name ="feat_type$org_num">
#};
    $feattype_menu .= qq{
   <OPTION VALUE=1 $cds_selected>CDS</option>
} if $has_cds;
    $feattype_menu .= qq{
   <OPTION VALUE=2 $genomic_selected>genomic</option>
};
    $feattype_menu .= "</select>";
    $message = "<span class='small alert'>No Coding Sequence in Genome</span>"
      unless $has_cds;

    return $html_dsg_info, $feattype_menu, $message, $chr_length, $org_num,
      $dsg->organism->name, $dsg->genomic_sequence_type_id;
}

sub gen_fasta {
    my %opts      = @_;
    my $dsgid     = $opts{dsgid};
    my $feat_type = $opts{feat_type};
    my $write_log = $opts{write_log} || 0;
    my ( $org_name, $title );
    ( $org_name, $title ) = gen_org_name(
        dsgid     => $dsgid,
        feat_type => $feat_type,
        write_log => $write_log
    );

    #we already have genomic sequences
    my $file;
    if ( $feat_type eq "genomic" ) {
        my $dsg = $coge->resultset('Genome')->find($dsgid);
        $file = $dsg->file_path;
    }
    else {
        $file = $FASTADIR . "/old_synmap/$dsgid-$feat_type.fasta";
    }
    my $res;
    if ($write_log) {
        CoGe::Accessory::Web::write_log( "#" x 20, $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Generating fasta file:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Generating $feat_type fasta sequence for " . $org_name,
            $cogeweb->logfile );
    }
    while ( -e "$file.running" ) {
        print STDERR "detecting $file.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $file ) {
        CoGe::Accessory::Web::write_log(
            "Fasta file for *" . $org_name . "* ($file) exists",
            $cogeweb->logfile )
          if $write_log;
        $res = 1;
    }
    else {
        system "touch $file.running"
          ;    #track that a blast anlaysis is running for this
        $res =
          generate_fasta( dsgid => $dsgid, file => $file, type => $feat_type )
          unless -r $file;
        system "/bin/rm $file.running" if -r "$file.running"; #remove track file
    }
    CoGe::Accessory::Web::write_log( "#" x (20) . "\n", $cogeweb->logfile )
      if $write_log;
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
        $org_name . " (v"
      . $dsg->version
      . ", dsgid"
      . $dsgid . ") "
      . $feat_type;
    $title =~ s/(`|')//g;

    if ($write_log) {
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Generating organism name:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( $title, $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "#" x (20) . "\n", $cogeweb->logfile );
    }
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

    $file = $FASTADIR . "/old_synmap/$file" unless $file =~ /$FASTADIR/;
    CoGe::Accessory::Web::write_log( "Creating fasta file:",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( $file, $cogeweb->logfile );

    open( OUT, ">$file" ) || die "Can't open $file for writing: $!";
    if ( $type eq "CDS" || $type eq "protein" ) {
        my $count = 1;
        my @feats =
          sort { $a->chromosome cmp $b->chromosome || $a->start <=> $b->start }
          $coge->resultset('Feature')->search(
            {
                feature_type_id => [ 3, 5, 8 ],
                genome_id       => $dsgid
            },
            {
                join => [ { dataset => 'dataset_connectors' } ],
                prefetch => ['feature_names']
            }
          );
        CoGe::Accessory::Web::write_log(
            "Getting sequence for "
              . scalar(@feats)
              . " features of types CDS, tRNA, and rRNA.",
            $cogeweb->logfile
        );
        foreach my $feat (@feats) {
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
            if ( $type eq "CDS" ) {
                my $seq = $feat->genomic_sequence( dsgid => $dsg->id );
                next unless $seq;

                #skip sequences that are only 'x' | 'n';
                next unless $seq =~ /[^x|n]/i;
                print OUT ">" . $title . "\n";
                print OUT $seq, "\n";
                $count++;
            }
            elsif ( $type eq "protein" ) {
                next unless $feat->feature_type_id == 3;
                my (@seqs) = $feat->protein_sequence( dsgid => $dsg->id );
                next unless scalar @seqs;
                next
                  if scalar @seqs > 1;   #didn't find the correct reading frame;
                next unless $seqs[0] =~ /[^x]/i;
                $title = ">" . $title . "\n";

                #		print OUT $title, join ($title, @seqs),"\n";
                print OUT $title, $seqs[0], "\n";
                $count++;
            }
        }
    }
    else {
        my @chr = sort $dsg->get_chromosomes;
        CoGe::Accessory::Web::write_log(
            "Getting sequence for "
              . scalar(@chr)
              . " chromosomes (genome sequence)",
            $cogeweb->logfile
        );
        $file = $dsg->file_path;

        #	foreach my $chr (@chr)
        #	  {
        #	    my $seq = $dsg->get_genomic_sequence(chr=>$chr);
        #	    next unless $seq;
        #	    print OUT ">".$chr."\n";
        #	    print OUT $seq,"\n";
        #	  }
    }
    close OUT;
    return 1 if -r $file;
    CoGe::Accessory::Web::write_log( "Error with fasta file creation",
        $cogeweb->logfile );
    return 0;
}

sub gen_blastdb {
    my %opts      = @_;
    my $dbname    = $opts{dbname};
    my $fasta     = $opts{fasta};
    my $org_name  = $opts{org_name};
    my $type      = $opts{type};
    my $write_log = $opts{write_log} || 0;
    my $blastdb   = "$BLASTDBDIR/$dbname";
    my $res       = 0;

    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile )
      if $write_log;
    CoGe::Accessory::Web::write_log( "Generating BlastDB file",
        $cogeweb->logfile );
    while ( -e "$blastdb.running" ) {

        print STDERR "detecting $blastdb.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $blastdb . ".nsq" || -r $blastdb . ".psq" ) {
        CoGe::Accessory::Web::write_log(
            "blastdb file for " . $org_name . " already exists",
            $cogeweb->logfile )
          if $write_log;
        $res = 1;
    }
    else {
        system "touch $blastdb.running"
          ;    #track that a blast anlaysis is running for this
        $res = generate_blast_db(
            fasta     => $fasta,
            blastdb   => $blastdb,
            org       => $org_name,
            write_log => $write_log,
            type      => $type
        );
        system "/bin/rm $blastdb.running"
          if -r "$blastdb.running";    #remove track file
    }
    CoGe::Accessory::Web::write_log( "blastdb file: $blastdb",
        $cogeweb->logfile )
      if $write_log;
    CoGe::Accessory::Web::write_log( "#" x (20) . "\n", $cogeweb->logfile )
      if $write_log;
    return $blastdb if $res;
    return 0;
}

sub generate_blast_db {
    my %opts  = @_;
    my $fasta = $opts{fasta};
    my $type  = $opts{type};

    #    exit if $fasta =~ /protein/;
    my $blastdb = $opts{blastdb};

    #    my $title= $opts{title};
    my $org       = $opts{org};
    my $write_log = $opts{write_log};
    my $command   = $FORMATDB;
    $command .= $type eq "protein" ? " -p T" : " -p F";
    $command .= " -i '$fasta'";
    $command .= " -t '$org'";
    $command .= " -n '$blastdb'";
    CoGe::Accessory::Web::write_log(
        "creating blastdb for *" . $org . "* ($blastdb): $command",
        $cogeweb->logfile )
      if $write_log;
    `$command`;

    if ( $type eq "protein" ) {
        return 1 if -r $blastdb . ".psq";
    }
    else {
        return 1 if -r "$blastdb.nsq";
    }
    CoGe::Accessory::Web::write_log(
        "error creating blastdb for $org ($blastdb)",
        $cogeweb->logfile )
      if $write_log;
    return 0;
}

sub generate_basefile {
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    return $cogeweb->basefilename;
}

sub run_blast {
    my %opts    = @_;
    my $fasta   = $opts{fasta};
    my $blastdb = $opts{blastdb};
    my $outfile = $opts{outfile};
    my $prog    = $opts{prog};
    my $gid1 = $opts{gid1};
    $prog = "1" unless $prog;
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $outfile || -r "$outfile.gz" ) {
        unless ( -s $outfile || -s "$outfile.gz" ) {
            CoGe::Accessory::Web::write_log(
                "WARNING: Blast output file ($outfile) contains no data!",
                $cogeweb->logfile );
            return 0;
        }
        CoGe::Accessory::Web::write_log( "blastfile $outfile already exists",
            $cogeweb->logfile );
        return 1;
    }
    my $pre_command =
      $ALGO_LOOKUP->{$prog}{algo};    #$prog =~ /tblastx/i ? $TBLASTX : $BLASTN;
    if ( $pre_command =~ /lastz/i ) {
        $pre_command .= " -i $fasta -d $blastdb -o $outfile";
    }
    elsif ( $pre_command =~ /last_wrapper/i ) {
        $pre_command .= "$gid1 $blastdb $fasta -o $outfile";
    }
    else {
        $pre_command .= " -out $outfile -query $fasta -db $blastdb";
    }
    my $x;
    system
      "touch $outfile.running"; #track that a blast anlaysis is running for this
    ( $x, $pre_command ) = CoGe::Accessory::Web::check_taint($pre_command);
    CoGe::Accessory::Web::write_log( "running:\n\t$pre_command",
        $cogeweb->logfile );
    `/usr/bin/nice --adjustment=10 $pre_command`;
    system "/bin/rm $outfile.running"
      if -r "$outfile.running";    #remove track file
    unless ( -s $outfile ) {
        CoGe::Accessory::Web::write_log(
"WARNING: Problem running $pre_command command.  Blast output file contains no data!",
            $cogeweb->logfile
        );
        return 0;
    }
    return 1 if -r $outfile;
}

######
#Target this program for integration into makeflow
######
sub blast2bed {
    my %opts     = @_;
    my $infile   = $opts{infile};
    my $outfile1 = $opts{outfile1};
    my $outfile2 = $opts{outfile2};
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
    CoGe::Accessory::Web::gunzip("$infile.gz") if -r "$infile.gz";
    ##target for makeflow
    #build command
    my $cmd =
      $BLAST2BED . " -infile $infile -outfile1 $outfile1 -outfile2 $outfile2";

    #log command to logfile for this SynMap run
    CoGe::Accessory::Web::write_log( "Creating bed files: $cmd",
        $cogeweb->logfile );

    #run command
    `$cmd`;
}

######
#Target this program for integration into makeflow
######
sub run_blast2raw {
    my %opts            = @_;
    my $blastfile       = $opts{blastfile};
    my $bedfile1        = $opts{bedfile1};
    my $bedfile2        = $opts{bedfile2};
    my $outfile         = $opts{outfile};
    my $tandem_distance = $opts{tandem_distance};
    my $cscore          = $opts{cscore};

    #set defaults if necessary
    $tandem_distance = 10 unless defined $tandem_distance;

#check to see if job is running.  If it is, wait.  Using a <filename>.running file to show that a program is running
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }

    #return output file if it already exists
    if (   ( -r $outfile && -s $outfile )
        || ( -r "$outfile.gz" && -s "$outfile.gz" ) )
    {
        CoGe::Accessory::Web::write_log(
"Filtered blast file found where tandem dups have been removed: $outfile",
            $cogeweb->logfile
        );
        return $outfile;
    }

    #decompress input files
    CoGe::Accessory::Web::gunzip("$blastfile.gz") if -r "$blastfile.gz";
    CoGe::Accessory::Web::gunzip("$bedfile1.gz")  if -r "$bedfile1.gz";
    CoGe::Accessory::Web::gunzip("$bedfile2.gz")  if -r "$bedfile2.gz";

    #target for makeflow
    #build command
    my $cmd = $BLAST2RAW
      . " $blastfile --localdups --qbed $bedfile1 --sbed $bedfile2 --tandem_Nmax $tandem_distance";
    $cmd .= " --cscore=$cscore" if $cscore;
    $cmd .= " > $outfile";

    #log proceedure
    CoGe::Accessory::Web::write_log( "finding and removing local duplications",
        $cogeweb->logfile );

    #creating .running file to show that job is running
    system "/usr/bin/touch $outfile.running";    #track that this is running

    #log that command is running
    CoGe::Accessory::Web::write_log( "running:\n\t$cmd", $cogeweb->logfile );

    #run command
    `$cmd`;

    #remove .running file so other programs know that it is complete
    system "/bin/rm $outfile.running"
      if -r "$outfile.running";                  #remove track file

    #return outfile
    return $outfile;
}

sub run_synteny_score {
    my %opts      = @_;
    my $blastfile = $opts{blastfile};
    my $bedfile1  = $opts{bedfile1};
    my $bedfile2  = $opts{bedfile2};
    my $outfile   = $opts{outfile};

    my $window_size      = 40;
    my $cutoff           = .1;
    my $scoring_function = "collinear";

    $outfile .=
      "_" . $window_size . "_" . $cutoff . "_" . $scoring_function . ".db";
    return $outfile if -r $outfile && -s $outfile;
    my $cmd = $SYNTENY_SCORE
      . " $blastfile --qbed $bedfile1 --sbed $bedfile2 --window $window_size --cutoff $cutoff --scoring $scoring_function --sqlite $outfile";

    CoGe::Accessory::Web::write_log( "Synteny Score:  running:\n\t$cmd",
        $cogeweb->logfile );
    system("$PYTHON26 $cmd &");
    return $outfile;

#python /opt/apache/CoGe/bin/quota-alignment/scripts/synteny_score.py 3068_8.CDS-CDS.blastn.blast.filtered --qbed 3068_8.blastn.blast.q.bed --sbed 3068_8.CDS-CDS.blastn.blast.s.bed --sqlite 3068_8.CDS-CDS.db --window $window_size --cutoff $cutoff
}

sub process_local_dups_file {
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};
    if (   ( -r $outfile && -s $outfile )
        || ( -r $outfile . ".gz" && -s $outfile . ".gz" ) )
    {
        CoGe::Accessory::Web::write_log(
            "Processed tandem duplicate file found: $outfile",
            $cogeweb->logfile );
        return $outfile;
    }
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    return unless -r $infile;
    CoGe::Accessory::Web::write_log(
"Adding coge links to tandem duplication file.  Infile $infile : Outfile $outfile",
        $cogeweb->logfile
    );
    $/ = "\n";
    open( IN,  $infile );
    open( OUT, ">$outfile" );
    print OUT "#",
      join( "\t",
        "FeatList_link", "GEvo_link", "FastaView_link",
        "chr||start||stop||name||strand||type||database_id||gene_order" ),
      "\n";

    while (<IN>) {
        chomp;
        next unless $_;
        my @line = split /\t/;
        my %fids;
        foreach (@line) {
            my @item = split /\|\|/;
            next unless $item[6];
            $fids{ $item[6] } = 1;
        }
        next unless keys %fids;
        my $featlist = $BASE_URL . "FeatList.pl?";
        map { $featlist .= "fid=$_;" } keys %fids;
        my $fastaview = $BASE_URL . "FastaView.pl?";
        map { $fastaview .= "fid=$_;" } keys %fids;
        my $gevo  = $BASE_URL . "GEvo.pl?";
        my $count = 1;

        foreach my $id ( keys %fids ) {
            $gevo .= "fid$count=$id;";
            $count++;
        }
        $gevo .= "num_seqs=" . scalar keys %fids;
        print OUT join( "\t", $featlist, $gevo, $fastaview, @line ), "\n";
    }
    close OUT;
    close IN;
    `/bin/rm $infile`;
    $infile =~ s/localdups/nolocaldups.bed/;
    `/bin/rm $infile`;
    return $outfile;
}

sub run_dag_tools {
    my %opts       = @_;
    my $query      = $opts{query};
    my $subject    = $opts{subject};
    my $blast      = $opts{blast};
    my $outfile    = $opts{outfile};
    my $feat_type1 = $opts{feat_type1};
    my $feat_type2 = $opts{feat_type2};
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $outfile || -r $outfile . ".gz" ) {
        CoGe::Accessory::Web::write_log(
            "run dag_tools: file $outfile already exists",
            $cogeweb->logfile );
        return 1;
    }
    CoGe::Accessory::Web::gunzip("$blast.gz") if -r "$blast.gz";
    unless ( -r $blast && -s $blast ) {
        CoGe::Accessory::Web::write_log(
"WARNING:   Cannot create input file for DAGChainer! Blast output file ($blast) contains no data!",
            $cogeweb->logfile
        );
        return 0;
    }
    my $query_dup_file   = $opts{query_dup_files};
    my $subject_dup_file = $opts{subject_dup_files};
    my $cmd = "$PYTHON $DAG_TOOL -q \"$query\" -s \"$subject\" -b $blast";
    $cmd .= " -c";    # if $feat_type1 eq "genomic" && $feat_type2 eq "genomic";
    $cmd .= " --query_dups $query_dup_file" if $query_dup_file;
    $cmd .= " --subject_dups $subject_dup_file" if $subject_dup_file;
    $cmd .= " > $outfile";
    system "/usr/bin/touch $outfile.running"
      ;               #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "run dag_tools: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm $outfile.running"
      if -r "$outfile.running";    #remove track file

    unless ( -s $outfile ) {
        CoGe::Accessory::Web::write_log(
            "WARNING: DAGChainer input file ($outfile) contains no data!",
            $cogeweb->logfile );
        return 0;
    }
    return 1 if -r $outfile;
}

sub run_tandem_finder    #FIXME mdb 2/28/13 - this doesn't get called ...?
{
    my %opts    = @_;
    my $infile  = $opts{infile};    #dag file produced by dat_tools.py
    my $outfile = $opts{outfile};
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    unless ( -r $infile && -s $infile ) {
        CoGe::Accessory::Web::write_log(
"WARNING:   Cannot run tandem finder! DAGChainer input file ($infile) contains no data!",
            $cogeweb->logfile
        );
        return 0;
    }
    if ( -r $outfile ) {
        CoGe::Accessory::Web::write_log(
            "run_tandem_filter: file $outfile already exists",
            $cogeweb->logfile );
        return 1;
    }
    my $cmd = "$PYTHON $TANDEM_FINDER -i $infile > $outfile";
    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "run_tandem_filter: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm $outfile.running"
      if -r "$outfile.running";    #remove track file
    return 1 if -r $outfile;
}

sub run_adjust_dagchainer_evals {
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};
    my $cvalue  = $opts{cvalue};
    $cvalue = 4 unless defined $cvalue;
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $outfile || -r $outfile . ".gz" ) {
        CoGe::Accessory::Web::write_log(
            "run_adjust_dagchainer_evals: file $outfile already exists",
            $cogeweb->logfile );
        return 1;
    }
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    unless ( -r $infile && -s $infile ) {
        CoGe::Accessory::Web::write_log(
"WARNING:   Cannot adjust dagchainer evals! DAGChainer input file ($infile) contains no data!",
            $cogeweb->logfile
        );
        return 0;
    }
    my $cmd = "$PYTHON $EVAL_ADJUST -c $cvalue $infile > $outfile";

#There is a parameter that can be passed into this to filter repetitive sequences more or less stringently:
# -c   2 gets rid of more stuff; 10 gets rid of less stuff; default is 4
#consider making this a parameter than can be adjusted from SynMap -- will need to actually play with this value to see how it works
#if implemented, this will require re-naming all the files to account for this parameter
#and updating the auto-SynMap link generator for redoing an analysis

    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log(
        "run_adjust_dagchainer_evals: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm $outfile.running" if -r "$outfile.running";
    ;      #remove track file
    return 1 if -r $outfile;

}

sub run_convert_to_gene_order {
    my %opts         = @_;
    my $infile       = $opts{infile};
    my $dsgid1       = $opts{dsgid1};
    my $dsgid2       = $opts{dsgid2};
    my $ft1          = $opts{ft1};
    my $ft2          = $opts{ft2};
    my $genomic_flag = 0;
    my %genomic_order;
    $/ = "\n";
    my $outfile = $infile . ".go";

    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $outfile || -r $outfile . ".gz" ) {
        CoGe::Accessory::Web::write_log(
            "run_convert_to_gene_order: file $outfile already exists",
            $cogeweb->logfile );
        return $outfile;
    }
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    unless ( -r $infile && -s $infile ) {
        CoGe::Accessory::Web::write_log(
"WARNING:   Cannot convert to gene order! DAGChainer input file ($infile) contains no data!",
            $cogeweb->logfile
        );
        return 0;
    }

    if ( $ft1 eq "genomic" || $ft2 eq "genomic" ) {

        #we need to process genomic names differently!
        $genomic_flag = 1;
        open( IN, $infile );
        while (<IN>) {

            #get positions;
            chomp;
            next if (/^#/);
            my @line  = split /\t/;
            my @item1 = split /\|\|/, $line[1];
            my @item2 = split /\|\|/, $line[5];
            $genomic_order{1}{ $line[0] }{ $line[1] }{start} = $line[2]
              if $ft1 eq "genomic";
            $genomic_order{2}{ $line[4] }{ $line[5] }{start} = $line[6]
              if $ft2 eq "genomic";
        }
        close IN;
        if ( $ft1 eq "genomic" ) {

            #sort positions
            foreach my $chr ( keys %{ $genomic_order{1} } ) {
                my $i = 1;
                foreach my $item ( sort { $a->{start} <=> $b->{start} }
                    values %{ $genomic_order{1}{$chr} } )
                {
                    $item->{order} = $i;
                    $i++;
                }
            }
        }
        if ( $ft2 eq "genomic" ) {

            #sort positions
            foreach my $chr ( keys %{ $genomic_order{2} } ) {
                my $i = 1;
                foreach my $item ( sort { $a->{start} <=> $b->{start} }
                    values %{ $genomic_order{2}{$chr} } )
                {
                    $item->{order} = $i;
                    $i++;
                }
            }
        }
    }

    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    open( OUT, ">$outfile" );
    open( IN,  $infile );
    while (<IN>) {
        chomp;
        if (/^#/) {
            print OUT $_, "\n";
            next;
        }
        my @line = split /\t/;

        if ( $ft1 eq "genomic" ) {
            $line[2] = $genomic_order{1}{ $line[0] }{ $line[1] }{order};
            $line[3] = $genomic_order{1}{ $line[0] }{ $line[1] }{order};
        }
        else {
            my @item1 = split /\|\|/, $line[1];
            $line[2] = $item1[7];
            $line[3] = $item1[7];
        }
        if ( $ft2 eq "genomic" ) {
            $line[6] = $genomic_order{2}{ $line[4] }{ $line[5] }{order};
            $line[7] = $genomic_order{2}{ $line[4] }{ $line[5] }{order};
        }
        else {
            my @item2 = split /\|\|/, $line[5];
            $line[6] = $item2[7];
            $line[7] = $item2[7];
        }
        print OUT join( "\t", @line ), "\n";
    }
    close IN;
    close OUT;

    CoGe::Accessory::Web::write_log(
        "running coversion to gene order for $infile",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log(
        "Completed conversion of gene order to file $outfile",
        $cogeweb->logfile );
    system "/bin/rm $outfile.running" if -r "$outfile.running";
    ;    #remove track filereturn $outfile;
    return $outfile;
}

sub replace_gene_order_with_genomic_positions {
    my %opts    = @_;
    my $file    = $opts{file};
    my $outfile = $opts{outfile};

    #must convert file's coordinates back to genomic
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if (   ( -r "$outfile" && -s "$outfile" )
        || ( -r "$outfile.gz" && -s "$outfile.gz" ) )
    {
        CoGe::Accessory::Web::write_log(
"  no conversion for $file back to genomic coordinates needed, convered file, $outfile,  exists",
            $cogeweb->logfile
        );
        return;
    }
    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log(
        "  converting $file back to genomic coordinates, $outfile",
        $cogeweb->logfile );

    #    `mv $file $file.orig`;
    $/ = "\n";    #just in case
    CoGe::Accessory::Web::gunzip( $file . ".gz" ) if -r $file . ".gz";
    open( IN,  "$file" );
    open( OUT, ">$outfile" );
    while (<IN>) {
        if (/^#/) { print OUT $_; next; }
        my @line  = split /\t/;
        my @item1 = split /\|\|/, $line[1];
        my @item2 = split /\|\|/, $line[5];
        my ( $start, $stop ) = ( $item1[1], $item1[2] );
        ( $start, $stop ) = ( $stop, $start ) if $item1[4] && $item1[4] =~ /-/;
        $line[2] = $start;
        $line[3] = $stop;
        ( $start, $stop ) = ( $item2[1], $item2[2] );
        ( $start, $stop ) = ( $stop, $start ) if $item2[4] && $item2[4] =~ /-/;
        $line[6] = $start;
        $line[7] = $stop;
        print OUT join( "\t", @line );
    }
    close IN;
    close OUT;
    system "/bin/rm $outfile.running" if -r "$outfile.running";
}

sub run_dagchainer {
    my %opts   = @_;
    my $infile = $opts{infile};
    my $D      = $opts{D};        #maximum distance allowed between two matches
    my $g =
      $opts{g}; #length of a gap (average distance expected between two syntenic genes)
    $g = floor( $D / 2 ) unless defined $g;
    my $A = $opts{A};    #Minium number of Aligned Pairs
    my $Dm =
      $opts{Dm}; #maximum distance between sytnenic blocks for merging syntenic blocks
    my $gm =
      $opts{gm}; #average distance between sytnenic blocks for merging syntenic blocks
    my $merge = $opts{merge};    #flag to use the merge function of this algo;
    my $self_self = $opts{self_self};    #flag that it is a self-self run

    unless ($merge)    #turn off merging options unless $merge is set to true
    {
        $Dm = 0;
        $gm = 0;
    }

    my $outfile = $infile;
    $outfile .= "_D$D"   if $D;
    $outfile .= "_g$g"   if $g;
    $outfile .= "_A$A"   if $A;
    $outfile .= "_Dm$Dm" if $Dm;
    $outfile .= "_gm$gm" if $gm;
    $outfile .= ".aligncoords";
    $outfile .= ".ma2.dag" if $Dm && $gm;
    my $running_file = $outfile;
    $running_file .= ".dag.merge" if $Dm && $gm;
    $running_file .= ".running";
    my $merged_file = $outfile . ".merge";
    my $return_file = $outfile;
    $return_file .= ".merge"
      if $Dm && $gm;    #created by $RUN_DAGCHAINER if merge option is true

    while ( -e "$running_file" ) {
        print STDERR "detecting $outfile.merge.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $return_file || -r $return_file . ".gz" ) {
        CoGe::Accessory::Web::write_log(
            "run dagchainer: file $return_file already exists",
            $cogeweb->logfile );
        return ( $outfile, $merged_file );
    }
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";

    unless ( -r $infile && -s $infile ) {
        CoGe::Accessory::Web::write_log(
"WARNING:   Cannot run DAGChainer! DAGChainer input file ($infile) contains no data!",
            $cogeweb->logfile
        );
        return 0;
    }
    my $cmd = "$PYTHON $RUN_DAGCHAINER -E 0.05 -i $infile";
    $cmd .= " -D $D"          if defined $D;
    $cmd .= " -g $g"          if defined $g;
    $cmd .= " -A $A"          if defined $A;
    $cmd .= " --Dm $Dm"       if $Dm;
    $cmd .= " --gm $gm"       if $gm;
    $cmd .= " --new_behavior" if $self_self;

    if ( $Dm && $gm ) {
        $cmd .= " --merge $outfile";
    }
    else {
        $cmd .= " > $outfile";
    }
    ###MERGING OF DIAGONALS FUNCTION
# --merge $outfile     #this will automatically cause the merging of diagonals to happen.  $outfile will be created and $outfile.merge (which is the inappropriately named merge of diagonals file.
# --gm  average distance between diagonals
# --Dm  max distance between diagonals
    ##Both of these parameters' default values is 4x -g and -D respectively.

    #

    system "/usr/bin/touch $running_file"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "run dagchainer: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;

    #    `mv $infile.aligncoords $outfile`;
    system "/bin/rm $running_file" if -r "$running_file";
    ;      #remove track file
    return ( $outfile, $merged_file );
}

sub run_quota_align_merge {
    my %opts     = @_;
    my $infile   = $opts{infile};
    my $max_dist = $opts{max_dist};
    my $returnfile =
      $infile . ".Dm" . $max_dist . ".ma1";    #ma stands for merge algo
    return $returnfile if -r $returnfile || -f $returnfile . ".gz";

    #convert to quota-align format
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    my $cmd = $CLUSTER_UTILS
      . " --format=dag --log_evalue $infile $infile.Dm$max_dist.qa";
    `$cmd`;
    CoGe::Accessory::Web::write_log(
        "Converting dag output to quota_align format: $cmd",
        $cogeweb->logfile );
    $cmd = $QUOTA_ALIGN . " --Dm=$max_dist --merge $infile.Dm$max_dist.qa";
    CoGe::Accessory::Web::write_log(
        "Running quota_align to merge diagonals:\n\t$cmd",
        $cogeweb->logfile );
    CoGe::Accessory::Web::gunzip("$infile.Dm$max_dist.qa.gz")
      if -r "$infile.Dm$max_dist.qa.gz";
    `$cmd`;
    CoGe::Accessory::Web::gunzip("$infile.Dm$max_dist.qa.merged.gz")
      if -r "$infile.Dm$max_dist.qa.merged.gz";

    if ( -r "$infile.Dm$max_dist.qa.merged" ) {
        my %data;
        $/ = "\n";
        open( IN, $infile );
        while (<IN>) {
            next if /^#/;
            my @line = split /\t/;
            $data{ join( "_", $line[0], $line[2], $line[4], $line[6] ) } = $_;
        }
        close IN;
        open( OUT, ">$returnfile" );
        open( IN,  "$infile.Dm$max_dist.qa.merged" );
        while (<IN>) {
            if (/^#/) {
                print OUT $_;
            }
            else {
                chomp;
                my @line = split /\t/;
                print OUT $data{
                    join( "_", $line[0], $line[1], $line[2], $line[3] ) };
            }
        }
        close IN;
        close OUT;
    }
    return "$returnfile";
}

sub run_quota_align_coverage {
    my %opts         = @_;
    my $infile       = $opts{infile};
    my $org1         = $opts{org1};           #ratio of org1
    my $org2         = $opts{org2};           #ratio of org2
    my $overlap_dist = $opts{overlap_dist};
    my $returnfile =
        $infile . ".qac"
      . $org1 . "."
      . $org2 . "."
      . $overlap_dist;                        #ma stands for merge algo
    return $returnfile if -r $returnfile || -r $returnfile . ".gz";
    CoGe::Accessory::Web::gunzip( $infile . ".gz" )    if -r $infile . ".gz";
    CoGe::Accessory::Web::gunzip( $infile . ".qa.gz" ) if -r $infile . ".qa.gz";

    #convert to quota-align format
    my $cov_cmd =
      $CLUSTER_UTILS . " --format=dag --log_evalue $infile $infile.qa";
    my $qa_cmd = $QUOTA_ALIGN
      . " --Nm=$overlap_dist --quota=$org1:$org2 $infile.qa > $returnfile.tmp";

    CoGe::Accessory::Web::write_log( "Convert command: $cov_cmd",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Quota Align command: $qa_cmd",
        $cogeweb->logfile );
    while ( -e "$infile.qa.running" ) {
        print STDERR "detecting $infile.qa.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r "$infile.qa" ) {
        CoGe::Accessory::Web::write_log(
"Dag output file already converted to quota_align input: $infile.qa",
            $cogeweb->logfile
        );
    }
    else {
        system "/usr/bin/touch $infile.qa.running";  #track that this is running
        CoGe::Accessory::Web::write_log(
            "Converting dag output to quota_align format.",
            $cogeweb->logfile );
        `$cov_cmd`;
        system "/bin/rm $infile.qa.running"
          if -r "$infile.qa.running";                #remove track file
    }
    while ( -e "$returnfile.running" ) {
        print STDERR "detecting $returnfile.running.  Waiting. . .\n";
        sleep 60;
    }

    if ( -r "$returnfile.tmp" || -r "$returnfile.tmp.gz" ) {
        CoGe::Accessory::Web::write_log(
"Quota_align syntenic coverage parameters already run: $returnfile.tmp"
        );
        CoGe::Accessory::Web::gunzip("$returnfile.tmp.gz")
          if -r "$returnfile.tmp.gz";

    }
    else {
        system "/usr/bin/touch $returnfile.running"; #track that this is running
        CoGe::Accessory::Web::write_log(
            "Running quota_align to find syntenic coverage.",
            $cogeweb->logfile );
        my $qa_output = `$qa_cmd`;
        system "/bin/rm $returnfile.running"
          if -r "$returnfile.running";               #remove track file
    }
    if ( -r "$returnfile.tmp" ) {
        my %data;
        $/ = "\n";
        open( IN, $infile );
        while (<IN>) {
            next if /^#/;
            my @line = split /\t/;
            $data{ join( "_", $line[0], $line[2], $line[4], $line[6] ) } = $_;
        }
        close IN;
        open( OUT, ">$returnfile" );
        open( IN,  "$returnfile.tmp" );
        while (<IN>) {
            if (/^#/) {
                print OUT $_;
            }
            else {
                chomp;
                my @line = split /\t/;
                print OUT $data{
                    join( "_", $line[0], $line[1], $line[2], $line[3] ) };
            }
        }
        close IN;
        close OUT;
    }
    return "$returnfile";
}

sub generate_grimm_input {
    my %opts   = @_;
    my $infile = $opts{infile};
    CoGe::Accessory::Web::gunzip( $infile . ".gz" )    if -r $infile . ".gz";
    CoGe::Accessory::Web::gunzip( $infile . ".qa.gz" ) if -r $infile . ".qa.gz";
    my $cmd = $CLUSTER_UTILS . " --format=dag --log_evalue $infile $infile.qa";
    CoGe::Accessory::Web::write_log( "\nGenerating input data for GRIMM",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log(
        "Converting dag output to quota_align format: $cmd",
        $cogeweb->logfile );
    `$cmd`;
    $cmd = $CLUSTER_UTILS . " --print_grimm $infile.qa";
    CoGe::Accessory::Web::write_log(
        "running  cluster_utils to generating grimm input:\n\t$cmd",
        $cogeweb->logfile );
    my $output;
    open( IN, "$cmd |" );

    while (<IN>) {
        $output .= $_;
    }
    close IN;
    my @seqs;
    foreach my $item ( split /\n>/, $output ) {
        $item =~ s/>//g;
        my ( $name, $seq ) = split /\n/, $item, 2;
        $seq =~ s/\n$//;
        push @seqs, $seq;
    }

    return \@seqs;
}

sub run_find_nearby {
    my %opts         = @_;
    my $infile       = $opts{infile};
    my $dag_all_file = $opts{dag_all_file};
    my $outfile      = $opts{outfile};
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r $outfile ) {
        CoGe::Accessory::Web::write_log(
            "run find_nearby: file $outfile already exists",
            $cogeweb->logfile );
        return 1;
    }
    my $cmd =
      "$PYTHON $FIND_NEARBY --diags=$infile --all=$dag_all_file > $outfile";
    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "run find_nearby: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm $outfile.running" if -r "$outfile.running";
    ;      #remove track file
    return 1 if -r $outfile;
}

sub gen_ks_db {
    my %opts   = @_;
    my $infile = $opts{infile};
    my ($outfile) = $infile =~ /^(.*?CDS-CDS)/;
    ($outfile) = $infile =~ /^(.*?protein-protein)/ unless $outfile;
    return unless $outfile;
    $outfile .= ".sqlite";
    CoGe::Accessory::Web::write_log( "Generating ks data.", $cogeweb->logfile );
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }

    unless ( -r $outfile ) {
        system "/usr/bin/touch $outfile.running";    #track that this is running
        CoGe::Accessory::Web::write_log( "initializing ks database $outfile",
            $cogeweb->logfile );
        my $create = qq{
CREATE TABLE ks_data
(
id INTEGER PRIMARY KEY,
fid1 integer,
fid2 integer,
dS varchar,
dN varchar,
dN_dS varchar,
protein_align_1,
protein_align_2,
DNA_align_1,
DNA_align_2
)
};

        my $dbh = DBI->connect( "dbi:SQLite:dbname=$outfile", "", "" );
        $dbh->do($create) if $create;
        $dbh->do('create INDEX fid1 ON ks_data (fid1)');
        $dbh->do('create INDEX fid2 ON ks_data (fid2)');
        $dbh->disconnect;
    }
    my $ksdata = get_ks_data( db_file => $outfile );
    $/ = "\n";
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    open( IN, $infile );
    my @data;
    while (<IN>) {

        next if /^#/;
        chomp;
        my @line  = split /\t/;
        my @item1 = split /\|\|/, $line[1];
        my @item2 = split /\|\|/, $line[5];
        unless ( $item1[6] && $item2[6] ) {
            warn "Line does not appear to contain coge feature ids:  $_\n";
            next;
        }
        next unless $item1[5] eq "CDS" && $item2[5] eq "CDS";
        next if $ksdata->{ $item1[6] }{ $item2[6] };
        push @data, [ $line[1], $line[5], $item1[6], $item2[6] ];
    }
    close IN;
    my $MAX_RUNS = $MAX_PROC;
    my $ports =
      initialize_nwalign_servers( start_port => 3000, procs => $MAX_RUNS );
    my $pm = new Parallel::ForkManager( $MAX_RUNS * 4 );
    my $i  = 0;

    my $start_time = new Benchmark;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$outfile", "", "" );

    $dbh->do("PRAGMA synchronous=OFF");
    $dbh->do("PRAGMA cache_size = 25000");
    $dbh->do("PRAGMA count_changes=OFF");
    $dbh->do("PRAGMA journal_mode=MEMORY");
    $dbh->do("PRAGMA temp_store=MEMORY");

    foreach my $item (@data) {
        $i++;
        $i = 0 if $i == $MAX_RUNS * 4;
        $pm->start and next;
        my ($fid1) = $item->[2] =~ /(^\d+$)/;
        my ($fid2) = $item->[3] =~ /(^\d+$)/;
        my ($feat1) = $coge->resultset('Feature')->find($fid1);
        my ($feat2) = $coge->resultset('Feature')->find($fid2);
        my $max_res;
        my $ks = new CoGe::Algos::KsCalc(config=>$CONFIG);
        $ks->nwalign_server_port( $ports->[$i] );
        $ks->feat1($feat1);
        $ks->feat2($feat2);

        #	for (1..5)
        #	  {
        my $res = $ks->KsCalc(config=>$CONFIG);    #send in port number?
        $max_res = $res unless $max_res;
        $max_res = $res
          if $res && $max_res && $res->{dS} && $max_res->{dS} && $res->{dS} < $max_res->{dS};

        #	  }
        unless ($max_res) {
            $max_res = {};
        }
        my ( $dS, $dN, $dNS ) = ( "NA", "NA", "NA" );
        if ( keys %$max_res ) {
            $dS  = $max_res->{dS}      if defined $max_res->{dS};
            $dN  = $max_res->{dN}      if defined $max_res->{dN};
            $dNS = $max_res->{'dN/dS'} if defined $max_res->{'dN/dS'};
        }

        #alignments
        my $dalign1 = $ks->dalign1;    #DNA sequence 1
        my $dalign2 = $ks->dalign2;    #DNA sequence 2
        my $palign1 = $ks->palign1;    #PROTEIN sequence 1
        my $palign2 = $ks->palign2;    #PROTEIN sequence 2

        my $insert = qq{
INSERT INTO ks_data (fid1, fid2, dS, dN, dN_dS, protein_align_1, protein_align_2, DNA_align_1, DNA_align_2) values ($fid1, $fid2, "$dS", "$dN", "$dNS", "$palign1", "$palign2", "$dalign1", "$dalign2")
};

        my $insert_success = 0;
	my $tries = 0;
        while ( !$insert_success ) {
            $insert_success = $dbh->do($insert);
            unless ($insert_success) {
	      $tries++;
                #				print STDERR $insert;
                sleep .1;
            }
	    if ($tries > 10)
	      {
		print STDERR "Tried sqlite insert $tries times without success.  Failing and existing.";
		last;
	      }
        }

        $pm->finish;
    }
    $pm->wait_all_children();
    $dbh->disconnect();
    my $finished_time = new Benchmark;
    my $completion_time = timestr( timediff( $finished_time, $start_time ) );
    say STDERR "Completed in: $completion_time";

    system "/bin/rm $outfile.running" if -r "$outfile.running";
    ;    #remove track file
    CoGe::Accessory::Web::write_log( "Completed generating ks data.",
        $cogeweb->logfile );
    return $outfile;
}

sub initialize_nwalign_servers {
    my %opts       = @_;
    my $start_port = $opts{start_port};
    my $procs      = $opts{procs};
    my @ports;
    use IO::Socket;

    for ( 1 .. $procs ) {
        unless (
            IO::Socket::INET->new(
                PeerAddr => 'localhost',
                PeerPort => $start_port
            )
          )
        {
            system("$NWALIGN --server $start_port &");
        }

        push @ports, $start_port;
        $start_port++;
    }
    sleep 1;    #let them get started;
    return \@ports;
}

sub get_ks_data {
    my %opts    = @_;
    my $db_file = $opts{db_file};
    my %ksdata;
    unless ( -r $db_file ) {
        return \%ksdata;
    }
    CoGe::Accessory::Web::write_log( "\tconnecting to ks database $db_file",
        $cogeweb->logfile );
    my $select = "select * from ks_data";
    my $dbh    = DBI->connect( "dbi:SQLite:dbname=$db_file", "", "" );
    my $sth    = $dbh->prepare($select);
    $sth->execute();
    CoGe::Accessory::Web::write_log(
        "\texecuting select all from ks database $db_file",
        $cogeweb->logfile );
    my $total   = 0;
    my $no_data = 0;

    while ( my $data = $sth->fetchrow_arrayref ) {
        $total++;
        if ( $data->[3] eq "" ) {
            $no_data++;

            #	    next; #uncomment to force recalculation of missing data
        }

        $ksdata{ $data->[1] }{ $data->[2] } =
          $data->[3]
          ? {
            dS      => $data->[3],
            dN      => $data->[4],
            'dN/dS' => $data->[5]
          }
          : {};    # unless $data->[3] eq "";
    }
    CoGe::Accessory::Web::write_log(
        "\tgathered data from ks database $db_file",
        $cogeweb->logfile );
    $sth->finish;
    undef $sth;
    $dbh->disconnect();
    CoGe::Accessory::Web::write_log(
        "\tdisconnecting from ks database $db_file",
        $cogeweb->logfile );
    return \%ksdata;
}

sub add_GEvo_links {
    my %opts   = @_;
    my $infile = $opts{infile};
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    $/ = "\n";
    return
      if ( -r "$infile.condensed" || -r "$infile.condensed.gz" );    #check this
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    open( IN,  $infile );
    open( OUT, ">$infile.tmp" );
    my %condensed;
    my %names;
    my $previously_generated = 0;

    while (<IN>) {
        chomp;
        if (/^#/) {
            print OUT $_, "\n";
            next;
        }
        if (/GEvo/) {
            $previously_generated = 1;
        }
        s/^\s+//;
        next unless $_;
        my @line  = split /\t/;
        my @feat1 = split /\|\|/, $line[1];
        my @feat2 = split /\|\|/, $line[5];
        my $link  = $BASE_URL . "GEvo.pl?";
        my ( $fid1, $fid2 );

        if ( $feat1[6] ) {
            $fid1 = $feat1[6];
            $link .= "fid1=" . $fid1;
        }
        else {
            my ($xmin) = sort ( $feat1[1], $feat1[2] );
            my $x = sprintf( "%.0f", $xmin + abs( $feat1[1] - $feat1[2] ) / 2 );
            $link .= "chr1=" . $feat1[0] . ";x1=" . $x;
        }
        if ( $feat2[6] ) {
            $fid2 = $feat2[6];
            $link .= ";fid2=" . $fid2;
        }
        else {
            my ($xmin) = sort ( $feat2[1], $feat2[2] );
            my $x = sprintf( "%.0f", $xmin + abs( $feat2[1] - $feat2[2] ) / 2 );
            $link .= ";chr2=" . $feat2[0] . ";x2=" . $x;
        }
        $link .= ";dsgid1=" . $dsgid1;
        $link .= ";dsgid2=" . $dsgid2;

        if ( $fid1 && $fid2 ) {
            $condensed{ $fid1 . "_" . $dsgid1 }{ $fid2 . "_" . $dsgid2 } = 1;
            $condensed{ $fid2 . "_" . $dsgid2 }{ $fid1 . "_" . $dsgid1 } = 1;
            $names{$fid1} = $feat1[3];
            $names{$fid2} = $feat2[3];
        }

#	accn1=".$feat1[3]."&fid1=".$feat1[6]."&accn2=".$feat2[3]."&fid2=".$feat2[6] if $feat1[3] && $feat1[6] && $feat2[3] && $feat2[6];
        print OUT $_;
        print OUT "\t", $link;
        print OUT "\n";
    }
    close IN;
    close OUT;
    if ($previously_generated) {
        `/bin/rm $infile.tmp` if -r "$infile.tmp";
    }
    else {
        my $cmd = "/bin/mv $infile.tmp $infile";
        `$cmd`;
    }
    if ( keys %condensed
        && !( -r "$infile.condensed" || -r "$infile.condensed.gz" ) )
    {

        #take into account transitivity
        foreach my $id1 ( keys %condensed ) {
            foreach my $id2 ( keys %{ $condensed{$id1} } ) {
                foreach my $id3 ( keys %{ $condensed{$id2} } ) {
                    next if $id1 eq $id2;
                    $condensed{$id1}{$id3} = 1;
                    $condensed{$id3}{$id1} = 1;
                }
            }
        }

        open( OUT, ">$infile.condensed" );
        print OUT join( "\t",
            qw(COUNT GEVO MASKED_GEVO FASTA_LINK GENE_LIST GENE_NAMES) ),
          "\n";
        my %seen;
        foreach my $id1 (
            sort {
                scalar( keys %{ $condensed{$b} } ) <=>
                  scalar( keys %{ $condensed{$a} } )
            } keys %condensed
          )
        {
            my ( $fid1, $dsgid1 ) = split /_/, $id1;
            next if $seen{$fid1};
            $seen{$fid1} = 1;
            my @names     = $names{$fid1};
            my $gevo_link = $BASE_URL . "GEvo.pl?fid1=$fid1;dsgid1=$dsgid1";
            my $fids      = "fid=$fid1";
            my $count     = 2;

            foreach my $id2 ( sort keys %{ $condensed{$id1} } ) {
                my ( $fid2, $dsgid2 ) = split /_/, $id2, 2;
                next if $fid1 == $fid2;
                $seen{$fid2} = 1;
                $gevo_link .= ";fid$count=$fid2;dsgid$count=$dsgid2";
                $fids      .= ",$fid2";
                push @names, $names{$fid2};
                $count++;
            }
            $count--;
            $gevo_link .= ";num_seqs=$count";
            my $gevo_link2 = $gevo_link;
            $gevo_link .= ";pad_gs=20000";
            for my $i ( 1 .. $count ) {
                $gevo_link2 .= ";mask$i=non-cds";
            }
            $gevo_link2 .= ";pad_gs=200000";
            $gevo_link2 .= ";autogo=1";
            my $fasta_link    = $BASE_URL . "FastaView.pl?$fids";
            my $featlist_link = $BASE_URL . "FeatList.pl?$fids";
            print OUT join( "\t",
                $count, $gevo_link, $gevo_link2, $fasta_link, $featlist_link,
                @names ),
              "\n";
        }
        close OUT;
    }
}

sub gen_ks_blocks_file {
    my %opts   = @_;
    my $infile = $opts{infile};
    CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
    my ($dbfile) = $infile =~ /^(.*?CDS-CDS)/;
    return unless $dbfile;
    my $outfile = $infile . ".ks";
    return $outfile if -r $outfile;
    my $ksdata = get_ks_data( db_file => $dbfile . '.sqlite' );
    $/ = "\n";
    open( IN,  $infile );
    open( OUT, ">" . $outfile );
    print OUT
      "#This file contains synonymous rate values in the first two columns:\n";
    my @block;
    my $block_title;

    while (<IN>) {
        if (/^#/) {
            my $output = process_block(
                ksdata => $ksdata,
                block  => \@block,
                header => $block_title
            ) if $block_title;
            print OUT $output if $output;
            @block       = ();
            $block_title = $_;

            #beginning of a block;
        }
        else {
            push @block, $_;
        }
    }
    close IN;
    my $output = process_block(
        ksdata => $ksdata,
        block  => \@block,
        header => $block_title
    ) if $block_title;
    print OUT $output;
    close OUT;
    return $outfile;
}

sub gen_svg_file {
    my %opts   = @_;
    my $infile = $opts{infile};
    CoGe::Accessory::Web::gunzip($infile);
    my $outfile = $opts{outfile};
    my $flip    = $opts{flip};
    my $xname   = $opts{xname};
    my $yname   = $opts{yname};

    $outfile = $infile unless $outfile;

#/opt/apache/CoGe/bin/SynMap/dotplot.py --dag_file 3068_11022.CDS-CDS.lastz.dag.go_c4_D20_g10_A5.aligncoords.qac1.1.40.gcoords.ks --output my
#    if (-r $infile.".svg")
    {
        CoGe::Accessory::Web::write_log(
            "generate svg dotplot: file exists: $infile.svg",
            $cogeweb->logfile );
    }

    #   else
    {
        my $cmd = $SVG_DOTPLOT . " --dag_file $infile";
        $cmd .= " --flip"           if $flip;
        $cmd .= " --xhead '$xname'" if $xname;
        $cmd .= " --yhead '$yname'" if $yname;
        $cmd .= " --output $infile";
        CoGe::Accessory::Web::write_log( "generate svg dotplot: $cmd",
            $cogeweb->logfile );
        `$cmd`;
    }
    return $infile . ".svg";
}

sub process_block {
    my %opts   = @_;
    my $ksdata = $opts{ksdata};
    my $block  = $opts{block};
    my $header = $opts{header};
    my $output;
    my @ks;
    my @kn;
    foreach my $item (@$block) {
        my @line = split /\t/,   $item;
        my @seq1 = split /\|\|/, $line[1];
        my @seq2 = split /\|\|/, $line[5];
        my $ks   = $ksdata->{ $seq1[6] }{ $seq2[6] };
        if ( defined $ks->{dS} ) {
            unshift @line, $ks->{dN};
            unshift @line, $ks->{dS};
            push @ks, $ks->{dS} unless $ks->{dS} eq "NA";
            push @kn, $ks->{dN} unless $ks->{dN} eq "NA";
        }
        else {
            unshift @line, "undef";
            unshift @line, "undef";
        }
        $output .= join "\t", @line;
    }
    my $mean_ks = 0;
    if ( scalar @ks ) {
        map { $mean_ks += $_ } @ks;
        $mean_ks = sprintf( "%.4f", $mean_ks / scalar @ks );
    }
    else {
        $mean_ks = "NA";
    }
    my $mean_kn = 0;
    if ( scalar @kn ) {
        map { $mean_kn += $_ } @kn;
        $mean_kn = sprintf( "%.4f", $mean_kn / scalar @kn );
    }
    else {
        $mean_kn = "NA";
    }
    chomp $header;
    $header .= "  Mean Ks:  $mean_ks\tMean Kn: $mean_kn\n";
    $header .= join( "\t",
        "#Ks",
        qw(Kn a<db_genome_id>_<chr> chr1||start1||stop1||name1||strand1||type1||db_feature_id1||percent_id1 start1 stop1 b<db_genome_id>_<chr> chr2||start2||stop2||name2||strand2||type2||db_feature_id2||percent_id2 start2 stop2 eval block_score GEVO_link)
    ) . "\n";
    return $header . $output;
}

sub add_reverse_match #this is for when there is a self-self comparison.  DAGchainer, for some reason, is only adding one diag.  For example, if chr_14 vs chr_10 have a diag, chr_10 vs chr_14 does not.
{
    my %opts   = @_;
    my $infile = $opts{infile};
    $/ = "\n";
    open( IN, $infile );
    my $stuff;
    my $skip = 0;
    while (<IN>) {
        chomp;
        s/^\s+//;
        $skip = 1
          if /GEvo\.pl/
        ; #GEvo links have been added, this file was generated on a previous run.  Skip!
        last if ($skip);
        next unless $_;
        my @line = split(/\s+/);
        if (/^#/) {
            my $chr1 = $line[2];
            my $chr2 = $line[4];
            $chr1 =~ s/^a//;
            $chr2 =~ s/^b//;
            next if $chr1 eq $chr2;
            $line[2] = "b" . $chr2;
            $line[4] = "a" . $chr1;
            $stuff .= join( " ", @line ) . "\n";
            next;
        }
        my $chr1 = $line[0];
        my $chr2 = $line[4];
        $chr1 =~ s/^a//;
        $chr2 =~ s/^b//;
        next if $chr1 eq $chr2;
        my @tmp1 = @line[ 1 .. 3 ];
        my @tmp2 = @line[ 5 .. 7 ];
        @line[ 1 .. 3 ] = @tmp2;
        @line[ 5 .. 7 ] = @tmp1;
        $line[0]        = "a" . $chr2;
        $line[4]        = "b" . $chr1;
        $stuff .= join( "\t", @line ) . "\n";
    }
    return if $skip;
    close IN;
    open( OUT, ">>$infile" );
    print OUT $stuff;
    close OUT;

}

sub generate_dotplot {
    my %opts    = @_;
    my $dag     = $opts{dag};
    my $coords  = $opts{coords};
    my $outfile = $opts{outfile};
    my $dsgid1  = $opts{dsgid1};
    my $dsgid2  = $opts{dsgid2};
    my $dagtype = $opts{dagtype};
    my $ks_db   = $opts{ks_db};
    my $ks_type = $opts{ks_type};
    my ($basename) = $coords =~ /([^\/]*aligncoords.*)/;    #.all.aligncoords/;
    my $regen_images = $opts{regen_images} =~ /true/i ? 1 : 0;
    my $width        = $opts{width} || 1000;
    my $assemble     = $opts{assemble};
    my $metric       = $opts{metric};
    my $relationship = $opts{relationship};
    my $min_chr_size = $opts{min_chr_size};
    my $color_type   = $opts{color_type};
    my $just_check =
      $opts{just_check};    #option to just check if the outfile already exists
    my $box_diags = $opts{box_diags};
    my $fid1      = $opts{fid1};     #fids for highlighting gene pair in dotplot
    my $fid2      = $opts{fid2};
    my $snsd =
      $opts{snsd};    #option for showing non syntenic (grey) dots in dotplot
    my $flip         = $opts{flip};            #flip axis
    my $clabel       = $opts{clabel};          #label chromosomes
    my $skip_rand    = $opts{skip_rand};
    my $color_scheme = $opts{color_scheme};    #color_scheme for dotplot
    my $chr_sort_order =
      $opts{chr_sort_order};    #sort order of chromosomes in dotplot
    my $codeml_min = $opts{codeml_min};
    my $codeml_max = $opts{codeml_max};
    my $log        = $opts{log};
    my $cmd        = $DOTPLOT;

    #add ks_db to dotplot command if requested
    $outfile .= ".spa$assemble"     if $assemble;
    $outfile .= ".gene"             if $metric =~ /gene/i;
    $outfile .= ".s"                if $relationship =~ /s/i;
    $outfile .= ".mcs$min_chr_size" if $min_chr_size;
    $outfile .= ".$fid1"            if $fid1;
    $outfile .= ".$fid2"            if $fid2;
    if ( $ks_db && -r $ks_db ) {
        $cmd     .= qq{ -ksdb $ks_db -kst $ks_type -log $log};
        $outfile .= ".$ks_type";
    }
    $outfile .= ".box"                if $box_diags;
    $outfile .= ".flip"               if $flip;
    $outfile .= ".c0"                 if $clabel eq 0;
    $outfile .= ".sr"                 if $skip_rand;
    $outfile .= ".cs$color_scheme"    if defined $color_scheme;
    $outfile .= ".cso$chr_sort_order" if defined $chr_sort_order;
    $outfile .= ".min$codeml_min"     if defined $codeml_min;
    $outfile .= ".max$codeml_max"     if defined $codeml_max;
    $outfile .= ".log"                if $log;

    #are non-syntenic dots being displayed
    if ($snsd) {
        $cmd .= qq{ -d $dag};
    }
    else {
        $outfile .= ".nsd";    #no syntenic dots, yes, nomicalture is confusing.
    }

    return $outfile if $just_check && -r "$outfile.html";

    CoGe::Accessory::Web::gunzip( $coords . ".gz" ) if -r "$coords.gz";
    CoGe::Accessory::Web::gunzip( $dag . ".gz" )    if -r "$dag.gz";
    $cmd .= qq{ -a $coords};
    $cmd .=
qq{ -b $outfile -l 'javascript:synteny_zoom("$dsgid1","$dsgid2","$basename",};
    $cmd .= $flip ? qq{"YCHR","XCHR"} : qq{"XCHR","YCHR"};
    $cmd .= qq{,"$ks_db"}            if $ks_db;
    $cmd .= qq{)' -dsg1 $dsgid1 -dsg2 $dsgid2 -w $width -lt 2};
    $cmd .= qq{ -assemble $assemble} if $assemble;
    $cmd .= qq{ -am $metric}         if $metric;
    $cmd .= qq{ -fb} if $relationship && $relationship =~ /s/;
    $cmd .= qq{ -mcs $min_chr_size}              if $min_chr_size;
    $cmd .= qq{ -cdt $color_type}                if $color_type;
    $cmd .= qq{ -bd 1}                           if $box_diags;
    $cmd .= qq{ -fid1 $fid1}                     if $fid1;
    $cmd .= qq{ -fid2 $fid2}                     if $fid2;
    $cmd .= qq{ -f 1}                            if $flip;
    $cmd .= qq{ -labels 0}                       if $clabel eq 0;
    $cmd .= qq{ -sr 1}                           if $skip_rand;
    $cmd .= qq{ -color_scheme $color_scheme}     if defined $color_scheme;
    $cmd .= qq{ -chr_sort_order $chr_sort_order} if defined $chr_sort_order;
    $cmd .= qq{ -min $codeml_min}                if defined $codeml_min;
    $cmd .= qq{ -max $codeml_max}                if defined $codeml_max;

#    $cmd .= qq{ -cf }.$ENV{COGE_HOME}. 'coge.conf'; #config file for getting defaults for coge server installation
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if ( -r "$outfile.png" && !$regen_images ) {
        CoGe::Accessory::Web::write_log(
            "generate dotplot: file $outfile already exists",
            $cogeweb->logfile );
        return $outfile;
    }
    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "generate dotplot: running\n\t$cmd",
        $cogeweb->logfile );
    system "/bin/rm $outfile.running" if -r "$outfile.running";
    ;      #remove track file
    `$cmd`;
    return ($outfile) if -r "$outfile.html";
}

sub get_query_link {
    my %url_options = @_;

    #  print STDERR Dumper \%url_options;

    my $dagchainer_D = $url_options{D};

 #    my $dagchainer_g = $url_options{g}; #depreciated -- will be a factor of -D
    my $dagchainer_A = $url_options{A};
    my $Dm           = $url_options{Dm};
    my $gm           = $url_options{gm};
    ($Dm) = $Dm =~ /(\d+)/;
    ($gm) = $gm =~ /(\d+)/;

#	my $repeat_filter_cvalue = $url_options{c}; #parameter to be passed to run_adjust_dagchainer_evals
    my $cscore =
      $url_options{csco}; #c-score for filtering low quality blast hits, fed to blast to raw
    my $dupdist =
      $url_options{tdd};    #tandem duplication distance, fed to blast to raw
    my $regen_images = $url_options{regen_images};
    my $email        = $url_options{email};
    my $job_title    = $url_options{jobtitle};
    my $width        = $url_options{width};
    my $basename     = $url_options{basename};
    my $blast        = $url_options{blast};

    my $feat_type1 = $url_options{feat_type1};
    my $feat_type2 = $url_options{feat_type2};

    my $dsgid1   = $url_options{dsgid1};
    my $dsgid2   = $url_options{dsgid2};
    my $ks_type  = $url_options{ks_type};
    my $assemble = $url_options{assemble} =~ /true/i ? 1 : 0;
    $assemble = 2 if $assemble && $url_options{show_non_syn} =~ /true/i;
    $assemble *= -1 if $assemble && $url_options{spa_ref_genome} < 0;
    my $axis_metric       = $url_options{axis_metric};
    my $axis_relationship = $url_options{axis_relationship};
    my $min_chr_size      = $url_options{min_chr_size};
    my $dagchainer_type   = $url_options{dagchainer_type};
    my $color_type        = $url_options{color_type};
    my $merge_algo = $url_options{merge_algo};    #is there a merging function?

    #options for finding syntenic depth coverage by quota align (Bao's algo)
    my $depth_algo        = $url_options{depth_algo};
    my $depth_org_1_ratio = $url_options{depth_org_1_ratio};
    my $depth_org_2_ratio = $url_options{depth_org_2_ratio};
    my $depth_overlap     = $url_options{depth_overlap};

    #fids that are passed in for highlighting the pair in the dotplot
    my $fid1 = $url_options{fid1};
    my $fid2 = $url_options{fid2};

    #will non-syntenic dots be shown?
    my $snsd = $url_options{show_non_syn_dots} =~ /true/i ? 1 : 0;
    my $algo_name = $ALGO_LOOKUP->{$blast}{displayname};

    #will the axis be flipped?
    my $flip = $url_options{flip};
    $flip = $flip =~ /true/i ? 1 : 0;

    #are axes labeled?
    my $clabel = $url_options{clabel};
    $clabel = $clabel =~ /true/i ? 1 : 0;

    #are random chr skipped
    my $skip_rand = $url_options{skip_rand};
    $skip_rand = $skip_rand =~ /true/i ? 1 : 0;

    #which color scheme for ks/kn dots?
    my $color_scheme = $url_options{color_scheme};

    #codeml min and max calues
    my $codeml_min = $url_options{codeml_min};
    $codeml_min = undef
      unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
    my $codeml_max = $url_options{codeml_max};
    $codeml_max = undef
      unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
    my $logks = $url_options{logks};
    $logks = $logks eq "true" ? 1 : 0;

    #how are the chromosomes to be sorted?
    my $chr_sort_order = $url_options{chr_sort_order};

    #draw a box around identified diagonals?
    my $box_diags = $url_options{box_diags};
    $box_diags = $box_diags eq "true" ? 1 : 0;
    my $synmap_link =
        $SERVER
      . "SynMap.pl?dsgid1=$dsgid1;dsgid2=$dsgid2"
      . ";D=$dagchainer_D;A=$dagchainer_A;w=$width;b=$blast;ft1=$feat_type1;"
      . "ft2=$feat_type2;regen_images=1;autogo=1";
    $synmap_link .= ";fid1=$fid1"   if $fid1;
    $synmap_link .= ";fid2=$fid2"   if $fid2;
    $synmap_link .= ";Dm=$Dm"       if defined $Dm;
    $synmap_link .= ";csco=$cscore" if $cscore;
    $synmap_link .= ";tdd=$dupdist" if defined $dupdist;
    $synmap_link .= ";gm=$gm"       if defined $gm;
    $synmap_link .= ";snsd=$snsd";

    $synmap_link .= ";bd=$box_diags"          if $box_diags;
    $synmap_link .= ";mcs=$min_chr_size"      if $min_chr_size;
    $synmap_link .= ";sp=$assemble"           if $assemble;
    $synmap_link .= ";ma=$merge_algo"         if $merge_algo;
    $synmap_link .= ";da=$depth_algo"         if $depth_algo;
    $synmap_link .= ";do1=$depth_org_1_ratio" if $depth_org_1_ratio;
    $synmap_link .= ";do2=$depth_org_2_ratio" if $depth_org_2_ratio;
    $synmap_link .= ";do=$depth_overlap"      if $depth_overlap;
    $synmap_link .= ";flip=1"                 if $flip;
    $synmap_link .= ";cs=$color_scheme";
    $synmap_link .= ";cmin=$codeml_min"
      if defined $codeml_min;   #$codeml_min=~/\d/ && $codeml_min=~/^\d*.?\d*$/;
    $synmap_link .= ";cmax=$codeml_max"
      if defined $codeml_max;   #$codeml_max=~/\d/ && $codeml_max=~/^\d*.?\d*$/;
    $synmap_link .= ";logks=$logks"        if defined $logks;
    $synmap_link .= ";cl=0"                if $clabel eq "0";
    $synmap_link .= ";sr=$skip_rand"       if defined $skip_rand;
    $synmap_link .= ";cso=$chr_sort_order" if $chr_sort_order;

    $email = 0 if check_address_validity($email) eq 'invalid';

    $feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
    $feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
    $feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
    $feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time

    $synmap_link .= ";dt=$dagchainer_type";

    if ($ks_type) {
        my $num;
        if    ( $ks_type eq "ks" )    { $num = 1; }
        elsif ( $ks_type eq "kn" )    { $num = 2; }
        elsif ( $ks_type eq "kn_ks" ) { $num = 3; }
        $synmap_link .= ";ks=$num";
    }
    $synmap_link .= ";am=g" if $axis_metric       && $axis_metric       =~ /g/i;
    $synmap_link .= ";ar=s" if $axis_relationship && $axis_relationship =~ /s/i;
    $synmap_link .= ";ct=$color_type" if $color_type;

    my ( $org_name1, $titleA ) = gen_org_name(
        dsgid     => $dsgid1,
        feat_type => $feat_type1,
        write_log => 0
    );
    my ( $org_name2, $titleB ) = gen_org_name(
        dsgid     => $dsgid2,
        feat_type => $feat_type2,
        write_log => 0
    );
    my $log_msg =
        "<span class=link onclick=window.open('OrganismView.pl?dsgid="
      . $dsgid1
      . "')>$org_name1</span> v. <span class=link"
      . "onclick=window.open('OrganismView.pl?dsgid=$dsgid2')>$org_name2</span>";

    $log_msg .= " Ks" if $ks_type;
    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    return CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $synmap_link,
        log_msg => $log_msg
    );
}

sub go {
    my %opts = @_;
    foreach my $k ( keys %opts ) {
        $opts{$k} =~ s/^\s+//;
        $opts{$k} =~ s/\s+$//;
    }
    my $dagchainer_D = $opts{D};
    my $tiny_link    = $opts{tiny_link};

    say STDERR "tiny_link is required for logging." unless defined($tiny_link);

    #print STDERR "SUB GO:";
    #    print STDERR Dumper \%opts;

    #    my $dagchainer_g = $opts{g}; #depreciated -- will be a factor of -D
    my $dagchainer_A = $opts{A};
    my $Dm           = $opts{Dm};
    my $gm           = $opts{gm};
    ($Dm) = $Dm =~ /(\d+)/;
    ($gm) = $gm =~ /(\d+)/;

#	my $repeat_filter_cvalue = $opts{c};              #parameter to be passed to run_adjust_dagchainer_evals
    my $cscore =
      $opts{csco}; #c-score for filtering low quality blast hits, fed to blast to raw
    my $dupdist = $opts{tdd};  #tandem duplication distance, fed to blast to raw
    $dupdist = 10 unless defined($dupdist);

    my $regen_images = $opts{regen_images};
    my $email        = $opts{email};
    my $job_title    = $opts{jobtitle};
    my $width        = $opts{width};
    my $basename     = $opts{basename};
    my $blast        = $opts{blast};

    my $feat_type1 = $opts{feat_type1};
    my $feat_type2 = $opts{feat_type2};

    my $dsgid1   = $opts{dsgid1};
    my $dsgid2   = $opts{dsgid2};
    my $ks_type  = $opts{ks_type};
    my $assemble = $opts{assemble} =~ /true/i ? 1 : 0;
    $assemble = 2 if $assemble && $opts{show_non_syn} =~ /true/i;
    $assemble *= -1 if $assemble && $opts{spa_ref_genome} < 0;
    my $axis_metric       = $opts{axis_metric};
    my $axis_relationship = $opts{axis_relationship};
    my $min_chr_size      = $opts{min_chr_size};
    my $dagchainer_type   = $opts{dagchainer_type};
    my $color_type        = $opts{color_type};
    my $merge_algo = $opts{merge_algo};    #is there a merging function?

    #options for finding syntenic depth coverage by quota align (Bao's algo)
    my $depth_algo        = $opts{depth_algo};
    my $depth_org_1_ratio = $opts{depth_org_1_ratio};
    my $depth_org_2_ratio = $opts{depth_org_2_ratio};
    my $depth_overlap     = $opts{depth_overlap};

    #fids that are passed in for highlighting the pair in the dotplot
    my $fid1 = $opts{fid1};
    my $fid2 = $opts{fid2};

    #will non-syntenic dots be shown?
    my $snsd = $opts{show_non_syn_dots} =~ /true/i ? 1 : 0;
    my $algo_name = $ALGO_LOOKUP->{$blast}{displayname};

    #will the axis be flipped?
    my $flip = $opts{flip};
    $flip = $flip =~ /true/i ? 1 : 0;

    #are axes labeled?
    my $clabel = $opts{clabel};
    $clabel = $clabel =~ /true/i ? 1 : 0;

    #are random chr skipped
    my $skip_rand = $opts{skip_rand};
    $skip_rand = $skip_rand =~ /true/i ? 1 : 0;

    #which color scheme for ks/kn dots?
    my $color_scheme = $opts{color_scheme};

    #codeml min and max calues
    my $codeml_min = $opts{codeml_min};
    $codeml_min = undef
      unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
    my $codeml_max = $opts{codeml_max};
    $codeml_max = undef
      unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
    my $logks = $opts{logks};
    $logks = $logks eq "true" ? 1 : 0;

    #draw a box around identified diagonals?
    my $box_diags = $opts{box_diags};
    $box_diags = $box_diags eq "true" ? 1 : 0;

    #how are the chromosomes to be sorted?
    my $chr_sort_order = $opts{chr_sort_order};

    #	$dagchainer_type = $dagchainer_type eq "true" ? "geneorder" : "distance";

    unless ( $dsgid1 && $dsgid2 ) {
        return "<span class=alert>You must select two genomes.</span>";
    }

    my ($dsg1) = $coge->resultset('Genome')->find($dsgid1);
    my ($dsg2) = $coge->resultset('Genome')->find($dsgid2);

    unless ( $dsg1 && $dsg2 ) {
        return
"<span class=alert>Problem generating dataset group objects for ids:  $dsgid1, $dsgid2.</span>";
    }

    unless ( $dsg1 && $dsg2 ) {
        return
"<span class=alert>Problem generating one of the genome objects for id1: $dsgid1 or id2: $dsgid2</span>";
    }

    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );
    $email = 0 if check_address_validity($email) eq 'invalid';

    $feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
    $feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
    $feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
    $feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time

    ##generate fasta files and blastdbs
    my $t0   = new Benchmark;
    my $pm   = new Parallel::ForkManager($MAX_PROC);
    my @dsgs = ( [ $dsgid1, $feat_type1 ] );
    push @dsgs, [ $dsgid2, $feat_type2 ]
      unless $dsgid1 == $dsgid2 && $feat_type1 eq $feat_type2;
    foreach my $item (@dsgs) {
        $pm->start and next;
        my $dsgid = $item->[0];

        my $feat_type = $item->[1];

        my ( $fasta, $org_name ) =
          gen_fasta( dsgid => $dsgid, feat_type => $feat_type, write_log => 1 );

#	gen_blastdb(dbname=>"$dsgid-$feat_type",fasta=>$fasta,org_name=>$org_name, write_log=>1) unless $ALGO_LOOKUP->{$blast}{filename}=~/lastz/; #don't need to do this for lastz
        $pm->finish;
    }
    $pm->wait_all_children();
    my ( $fasta1, $org_name1, $title1 );
    my ( $fasta2, $org_name2, $title2 );
    ( $fasta1, $org_name1, $title1 ) =
      gen_fasta( dsgid => $dsgid1, feat_type => $feat_type1 );
    ( $fasta2, $org_name2, $title2 ) =
      gen_fasta( dsgid => $dsgid2, feat_type => $feat_type2 );
    (
        $dsgid1,     $dsg1,              $org_name1,  $fasta1,
        $feat_type1, $depth_org_1_ratio, $dsgid2,     $dsg2,
        $org_name2,  $fasta2,            $feat_type2, $depth_org_2_ratio
      )
      = (
        $dsgid2,     $dsg2,              $org_name2,  $fasta2,
        $feat_type2, $depth_org_2_ratio, $dsgid1,     $dsg1,
        $org_name1,  $fasta1,            $feat_type1, $depth_org_1_ratio
      ) if ( $org_name2 lt $org_name1 );
    unless ( $fasta1 && $fasta2 ) {
        my $log = $cogeweb->logfile;
        $log =~ s/$DIR/$URL/;
        return
"<span class=alert>Something went wrong generating the fasta files: <a href=$log>log file</a></span>";
    }
    else {
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Fasta creation passed final check.",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "#" x (20) . "\n", $cogeweb->logfile );
    }

    my ($blastdb);
    if ( $ALGO_LOOKUP->{$blast}{formatdb} ) {
        $blastdb = gen_blastdb(
            dbname   => "$dsgid2-$feat_type2",
            fasta    => $fasta2,
            org_name => $org_name2,
            type     => $feat_type2
        );
    }
    else {
        $blastdb = $fasta2;
    }

    #need to convert the blastdb to a fasta file if the algo used is blastz

    unless ($blastdb) {
        my $log = $cogeweb->logfile;
        $log =~ s/$DIR/$URL/;
        return
"<span class=alert>Something went wrong generating the blastdb files: <a href=$log>log file</a></span>";
    }
    else {
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "BlastDB creation passed final check.",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "#" x (20) . "\n", $cogeweb->logfile );
    }
    my $html;

#need to blast each org against itself for finding local dups, then to one another
    my $tmp1 = $org_name1;

    my $tmp2 = $org_name2;
    foreach my $tmp ( $tmp1, $tmp2 ) {
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
    my $orgkey1  = $title1;
    my $orgkey2  = $title2;
    my %org_dirs = (
        $orgkey1 . "_"
          . $orgkey2 => {
            fasta    => $fasta1,
            db       => $blastdb,
            basename => $dsgid1 . "_"
              . $dsgid2
              . ".$feat_type1-$feat_type2."
              . $ALGO_LOOKUP->{$blast}{filename},
            dir => $DIAGSDIR . "/" . $tmp1 . "/" . $tmp2,
          },
    );
    foreach my $org_dir ( keys %org_dirs ) {
        my $outfile = $org_dirs{$org_dir}{dir};
        mkpath( $outfile, 0, 0777 ) unless -d $outfile;
        warn "didn't create path $outfile: $!" unless -d $outfile;
        $outfile .= "/" . $org_dirs{$org_dir}{basename};
        $org_dirs{$org_dir}{blastfile} = $outfile;    #.".blast";
    }

    #blast! use Parallel::ForkManager
    foreach my $key ( keys %org_dirs ) {
        $pm->start and next;
        my $fasta   = $org_dirs{$key}{fasta};
        my $db      = $org_dirs{$key}{db};
        my $outfile = $org_dirs{$key}{blastfile};
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Running genome comparison",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Running " . $ALGO_LOOKUP->{$blast}{displayname},
            $cogeweb->logfile );
        run_blast(
            fasta   => $fasta,
            blastdb => $db,
            outfile => $outfile,
            prog    => $blast,
		  gid1=>$dsgid1
        );    # unless -r $outfile;
        $pm->finish;
    }
    $pm->wait_all_children();

    #check blast runs for problems;  Not forked in order to keep variables
    my $problem = 0;
    foreach my $key ( keys %org_dirs ) {
        my $fasta     = $org_dirs{$key}{fasta};
        my $db        = $org_dirs{$key}{db};
        my $outfile   = $org_dirs{$key}{blastfile};
        my $blast_run = run_blast(
            fasta   => $fasta,
            blastdb => $db,
            outfile => $outfile,
            prog    => $blast,
            ft1     => $feat_type1,
            ft2     => $feat_type2
        );
        $problem = 1 unless $blast_run;
    }
    CoGe::Accessory::Web::write_log( "Completed blast run(s)",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    my $t1 = new Benchmark;
    my $blast_time = timestr( timediff( $t1, $t0 ) );

    #local dup finder
    my $t2            = new Benchmark;
    my $raw_blastfile = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{blastfile};
    my $bedfile1      = $raw_blastfile . ".q.bed";
    my $bedfile2      = $raw_blastfile . ".s.bed";
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Creating .bed files", $cogeweb->logfile );
    blast2bed(
        infile   => $raw_blastfile,
        outfile1 => $bedfile1,
        outfile2 => $bedfile2
    );
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    my $filtered_blastfile = $raw_blastfile;
    $filtered_blastfile .= ".new" if $raw_blastfile =~ /genomic/;
    $filtered_blastfile .= ".tdd$dupdist";
    $filtered_blastfile .= ".cs$cscore" if $cscore;
    $filtered_blastfile .= ".filtered";
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Filtering results of tandem duplicates",
        $cogeweb->logfile );
    run_blast2raw(
        blastfile       => $raw_blastfile,
        bedfile1        => $bedfile1,
        bedfile2        => $bedfile2,
        outfile         => $filtered_blastfile,
        tandem_distance => $dupdist,
        cscore          => $cscore
    );
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    $filtered_blastfile = $raw_blastfile
      unless ( -r $filtered_blastfile && -s $filtered_blastfile )
      || ( -r "$filtered_blastfile.gz" && -s "$filtered_blastfile.gz" );

#    my $synteny_score_db = run_synteny_score (blastfile=>$filtered_blastfile, bedfile1=>$bedfile1, bedfile2=>$bedfile2, outfile=>$org_dirs{$orgkey1."_".$orgkey2}{dir}."/".$dsgid1."_".$dsgid2.".$feat_type1-$feat_type2"); #needed to comment out as the bed files and blast files have changed in SynFind
    my $local_dup_time = timestr( timediff( $t2, $t1 ) );

#prepare dag for synteny analysis
#	my $dag_file12      = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{dir} . "/" . $org_dirs{ $orgkey1 . "_" . $orgkey2 }{basename} . ".dag";
    my $dag_file12     = $filtered_blastfile . ".dag";
    my $dag_file12_all = $dag_file12 . ".all";
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log(
        "Converting blast file to dagchainer input file",
        $cogeweb->logfile );
    $problem = 1
      unless run_dag_tools(
        query      => "a" . $dsgid1,
        subject    => "b" . $dsgid2,
        blast      => $filtered_blastfile,
        outfile    => $dag_file12_all,
        feat_type1 => $feat_type1,
        feat_type2 => $feat_type2
      );
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    my $t2_5 = new Benchmark;
    my $dag_tool_time = timestr( timediff( $t2_5, $t2 ) );

    #is this an ordered gene run?
    my $dag_file12_all_geneorder;
    if ( $dagchainer_type eq "geneorder" ) {
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Converting dagchainer input into gene order coordinates",
            $cogeweb->logfile );
        $dag_file12_all_geneorder = run_convert_to_gene_order(
            infile => $dag_file12_all,
            dsgid1 => $dsgid1,
            dsgid2 => $dsgid2,
            ft1    => $feat_type1,
            ft2    => $feat_type2
        );
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    }
    my $t3 = new Benchmark;
    my $convert_to_gene_order_time = timestr( timediff( $t3, $t2_5 ) );
    my $all_file =
        $dagchainer_type eq "geneorder"
      ? $dag_file12_all_geneorder
      : $dag_file12_all;
    $dag_file12 .= ".go" if $dagchainer_type eq "geneorder";

#B Pedersen's program for automatically adjusting the evals in the dag file to remove bias from local gene duplicates and transposons
#	$dag_file12 .= "_c" . $repeat_filter_cvalue;
#	CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
#	CoGe::Accessory::Web::write_log( "Adjusting evalue of blast hits to correct for repeat sequences", $cogeweb->logfile );
#	run_adjust_dagchainer_evals( infile => $all_file, outfile => $dag_file12, cvalue => $repeat_filter_cvalue );
#	CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
#	CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    my $t3_5 = new Benchmark;

#this step will fail if the dag_file_all is larger than the system memory limit.  If this file does not exist, let's send a warning to the log file and continue with the analysis using the dag_all file
    unless ( ( -r $dag_file12 && -s $dag_file12 )
        || ( -r $dag_file12 . ".gz" && -s $dag_file12 . ".gz" ) )
    {
        $dag_file12 = $all_file;

#		CoGe::Accessory::Web::write_log( "WARNING:  sub run_adjust_dagchainer_evals failed.  Perhaps due to Out of Memory error.  Proceeding without this step!", $cogeweb->logfile );
    }
    my $run_adjust_eval_time = timestr( timediff( $t3_5, $t3 ) );
    #############

    #run dagchainer
    my $dag_merge = 1
      if $merge_algo == 2;    #this is for using dagchainer's merge function;
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Running DagChainer", $cogeweb->logfile );
    my $self_self = 0;
    $self_self = 1 if $dsgid1 eq $dsgid2;
    my ( $dagchainer_file, $merged_dagchainer_file ) = run_dagchainer(
        infile    => $dag_file12,
        D         => $dagchainer_D,
        A         => $dagchainer_A,
        type      => $dagchainer_type,
        Dm        => $Dm,
        gm        => $gm,
        merge     => $dag_merge,
        self_self => $self_self
    );
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    CoGe::Accessory::Web::write_log( "Completed dagchainer run",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    my $t4 = new Benchmark;
    my $run_dagchainer_time = timestr( timediff( $t4, $t3_5 ) );
    my (
        $find_nearby_time, $gen_ks_db_time,
        $dotplot_time,     $add_gevo_links_time
    );
    my $final_results_files;

    if ( -r $dagchainer_file || -r $dagchainer_file . ".gz" ) {
        if ( $merge_algo == 1 )  #id 1 is to specify quota align as a merge algo
        {
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Merging Syntenic Blocks",
                $cogeweb->logfile );
            $merged_dagchainer_file = run_quota_align_merge(
                infile   => $dagchainer_file,
                max_dist => $Dm
            );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        }
        my $post_dagchainer_file =
          -r $merged_dagchainer_file || -r $merged_dagchainer_file . ".gz"
          ? $merged_dagchainer_file
          : $dagchainer_file;  #temp file name for the final post-processed data
        my $post_dagchainer_file_w_nearby = $post_dagchainer_file;
        $post_dagchainer_file_w_nearby =~ s/aligncoords/all\.aligncoords/;

        #add pairs that were skipped by dagchainer
        $post_dagchainer_file_w_nearby = $post_dagchainer_file;

# 	run_find_nearby(infile=>$post_dagchainer_file, dag_all_file=>$all_file, outfile=>$post_dagchainer_file_w_nearby); #program is not working correctly.

        my $t5 = new Benchmark;
        $find_nearby_time = timestr( timediff( $t5, $t4 ) );

        my $quota_align_coverage;
        my $grimm_stuff;
        if ( $depth_algo == 1 )    #id 1 is to specify quota align
        {
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Running Quota Align",
                $cogeweb->logfile );
            $quota_align_coverage = run_quota_align_coverage(
                infile       => $post_dagchainer_file_w_nearby,
                org1         => $depth_org_1_ratio,
                org2         => $depth_org_2_ratio,
                overlap_dist => $depth_overlap
            );
            $grimm_stuff =
              generate_grimm_input( infile => $quota_align_coverage );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        }
        my $final_dagchainer_file =
          $quota_align_coverage
          && ( -r $quota_align_coverage || -r $quota_align_coverage . ".gz" )
          ? $quota_align_coverage
          : $post_dagchainer_file_w_nearby;

        #convert to genomic coordinates if gene order was used
        if ( $dagchainer_type eq "geneorder" ) {
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log(
                "Converting gene order coordinates back to genomic coordinates",
                $cogeweb->logfile
            );
            replace_gene_order_with_genomic_positions(
                file    => $final_dagchainer_file,
                outfile => $final_dagchainer_file . ".gcoords"
            );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

            $final_dagchainer_file = $final_dagchainer_file . ".gcoords";
        }

        #generate dotplot images
        my $org1_length = 0;
        my $org2_length = 0;
        my $chr1_count  = 0;
        my $chr2_count  = 0;
        foreach my $gs ( $dsg1->genomic_sequences ) {
            $chr1_count++;
            $org1_length += $gs->sequence_length;
        }
        foreach my $gs ( $dsg2->genomic_sequences ) {
            $chr2_count++;
            $org2_length += $gs->sequence_length;
        }

        unless ($width) {
            my $max_chr = $chr1_count;
            $max_chr = $chr2_count if $chr2_count > $chr1_count;
            $width   = int( $max_chr * 100 );
            $width   = 1000 if $width > 1000;
            $width   = 500 if $width < 500;
        }
        my $qlead = "a";
        my $slead = "b";
        my $out   = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{dir} . "/html/";
        mkpath( $out, 0, 0777 ) unless -d $out;
        $out .= "master_";
        my ($base) = $final_dagchainer_file =~ /([^\/]*$)/;
        $out .= $base;
        $out .= "_ct$color_type" if defined $color_type;
        $out .= ".w$width";

        #deactivation ks calculations due to how slow it is
        my $ks_db;
        my $ks_blocks_file;
        my $svg_file;
        if ($ks_type) {
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log(
                "Running CodeML for synonymous/nonsynonmous rate calculations",
                $cogeweb->logfile
            );
            $ks_db = gen_ks_db( infile => $final_dagchainer_file );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Generate Ks Blocks File",
                $cogeweb->logfile );
            $ks_blocks_file =
              gen_ks_blocks_file( infile => $final_dagchainer_file );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Generate SVG File",
                $cogeweb->logfile );
            $svg_file = gen_svg_file(
                infile => $ks_blocks_file,
                flip   => $flip,
                xname  => $org_name1,
                yname  => $org_name2
            );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

        }

        my $t6 = new Benchmark;
        $gen_ks_db_time = timestr( timediff( $t6, $t5 ) );
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Generating dotplot",
            $cogeweb->logfile );
        $out = generate_dotplot(
            dag          => $dag_file12_all,
            coords       => $final_dagchainer_file,
            outfile      => "$out",
            regen_images => $regen_images,
            dsgid1       => $dsgid1,
            dsgid2       => $dsgid2,
            width        => $width,
            dagtype      => $dagchainer_type,
            ks_db        => $ks_db,
            ks_type      => $ks_type,
            assemble     => $assemble,
            metric       => $axis_metric,
            ,
            relationship   => $axis_relationship,
            min_chr_size   => $min_chr_size,
            color_type     => $color_type,
            box_diags      => $box_diags,
            fid1           => $fid1,
            fid2           => $fid2,
            snsd           => $snsd,
            flip           => $flip,
            clabel         => $clabel,
            skip_rand      => $skip_rand,
            color_scheme   => $color_scheme,
            chr_sort_order => $chr_sort_order,
            codeml_min     => $codeml_min,
            codeml_max     => $codeml_max,
            log            => $logks
        );
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

        my $hist     = $out . ".hist.png";
        my $spa_file = $out
          . ".spa_info.txt"
          ; #this would be generated by the DOTPLOT program is Syntenic path assembly was requested
        my $t7 = new Benchmark;
        $dotplot_time = timestr( timediff( $t7, $t6 ) );
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Adding GEvo links to final output files",
            $cogeweb->logfile );
        add_GEvo_links(
            infile => $final_dagchainer_file,
            dsgid1 => $dsgid1,
            dsgid2 => $dsgid2
        );
        CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        my $t8 = new Benchmark;
        $add_gevo_links_time = timestr( timediff( $t8, $t7 ) );

        if ( -r "$out.html" ) {
            $html .= qq{
<div class="ui-widget-content ui-corner-all" id="synmap_zoom_box" style="float:left">
 Zoomed SynMap:
 <table class=small>
 <tr>
 <td>Image Width
 <td><input class="backbox" type=text name=zoom_width id=zoom_width size=6 value="400">
<tr>
 <td>Ks, Kn, Kn/Ks cutoffs:
 <td>Min: <input class="backbox" type=text name=zoom_min id=zoom_min size=6 value="};
            $html .= $codeml_min if defined $codeml_min;
            $html .= qq{">
 <td>Max: <input class="backbox" type=text name=zoom_max id=zoom_max size=6 value="};
            $html .= $codeml_max if defined $codeml_max;
            $html .= qq{">
 </table>
</div>
<div style="clear: both;"> </div>
 };

            $/ = "\n";
            open( IN, "$out.html" )
              || warn "problem opening $out.html for reading\n";
            $axis_metric = $axis_metric =~ /g/ ? "genes" : "nucleotides";
            $html .=
              "<span class='small'>Axis metrics are in $axis_metric</span><br>";

            #add version of genome to organism names
            $org_name1 .= " (v" . $dsg1->version . ")";
            $org_name2 .= " (v" . $dsg2->version . ")";
            $html .=
              $flip
              ? "<span class='species small'>y-axis organism: $org_name1</span><br>"
              : "<span class='species small'>y-axis organism: $org_name2</span><br>";
            $html .= qq{<table><tr><td valign=top>};
            my $y_lab = "$out.y.png";
            my $x_lab = "$out.x.png";
            $out =~ s/$DIR/$URL/;

            if ( -r $y_lab ) {
                $y_lab =~ s/$DIR/$URL/;
                $html .= qq
		  {
<div><img src="$y_lab"></div><td>
};

            }

            $/ = "\n";
            my $tmp;
            while (<IN>) {
                next if /<\/?html>/;
                $tmp .= $_;
            }
            close IN;
            $tmp =~ s/master.*\.png/$out.png/;
            warn "$out.html did not parse correctly\n" unless $tmp =~ /map/i;
            $html .= $tmp;

            #check for x-axis label
            if ( -r $x_lab ) {
                $x_lab =~ s/$DIR/$URL/;
                $html .= qq
		  {
<br><img src="$x_lab">
};
            }

            $html .=
              $flip
              ? qq{<br><span class="species small">x-axis: $org_name2</span></table>}
              : qq{<br><span class="species small">x-axis: $org_name1</span></table>};
            $html .=
              "<span class='small'>Axis metrics are in $axis_metric</span><br>";
            $html .=
              "<span class='small'>Algorithm:  " . $algo_name . "</span><br>";

            $html .=
"<br><span class='small link' onclick=window.open('$out.png')>Image File</span><br>";
            $html .=
"<div class='small link ui-widget-content ui-corner-all' style='float:left' onclick=window.open('$out.hist.png')>Histogram of $ks_type values.<br><img src='$out.hist.png'></div><div style='clear: both;'> </div>"
              if -r $hist;

            my $log = $cogeweb->logfile;
            $log =~ s/$DIR/$URL/;
            $html .=
"<span class='small link' onclick=window.open('$log')>Analysis Log (id: $basename)</span><br>";

            $html .= "Links and Downloads:";
            $html .= qq{<table class="small ui-widget-content ui-corner-all">};
            $html .=
              qq{<TR valign=top><td>Homolog search<td>Diagonals<td>Results};
            $html .= qq{<tr valign=top><td>};

    #	    $html .= qq{<span class='small link' onclick=window.open('')></span>};
            $fasta1 =~ s/$DIR/$URL/;
            $html .=
qq{<span class='small link' onclick=window.open('$fasta1')>Fasta file for $org_name1: $feat_type1</span><br>};
            $fasta2 =~ s/$DIR/$URL/;
            $html .=
qq{<span class='small link' onclick=window.open('$fasta2')>Fasta file for $org_name2: $feat_type2</span><br>};
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Processing Tandem Duplicate File",
                $cogeweb->logfile );
            my $org1_localdups = process_local_dups_file(
                infile  => $raw_blastfile . ".q.localdups",
                outfile => $raw_blastfile . ".q.tandems"
            );
            my $org2_localdups = process_local_dups_file(
                infile  => $raw_blastfile . ".s.localdups",
                outfile => $raw_blastfile . ".s.tandems"
            );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

            my $final_dagchainer_file_condensed =
              $final_dagchainer_file . ".condensed";
            my $qa_file = $merged_dagchainer_file;
            $qa_file =~ s/\.ma\d$/\.qa/;
            my $qa_merged_file  = $qa_file . ".merged";
            my $qa_coverage_tmp = $quota_align_coverage . ".tmp"
              if $quota_align_coverage;
            my $qa_coverage_qa = $quota_align_coverage . ".qa"
              if $quota_align_coverage;

            my $file_list = [
                \$raw_blastfile,          \$filtered_blastfile,
                \$bedfile1,               \$bedfile2,
                \$org1_localdups,         \$org2_localdups,
                \$dag_file12_all,         \$dag_file12_all_geneorder,
                \$dag_file12,             \$dagchainer_file,
                \$final_dagchainer_file,  \$final_dagchainer_file_condensed,
                \$merged_dagchainer_file, \$qa_file,
                \$qa_merged_file,         \$ks_blocks_file,
                \$quota_align_coverage,   \$qa_coverage_qa
            ];    #, \$spa_file];

            foreach my $item (@$file_list) {
                $pm->start and next;

                #				$$item = CoGe::Accessory::Web::gzip($$item);
                $pm->finish;
            }
            $pm->wait_all_children();

            foreach my $item (@$file_list) {

                #				$$item = CoGe::Accessory::Web::gzip($$item);
            }
            $raw_blastfile =~ s/$DIR/$URL/;
            $html .=
"<span class='link small' onclick=window.open('$raw_blastfile')>Unfiltered $algo_name results</span><br>";

            $filtered_blastfile =~ s/$DIR/$URL/;
            $html .=
"<span class='link small' onclick=window.open('$filtered_blastfile')>Filtered $algo_name results (no tandem duplicates)</span><br>";
            $org1_localdups =~ s/$DIR/$URL/;
            $html .=
"<span class='link small' onclick=window.open('$org1_localdups')>Tandem Duplicates for $org_name1</span><br>";
            $org2_localdups =~ s/$DIR/$URL/;
            $html .=
"<span class='link small' onclick=window.open('$org2_localdups')>Tandem Duplicates for $org_name2</span><br>";

            $html .= "<td>";
            $dag_file12_all =~ s/$DIR/$URL/;
            $html .=
qq{<span class='small link' onclick=window.open('$dag_file12_all')>DAGChainer Initial Input file</span><br>};
            if ( $dag_file12_all_geneorder && -r $dag_file12_all_geneorder ) {
                $dag_file12_all_geneorder =~ s/$DIR/$URL/;
                $html .=
qq{<span class='small link' onclick=window.open('$dag_file12_all_geneorder')>DAGChainer Input file converted to gene order</span><br>};
            }
            $dag_file12 =~ s/$DIR/$URL/;
            $html .=
qq{<span class='small link' onclick=window.open('$dag_file12')>DAGChainer Input file post repetivive matches filter</span><br>};

            $html .= "<td>";

            $dagchainer_file =~ s/$DIR/$URL/;
            $html .=
qq{<span class='small link' onclick=window.open('$dagchainer_file')>DAGChainer output</span>};
            if ( -r $merged_dagchainer_file ) {
                $merged_dagchainer_file =~ s/$DIR/$URL/;
                $html .=
qq{<br><span class='small link' onclick=window.open('$merged_dagchainer_file')>Merged DAGChainer output</span>};
            }

#	    $final_dagchainer_file =~ s/$DIR/$URL/;
#	    if ($final_dagchainer_file=~/gcoords/)
#	      {
#		$post_dagchainer_file_w_nearby =~ s/$DIR/$URL/;
#		$html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results with nearby genes added</span><br>};
#		$html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results converted back to genomic coordinates</span>};
#	      }
#	    else
#	      {
#		$post_dagchainer_file_w_nearby =~ s/$DIR/$URL/;
#		$html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results with nearby genes added</span>};
#	      }
            if ( $quota_align_coverage && -r $quota_align_coverage ) {
                $quota_align_coverage =~ s/$DIR/$URL/;
                $html .=
qq{<br><span class='small link' onclick=window.open('$quota_align_coverage')>Quota Alignment Results</span>};
            }

            if ( $final_dagchainer_file =~ /gcoords/ ) {

                #		print STDERR $final_dagchainer_file,"\n";
                #		my $tmp= $final_dagchainer_file;
                $final_dagchainer_file =~ s/$DIR/$URL/;

                #		$tmp =~ s/\.gcoords//;

#		$html .= qq{<br><span class='small link' onclick=window.open('$tmp')>DAGChainer output in gene coordinates</span>};
                $html .=
qq{<br><span class='small link' onclick=window.open('$final_dagchainer_file')>DAGChainer output in genomic coordinates</span>};

#		$html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results converted back to genomic coordinates</span>};
            }

            if ($ks_blocks_file) {
                $ks_blocks_file =~ s/$DIR/$URL/;
                $html .=
"<br><span class='small link' onclick=window.open('$ks_blocks_file') target=_new>Results with synonymous/nonsynonymous rate values</span>";
            }
            my $final_file = $final_dagchainer_file;
            $final_file =~ s/$DIR/$URL/;
            $html .=
"<br><span class='small link' onclick=window.open('$final_file')>Final syntenic gene-set output with GEvo links</span>";
            if ( -r $final_dagchainer_file_condensed ) {
                $final_dagchainer_file_condensed =~ s/$DIR/$URL/;
                $html .=
"<br><span class='small link' onclick=window.open('$final_dagchainer_file_condensed')>Condensed syntelog file with GEvo links</span>";
            }
            if ( $svg_file && -r $svg_file ) {
                $svg_file =~ s/$DIR/$URL/;
                $html .=
"<br><span class='small link' onclick=window.open('$svg_file')>SVG Version of Syntenic Dotplot</span>";
            }
            if ( $spa_file && -r $spa_file ) {
                $spa_file =~ s/$DIR/$URL/;
                $html .=
"<br><span class='small link' onclick=window.open('$spa_file')>Syntenic Path Assembly mapping</span>";
            }

            $html .= "<tr><td>";
            my $conffile = $ENV{COGE_HOME} . 'coge.conf';
            $dagchainer_file =~ s/^$URL/$DIR/;
            $html .= "<br>"
              . qq{<span class="small link" id="" onClick="window.open('bin/SynMap/order_contigs_to_chromosome.pl?f=$dagchainer_file&cf=$conffile;l=$tiny_link');" >Generate Pseudo-Assembled Genomic Sequence</span>}
              if $assemble;
            $html .= qq{</table>};

            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Link to Regenerate Analysis",
                $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "$tiny_link", $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

            $html .=
"<a href='$tiny_link' class='ui-button ui-corner-all' style='color: #000000' target=_new_synmap>Regenerate this analysis: $tiny_link</a>";
            if ($ks_type) {
                $html .=
qq{<span  class='ui-button ui-corner-all' onclick="window.open('SynSub.pl?dsgid1=$dsgid1;dsgid2=$dsgid2')">Generate Substitution Matrix of Syntelogs</span>};
            }
            if ($grimm_stuff) {
                my $seq1 = ">$org_name1||" . $grimm_stuff->[0];
                $seq1 =~ s/\n/\|\|/g;
                my $seq2 = ">$org_name2||" . $grimm_stuff->[1];
                $seq2 =~ s/\n/\|\|/g;
                $html .= qq
		  {
<br>
<span class="ui-button ui-corner-all" id = "grimm_link" onclick="post_to_grimm('$seq1','$seq2')" > Rearrangement Analysis</span> <a class="small" href=http://grimm.ucsd.edu/GRIMM/index.html target=_new>(Powered by GRIMM!)</a>
};
            }
            $html .= "<br>";

        }
    }
    else {
        $problem = 1;
        my $log = $cogeweb->logfile;
        $log =~ s/$DIR/$URL/;
        $html .=
"<span class='small link' onclick=window.open('$log')>Analysis Log (id: $basename)</span><br>";

    }
    ##print out all the datafiles created

    $html .= "<br>";

    if ($problem) {
        $html .= qq{
 <span class=alert>There was a problem running your analysis.  Please check the log file for details.</span><br>
 	  };
    }
    email_results(
        email    => $email,
        html     => $html,
        org1     => $org_name1,
        org2     => $org_name2,
        jobtitle => $job_title,
        link     => $tiny_link
    ) if $email;
    my $benchmarks = qq{
$org_name1 versus $org_name2
Benchmarks:
$algo_name:                    $blast_time
Find Local Dups:          $local_dup_time
Dag tools time:           $dag_tool_time
Convert Gene Order:       $convert_to_gene_order_time
Adjust eval time          $run_adjust_eval_time
DAGChainer:               $run_dagchainer_time
find nearby:              $find_nearby_time
Ks calculations:          $gen_ks_db_time
Dotplot:                  $dotplot_time
GEvo links:               $add_gevo_links_time
};

    #    print STDERR $benchmarks;
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Benchmark", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( $benchmarks, $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    $html =~ s/<script src="\/CoGe\/js\/jquery-1.3.2.js"><\/script>//g
      ; #need to remove this from the output from dotplot -- otherwise it over-loads the stuff in the web-page already. This can mess up other loaded js such as tablesoter
    return $html;
}

sub get_previous_analyses {
    my %opts = @_;
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    return unless $oid1 && $oid2;
    my ($org1) = $coge->resultset('Organism')->find($oid1);
    my ($org2) = $coge->resultset('Organism')->find($oid2);
    return
      if ( $USER->user_name =~ /public/i
        && ( $org1->restricted || $org2->restricted ) );
    my ($org_name1) = $org1->name;
    my ($org_name2) = $org2->name;
    ( $oid1, $org_name1, $oid2, $org_name2 ) =
      ( $oid2, $org_name2, $oid1, $org_name1 )
      if ( $org_name2 lt $org_name1 );

    my $tmp1 = $org_name1;
    my $tmp2 = $org_name2;
    foreach my $tmp ( $tmp1, $tmp2 ) {
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

    my $dir = $tmp1 . "/" . $tmp2;
    $dir = "$DIAGSDIR/" . $dir;
    my $sqlite = 0;
    my @items;
    if ( -d $dir ) {
        opendir( DIR, $dir );
        while ( my $file = readdir(DIR) ) {
            $sqlite = 1 if $file =~ /sqlite$/;
            next unless $file =~ /\.aligncoords/;    #$/\.merge$/;
            my ( $D, $g, $A ) = $file =~ /D(\d+)_g(\d+)_A(\d+)/;
            my ($Dm) = $file =~ /Dm(\d+)/;
            my ($gm) = $file =~ /gm(\d+)/;
            my ($ma) = $file =~ /ma(\d+)/;
            $Dm = " " unless defined $Dm;
            $gm = " " unless defined $gm;
            $ma = 0   unless $ma;
            my $merge_algo;
            $merge_algo = "DAGChainer" if $ma && $ma == 2;

            if ( $ma && $ma == 1 ) {
                $merge_algo = "Quota Align";
                $gm         = " ";
            }
            unless ($ma) {
                $merge_algo = "--none--";
                $gm         = " ";
                $Dm         = " ";
            }

            #	    $Dm = 0 unless $Dm;
            #	    $gm = 0 unless $gm;
            next unless ( $D && $g && $A );

            my ($blast) = $file =~
              /^[^\.]+\.[^\.]+\.([^\.]+)/;    #/blastn/ ? "BlastN" : "TBlastX";
            my $select_val;
            foreach my $item ( values %$ALGO_LOOKUP ) {
                if ( $item->{filename} eq $blast ) {
                    $blast      = $item->{displayname};
                    $select_val = $item->{html_select_val};
                }
            }
            my ( $dsgid1, $dsgid2, $type1, $type2 ) =
              $file =~ /^(\d+)_(\d+)\.(\w+)-(\w+)/;
            $type1 = "CDS" if $type1 eq "protein";
            $type2 = "CDS" if $type2 eq "protein";

            #			print STDERR $file,"\n";
            #			my ($repeat_filter) = $file =~ /_c(\d+)/;
            next unless ( $dsgid1 && $dsgid2 && $type1 && $type2 );
            my ($dupdist) = $file =~ /tdd(\d+)/;
            my %data = (

                #									 repeat_filter => $repeat_filter,
                tdd        => $dupdist,
                D          => $D,
                g          => $g,
                A          => $A,
                Dm         => $Dm,
                gm         => $gm,
                ma         => $ma,
                merge_algo => $merge_algo,
                blast      => $blast,
                dsgid1     => $dsgid1,
                dsgid2     => $dsgid2,
                select_val => $select_val
            );
            my $geneorder = $file =~ /\.go/;
            my $dsg1 = $coge->resultset('Genome')->find($dsgid1);
            next unless $dsg1;
            next
              if ( $dsg1->restricted && !$USER->has_access_to_genome($dsg1) );
            my ($ds1) = $dsg1->datasets;
            my $dsg2 = $coge->resultset('Genome')->find($dsgid2);
            next unless $dsg2;
            next
              if ( $dsg2->restricted && !$USER->has_access_to_genome($dsg2) );
            my ($ds2) = $dsg2->datasets;
            $data{dsg1} = $dsg1;
            $data{dsg2} = $dsg2;
            $data{ds1}  = $ds1;
            $data{ds2}  = $ds2;
            my $genome1;
            $genome1 .= $dsg1->name if $dsg1->name;
            $genome1 .= ": "        if $genome1;
            $genome1 .= $ds1->data_source->name;
            my $genome2;
            $genome2 .= $dsg2->name if $dsg2->name;
            $genome2 .= ": "        if $genome2;
            $genome2 .= $ds2->data_source->name;
            $data{genome1}    = $genome1;
            $data{genome2}    = $genome2;
            $data{type_name1} = $type1;
            $data{type_name2} = $type2;
            $type1 = $type1 eq "CDS" ? 1 : 2;
            $type2 = $type2 eq "CDS" ? 1 : 2;
            $data{type1}   = $type1;
            $data{type2}   = $type2;
            $data{dagtype} = $geneorder ? "Ordered genes" : "Distance";
            push @items, \%data;
        }
        closedir(DIR);
    }
    return unless @items;
    my $size = scalar @items;
    $size = 8 if $size > 8;
    my $html;
    my $prev_table = qq{<table id=prev_table class="small resultborder">};
    $prev_table .= qq{<THEAD><TR><TH>}
      . join( "<TH>",
        qw(Org1 Genome1 Ver1 Genome%20Type1 Sequence%20Type1 Org2 Genome2 Ver2 Genome%20Type2 Sequence%20type2 Algo Dist%20Type Dup%20Dist Ave%20Dist(g) Max%20Dist(D) Min%20Pairs(A))
      ) . "</THEAD><TBODY>\n";
    my %seen;

    foreach my $item (
        sort { $b->{dsgid1} <=> $a->{dsgid1} || $b->{dsgid2} <=> $a->{dsgid2} }
        @items )
    {
        my $val = join( "_",
            $item->{g},          $item->{D},       $item->{A},
            $oid1,               $item->{dsgid1},  $item->{type1},
            $oid2,               $item->{dsgid2},  $item->{type2},
            $item->{select_val}, $item->{dagtype}, $item->{tdd} );
        next if $seen{$val};
        $seen{$val} = 1;
        $prev_table .=
          qq{<TR class=feat onclick="update_params('$val')" align=center><td>};
        my $ver1 = $item->{dsg1}->version;
        $ver1 = "0" . $ver1 if $ver1 =~ /^\./;
        my $ver2 = $item->{dsg2}->version;
        $ver2 = "0" . $ver2 if $ver2 =~ /^\./;
        $prev_table .= join( "<td>",
            $item->{dsg1}->organism->name, $item->{genome1},
            $ver1,                         $item->{dsg1}->type->name,
            $item->{type_name1},           $item->{dsg2}->organism->name,
            $item->{genome2},              $ver2,
            $item->{dsg2}->type->name,     $item->{type_name2},
            $item->{blast},                $item->{dagtype},
            $item->{tdd},                  $item->{g},
            $item->{D},                    $item->{A} )
          . "\n";
    }
    $prev_table .= qq{</TBODY></table>};
    $html .= $prev_table;
    $html .=
"<br><span class=small>Synonymous substitution rates previously calculated</span>"
      if $sqlite;
    return "$html";
}

sub get_dsg_info {
    my $dsg       = shift;
    my $length    = 0;
    my $chr_count = 0;

    $length = $dsg->length;    #$rs->first->get_column('total_length');
    $chr_count = $dsg->genomic_sequences->count();
    my ( $gc, $at, $n, $x ) = ( 0, 0, 0, 0 );
    if ( $chr_count < 100 && $length < 50000000 ) {
        ( $gc, $at, $n, $x ) = get_dsg_gc( dsg => $dsg );
    }
    my ( $plasmid, $contig, $scaffold ) = get_chr_types( dsg => $dsg );
    return $gc, $at, $n, $x, $length, $chr_count, $plasmid, $contig, $scaffold;
}

sub get_dsg_gc {
    my %opts  = @_;
    my $dsg   = $opts{dsg};
    my $dsgid = $opts{dsgid};
    my $text  = $opts{text};
    $dsg = $coge->resultset('Genome')->find($dsgid) if $dsgid;
    my ( $gc, $at, $n, $x ) = $dsg->percent_gc;
    $gc *= 100 if ( defined $gc );
    $at *= 100 if ( defined $at );
    $n  *= 100 if ( defined $n );
    $x  *= 100 if ( defined $x );

    if ($text) {
        return "GC: $gc%, AT: $at%, N: $n%, X: $x%";
    }
    else {
        return ( $gc, $at, $n, $x );
    }
}

sub get_chr_types {
    my %opts  = @_;
    my $dsg   = $opts{dsg};
    my $dsgid = $opts{dsgid};
    $dsg = $coge->resultset('Genome')->find($dsgid) if $dsgid;
    my $plasmid  = 0;
    my $contig   = 0;
    my $scaffold = 0;
    my @gs       = $dsg->genomic_sequences;
    if ( @gs > 100 ) {
        return ( 0, 1, 0 );
    }
    foreach my $chr ( map { $_->chromosome } @gs ) {
        $plasmid  = 1 if !$plasmid  && $chr =~ /plasmid/i;
        $contig   = 1 if !$contig   && $chr =~ /contig/i;
        $scaffold = 1 if !$scaffold && $chr =~ /scaffold/i;
    }
    return ( $plasmid, $contig, $scaffold );
}

sub get_pair_info {
    my @anno;
    foreach my $fid (@_) {
        unless ( $fid =~ /^\d+$/ ) {
            push @anno, $fid . "<br>genomic";
            next;
        }
        my $feat = $coge->resultset('Feature')->find($fid);

#		my $anno     = "Name: " . join( ", ", map { "<a class=\"data link\" href=\"$URL/FeatView.pl?accn=" . $_ . "\" target=_new>" . $_ . "</a>" } $feat->names );
#		my $location = "Chr " . $feat->chromosome . " ";
#		$location .= commify( $feat->start ) . " - " . commify( $feat->stop );

        #	$location .=" (".$feat->strand.")";
        #		push @anno, $location . "<br>" . $anno;
        push @anno, $feat->annotation_pretty_print_html;
    }
    return unless @anno;
    my $output =
        "<table class=small valign=top>"
      . join( "\n", ( map { "<tr><td>" . $_ . "</td></tr>" } @anno ) )
      . "</table>";
    my $URL = $P->{URL};
    $output =~ s/window\.open\('(.*?)'\)/window.open('$URL$1')/g;
    return $output;
}

sub check_address_validity {
    my $address = shift;
    return 'valid' unless $address;
    my $validity =
      $address =~
/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/
      ? 'valid'
      : 'invalid';
    return $validity;
}

sub email_results {
    my %opts          = @_;
    my $email_address = $opts{email};

    #my $html 			= $opts{html};
    my $org_name1 = $opts{org1};
    my $org_name2 = $opts{org2};
    my $job_title = $opts{jobtitle};
    my $link      = $opts{link};

    my $subject = "CoGe's SynMap Results are ready";
    $subject .= ": $job_title" if $job_title;
    my $name = $USER->display_name;
    my $body = qq{Dear $name,

Thank you for using SynMap! The results from your latest analysis between $org_name1 and $org_name2 are ready.  To view your results, follow this link and press "Generate SynMap":

$link

Thank you for using the CoGe Software Package.

- The CoGe Team
};

    send_email(
        {
            from    => 'CoGe <coge_results@genomevolution.org>',
            to      => $email_address,
            subject => $subject,
            body    => $body
        }
    );
}

sub get_dotplot {
    my %opts         = @_;
    my $url          = $opts{url};
    my $loc          = $opts{loc};
    my $flip         = $opts{flip} eq "true" ? 1 : 0;
    my $regen        = $opts{regen_images} eq "true" ? 1 : 0;
    my $width        = $opts{width};
    my $ksdb         = $opts{ksdb};
    my $kstype       = $opts{kstype};
    my $metric       = $opts{am};                             #axis metrix
    my $relationship = $opts{ar};                             #axis relationship
    my $max          = $opts{max};
    my $min          = $opts{min};
    my $color_type   = $opts{ct};
    my $box_diags    = $opts{bd};
    my $color_scheme = $opts{color_scheme};
    my $fid1         = $opts{fid1};
    my $fid2         = $opts{fid2};

    #print STDERR Dumper \%opts;
    $box_diags = $box_diags eq "true" ? 1 : 0;

    # base=8_8.CDS-CDS.blastn.dag_geneorder_D60_g30_A5;

    $url = $P->{SERVER} . "run_dotplot_old.pl?" . $url;
    $url .= ";flip=$flip"       if $flip;
    $url .= ";regen=$regen"     if $regen;
    $url .= ";width=$width"     if $width;
    $url .= ";ksdb=$ksdb"       if $ksdb;
    $url .= ";kstype=$kstype"   if $kstype;
    $url .= ";log=1"            if $kstype;
    $url .= ";min=$min"         if defined $min;
    $url .= ";max=$max"         if defined $max;
    $url .= ";am=$metric"       if defined $metric;
    $url .= ";ar=$relationship" if defined $relationship;
    $url .= ";ct=$color_type"   if $color_type;
    $url .= ";bd=$box_diags"    if $box_diags;
    $url .= ";cs=$color_scheme" if defined $color_scheme;
    $url .= ";fid1=$fid1"       if defined $fid1 && $fid1 =~ /^\d+$/;
    $url .= ";fid2=$fid2"       if defined $fid2 && $fid2 =~ /^\d+$/;
    my $content = LWP::Simple::get($url);

    unless ($content) {
        return "Unable to get image for dotplot: $url";
    }
    ($url) = $content =~ /url=(.*?)"/is;
    my $png = $url;
    $png =~ s/html$/png/;
    $png =~ s/$URL/$DIR/;
    my $img = GD::Image->new($png);
    my ( $w, $h ) = $img->getBounds();
    $w += 600;
    $h += 250;

    if ($loc) {
        return ( $url, $loc, $w, $h );
    }
    my $html =
qq{<iframe src=$url frameborder=0 width=$w height=$h scrolling=no></iframe>};
    return $html;
}
