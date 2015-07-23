#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

use CoGeX;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Workflow;
use CoGe::Accessory::Web qw(url_for);
use CoGe::Accessory::Utils qw( commify );
use CoGe::Builder::Tools::SynMap qw(generate_pseudo_assembly);
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use DBIxProfiler;
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
use HTML::Template;
use JSON::XS;
use LWP::UserAgent;
use Parallel::ForkManager;
use GD;
use File::Path;
use File::Spec::Functions;
use Mail::Mailer;
use Benchmark;
use DBI;
use POSIX;
use Sort::Versions;

our (
    $config,       $DEBUG,         $DIR,            $URL,
    $SERVER,       $USER,          $FORM,           $coge,
    $cogeweb,      $PAGE_NAME,     $FORMATDB,       $BLAST,
    $TBLASTX,      $BLASTN,        $BLASTP,         $LASTZ,
    $LAST,         $DATADIR,       $FASTADIR,       $BLASTDBDIR,
    $DIAGSDIR,     $MAX_PROC,      $DAG_TOOL,       $PYTHON,
    $PYTHON26,     $TANDEM_FINDER, $RUN_DAGCHAINER, $EVAL_ADJUST,
    $FIND_NEARBY,  $DOTPLOT,       $SVG_DOTPLOT,    $NWALIGN,
    $QUOTA_ALIGN,  $CLUSTER_UTILS, $BLAST2RAW,      $BASE_URL,
    $BLAST2BED,    $SYNTENY_SCORE, $TEMPDIR,        $TEMPURL,
    $ALGO_LOOKUP,  $GZIP,          $GUNZIP,         %FUNCTIONS,
    $JEX,          $GENE_ORDER,    $PAGE_TITLE,     $KSCALC,
    $GEN_FASTA,    $RUN_ALIGNMENT, $RUN_COVERAGE,   $GEVO_LINKS,
    $PROCESS_DUPS, $DOTPLOT_DOTS,  $SEQUENCE_SIZE_LIMIT
);

$DEBUG = 0;
$|     = 1;    # turn off buffering

# Limit the maximum genome size for genomic-genomic
$SEQUENCE_SIZE_LIMIT = 50_000_000;

$FORM       = new CGI;
$PAGE_TITLE = "SynMap";
$PAGE_NAME  = "$PAGE_TITLE.pl";

( $coge, $USER, $config ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE,
);

$JEX = CoGe::Accessory::Jex->new( host => $config->{JOBSERVER}, port => $config->{JOBPORT} );

$ENV{PATH} = join ":",
  (
    $config->{COGEDIR}, $config->{BINDIR}, $config->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin"
  );
$ENV{BLASTDB}    = $config->{BLASTDB};
$ENV{BLASTMAT}   = $config->{BLASTMATRIX};
$ENV{PYTHONPATH} = "/opt/apache/CoGe/bin/dagchainer_bp";

$BASE_URL = $config->{SERVER};
$DIR      = $config->{COGEDIR};
$URL      = $config->{URL};
$SERVER   = $config->{SERVER};
$TEMPDIR  = $config->{TEMPDIR} . "SynMap";
$TEMPURL  = $config->{TEMPURL} . "SynMap";
$FORMATDB = $config->{FORMATDB};
$MAX_PROC = $config->{MAX_PROC};
$BLAST    = $config->{BLAST} . " -a " . $MAX_PROC . " -K 80 -m 8 -e 0.0001";
my $blast_options = " -num_threads $MAX_PROC -evalue 0.0001 -outfmt 6";
$TBLASTX = $config->{TBLASTX} . $blast_options;
$BLASTN  = $config->{BLASTN} . $blast_options;
$BLASTP  = $config->{BLASTP} . $blast_options;
$LASTZ =
    $config->{PYTHON} . " "
  . $config->{MULTI_LASTZ}
  . " -A $MAX_PROC --path="
  . $config->{LASTZ};
$LAST =
    $config->{MULTI_LAST}
  . " -a $MAX_PROC --path="
  . $config->{LAST_PATH};
# mdb removed 9/20/13 issue 213
#  . " --dbpath="
#  . $config->{LASTDB};
$GZIP          = $config->{GZIP};
$GUNZIP        = $config->{GUNZIP};
$KSCALC        = $config->{KSCALC};
$GEN_FASTA     = $config->{GEN_FASTA};
$RUN_ALIGNMENT = $config->{RUN_ALIGNMENT};
$RUN_COVERAGE  = $config->{RUN_COVERAGE};
$PROCESS_DUPS  = $config->{PROCESS_DUPS};
$GEVO_LINKS =  $config->{GEVO_LINKS};
$DOTPLOT_DOTS = $config->{DOTPLOT_DOTS};

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

$DATADIR  = $config->{DATADIR};
$DIAGSDIR = $config->{DIAGSDIR};
$FASTADIR = $config->{FASTADIR};

mkpath( $FASTADIR,    0, 0777 );
mkpath( $DIAGSDIR,    0, 0777 );    # mdb added 7/9/12
mkpath( $config->{LASTDB}, 0, 0777 );    # mdb added 7/9/12
$BLASTDBDIR = $config->{BLASTDB};

$PYTHON        = $config->{PYTHON};                         #this was for python2.5
$PYTHON26      = $config->{PYTHON};
$DAG_TOOL      = 'nice ' . $config->{DAG_TOOL};
$BLAST2BED     = $config->{BLAST2BED};
$GENE_ORDER    = $DIR . "/bin/SynMap/gene_order.py";
$TANDEM_FINDER = $config->{TANDEM_FINDER}
  . " -d 5 -s -r"
  ; #-d option is the distance (in genes) between dups -- not sure if the -s and -r options are needed -- they create dups files based on the input file name

#$RUN_DAGHAINER = $DIR."/bin/dagchainer/DAGCHAINER/run_DAG_chainer.pl -E 0.05 -s";
$RUN_DAGCHAINER = 'nice ' . $PYTHON26 . " " . $config->{DAGCHAINER};
$EVAL_ADJUST    = $config->{EVALUE_ADJUST};

$FIND_NEARBY = $config->{FIND_NEARBY}
  . " -d 20"
  ; #the parameter here is for nucleotide distances -- will need to make dynamic when gene order is selected -- 5 perhaps?

#programs to run Haibao Tang's quota_align program for merging diagonals and mapping coverage
$QUOTA_ALIGN   = 'nice ' . $config->{QUOTA_ALIGN};     #the program
$CLUSTER_UTILS = $config->{CLUSTER_UTILS};   #convert dag output to quota_align input
$BLAST2RAW     = $config->{BLAST2RAW};       #find local duplicates
$SYNTENY_SCORE = $config->{SYNTENY_SCORE};

$DOTPLOT     = $config->{DOTPLOT} . " -cf " . $config->{_CONFIG_PATH};
$SVG_DOTPLOT = $config->{SVG_DOTPLOT};

#$CONVERT_TO_GENE_ORDER = $DIR."/bin/SynMap/convert_to_gene_order.pl";
#$NWALIGN = $DIR."/bin/nwalign-0.3.0/bin/nwalign";
$NWALIGN = $config->{NWALIGN};

my %ajax = CoGe::Accessory::Web::ajax_func();

#$ajax{read_log}=\&read_log_test;
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
    #read_log               => \&CoGe::Accessory::Web::read_log,
    get_results            => \&get_results,
    generate_assembly      => \&generate_assembly,
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

    #   print $FORM->header; print gen_html();
}

################################################################################
# Web functions
################################################################################

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
      HTML::Template->new( filename => $config->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( PAGE_TITLE => 'SynMap',
                      TITLE      => 'SynMap: Whole Genome Synteny Analysis',
                      HEAD       => qq{},
                      USER       => $USER->display_name || '' );

    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";

    #$template->param(ADJUST_BOX=>1);
    $template->param( BODY       => $body );
    $template->param( HOME       => $config->{SERVER},
                      HELP       => 'SynMap',
                      WIKI_URL   => $config->{WIKI_URL} || '',
                      ADMIN_ONLY => $USER->is_admin,
                      CAS_URL    => $config->{CAS_URL} || '' );
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $form = shift || $FORM;
    my $template = HTML::Template->new( filename => $config->{TMPLDIR} . 'SynMap.tmpl' );

    $template->param( MAIN => 1 );
    #$template->param( EMAIL       => $USER->email )  if $USER->email;

    my $master_width = $FORM->param('w') || 0;
    $template->param( MWIDTH => $master_width );

    #set search algorithm on web-page
    if ( defined( $FORM->param('b') ) ) {
        $template->param( $ALGO_LOOKUP->{ $FORM->param('b') }{opt} => "selected" );
    }
    else {
        $template->param( $ALGO_LOOKUP->{6}{opt} => "selected" );
    }
    my ( $D, $A, $Dm, $gm, $dt, $dupdist, $cscore );
    $D  = $FORM->param('D');
    $A  = $FORM->param('A');
    $Dm = $FORM->param('Dm');
    $gm = $FORM->param('gm');
    $gm //= 40; #/
    $dt     = $FORM->param('dt');
    $cscore = $FORM->param('csco');
    $cscore //= 0; #/
    $dupdist = $FORM->param('tdd');

#   $cvalue = $FORM->param('c');       #different c value than the one for cytology.  But if you get that, you probably shouldn't be reading this code

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

    #   $cvalue = 4 unless defined $cvalue;
    #   $template->param( 'CVALUE'                      => $cvalue );
    $dupdist = 10 unless defined $dupdist;
    $template->param( 'DUPDIST' => $dupdist );
    $template->param( 'CSCORE'  => $cscore );
    $template->param(
        'DISPLAY_DAGCHAINER_SETTINGS' => $display_dagchainer_settings );
    $template->param( 'MIN_CHR_SIZE' => $FORM->param('mcs') )
      if $FORM->param('mcs');

    #will the program automatically run?
    my $autogo = $FORM->param('autogo');
    $autogo = 0 unless defined $autogo;
    $template->param( AUTOGO => $autogo );

#if the page is loading with genomes, there will be a check for whether the genome is rest
#populate organism menus
    my $error = 0;

    for ( my $i = 1 ; $i <= 2 ; $i++ ) {
        my $dsgid = 0;
        $dsgid = $form->param( 'dsgid' . $i )
          if $form->param( 'dsgid' . $i );    #old method for specifying genome
        $dsgid = $form->param( 'gid' . $i )
          if $form->param( 'gid' . $i );      #new method for specifying genome
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

        my ($dsg) = $coge->resultset('Genome')->find($dsgid);

        if($dsgid > 0 and !$USER->has_access_to_genome($dsg)) {
            $error = 1;
        }
    }

    if ($error) {
        $template->param("error" => 'The genome was not found or is restricted.');
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
    $template->param(beta => 1) if $FORM->param("beta");

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

    $name = "Search" unless $name;
    $desc = "Search" unless $desc;
    my ($dsg) = $coge->resultset('Genome')->find($dsgid);

    my $template = HTML::Template->new( filename => $config->{TMPLDIR} . 'partials/organism_menu.tmpl' );
    $template->param(
        ORG_MENU => 1,
        NUM      => $num,
        SEARCH   => $name,
    );

    if ($dsg and $USER->has_access_to_genome($dsg)) {
        my $org = $dsg->organism;
        $oid = $org->id;

        my ( $dsg_info, $feattype_menu, $message ) = get_genome_info(
            dsgid    => $dsgid,
            org_num  => $num,
            feattype => $feattype_param
        );

        $template->param(
            DSG_INFO       => $dsg_info,
            FEATTYPE_MENU  => $feattype_menu,
            GENOME_MESSAGE => $message
        );
    }
    else {
        $oid = 0;
        $dsgid = 0;
    }

    $template->param( 'ORG_LIST' => get_orgs( search => $name, i => $num, oid => $oid ) );

    my ($dsg_menu) = gen_dsg_menu( oid => $oid, dsgid => $dsgid, num => $num );
    $template->param( DSG_MENU => $dsg_menu );

    return $template->output;
}

sub gen_dsg_menu {
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $num   = $opts{num};
    my $dsgid = $opts{dsgid};
    my @dsg_menu;
    my $message;
    my $org_name;

    foreach my $dsg (
        $coge->resultset('Genome')->search(
            { organism_id => $oid },
            {
                prefetch => ['genomic_sequence_type'],
                join     => ['genomic_sequence_type']
            }
        )
      )
    {
        my $name;
        my $has_cds = 0;

        if ( $dsg->restricted && !$USER->has_access_to_genome($dsg) ) {
            next unless $dsgid && $dsg->id == $dsgid;
            $name = "Restricted";
        }
        elsif ($dsg->deleted) {
            if ($dsgid && $dsgid == $dsg->id) {
                $name = "DELETED: ".$dsg->type->name . " (v" . $dsg->version . ",id" . $dsg->id . ")";
            }
            else {
                next;
            }
        }
        else {
            $name .= $dsg->name . ": " if $dsg->name;
            $name .= $dsg->type->name . " (v" . $dsg->version . ",id" . $dsg->id . ")";
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

    return ( qq{<span id="dsgid$num" class="hidden"></span>}, '') unless (@dsg_menu);

    #my $dsg_menu = qq{<select id="dsgid$num" onChange="\$('#dsg_info$num').html('<div class=dna_small class="loading" class="small">loading. . .</div>'); get_genome_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],[handle_dsg_info])">};
    my $dsg_menu =
        qq{<span class="coge-padded-top">} .
        qq{<span class="small text">Genomes: </span>} .
        qq{<select id="dsgid$num" style="max-width:400px;" onChange="get_genome_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],[handle_dsg_info])">} ;

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
        $dsg_menu .= qq{<OPTION VALUE=$numt $selected>$name</option>};
    }
    $dsg_menu .= "</select>";
    $dsg_menu .= "</span>";

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
    my $search = $opts{search};
    my $oid  = $opts{oid};
    my $i    = $opts{i};

    #get rid of trailing white-space
    $search =~ s/^\s+//g if $search;
    $search =~ s/\s+$//g if $search;
    $search = "" if $search && $search =~ /Search/; #need to clear to get full org count

    my @organisms;
    my $org_count;

    # Create terms for search
    my @terms = split /\s+/, $search if defined $search;

    if (scalar @terms or $oid)  {
        my @constraints = map {
            -or => [{ name => {like => qq{%$_%}}},
                    { description => {like => qq{%$_%}}}]
        } @terms;

        @organisms = $coge->resultset("Organism")->search({
            -or => [
                -and => \@constraints,
                { organism_id => $oid },
            ]
        });
    } else {
        $org_count = $coge->resultset("Organism")->count;
    }

    my @opts;
    foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @organisms ) {
        my $option = "<OPTION value=\"" . $item->id . "\"";
        $option .= " selected" if $oid && $oid == $item->id;
        $option .= ">" . $item->name . " (id" . $item->id . ")</OPTION>";
        push @opts, $option;
    }

    unless ( @opts && @organisms) {
        return qq{<span name="org_id$i" id="org_id$i"></span>};
    }

    $org_count = scalar @opts unless $org_count;
    my $html;
    $html .= qq{<span class="small info">Organisms: (}
      . $org_count
      . qq{)</span>\n<BR>\n};


    $html .= qq{<SELECT id="org_id$i" SIZE="5" MULTIPLE onChange="get_genome_info_chain($i)" class="coge-fill-width">\n}
          . join( "\n", @opts )
          . "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/ unless $oid;
    return $html;
}

sub get_genome_info {
    my %opts     = @_;
    my $dsgid    = $opts{dsgid};
    my $org_num  = $opts{org_num};
    my $feattype = $opts{feattype};
    $feattype = 1 unless defined $feattype;
    return ("<div class='small note indent'>No matching results found</div>", " ", " ", '', $org_num, '', '') unless ($dsgid && $dsgid =~ /\d+/);

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
                    qq{<span class="link" onclick="}
                  . qq{search_bar('$_', '#org_name$org_num'); timing('org_name$org_num')">$_</span>}
              } split /\s*;\s*/,
            $org->description
        );
    }
    
    my $i = 0;
    my (
        $percent_gc, $percent_at, $percent_n, $percent_x, $chr_length,
        $chr_count,  $plasmid,    $contig,    $scaffold
    ) = get_dsg_info($dsg);
    my ($ds) = $dsg->datasets;
    my $link = $ds->data_source->link;
    $link = $BASE_URL unless $link;
    $link = "http://" . $link unless $link && $link =~ /^http/;
    $html_dsg_info .= qq{<table class="xsmall" style="margin-left:1em; border-top:1px solid lightgray;">};
    $html_dsg_info .= qq{<tr><td>Description:</td><td class="link" onclick=window.open('GenomeInfo.pl?gid=$dsgid')>}.$dsg->info.qq{</td></tr>};
    #$html_dsg_info .= qq{<tr><td>Organism:</td><td>$orgname</td></tr>}; # redundant
    $html_dsg_info .= qq{<tr><td>Taxonomy:</td><td>$org_desc</td></tr>};
    $html_dsg_info .= "<tr><td>Name: <td>" . $dsg->name if $dsg->name;
    $html_dsg_info .= "<tr><td>Description: <td>" . $dsg->description if $dsg->description;
    $html_dsg_info .= "<tr><td>Source:  <td><a href=" . $link . " target=_new>" . $ds->data_source->name . "</a>";

    #$html_dsg_info .= $dsg->chr_info(summary=>1);
    $html_dsg_info .= "<tr><td>Dataset: <td>" . $ds->name;
    $html_dsg_info .= ": " . $ds->description if $ds->description;
    $html_dsg_info .= "<tr><td>Chromosomes: <td>" . commify($chr_count);
    if ( $percent_gc > 0 ) {
        $html_dsg_info .= "<tr><td>DNA content: <td>GC: $percent_gc%, AT: $percent_at%, N: $percent_n%, X: $percent_x%";
    }
    else {
        $html_dsg_info .= qq{<tr><td>DNA content: <td id='gc_content$org_num' class='link' onclick="get_gc($dsgid, 'gc_content$org_num')">Click to retrieve};
    }
    $html_dsg_info .= "<tr><td>Total length: <td>" . commify($chr_length);
    $html_dsg_info .= "<tr><td>Contains plasmid" if $plasmid;
    $html_dsg_info .= "<tr><td>Contains contigs" if $contig;
    $html_dsg_info .= "<tr><td>Contains scaffolds" if $scaffold;
    $html_dsg_info .= qq{<tr class="alert"><td>Restricted:</td><td>Yes} if $dsg->restricted;
    $html_dsg_info .= "</table>";
    if ( $dsg->restricted && !$USER->has_access_to_genome($dsg) ) {
        $html_dsg_info = "Restricted";
    }
    if ($dsg->deleted)
      {
    $html_dsg_info = "<span class='alert'>This genome has been deleted and cannot be used in this analysis.</span>  <a href='GenomeInfo.pl?gid=$dsgid' target=_new>More information</a>.";
      }

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

    my $feattype_menu = qq{<select id="feat_type$org_num" name ="feat_type$org_num">#};
    $feattype_menu .= qq{<OPTION VALUE=1 $cds_selected>CDS</option>} if $has_cds;
    $feattype_menu .= qq{<OPTION VALUE=2 $genomic_selected>genomic</option>};
    $feattype_menu .= "</select>";
    $message = "<span class='small alert'>No Coding Sequence in Genome</span>"
      unless $has_cds;

    return $html_dsg_info, $feattype_menu, $message, $chr_length, $org_num,
      $dsg->organism->name, $dsg->genomic_sequence_type_id;
}

sub get_previous_analyses {

    #FIXME:  THis whole sub needs updating or removal!  Lyons 6/12/13
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

            #       $Dm = 0 unless $Dm;
            #       $gm = 0 unless $gm;
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

            #           print STDERR $file,"\n";
            #           my ($repeat_filter) = $file =~ /_c(\d+)/;
            next unless ( $dsgid1 && $dsgid2 && $type1 && $type2 );
            my ($dupdist) = $file =~ /tdd(\d+)/;
            my %data = (

           #                                    repeat_filter => $repeat_filter,
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
            my $genome1 = $coge->resultset('Genome')->find($dsgid1);
            next unless $genome1;
            next
              if ( $genome1->restricted && !$USER->has_access_to_genome($genome1) );
            my ($ds1) = $genome1->datasets;
            my $genome2 = $coge->resultset('Genome')->find($dsgid2);
            next unless $genome2;
            next
              if ( $genome2->restricted && !$USER->has_access_to_genome($genome2) );
            my ($ds2) = $genome2->datasets;
            $data{dsg1} = $genome1;
            $data{dsg2} = $genome2;
            $data{ds1}  = $ds1;
            $data{ds2}  = $ds2;
            
            $genome1 .= $genome1->name if $genome1->name;
            $genome1 .= ": "        if $genome1;
            $genome1 .= $ds1->data_source->name;
            $genome2 .= $genome2->name if $genome2->name;
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
    $chr_count = $dsg->chromosome_count;
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
    $gc *= 100;
    $at *= 100;
    $n  *= 100;
    $x  *= 100;

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
    my @chromosomes = $dsg->chromosomes;
    if ( @chromosomes > 100 ) {
        return ( 0, 1, 0 );
    }
    foreach my $chr ( @chromosomes ) {
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

#       my $anno     = "Name: " . join( ", ", map { "<a class=\"data link\" href=\"$URL/FeatView.pl?accn=" . $_ . "\" target=_new>" . $_ . "</a>" } $feat->names );
#       my $location = "Chr " . $feat->chromosome . " ";
#       $location .= commify( $feat->start ) . " - " . commify( $feat->stop );

        #   $location .=" (".$feat->strand.")";
        #       push @anno, $location . "<br>" . $anno;
        push @anno, $feat->annotation_pretty_print_html;
    }
    return unless @anno;
    my $output =
        "<table class=small valign=top>"
      . join( "\n", ( map { "<tr><td>" . $_ . "</td></tr>" } @anno ) )
      . "</table>";
    my $URL = $config->{URL};
    $output =~ s/window\.open\('(.*?)'\)/window.open('$URL$1')/g;
    return $output;
}

sub get_query_link {
    my %url_options = @_;
    my $dagchainer_D = $url_options{D};

 #  my $dagchainer_g = $url_options{g}; #depreciated -- will be a factor of -D
    my $dagchainer_A = $url_options{A};
    my $Dm           = $url_options{Dm};
    my $gm           = $url_options{gm};
    ($Dm) = $Dm =~ /(\d+)/;
    ($gm) = $gm =~ /(\d+)/;

#   my $repeat_filter_cvalue = $url_options{c}; #parameter to be passed to run_adjust_dagchainer_evals
    my $cscore  = $url_options{csco}; #c-score for filtering low quality blast hits, fed to blast to raw
    my $dupdist = $url_options{tdd}; #tandem duplication distance, fed to blast to raw
    my $regen_images = $url_options{regen_images};
    my $email        = $url_options{email};
    my $job_title    = $url_options{jobtitle};
    my $width        = $url_options{width};
    my $basename     = $url_options{basename};
    my $blast        = $url_options{blast};

    my $feat_type1 = $url_options{feat_type1};
    my $feat_type2 = $url_options{feat_type2};

    my $dsgid1 = $url_options{dsgid1};
    my $dsgid2 = $url_options{dsgid2};

    unless($dsgid1 and $dsgid2) {
        return encode_json({
            error => "Missing a genome id."
        });
    }

    my $ks_type  = $url_options{ks_type};
    my $assemble = $url_options{assemble} =~ /true/i ? 1 : 0;
    $assemble = 2 if $assemble && $url_options{show_non_syn} =~ /true/i;
    $assemble *= -1 if $assemble && $url_options{spa_ref_genome} < 0;
    my $axis_metric       = $url_options{axis_metric};
    my $axis_relationship = $url_options{axis_relationship};
    my $min_chr_size      = $url_options{min_chr_size};
    my $dagchainer_type   = $url_options{dagchainer_type};
    my $color_type        = $url_options{color_type};
    my $merge_algo        = $url_options{merge_algo};    #is there a merging function?

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

    # Sort by genome id
    (
        $dsgid1, $org_name1, $feat_type1, $depth_org_1_ratio,
        $dsgid2, $org_name2, $feat_type2, $depth_org_2_ratio
    ) = (
        $dsgid2, $org_name2, $feat_type2, $depth_org_2_ratio,
        $dsgid1, $org_name1, $feat_type1, $depth_org_1_ratio
    ) if ( $dsgid2 lt $dsgid1 );

    my $synmap_link =
        $SERVER
      . "SynMap.pl?dsgid1=$dsgid1;dsgid2=$dsgid2"
      . ";D=$dagchainer_D;A=$dagchainer_A;w=$width;b=$blast;ft1=$feat_type1;"
      . "ft2=$feat_type2;autogo=1";

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

    my $tiny_link = CoGe::Accessory::Web::get_tiny_link(url => $synmap_link);

    return $tiny_link;
}

sub generate_basefile {
    $cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
    return $cogeweb->basefilename;
}

sub go {
    my %opts = @_;
    foreach my $k ( keys %opts ) {
        $opts{$k} =~ s/^\s+//;
        $opts{$k} =~ s/\s+$//;
    }
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};

    return encode_json({
        success =>  JSON::false,
        error => "You must select two genomes."
    }) unless ($dsgid1 && $dsgid2);

    my ($genome1) = $coge->resultset('Genome')->find($dsgid1);
    my ($genome2) = $coge->resultset('Genome')->find($dsgid2);

    return encode_json({
        success =>  JSON::false,
        error => "The Genome $dsgid1 could not be found."
    }) unless $genome1;

    return encode_json({
        success =>  JSON::false,
        error => "The Genome $dsgid2 could not be found."
    }) unless $genome2;

    return encode_json({
        success =>  JSON::false,
        error => "Genome $dsgid1 primary data is missing."
    }) unless -r $genome1->file_path;

    return encode_json({
        success =>  JSON::false,
        error => "Genome $dsgid2 primary data is missing."
    }) unless -r $genome2->file_path;

    my ( $dir1, $dir2 ) = sort ( $dsgid1, $dsgid2 );
    ############################################################################
    # Fetch organism name and title
    ############################################################################
    my $feat_type1 = $opts{feat_type1};
    my $feat_type2 = $opts{feat_type2};

    # Block large genomic-genomic jobs from running
    if (($feat_type1 == 2 && $genome1->length > $SEQUENCE_SIZE_LIMIT && !$genome1->type->name =~/hard/i) &&
        ($feat_type2 == 2 && $genome2->length > $SEQUENCE_SIZE_LIMIT && !$genome2->type->name =~/hard/i)) {
         
        return encode_json({
            success => JSON::false,
            error => "The analysis was blocked: " .
                     "a comparison of two unmasked and unannotated genomes larger than 50Mb requires many days to weeks to finish. " .
                     "Please use at least one annotated genome in the analysis."
        });
    }

    my $basename = $opts{basename};
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    my ( $org_name1, $title1 ) = gen_org_name(
        dsgid     => $dsgid1,
        feat_type => $feat_type1,
        write_log => 1
    );

    my ( $org_name2, $title2 ) = gen_org_name(
        dsgid     => $dsgid2,
        feat_type => $feat_type2,
        write_log => 1
    );

    my $ks_type = $opts{ks_type};
    ############################################################################
    # Initialize Jobs
    ############################################################################

    my $tiny_link = get_query_link(@_);

    my $log_msg =
        #"<a href='OrganismView.pl?dsgid=$dsgid1' target='_blank'>$org_name1</a> v. <a href='OrganismView.pl?dsgid=$dsgid2' target='_blank'>$org_name2</a>";
        "$org_name1 v. $org_name2";
    $log_msg .= " Ks" if $ks_type;

    say STDERR "tiny_link is required for logging." unless defined($tiny_link);

    my ($tiny_id) = $tiny_link =~ /\/(\w+)$/;
    my $workflow_name = "synmap-$tiny_id";

    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Creating Workflow", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Link to Regenerate Analysis",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "$tiny_link", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "",           $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Created Workflow: synmap-$workflow_name",
        $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    ############################################################################
    # Parameters
    ############################################################################

    # blast options
    my $blast = $opts{blast};

    # blast2bed options

    # dagchainer options
    #my $dagchainer_g = $opts{g}; #depreciated -- will be a factor of -D
    my $dagchainer_D = $opts{D};
    my $dagchainer_A = $opts{A};
    my $Dm           = $opts{Dm};
    my $gm           = $opts{gm};
    ($Dm) = $Dm =~ /(\d+)/;
    ($gm) = $gm =~ /(\d+)/;

    #$dagchainer_type = $dagchainer_type eq "true" ? "geneorder" : "distance";

    #my $repeat_filter_cvalue = $opts{c};
    ##parameter to be passed to run_adjust_dagchainer_evals

    #c-score for filtering low quality blast hits, fed to blast to raw
    my $cscore = $opts{csco};

    #tandem duplication distance, fed to blast to raw
    my $dupdist = defined( $opts{tdd} ) ? $opts{tdd} : 10;

    # dotplot options
    my $regen = $opts{regen_images};
    my $regen_images      = ( $regen && $regen eq "true" ) ? 1 : 0;
    my $job_title         = $opts{jobtitle};
    my $width             = $opts{width};
    my $axis_metric       = $opts{axis_metric};
    my $axis_relationship = $opts{axis_relationship};
    my $min_chr_size      = $opts{min_chr_size};
    my $dagchainer_type   = $opts{dagchainer_type};
    my $color_type        = $opts{color_type};
    my $merge_algo = $opts{merge_algo};    #is there a merging function?
                                           #will non-syntenic dots be shown?
    my $snsd      = $opts{show_non_syn_dots} =~ /true/i ? 1 : 0;
    my $algo_name = $ALGO_LOOKUP->{$blast}{displayname};

    #will the axis be flipped?
    my $flip = $opts{flip} =~ /true/i ? 1 : 0;

    #are axes labeled?
    my $clabel = $opts{clabel} =~ /true/i ? 1 : 0;

    #are random chr skipped
    my $skip_rand = $opts{skip_rand} =~ /true/i ? 1 : 0;

    #which color scheme for ks/kn dots?
    my $color_scheme = $opts{color_scheme};

    #fids that are passed in for highlighting the pair in the dotplot
    my $fid1 = $opts{fid1};
    my $fid2 = $opts{fid2};

    #draw a box around identified diagonals?
    my $box_diags = $opts{box_diags};
    $box_diags = $box_diags eq "true" ? 1 : 0;

    #how are the chromosomes to be sorted?
    my $chr_sort_order = $opts{chr_sort_order};

    #codeml min and max calues
    my $codeml_min = $opts{codeml_min};
    $codeml_min = undef
      unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
    my $codeml_max = $opts{codeml_max};
    $codeml_max = undef
      unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
    my $logks = $opts{logks};
    $logks = $logks eq "true" ? 1 : 0;

    my $email = 0 if check_address_validity( $opts{email} ) eq 'invalid';

    my $assemble = $opts{assemble} =~ /true/i ? 1 : 0;
    $assemble = 2 if $assemble && $opts{show_non_syn} =~ /true/i;
    $assemble *= -1 if $assemble && $opts{spa_ref_genome} < 0;

    #options for finding syntenic depth coverage by quota align (Bao's algo)
    my $depth_algo        = $opts{depth_algo};
    my $depth_org_1_ratio = $opts{depth_org_1_ratio};
    my $depth_org_2_ratio = $opts{depth_org_2_ratio};
    my $depth_overlap     = $opts{depth_overlap};

    $feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
    $feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
    $feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
    $feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time


    # Sort by genome id
    (
        $dsgid1,     $genome1,              $org_name1,
        $feat_type1, $depth_org_1_ratio, $dsgid2,     $genome2,
        $org_name2,  $feat_type2, $depth_org_2_ratio
      )
      = (
        $dsgid2,     $genome2,              $org_name2,
        $feat_type2, $depth_org_2_ratio, $dsgid1,     $genome1,
        $org_name1,  $feat_type1, $depth_org_1_ratio
      ) if ( $dsgid2 lt $dsgid1 );


    ############################################################################
    # Generate Fasta files
    ############################################################################
    my ( $fasta1, $fasta2 );
    my $workflow = undef;
    my $status   = undef;

    #my @dsgs = ([$dsgid1, $feat_type1]);
    #push @dsgs, [$dsgid2, $feat_type2]
    #  unless $dsgid1 == $dsgid2 && $feat_type1 eq $feat_type2;
    #
    #    foreach my $item (@dsgs) {
    #        my $dsgid = $item->[0];
    #        my $feat_type = $item->[1];
    #
    #        #TODO: Schedule fasta generation only if feat_type not genomic
    #        #if ($feat_type eq "genomic") {
    #        #    my $genome = $coge->resultset('Genome')->find($gid);
    #        #    $file = $genome->file_path;
    #        #} else {
    #        #    $file = $FASTADIR . "/$gid-$feat_type.fasta";
    #        #}
    #    }
    $workflow = $JEX->create_workflow(
        name     => $workflow_name,
        logfile  => $cogeweb->logfile,
    );

    if ( $feat_type1 eq "genomic" ) {
        $fasta1 = $genome1->file_path;

        CoGe::Accessory::Web::write_log( "Fetched fasta file for:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( " " x (2) . $org_name1,
            $cogeweb->logfile );
    }
    else {
        my @fasta1args = ();
        $fasta1 = $FASTADIR . "/$dsgid1-$feat_type1.fasta";
        push @fasta1args, [ "--config",       $config->{_CONFIG_PATH},     0 ];
        push @fasta1args, [ "--genome_id",    $dsgid1,     1 ];
        push @fasta1args, [ "--feature_type", $feat_type1, 1 ];
        push @fasta1args, [ "--fasta",        $fasta1,     1 ];

        $workflow->add_job({
            cmd         => $GEN_FASTA,
            script      => undef,
            args        => \@fasta1args,
            inputs      => undef,
            outputs     => [$fasta1],
            description => "Generating fasta file...",
        });

        CoGe::Accessory::Web::write_log( "Added fasta file generation for:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( " " x (2) . $org_name1,
            $cogeweb->logfile );
    }

    if ( $feat_type2 eq "genomic" ) {
        $fasta2 = $genome2->file_path;

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Fetched fasta file for:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( " " x (2) . $org_name2,
            $cogeweb->logfile );
    }
    else {
        $fasta2 = $FASTADIR . "/$dsgid2-$feat_type2.fasta";
        my @fasta2args = ();
        push @fasta2args, [ "--config",       $config->{_CONFIG_PATH},     0 ];
        push @fasta2args, [ "--genome_id",    $dsgid2,     1 ];
        push @fasta2args, [ "--feature_type", $feat_type2, 1 ];
        push @fasta2args, [ "--fasta",        $fasta2,     1 ];

        $workflow->add_job({
            cmd         => $GEN_FASTA,
            script      => undef,
            args        => \@fasta2args,
            inputs      => undef,
            outputs     => [$fasta2],
            description => "Generating fasta file...",
        });

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Added fasta file generation for:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( " " x (2) . $org_name2,
            $cogeweb->logfile );
    }

    ############################################################################
    # Generate blastdb files
    ############################################################################
    my ( $blastdb, @blastdb_files );

    if ( $ALGO_LOOKUP->{$blast}{formatdb} ) {
        my $write_log = 0;
        my $basename  = "$BLASTDBDIR/$dsgid2-$feat_type2";

        my @blastdbargs = ();
        push @blastdbargs, [ '-p', $feat_type2 eq "protein" ? "T" : "F", 1 ];
        push @blastdbargs, [ '-i', $fasta2, 1 ];
        push @blastdbargs, [ '-t', '"' . $org_name2 . '"', 1 ];
        push @blastdbargs, [ '-n', $basename, 1 ];

        $blastdb = $basename;
        $basename .= $feat_type2 eq "protein" ? ".p" : ".n";

        push @blastdb_files, "$basename" . "sq";
        push @blastdb_files, "$basename" . "in";
        push @blastdb_files, "$basename" . "hr";

        $workflow->add_job({
            cmd         => $FORMATDB,
            script      => undef,
            args        => \@blastdbargs,
            inputs      => [$fasta2],
            outputs     => \@blastdb_files,
            description => "Generating BlastDB...",
        });

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Added BlastDB generation",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( $blastdb, $cogeweb->logfile );
    }
    else {
        $blastdb = $fasta2;
        push @blastdb_files, $blastdb;
    }
    my ( $orgkey1, $orgkey2 ) = ( $title1, $title2 );
    my %org_dirs = (
        $orgkey1 . "_"
          . $orgkey2 => {
            fasta    => $fasta1,
            db       => $blastdb,
            basename => $dsgid1 . "_"
              . $dsgid2
              . ".$feat_type1-$feat_type2."
              . $ALGO_LOOKUP->{$blast}{filename},
            dir => "$DIAGSDIR/$dir1/$dir2",
          },
    );

    foreach my $org_dir ( keys %org_dirs ) {
        my $outfile = $org_dirs{$org_dir}{dir};
        mkpath( $outfile, 0, 0777 ) unless -d $outfile;
        warn "didn't create path $outfile: $!" unless -d $outfile;
        $outfile .= "/" . $org_dirs{$org_dir}{basename};
        $org_dirs{$org_dir}{blastfile} = $outfile;    #.".blast";
    }

    ############################################################################
    # Run Blast
    ############################################################################

    #check blast runs for problems;  Not forked in order to keep variables
    my $raw_blastfile = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{blastfile};

    foreach my $key ( keys %org_dirs ) {
        my $cmd = 'nice ' . $ALGO_LOOKUP->{$blast}{algo};    #$prog =~ /tblastx/i ? $TBLASTX : $BLASTN;

        my $fasta   = $org_dirs{$key}{fasta};
        my $db      = $org_dirs{$key}{db};
        my $outfile = $org_dirs{$key}{blastfile};
        my @blastargs;

        if ( $cmd =~ /lastz/i ) {
            push @blastargs, [ "-i", $fasta,   0 ];
            push @blastargs, [ "-d", $db,      0 ];
            push @blastargs, [ "-o", $outfile, 1 ];
        }
        elsif ( $cmd =~ /last_wrapper/i ) {
            # mdb added 9/20/13 issue 213
            my $dbpath = $config->{LASTDB} . '/' . $dsgid2;
            mkpath($dbpath, 0, 0777);
            push @blastargs, [ "--dbpath", $dbpath, 0 ];

            push @blastargs, [ "",   $db,      0 ];
            push @blastargs, [ "",   $fasta,   0 ];
            push @blastargs, [ "-o", $outfile, 1 ];
        }
        else {
            push @blastargs, [ "-out",   $outfile, 1 ];
            push @blastargs, [ "-query", $fasta,   0 ];
            push @blastargs, [ "-db",    $db,      0 ];
        }

        ( undef, $cmd ) = CoGe::Accessory::Web::check_taint($cmd);
        push @blastdb_files, $fasta;
        $workflow->add_job({
            cmd         => $cmd,
            script      => undef,
            args        => \@blastargs,
            inputs      => \@blastdb_files,
            outputs     => [$outfile],
            description => "Running genome comparison...",
        });
    }

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log(
        "Added genome comparison (algorithm: "
          . $ALGO_LOOKUP->{$blast}{displayname} . ")",
        $cogeweb->logfile
    );

    ###########################################################################
    # Converting blast to bed and finding local duplications
    ###########################################################################
    #
    # NOTES: The blast2bed program will take the input rawblast file and
    # filter it and creating a new rawblast and moving the old to
    # rawblast.orig
    #
    my $query_bed   = $raw_blastfile . ".q.bed";
    my $subject_bed = $raw_blastfile . ".s.bed";

    my @blastargs = ();
    push @blastargs, [ '-infile',   $raw_blastfile, 1 ];
    push @blastargs, [ '-outfile1', $query_bed,     1 ];
    push @blastargs, [ '-outfile2', $subject_bed,   1 ];

    my @bedoutputs = ();
    push @bedoutputs, $query_bed;
    push @bedoutputs, $subject_bed;
    push @bedoutputs, $raw_blastfile if ( $raw_blastfile =~ /genomic/ );
#    push @bedoutputs, "$raw_blastfile.orig" if ( $raw_blastfile =~ /genomic/ );
    push @bedoutputs, "$raw_blastfile.new" if ( $raw_blastfile =~ /genomic/ );
    $workflow->add_job({
        cmd         => $BLAST2BED,
        script      => undef,
        args        => \@blastargs,
        inputs      => [$raw_blastfile],
        outputs     => \@bedoutputs,
        description => "Creating .bed files...",
    });

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Added .bed files creation",
        $cogeweb->logfile );

    ###########################################################################
    # Converting blast to raw and finding local duplications
    ###########################################################################
    my $qlocaldups = $raw_blastfile.".q.localdups";
    my $slocaldups = $raw_blastfile.".s.localdups";
    $raw_blastfile .= ".new" if $raw_blastfile =~ /genomic/;
    my $filtered_blastfile = $raw_blastfile;
#    $filtered_blastfile .=".new" if $raw_blastfile =~ /genomic/;
    $filtered_blastfile .= ".tdd$dupdist";
    $filtered_blastfile .= ".cs$cscore" if $cscore < 1;
    $filtered_blastfile .= ".filtered";

    my @rawargs = ();
    push @rawargs, [ "", $raw_blastfile, 1 ];
    push @rawargs, [ "--localdups", "", 1 ];
    push @rawargs, [ "--qbed",        $query_bed,          1 ];
    push @rawargs, [ "--sbed",        $subject_bed,        1 ];
    push @rawargs, [ "--tandem_Nmax", $dupdist,            1 ];
    push @rawargs, [ "--cscore",      $cscore,             1 ] if $cscore < 1;
    push @rawargs, [ ">",             $filtered_blastfile, 1 ];

    my @rawoutputs = ();
    push @rawoutputs, $filtered_blastfile;
    push @rawoutputs, $qlocaldups;
    push @rawoutputs, $slocaldups;

    $workflow->add_job({
        cmd         => $BLAST2RAW,
        script      => undef,
        args        => \@rawargs,
        inputs      => [ $raw_blastfile, $query_bed, $subject_bed ],
        outputs     => \@rawoutputs,
        description => "Filtering tandem dups...",
    });

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    my $msg = "Added Filtering results of tandem ";
    $msg .= "duplicates (tandem duplication distance: $dupdist";
    $msg .= ", c-score: $cscore" if $cscore;
    $msg .= ")";

    CoGe::Accessory::Web::write_log( $msg, $cogeweb->logfile );

 #TODO: This feature is currently disabled
 #needed to comment out as the bed files and blast files have changed in SynFind
 #   my $synteny_score_db = run_synteny_score(
 #       blastfile => $filtered_blastfile,
 #       bedfile1  => $bedfile1,
 #       bedfile2  => $bedfile2,
 #       outfile   => $org_dirs{$orgkey1 . "_" . $orgkey2}{dir} . "/"
 #         . $dsgid1 . "_"
 #         . $dsgid2
 #         . ".$feat_type1-$feat_type2");

    ############################################################################
    # Run dag tools - Prepares DAG for syntenic analysis.
    ############################################################################
    my $dag_file12       = $filtered_blastfile . ".dag";
    my $dag_file12_all   = $dag_file12 . ".all";
    my $query_dup_file   = $opts{query_dup_files};
    my $subject_dup_file = $opts{subject_dup_files};
    my $query            = "a" . $dsgid1;
    my $subject          = "b" . $dsgid2;

    my @dagtoolargs = ();
    push @dagtoolargs, [ '-q', $query,              1 ];
    push @dagtoolargs, [ '-s', $subject,            1 ];
    push @dagtoolargs, [ '-b', $filtered_blastfile, 1 ];
    push @dagtoolargs, [ '-c', "", 1 ];

    push @dagtoolargs, [ '--query_dups', $query_dup_file, 1 ]
      if $query_dup_file;
    push @dagtoolargs, [ '--subject_dups', $subject_dup_file, 1 ]
      if $subject_dup_file;
    push @dagtoolargs, [ '>', $dag_file12_all, 1 ];

    $workflow->add_job({
        cmd         => $DAG_TOOL,
        script      => undef,
        args        => \@dagtoolargs,
        inputs      => [$filtered_blastfile],
        outputs     => [$dag_file12_all],
        description => "Formatting for DagChainer...",
    });

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log(
        "Added convertion of blast file to dagchainer input file",
        $cogeweb->logfile );

    ############################################################################
    # Convert to gene order
    ############################################################################
    my $dag_file12_all_geneorder = "$dag_file12_all.go";
    my $all_file;

    if ( $dagchainer_type eq "geneorder" ) {

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Added convertion of dagchainer input into gene order coordinates",
            $cogeweb->logfile
        );

        my @geneorderargs = ();
        push @geneorderargs, [ "",           $dag_file12_all,           1 ];
        push @geneorderargs, [ "",           $dag_file12_all_geneorder, 1 ];
        push @geneorderargs, [ "--gid1",     $dsgid1,                   1 ];
        push @geneorderargs, [ "--gid2",     $dsgid2,                   1 ];
        push @geneorderargs, [ "--feature1", $feat_type1,               1 ];
        push @geneorderargs, [ "--feature2", $feat_type2,               1 ];

        $workflow->add_job({
            cmd         => $GENE_ORDER,
            script      => undef,
            args        => \@geneorderargs,
            inputs      => [$dag_file12_all],
            outputs     => [$dag_file12_all_geneorder],
            description => "Converting to genomic order...",
        });

        $all_file = $dag_file12_all_geneorder;
        $dag_file12 .= ".go";
    }
    else {
        $all_file = $dag_file12_all;
    }

 #B Pedersen's program for automatically adjusting the evals in the dag file to
 # remove bias from local gene duplicates and transposons
 #   $dag_file12 .= "_c" . $repeat_filter_cvalue;
 #   CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
 #   CoGe::Accessory::Web::write_log( "Adjusting evalue of blast hits to correct
 #   for repeat sequences", $cogeweb->logfile );
 #   run_adjust_dagchainer_evals( infile => $all_file, outfile => $dag_file12,
 #   cvalue => $repeat_filter_cvalue );
 #   CoGe::Accessory::Web::write_log( "#" x (20), $cogeweb->logfile );
 #   CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    # This step will fail if the dag_file_all is larger than the system memory
    # limit. If this file does not exist, let's send a warning to the log file
    # and continue with the analysis using the dag_all file
    unless ( ( -r $dag_file12 && -s $dag_file12 )
        || ( -r $dag_file12 . ".gz" && -s $dag_file12 . ".gz" ) )
    {
        $dag_file12 = $all_file;
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "WARNING: sub run_adjust_dagchainer_evals failed. "
              . "Perhaps due to Out of Memory error. "
              . "Proceeding without this step!",
            $cogeweb->logfile
        );
    }

    ############################################################################
    # Run dagchainer
    ############################################################################
    my ( $dagchainer_file, $merged_dagchainer_file );

    #this is for using dagchainer's merge function
    my $dag_merge_enabled = ( $merge_algo == 2 ) ? 1 : 0;
    my $self_comparision = ( $dsgid1 eq $dsgid2 ) ? 1 : 0;

    #length of a gap (average distance expected between two syntenic genes)
    my $gap = defined( $opts{g} ) ? $opts{g} : floor( $dagchainer_D / 2 );
    $gap = 1 if $gap < 1;
    $dagchainer_file = $dag_file12;
    $dagchainer_file .= "_D$dagchainer_D" if defined $dagchainer_D;
    $dagchainer_file .= "_g$gap"          if defined $gap;
    $dagchainer_file .= "_A$dagchainer_A" if defined $dagchainer_A;
    $dagchainer_file .= "_Dm$Dm"          if $dag_merge_enabled;
    $dagchainer_file .= "_gm$gm"          if $dag_merge_enabled;
    $dagchainer_file .= ".aligncoords";
    $dagchainer_file .= ".ma2.dag"        if $dag_merge_enabled;

    my @dagargs = ();
    push @dagargs, [ "-E", "0.05", 1 ];
    push @dagargs, [ "-i",   $dag_file12,   1 ];
    push @dagargs, [ "-D",   $dagchainer_D, 1 ] if defined $dagchainer_D;
    push @dagargs, [ "-g",   $gap,          1 ] if defined $gap;
    push @dagargs, [ "-A",   $dagchainer_A, 1 ] if defined $dagchainer_A;
    push @dagargs, [ "--Dm", $Dm,           1 ] if $dag_merge_enabled;
    push @dagargs, [ "--m",  $gm,           1 ] if $dag_merge_enabled;
    push @dagargs, [ "--new_behavior", "", 1 ] if $self_comparision;

    ###MERGING OF DIAGONALS FUNCTION
    # --merge $outfile
    # this will automatically cause the merging of diagonals to happen.
    # $outfile will be created and $outfile.merge (which is the inappropriately
    # named merge of diagonals file.
    # --gm  average distance between diagonals
    # --Dm  max distance between diagonals
    ##Both of these parameters' default values is 4x -g and -D respectively.
    my $post_dagchainer_file;

    if ($dag_merge_enabled) {
        print_debug( msg => "DAG-Merge is disabled.", enabled => $DEBUG );
        $merged_dagchainer_file = "$dagchainer_file.merged";
        push @dagargs, [ "--merge", $merged_dagchainer_file, 1 ];

        $workflow->add_job({
            cmd         => $RUN_DAGCHAINER,
            script      => undef,
            args        => \@dagargs,
            inputs      => [$dag_file12],
            outputs     => [$merged_dagchainer_file],
            description => "Running DAGChainer (with merge)...",
        });
        $post_dagchainer_file = $merged_dagchainer_file;
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Added DagChainer (merge: enabled,"
              . " Maximum distance between two blocks: $Dm genes, "
              . " Average distance expected between syntenic blocks: $gm genes)",
            $cogeweb->logfile
        );
    }
    else {
        push @dagargs, [ ">", $dagchainer_file, 1 ];

        $workflow->add_job({
            cmd         => $RUN_DAGCHAINER,
            script      => undef,
            args        => \@dagargs,
            inputs      => [$dag_file12],
            outputs     => [$dagchainer_file],
            description => "Running DAGChainer...",
        });

        $post_dagchainer_file = $dagchainer_file;
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Added DagChainer (merge: disabled)",
            $cogeweb->logfile );
    }

    ############################################################################
    # Run quota align merge
    ############################################################################
    my $final_results_files;

    #id 1 is to specify quota align as a merge algo
    if ( $merge_algo == 1 ) {
        $merged_dagchainer_file = "$dagchainer_file.Dm$Dm.ma1";

        my @mergeargs = ();
        push @mergeargs, [ '--config',       $config->{_CONFIG_PATH},                 0 ];
        push @mergeargs, [ '--infile',       $dagchainer_file,        1 ];
        push @mergeargs, [ '--outfile',      $merged_dagchainer_file, 1 ];
        push @mergeargs, [ '--max_distance', $Dm,                     1 ];

        $workflow->add_job({
            cmd         => $RUN_ALIGNMENT,
            script      => undef,
            args        => \@mergeargs,
            inputs      => [$dagchainer_file],
            outputs     => [$merged_dagchainer_file],
            description => "Merging Syntenic Blocks...",
        });
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Added Merge Syntenic Blocks"
              . " (Maximum distance between two blocks: $Dm genes)",
            $cogeweb->logfile
        );

        $post_dagchainer_file = $merged_dagchainer_file;
    }
    my $post_dagchainer_file_w_nearby = $post_dagchainer_file;
    $post_dagchainer_file_w_nearby =~ s/aligncoords/all\.aligncoords/;

    print_debug(
        msg     => "Post dag chainer file: $post_dagchainer_file",
        enabled => $DEBUG
    );

    #add pairs that were skipped by dagchainer
    $post_dagchainer_file_w_nearby = $post_dagchainer_file;

    #TODO: Run find nearby is currently disabled
    #run_find_nearby(
    #    infile       => $post_dagchainer_file,
    #    dag_all_file => $all_file,
    #    outfile      => $post_dagchainer_file_w_nearby
    #);    #program is not working correctly.

    ############################################################################
    # Run quota align coverage
    ############################################################################
    my ( $quota_align_coverage, $grimm_stuff, $final_dagchainer_file );

    if ( $depth_algo == 1 )    #id 1 is to specify quota align
    {
        $quota_align_coverage = $post_dagchainer_file_w_nearby;
        $quota_align_coverage .= ".qac" . $depth_org_1_ratio . ".";
        $quota_align_coverage .= $depth_org_2_ratio . "." . $depth_overlap;

        print_debug( msg => $post_dagchainer_file_w_nearby, enabled => $DEBUG );

        my @depthargs = ();
        push @depthargs, [ '--config',  $config->{_CONFIG_PATH},        0 ];
        push @depthargs, [ '--infile',  $post_dagchainer_file_w_nearby, 1 ];
        push @depthargs, [ '--outfile', $quota_align_coverage,          1 ];
        push @depthargs, [ '--depth_ratio_org1', $depth_org_1_ratio, 1 ];
        push @depthargs, [ '--depth_ratio_org2', $depth_org_2_ratio, 1 ];
        push @depthargs, [ '--depth_overlap',    $depth_overlap,     1 ];

        $workflow->add_job({
            cmd         => $RUN_COVERAGE,
            script      => undef,
            args        => \@depthargs,
            inputs      => [$post_dagchainer_file_w_nearby],
            outputs     => [$quota_align_coverage],
            description => "Calculating Syntenic depth...",
        });

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
            "Added Syntenic Depth "
              . "(ratio: $depth_org_1_ratio/$depth_org_2_ratio, "
              . "overlap: $depth_overlap)",
            $cogeweb->logfile
        );
        CoGe::Accessory::Web::write_log( "$org_name1 -vs- $org_name2",
            $cogeweb->logfile );
        $final_dagchainer_file = $quota_align_coverage;
    }
    else {
        $final_dagchainer_file = $post_dagchainer_file_w_nearby;
    }

    print_debug(
        msg     => "Final dag chainer file: $final_dagchainer_file",
        enabled => $DEBUG
    );

    if ( $dagchainer_type eq "geneorder" ) {
        my @positionargs = ();
        push @positionargs, [ '', $final_dagchainer_file, 1 ];
        push @positionargs, [ '', "$final_dagchainer_file.gcoords", 1 ];
        push @positionargs, [ "--positional", '', 1 ];

        $workflow->add_job({
            cmd         => $GENE_ORDER,
            script      => undef,
            args        => \@positionargs,
            inputs      => [$final_dagchainer_file],
            outputs     => ["$final_dagchainer_file.gcoords"],
            description => "Converting to genomic coordinates...",
        });

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log(
"Added conversion gene order coordinates back to genomic coordinates",
            $cogeweb->logfile
        );

        $final_dagchainer_file = $final_dagchainer_file . ".gcoords";
    }

    #generate dotplot images
    unless ($width) {
    	my $chr1_count = $genome1->chromosome_count;
    	my $chr2_count = $genome2->chromosome_count;
        my $max_chr = $chr1_count;
        $max_chr = $chr2_count if $chr2_count > $chr1_count;
        $width   = int( $max_chr * 100 );
        $width   = 1000 if $width > 1000;
        $width   = 500 if $width < 500;
    }

    ############################################################################
    # Create html output directory
    ############################################################################
    my ( $qlead, $slead ) = ( "a", "b" );
    my $out = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{dir} . "/html/";
    mkpath( $out, 0, 0777 ) unless -d $out;
    $out .= "master_";
    my ($base) = $final_dagchainer_file =~ /([^\/]*$)/;
    $out .= $base;

    my $json_basename = "$out";
    $out .= "_ct$color_type" if defined $color_type;
    $out .= ".w$width";

    ############################################################################
    # KS Calculations
    ############################################################################
    my ( $ks_db, $ks_blocks_file, $svg_file );
    my $check_ks = $final_dagchainer_file =~ /^(.*?CDS-CDS)/;
    $check_ks = $final_dagchainer_file =~ /^(.*?protein-protein)/
        unless $check_ks;

    if ($ks_type and $check_ks) {
            $ks_db          = "$final_dagchainer_file.sqlite";
            $ks_blocks_file = "$final_dagchainer_file.ks";

            my @ksargs = ();
            push @ksargs, [ '--config',    $config->{_CONFIG_PATH}, 0 ];
            push @ksargs, [ '--infile',    $final_dagchainer_file,  1 ];
            push @ksargs, [ '--dbfile',    $ks_db,                  1 ];
            push @ksargs, [ '--blockfile', $ks_blocks_file,         1 ];

            $workflow->add_job({
                cmd         => $KSCALC,
                script      => undef,
                args        => \@ksargs,
                inputs      => [$final_dagchainer_file],
                outputs     => [ $ks_blocks_file, $ks_db ],
                description => "Calculating synonymous changes (slow)...",
            });

            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
            CoGe::Accessory::Web::write_log(
"Added ($ks_type) calculation of syntenic CDS pairs and color dots",
                $cogeweb->logfile
            );

            ####################################################################
            # Generate svg dotplot
            ####################################################################

            my @svgargs = ();
            push @svgargs, [ '--dag_file', $ks_blocks_file, 1 ];
            push @svgargs, [ '--flip', "", 1 ] if $flip;
            push @svgargs, [ '--xhead', '"' . $org_name1 . '"', 1 ]
              if $org_name1;
            push @svgargs, [ '--yhead', '"' . $org_name2 . '"', 1 ]
              if $org_name2;
            push @svgargs, [ '--output', $ks_blocks_file, 1 ];

            $svg_file = $ks_blocks_file . ".svg";
            $workflow->add_job({
                cmd         => $SVG_DOTPLOT,
                script      => undef,
                args        => \@svgargs,
                inputs      => [$ks_blocks_file],
                outputs     => [$svg_file],
                description => "Generating svg image...",
            });

            CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
            CoGe::Accessory::Web::write_log( "Added generation of svg dotplot",
                $cogeweb->logfile );
    } else {
        $ks_type = undef;

        ####################################################################
        # Generate svg dotplot
        ####################################################################

        my @svgargs = ();
        push @svgargs, [ '--dag_file', $final_dagchainer_file, 1 ];
        push @svgargs, [ '--flip', "", 1 ] if $flip;
        push @svgargs, [ '--xhead', '"' . $org_name1 . '"', 1 ]
            if $org_name1;
        push @svgargs, [ '--yhead', '"' . $org_name2 . '"', 1 ]
            if $org_name2;
        push @svgargs, [ '--output', $final_dagchainer_file, 1 ];
        push @svgargs, [ '--no-ks', "", 1 ];

        $svg_file = $final_dagchainer_file . ".svg";

        $workflow->add_job({
            cmd         => $SVG_DOTPLOT,
            script      => undef,
            args        => \@svgargs,
            inputs      => [$final_dagchainer_file],
            outputs     => [$svg_file],
            description => "Generating svg image...",
        });

        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "Added generation of svg dotplot",
            $cogeweb->logfile );

    }

    ############################################################################
    # Generate dot plot
    ############################################################################
    my @plotargs   = ();
    my @plotinputs = ();

    ($basename) =
      $final_dagchainer_file =~ /([^\/]*aligncoords.*)/;    #.all.aligncoords/;
    $width = 1000 unless defined($width);

    my $dotfile = "$out";
    $dotfile .= ".spa$assemble"     if $assemble;
    $dotfile .= ".gene"             if $axis_metric =~ /gene/i;
    $dotfile .= ".s"                if $axis_relationship =~ /s/i;
    $dotfile .= ".mcs$min_chr_size" if $min_chr_size;
    $dotfile .= ".$fid1"            if $fid1;
    $dotfile .= ".$fid2"            if $fid2;

    #add ks_db to dotplot command if requested
    if ($ks_type) {
        $dotfile .= ".$ks_type";

        push @plotargs, [ '-ksdb', $ks_db,   1 ];
        push @plotargs, [ '-kst',  $ks_type, 1 ];
        push @plotargs, [ '-log',  $logks,   1 ];
        push @plotinputs, $ks_db;
    }
    $dotfile .= ".box"                if $box_diags;
    $dotfile .= ".flip"               if $flip;
    $dotfile .= ".c0"                 if $clabel eq 0;
    $dotfile .= ".sr"                 if $skip_rand;
    $dotfile .= ".cs$color_scheme"    if defined $color_scheme;
    $dotfile .= ".cso$chr_sort_order" if defined $chr_sort_order;
    $dotfile .= ".min$codeml_min"     if defined $codeml_min;
    $dotfile .= ".max$codeml_max"     if defined $codeml_max;
    $dotfile .= ".log"                if $logks;

    #are non-syntenic dots being displayed
    if ($snsd) {
        push @plotargs, [ '-d', $dag_file12_all, 0 ];
        push @plotinputs, $dag_file12_all;
    }
    else {
        $dotfile .= ".nsd";    #no syntenic dots, yes, nomicalture is confusing.
    }

    push @plotargs, [ '-a', $final_dagchainer_file, 1 ];
    push @plotargs, [ '-b', $dotfile, 1 ];

    my $jsoption = "";
    $jsoption .= qq{'javascript:synteny_zoom(};
    $jsoption .= qq{"$dsgid1",};
    $jsoption .= qq{"$dsgid2",};
    $jsoption .= qq{"$basename",};
    $jsoption .= $flip ? qq{"YCHR","XCHR"} : qq{"XCHR","YCHR"};
    $jsoption .= qq{,"$ks_db"} if $ks_db;
    $jsoption .= qq{)'};

    push @plotargs, [ '-l',    $jsoption, 0 ];
    push @plotargs, [ '-dsg1', $dsgid1,   1 ];
    push @plotargs, [ '-dsg2', $dsgid2,   1 ];
    push @plotargs, [ '-w',    $width,    1 ];
    push @plotargs, [ '-lt', 2, 1 ];
    push @plotargs, [ '-assemble', $assemble,    1 ] if $assemble;
    push @plotargs, [ '-am',       $axis_metric, 1 ] if $axis_metric;
    push @plotargs, [ '-fb', '', 1 ]
      if $axis_relationship && $axis_relationship =~ /s/;
    push @plotargs, [ '-mcs', $min_chr_size, 1 ] if $min_chr_size;
    push @plotargs, [ '-cdt', $color_type,   1 ] if $color_type;
    push @plotargs, [ '-bd', 1, 1 ] if $box_diags;
    push @plotargs, [ '-fid1', $fid1, 1 ] if $fid1;
    push @plotargs, [ '-fid2', $fid2, 1 ] if $fid2;
    push @plotargs, [ '-f',      1, 1 ] if $flip;
    push @plotargs, [ '-labels', 0, 1 ] if $clabel eq 0;
    push @plotargs, [ '-sr',     1, 1 ] if $skip_rand;
    push @plotargs, [ '-color_scheme', $color_scheme, 1 ]
      if defined $color_scheme;
    push @plotargs, [ '-chr_sort_order', $chr_sort_order, 1 ]
      if defined $chr_sort_order;
    push @plotargs, [ '-min', $codeml_min, 1 ] if defined $codeml_min;
    push @plotargs, [ '-max', $codeml_max, 1 ] if defined $codeml_max;

    push @plotinputs, $final_dagchainer_file;

    my $hist = $dotfile . ".hist.png";

    # this would be generated by the DOTPLOT program is Syntenic path assembly
    # was requested
    my $spa_file = $dotfile . ".spa_info.txt";
    my @plotoutputs = ( "$dotfile.html", "$dotfile.png" );

    push @plotoutputs, ("$dotfile.x.png", "$dotfile.y.png") if $clabel eq "1";
    push @plotoutputs, $hist     if $ks_db;
    push @plotoutputs, $spa_file if $assemble;

    $workflow->add_job({
        cmd         => $DOTPLOT,
        script      => undef,
        args        => \@plotargs,
        inputs      => \@plotinputs,
        outputs     => \@plotoutputs,
        overwrite   => $regen_images,
        description => "Generating images...",
    });

#    my $dot_args = [
#        [ '-cf', $config->{_CONFIG_PATH}, 0 ],
#        [ '-genome1', $dsgid1, 0 ],
#        [ '-genome2', $dsgid2, 0 ],
#        [ '-a', $final_dagchainer_file, 1 ],
#        [ '-b', $json_basename, 1 ],
#    ];
#
#    push @$dot_args, [ '-ksdb', $ks_db,   1 ] if $ks_db;
#
#    my $dot_inputs = [
#        $final_dagchainer_file,
#    ];
#
#    my $dot_outputs = [
#        "$json_basename.json",
#    ];
#
#    if ($dag_file12_all) {
#        push @$dot_args, [ '-d', $dag_file12_all, 0 ];
#        push @$dot_inputs, $dag_file12_all;
#        push @$dot_outputs, "$json_basename.all.json";
#    }
#
#    if ($ks_db) {
#        push @$dot_inputs, $ks_db if $ks_db;
#        push @$dot_outputs, "$json_basename.datasets.json";
#    }
#
#    $workflow->add_job(
#        cmd         => $DOTPLOT_DOTS,
#        script      => undef,
#        args        => $dot_args,
#        inputs      => $dot_inputs,
#        outputs     => $dot_outputs,
#        description => "Generating dotplot dots...",
#    );
#
#    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
#    CoGe::Accessory::Web::write_log( "Added dotplot generation",
#        $cogeweb->logfile );

    ############################################################################
    # Post Processing
    ############################################################################

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Final Post Processing",
        $cogeweb->logfile);

    my $subject_dup_args = [
        ['--config', $config->{_CONFIG_PATH},                         0 ],
        ['--infile', $slocaldups, 1],#$raw_blastfile . ".s.localdups", 1 ],
        ['--outfile', $raw_blastfile . ".s.tandems",  1 ],
    ];

    $workflow->add_job({
        cmd         => $PROCESS_DUPS,
        script      => undef,
        args        => $subject_dup_args,
        inputs      => [$slocaldups], #[$raw_blastfile . ".s.localdups"],
        outputs     => [$raw_blastfile . ".s.tandems"],
        description => "Processing Subject Tandem Duplicate File...",
    });

    my $query_dup_args = [
        ['--config',  $config->{_CONFIG_PATH},                         0 ],
        ['--infile',  $qlocaldups,1], #$raw_blastfile . ".q.localdups", 1 ],
        ['--outfile', $raw_blastfile . ".q.tandems",  1 ],
    ];

    $workflow->add_job({
        cmd         => $PROCESS_DUPS,
        script      => undef,
        args        => $query_dup_args,
        inputs      => [$qlocaldups], #[$raw_blastfile . ".q.localdups"],
        outputs     => [$raw_blastfile . ".q.tandems"],
        description => "Processing Query Tandem Duplicate File...",
    });
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Added Processing of Tandem Duplicate Files",
        $cogeweb->logfile );

    my $condensed = "$final_dagchainer_file.condensed";

    my $link_args = [
        ['--config', $config->{_CONFIG_PATH}, 0],
        ['--infile', $final_dagchainer_file, 1],
        ['--dsgid1', $dsgid1, 1],
        ['--dsgid2', $dsgid2, 1],
        ['--outfile', $condensed, 1],
    ];

    $workflow->add_job({
        cmd         => $GEVO_LINKS,
        script      => undef,
        args        => $link_args,
        inputs      => [$final_dagchainer_file],
        outputs      => [$condensed],
        description => "Generating GEvo links...",
    });

    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Added GEvo links generation", $cogeweb->logfile);

    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Running Workflow", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "#" x (25), $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );

    my $response = $JEX->submit_workflow($workflow);
    
    if ($response and $response->{id}) {
        my $log = CoGe::Accessory::Web::log_history(
            db          => $coge,
            user_id     => $USER->id,
            description => $log_msg,
            page        => $PAGE_TITLE,
            link        => $tiny_link,
            parent_id   => $response->{id},
            parent_type => 7 #FIXME magic number
        );
    
        return encode_json({
            link     => $tiny_link,
            request  => "jex/synmap/status/" . $response->{id},
            status   => $response->{status},
            success  => $JEX->is_successful($response) ? JSON::true : JSON::false
        });
    }
    else {
        return encode_json({
            success  => JSON::false
        });
    }
}

sub get_results {
    my %opts = @_;
    foreach my $k ( keys %opts ) {
        $opts{$k} =~ s/^\s+//;
        $opts{$k} =~ s/\s+$//;
    }
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};

    return encode_json({
        error => "You must select two genomes."
    }) unless ( $dsgid1 && $dsgid2 );

    my ($genome1) = $coge->resultset('Genome')->find($dsgid1);
    my ($genome2) = $coge->resultset('Genome')->find($dsgid2);

    return encode_json({
        error => "Problem generating dataset group objects for ids:  $dsgid1, $dsgid2."
    }) unless ( $genome1 && $genome2 );

    my ( $dir1, $dir2 ) = sort ( $dsgid1, $dsgid2 );

    ############################################################################
    # Initialize Job info
    ############################################################################
    my $tiny_link = get_query_link(@_);

    say STDERR "tiny_link is required for logging." unless defined($tiny_link);

#    my ($tiny_id) = $tiny_link =~ /\/(\w+)$/;
#    my $workflow_name .= "-$tiny_id";

    my $basename = $opts{basename};
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR
    );

    ############################################################################
    # Parameters
    ############################################################################

    # blast options
    my $blast = $opts{blast};

    # blast2bed options

    # dagchainer options
    #my $dagchainer_g = $opts{g}; #depreciated -- will be a factor of -D
    my $dagchainer_D = $opts{D};
    my $dagchainer_A = $opts{A};
    my $Dm           = $opts{Dm};
    my $gm           = $opts{gm};
    ($Dm) = $Dm =~ /(\d+)/;
    ($gm) = $gm =~ /(\d+)/;

    #$dagchainer_type = $dagchainer_type eq "true" ? "geneorder" : "distance";

#   my $repeat_filter_cvalue = $opts{c};              #parameter to be passed to run_adjust_dagchainer_evals

    #c-score for filtering low quality blast hits, fed to blast to raw
    my $cscore = $opts{csco};

    #tandem duplication distance, fed to blast to raw
    my $dupdist = defined( $opts{tdd} ) ? $opts{tdd} : 10;

    # dotplot options
    my $regen = $opts{regen_images};
    my $regen_images      = ( $regen and $regen eq "true" ) ? 1 : 0;
    my $job_title         = $opts{jobtitle};
    my $width             = $opts{width};
    my $axis_metric       = $opts{axis_metric};
    my $axis_relationship = $opts{axis_relationship};
    my $min_chr_size      = $opts{min_chr_size};
    my $dagchainer_type   = $opts{dagchainer_type};
    my $color_type        = $opts{color_type};
    my $merge_algo = $opts{merge_algo};    #is there a merging function?
                                           #will non-syntenic dots be shown?
    my $snsd      = $opts{show_non_syn_dots} =~ /true/i ? 1 : 0;
    my $algo_name = $ALGO_LOOKUP->{$blast}{displayname};

    #will the axis be flipped?
    my $flip = $opts{flip} =~ /true/i ? 1 : 0;

    #are axes labeled?
    my $clabel = $opts{clabel} =~ /true/i ? 1 : 0;

    #are random chr skipped
    my $skip_rand = $opts{skip_rand} =~ /true/i ? 1 : 0;

    #which color scheme for ks/kn dots?
    my $color_scheme = $opts{color_scheme};

    #fids that are passed in for highlighting the pair in the dotplot
    my $fid1 = $opts{fid1};
    my $fid2 = $opts{fid2};

    #draw a box around identified diagonals?
    my $box_diags = $opts{box_diags};
    $box_diags = $box_diags eq "true" ? 1 : 0;

    #how are the chromosomes to be sorted?
    my $chr_sort_order = $opts{chr_sort_order};

    #codeml min and max calues
    my $codeml_min = $opts{codeml_min};
    $codeml_min = undef
      unless $codeml_min =~ /\d/ && $codeml_min =~ /^-?\d*.?\d*$/;
    my $codeml_max = $opts{codeml_max};
    $codeml_max = undef
      unless $codeml_max =~ /\d/ && $codeml_max =~ /^-?\d*.?\d*$/;
    my $logks = $opts{logks};
    $logks = $logks eq "true" ? 1 : 0;

    my $email = 0 if check_address_validity( $opts{email} ) eq 'invalid';

    my $feat_type1 = $opts{feat_type1};
    my $feat_type2 = $opts{feat_type2};

    my $ks_type = $opts{ks_type};
    my $assemble = $opts{assemble} =~ /true/i ? 1 : 0;
    $assemble = 2 if $assemble && $opts{show_non_syn} =~ /true/i;
    $assemble *= -1 if $assemble && $opts{spa_ref_genome} < 0;

    #options for finding syntenic depth coverage by quota align (Bao's algo)
    my $depth_algo        = $opts{depth_algo};
    my $depth_org_1_ratio = $opts{depth_org_1_ratio};
    my $depth_org_2_ratio = $opts{depth_org_2_ratio};
    my $depth_overlap     = $opts{depth_overlap};

    $feat_type1 = $feat_type1 == 2 ? "genomic" : "CDS";
    $feat_type2 = $feat_type2 == 2 ? "genomic" : "CDS";
    $feat_type1 = "protein" if $blast == 5 && $feat_type1 eq "CDS"; #blastp time
    $feat_type2 = "protein" if $blast == 5 && $feat_type2 eq "CDS"; #blastp time

    ############################################################################
    # Fetch organism name and title
    ############################################################################
    my ( $org_name1, $title1 ) = gen_org_name(
        dsgid     => $dsgid1,
        feat_type => $feat_type1,
    );

    my ( $org_name2, $title2 ) = gen_org_name(
        dsgid     => $dsgid2,
        feat_type => $feat_type2,
    );

    ############################################################################
    # Generate Fasta files
    ############################################################################
    my ( $fasta1, $fasta2 );
    my $workflow = undef;
    my $status   = undef;

    if ( $feat_type1 eq "genomic" ) {
        my $genome = $coge->resultset('Genome')->find($dsgid1);
        $fasta1 = $genome->file_path;
    }
    else {
        $fasta1 = $FASTADIR . "/$dsgid1-$feat_type1.fasta";
    }

    if ( $feat_type2 eq "genomic" ) {
        my $genome = $coge->resultset('Genome')->find($dsgid2);
        $fasta2 = $genome->file_path;
    }
    else {
        $fasta2 = $FASTADIR . "/$dsgid2-$feat_type2.fasta";
    }

    # Sort by genome id
    (
        $dsgid1,     $genome1,              $org_name1,  $fasta1,
        $feat_type1, $depth_org_1_ratio, $dsgid2,     $genome2,
        $org_name2,  $fasta2,            $feat_type2, $depth_org_2_ratio
      )
      = (
        $dsgid2,     $genome2,              $org_name2,  $fasta2,
        $feat_type2, $depth_org_2_ratio, $dsgid1,     $genome1,
        $org_name1,  $fasta1,            $feat_type1, $depth_org_1_ratio
      ) if ( $dsgid2 lt $dsgid1 );

    ############################################################################
    # Generate blastdb files
    ############################################################################
    my ( $blastdb, @blastdb_files );

    if ( $ALGO_LOOKUP->{$blast}{formatdb} ) {
        my $write_log = 0;
        my $basename  = "$BLASTDBDIR/$dsgid2-$feat_type2";

        $blastdb = $basename;
        $basename .= $feat_type2 eq "protein" ? ".p" : ".n";
    }
    else {
        $blastdb = $fasta2;
    }

    my ( $html, $warn );

    my ( $orgkey1, $orgkey2 ) = ( $title1, $title2 );
    my %org_dirs = (
        $orgkey1 . "_"
          . $orgkey2 => {
            fasta    => $fasta1,
            db       => $blastdb,
            basename => $dsgid1 . "_"
              . $dsgid2
              . ".$feat_type1-$feat_type2."
              . $ALGO_LOOKUP->{$blast}{filename},
            dir => "$DIAGSDIR/$dir1/$dir2",
          },
    );

    foreach my $org_dir ( keys %org_dirs ) {
        my $outfile = $org_dirs{$org_dir}{dir};
        $outfile .= "/" . $org_dirs{$org_dir}{basename};
        $org_dirs{$org_dir}{blastfile} = $outfile;    #.".blast";
    }

    ############################################################################
    # Run Blast
    ############################################################################
    my $raw_blastfile = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{blastfile};
    $raw_blastfile .=".new" if $raw_blastfile =~ /genomic/;

    ###########################################################################
    # Converting blast to bed and finding local duplications
    ###########################################################################
    #
    # NOTES: The blast2bed program will take the input rawblast file and
    # filter it and creating a new rawblast and moving the old to
    # rawblast.orig
    #
    my $query_bed   = $raw_blastfile . ".q.bed";
    my $subject_bed = $raw_blastfile . ".s.bed";

    ###########################################################################
    # Converting blast to raw and finding local duplications
    ###########################################################################
    my $filtered_blastfile = $raw_blastfile;
    $filtered_blastfile .= ".tdd$dupdist";
    $filtered_blastfile .= ".cs$cscore" if $cscore < 1;
    $filtered_blastfile .= ".filtered";

    if ( $cscore == 1 ) {
        $warn = 'Please choose a cscore less than 1 (cscore defaulted to 0).';
    }

    ############################################################################
    # Run dag tools - Prepares DAG for syntenic analysis.
    ############################################################################
    my $dag_file12       = $filtered_blastfile . ".dag";
    my $dag_file12_all   = $dag_file12 . ".all";
    my $query_dup_file   = $opts{query_dup_files};
    my $subject_dup_file = $opts{subject_dup_files};
    my $query            = "a" . $dsgid1;
    my $subject          = "b" . $dsgid2;

    ############################################################################
    # Convert to gene order
    ############################################################################
    my $dag_file12_all_geneorder = "$dag_file12_all.go";
    my $all_file;

    if ( $dagchainer_type eq "geneorder" ) {
        $all_file = $dag_file12_all_geneorder;
        $dag_file12 .= ".go";
    }
    else {
        $all_file = $dag_file12_all;
    }

    #FIXME: This is currently the output produced from the go function.
    $dag_file12 = $all_file;

    ############################################################################
    # Run dagchainer
    ############################################################################
    my ( $dagchainer_file, $merged_dagchainer_file );

    #this is for using dagchainer's merge function
    my $dag_merge_enabled = ( $merge_algo == 2 ) ? 1 : 0;
    my $self_comparision = ( $dsgid1 eq $dsgid2 ) ? 1 : 0;

    #length of a gap (average distance expected between two syntenic genes)
    my $gap = defined( $opts{g} ) ? $opts{g} : floor( $dagchainer_D / 2 );
    $gap = 1 if $gap < 1;

    $dagchainer_file = $dag_file12;
    $dagchainer_file .= "_D$dagchainer_D" if defined $dagchainer_D;
    $dagchainer_file .= "_g$gap"          if defined $gap;
    $dagchainer_file .= "_A$dagchainer_A" if defined $dagchainer_A;
    $dagchainer_file .= "_Dm$Dm"          if $dag_merge_enabled;
    $dagchainer_file .= "_gm$gm"          if $dag_merge_enabled;
    $dagchainer_file .= ".aligncoords";
    $dagchainer_file .= ".ma2.dag"        if $dag_merge_enabled;

    my $post_dagchainer_file;

    if ($dag_merge_enabled) {
        $merged_dagchainer_file = "$dagchainer_file.merged";
        $post_dagchainer_file   = $merged_dagchainer_file;
    }
    else {
        $post_dagchainer_file = $dagchainer_file;
    }

    ############################################################################
    # Run quota align merge
    ############################################################################
    my (
        $find_nearby_time,    $gen_ks_db_time, $dotplot_time,
        $add_gevo_links_time, $final_results_files
    );

    #id 1 is to specify quota align as a merge algo
    if ( $merge_algo == 1 ) {
        $merged_dagchainer_file = "$dagchainer_file.Dm$Dm.ma1";
        $post_dagchainer_file   = $merged_dagchainer_file;
    }

    my $post_dagchainer_file_w_nearby = $post_dagchainer_file;
    $post_dagchainer_file_w_nearby =~ s/aligncoords/all\.aligncoords/;

    #add pairs that were skipped by dagchainer
    $post_dagchainer_file_w_nearby = $post_dagchainer_file;
    ############################################################################
    # Run quota align coverage
    ############################################################################
    my ( $quota_align_coverage, $grimm_stuff, $final_dagchainer_file );

    if ( $depth_algo == 1 )    #id 1 is to specify quota align
    {
        $quota_align_coverage = $post_dagchainer_file_w_nearby;
        $quota_align_coverage .= ".qac" . $depth_org_1_ratio . ".";
        $quota_align_coverage .= $depth_org_2_ratio . "." . $depth_overlap;
        $final_dagchainer_file = $quota_align_coverage;
    }
    else {
        $final_dagchainer_file = $post_dagchainer_file_w_nearby;
    }

    if ( $dagchainer_type eq "geneorder" ) {
        $final_dagchainer_file = $final_dagchainer_file . ".gcoords";
    }

    #generate dotplot images
    unless ($width) {
    	my $chr1_count = $genome1->chromosome_count;
    	my $chr2_count = $genome2->chromosome_count;
        my $max_chr = $chr1_count;
        $max_chr = $chr2_count if $chr2_count > $chr1_count;
        $width   = int( $max_chr * 100 );
        $width   = 1000 if $width > 1000;
        $width   = 500 if $width < 500;
    }

    ############################################################################
    # Create html output directory
    ############################################################################
    my ( $qlead, $slead ) = ( "a", "b" );
    my $out = $org_dirs{ $orgkey1 . "_" . $orgkey2 }{dir} . "/html/";
    mkpath( $out, 0, 0777 ) unless -d $out;
    $out .= "master_";
    my ($base) = $final_dagchainer_file =~ /([^\/]*$)/;
    $out .= $base;
    my $json_basename = "$out";

    $out .= "_ct$color_type" if defined $color_type;
    $out .= ".w$width";

    ############################################################################
    # KS Calculations (Slow and needs to be optimized)
    ############################################################################
    my ( $ks_db, $ks_blocks_file, $svg_file );

    if ($ks_type) {
        my $check_ks = $final_dagchainer_file =~ /^(.*?CDS-CDS)/;
        $check_ks = $final_dagchainer_file =~ /^(.*?protein-protein)/
          unless $check_ks;

        if ($check_ks) {
            $ks_db          = "$final_dagchainer_file.sqlite";
            $ks_blocks_file = "$final_dagchainer_file.ks";
            $svg_file = $ks_blocks_file . ".svg";
        }
        else {
            $warn = "Unable to calculate Ks or Kn values due to at least"
              . " one genome lacking CDS features.";
            $ks_type = undef;
        }
    } else {
        $svg_file = "$final_dagchainer_file.svg";
    }

    ############################################################################
    # Generate dot plot
    ############################################################################
    ($basename) =
      $final_dagchainer_file =~ /([^\/]*aligncoords.*)/;    #.all.aligncoords/;
    $width = 1000 unless defined($width);

    my $dotfile = "$out";
    $dotfile .= ".spa$assemble"       if $assemble;
    $dotfile .= ".gene"               if $axis_metric =~ /gene/i;
    $dotfile .= ".s"                  if $axis_relationship =~ /s/i;
    $dotfile .= ".mcs$min_chr_size"   if $min_chr_size;
    $dotfile .= ".$fid1"              if $fid1;
    $dotfile .= ".$fid2"              if $fid2;
    $dotfile .= ".$ks_type"           if $ks_type;
    $dotfile .= ".box"                if $box_diags;
    $dotfile .= ".flip"               if $flip;
    $dotfile .= ".c0"                 if $clabel eq 0;
    $dotfile .= ".sr"                 if $skip_rand;
    $dotfile .= ".cs$color_scheme"    if defined $color_scheme;
    $dotfile .= ".cso$chr_sort_order" if defined $chr_sort_order;
    $dotfile .= ".min$codeml_min"     if defined $codeml_min;
    $dotfile .= ".max$codeml_max"     if defined $codeml_max;
    $dotfile .= ".log"                if $logks;

    #no syntenic dots, yes, nomicalture is confusing.
    $dotfile .= ".nsd" unless $snsd;

    my $hist = $dotfile . ".hist.png";

    # this would be generated by the DOTPLOT program is Syntenic path assembly
    # was requested
    my $spa_file = $dotfile . ".spa_info.txt";
    my $json_file = "$json_basename.json";
    my $all_json_file = "$json_basename.all.json";
    my $hist_json_file = "$json_basename.datasets.json";

    $out = $dotfile;

    ############################################################################
    # Generate html
    ############################################################################
    my $results =
      HTML::Template->new(
        filename => $config->{TMPLDIR} . 'partials/synmap_results.tmpl' );

    my ($x_label, $y_label);

    if ($clabel) {
        $y_label = "$out.y.png";
        $x_label = "$out.x.png";

        $warn .= qq{Unable to display the y-axis.} unless -r $y_label;
        $warn .= qq{Unable to display the x-axis.} unless -r $x_label;
    }

    my $problem;

    if ( -r "$out.html" ) {

        #Dotplot
        $/ = "\n";
        open( IN, "$out.html" )
          || warn "problem opening $out.html for reading\n";
        $axis_metric = $axis_metric =~ /g/ ? "genes" : "nucleotides";
        $html .=
          "<span class='small'>Axis metrics are in $axis_metric</span><br>";

        #add version of genome to organism names
        $org_name1 .= " (v" . $genome1->version . ")";
        $org_name2 .= " (v" . $genome2->version . ")";

        my $out_url = $out;
        $out_url =~ s/$DIR/$URL/;
        $y_label =~ s/$DIR/$URL/;
        $x_label =~ s/$DIR/$URL/;

        $/ = "\n";
        my $tmp;
        while (<IN>) {
            next if /<\/?html>/;
            $tmp .= $_;
        }
        close IN;
        $tmp =~ s/master.*\.png/$out_url.png/;
        warn "$out_url.html did not parse correctly\n" unless $tmp =~ /map/i;
        $html .= $tmp;

        #Synteny Zoom
        $results->param( codeml_min  => $codeml_min );
        $results->param( codeml_max  => $codeml_max );
        $results->param( axis_metric => $axis_metric );
        $results->param( ylabel      => $y_label );
        $results->param( xlabel      => $x_label );

        if ($flip) {
            $results->param( yorg_name => $org_name1 );
            $results->param( xorg_name => $org_name2 );
        }
        else {
            $results->param( yorg_name => $org_name2 );
            $results->param( xorg_name => $org_name1 );
        }

        $results->param( dotplot   => $tmp );
        $results->param( algorithm => $algo_name );

        if ($hist and $ks_type) {
            if (-r $hist and -s $hist) {
                $results->param( histogram => $out_url . '.hist.png' );
                $results->param( ks_type => $ks_type );
            } else {
                $warn = qq{The histogram was not generated no ks or kn data found.};
            }
        }

        my $final_dagchainer_file_condensed =
          $final_dagchainer_file . ".condensed";
        my $qa_file = $merged_dagchainer_file;
        $qa_file =~ s/\.ma\d$/\.qa/ if $qa_file;
        my $qa_merged_file = $qa_file . ".merged" if $qa_file;
        my $qa_coverage_tmp = $quota_align_coverage . ".tmp"
          if $quota_align_coverage;
        my $qa_coverage_qa = $quota_align_coverage . ".qa"
          if $quota_align_coverage;

        ########################################################################
        # Compress Results
        ########################################################################

        #my $file_list = [
        #    \$raw_blastfile,          \$filtered_blastfile,
        #    \$query_bed,              \$subject_bed,
        #    \$org1_localdups,         \$org2_localdups,
        #    \$dag_file12_all,         \$dag_file12_all_geneorder,
        #    \$dag_file12,             \$dagchainer_file,
        #    \$final_dagchainer_file,  \$final_dagchainer_file_condensed,
        #    \$merged_dagchainer_file, \$qa_file,
        #    \$qa_merged_file,         \$ks_blocks_file,
        #    \$quota_align_coverage,   \$qa_coverage_qa
        #];    #, \$spa_file];
        #my $pm = new Parallel::ForkManager($MAX_PROC);

        #foreach my $item (@$file_list) {
        #    $pm->start and next;

        #    #$$item = CoGe::Accessory::Web::gzip($$item);
        #    $pm->finish;
        #}
        #$pm->wait_all_children();

        #foreach my $item (@$file_list) {

        #    #$$item = CoGe::Accessory::Web::gzip($$item);
        #}

        #######################################################################
        # General
        #######################################################################

        my $log_url = _filename_to_link(
            file     => $cogeweb->logfile,
            msg      => qq{Analysis Log},
            required => 1,
        );

        my $image_url = _filename_to_link(
            file     => "$out.png",
            msg      => qq{Image File},
            required => 1,
        );

        #######################################################################
        # Homologs
        #######################################################################
# mdb removed 9/20/13 issue 77
#        my $fasta1_url = _filename_to_link(
#            file => $fasta1,
#            msg  => qq{Fasta file for $org_name1: $feat_type1}
#        );
#        my $fasta2_url = _filename_to_link(
#            file => $fasta2,
#            msg  => qq{Fasta file for $org_name2: $feat_type2}
#        );

        # mdb added 9/20/13 issue 77
        print_debug(msg => "$feat_type1 $feat_type2");

        my $sequence_url = "api/v1/legacy/sequence"; #"services/JBrowse/service.pl/sequence"; # mdb changed 2/5/15, COGE-289

        my $fasta1_url = _filename_to_link(
            url  => ( $feat_type1 eq "genomic" ? "$sequence_url/$dsgid1" : undef ),
            file => ( $feat_type1 eq "genomic" ? undef : $FASTADIR . "/$dsgid1-$feat_type1.fasta" ),
            msg  => qq{Fasta file for $org_name1: $feat_type1}
        );
        my $fasta2_url = _filename_to_link(
            url  => ( $feat_type2 eq "genomic" ? "$sequence_url/$dsgid2" : undef ),
            file => ( $feat_type2 eq "genomic" ? undef : $FASTADIR . "/$dsgid2-$feat_type2.fasta" ),
            msg  => qq{Fasta file for $org_name2: $feat_type2}
        );

        my $raw_blast_url = _filename_to_link(
            file => $raw_blastfile,
            msg  => qq{Unfiltered $algo_name results}
        );

        my $filtered_blast_url = _filename_to_link(
            file => $filtered_blastfile,
            msg  => qq{Filtered $algo_name results (no tandem duplicates)}
          ),

          my $tandem_dups1_url = _filename_to_link(
            file => $raw_blastfile . '.q.tandems',
            msg  => qq{Tandem Duplicates for $org_name1},
          );

        my $tandem_dups2_url = _filename_to_link(
            file => $raw_blastfile . '.s.tandems',
            msg  => qq{Tandem Duplicates for $org_name2},
        );

        #######################################################################
        # Diagonals
        #######################################################################
        my $dagchainer_input_url = _filename_to_link(
            file => $dag_file12_all,
            msg  => qq{DAGChainer Initial Input file}
        );

        my $geneorder_url = _filename_to_link(
            file => $dag_file12_all_geneorder,
            msg  => qq{DAGChainer Input file converted to gene order}
        );

        my $dagchainer12_url = _filename_to_link(
            file => $dag_file12,
            msg  => qq{DAGChainer Input file post repetitve matches filtered}
        );

        #######################################################################
        # Results
        #######################################################################
        my $dagchainer_url = _filename_to_link(
            file => $dagchainer_file,
            msg  => qq{DAGChainer Output}
        );

        my $merged_dag_url = _filename_to_link(
            file => $merged_dagchainer_file,
            msg  => qq{Merged DAGChainer output}
        );

        #merged dagchainer output is not specified in results.  This hack gets it there.
        if ($merged_dagchainer_file and -r $merged_dagchainer_file) {
            $dagchainer_url .= $merged_dag_url;
        };

        my $quota_align_coverage_url = _filename_to_link(
            file => $quota_align_coverage,
            msg  => qq{Quota Alignment Results}
        );

        my $final_result_url = _filename_to_link(
            file => $final_dagchainer_file,
            msg  => qq{DAGChainer output in genomic coordinates}
        );

        my $ks_blocks_url = _filename_to_link(
            file     => $ks_blocks_file,
            msg      => qq{Results with synonymous/non-synonymous rate values},
            required => $ks_type
        );

        my $final_url = _filename_to_link(
            file => $final_dagchainer_file,
            msg  => qq{Final syntenic gene-set output with GEvo links}
        );

        my $final_condensed_url = _filename_to_link(
            file => $final_dagchainer_file_condensed,
            msg  => qq{Condensed syntelog file with GEvo links}
        );

        my $svg_url = _filename_to_link(
            file => $svg_file,
            msg  => qq{SVG Version of Syntenic Dotplot},
        );

        my $spa_url = _filename_to_link(
            file => $spa_file,
            msg  => qq{Syntenic Path Assembly mapping},
        );

        $spa_url = "" unless $spa_url;

        my $spa_result = "";

        if ($spa_url and $assemble) {
	    #added by EHL: 6/17/15
            #fixing problem where 
	    my $genome1_chr_count = $genome1->chromosome_count();
	    my $genome2_chr_count = $genome2->chromosome_count();
	    my $flip =0;
	    if ($genome1_chr_count >= $genome2_chr_count && $assemble > 0) {
		$flip = 1;
	    }

	    if ($genome1_chr_count <= $genome2_chr_count && $assemble < 0) {
		$flip = 1;
	    }
            $spa_result = $spa_url
                . qq{<a href="#" onclick="coge.synmap.submit_assembly(window.event, '$dagchainer_file', '$dsgid1', '$dsgid2', '$flip');">}
                . qq{Generate Pseudo-Assembled Genomic Sequence}
                . qq{</a>};
        }

        my $json_url = _filename_to_link(
            file => $json_file,
            msg  => qq{Dotplot JSON},
        );

        $dagchainer_file =~ s/^$URL/$DIR/;


        my $rows = [
            {
                general  => 'General',
                homolog  => 'Homolog search',
                diagonal => 'Diagonals',
                result   => 'Results',
            },
            {
                general  => $log_url,
                homolog  => $fasta1_url,
                diagonal => $dagchainer_input_url,
                result   => $dagchainer_url,
            },
            {
                general  => $image_url,
                homolog  => $fasta2_url,
                diagonal => $geneorder_url,
                result   => $final_result_url,
            },
            {
                general  => undef,
                homolog  => $raw_blast_url,
                diagonal => $dagchainer12_url,
                result   => $final_url,
            },
            {
                general  => undef,
                homolog  => $filtered_blast_url,
                diagonal => undef,
                result   => $final_condensed_url,
            },
            {
                general  => undef,
                homolog  => $tandem_dups1_url,
                diagonal => undef,
                result   => $quota_align_coverage_url,
            },
            {
                general  => undef,
                homolog  => $tandem_dups2_url,
                diagonal => undef,
                result   => $ks_blocks_url,
            },
            {
                general  => undef,
                homolog  => undef,
                diagonal => undef,
                result   => $spa_result
,
            },
            {
                general  => undef,
                homolog  => undef,
                diagonal => undef,
                result   => $svg_url,
            },
            {
                general  => undef,
                homolog  => undef,
                diagonal => undef,
                result   => $json_url,
            },
        ];

        $results->param( files => $rows );

        ########################################################################
        # Regenerate Analysis Link - HTML
        ########################################################################

        $results->param( link => $tiny_link );
        if ($ks_type) {
	        my ($ks_file) = $ks_db =~ /([^\/]*)$/;
            my $link = "SynSub.pl?dsgid1=$dsgid1;dsgid2=$dsgid2;file=$ks_file";
            $results->param( synsub => $link );
        }

        if ($grimm_stuff) {
            my $seq1 = ">$org_name1||" . $grimm_stuff->[0];
            $seq1 =~ s/\n/\|\|/g;
            my $seq2 = ">$org_name2||" . $grimm_stuff->[1];
            $seq2 =~ s/\n/\|\|/g;
            $html .= qq{
            <br>
            <span class="ui-button ui-corner-all" id = "grimm_link" onclick="post_to_grimm('$seq1','$seq2')" > Rearrangement Analysis</span> <a class="small" href=http://grimm.ucsd.edu/GRIMM/index.html target=_new>(Powered by GRIMM!)</a>
            };

            my $grimm_data = {
                seq1       => $seq1,
                seq2       => $seq2,
                grimm_link => qq{http://grimm.ucsd.edu/GRIMM/index.html},
            };

            $results->param( grimm => $grimm_data );
        }
        $html .= "<br>";
    }
    else {
        return encode_json({
            error => "The output $out.html could not be found."
        });
    }

    #unless(-r $json_file and -s $json_file ) {
    #    return encode_json({
    #        error => "The json file could not be found."
    #    });
    #}

    #if($ks_type and -s $hist_json_file == 0) {
    #    return encode_json({
    #        error => "The histogram json file could not be found."
    #    });
    #}

    my $log = $cogeweb->logfile;
    $log =~ s/$DIR/$URL/;
    $json_file =~ s/$DIR/$URL/;
    $all_json_file =~ s/$DIR/$URL/;
    $hist_json_file =~ s/$DIR/$URL/;

    $results->param( error    => $problem ) if $problem;
    $results->param( warning  => $warn )    if $warn;
    $results->param( log      => $log );
    #$results->param( json     => $json_file );
    #$results->param( allpairs => $all_json_file );
    #$results->param( hist     => $hist_json_file );
    $results->param( beta     => 1) if $opts{beta};

    ##print out all the datafiles created
    $html .= "<br>";

    $html .=
qq{<span id="clear" style="font-size: 0.8em" class="ui-button ui-corner-all"
        onClick="\$('#results').hide(); \$(this).hide(); \$('#intro').fadeIn();" >Clear Results</span>};

    $html .= qq{</div>};

    $warn = qq{There was a problem running your analysis.}
      . qq{ Please check the log file for details};

    ############################################################################
    # Email results, output benchmark and return results
    ############################################################################

    email_results(
        email    => $email,
        html     => $html,
        org1     => $org_name1,
        org2     => $org_name2,
        jobtitle => $job_title,
        link     => $tiny_link
    ) if $email;

    # Need to remove this from the output from dotplot
    # -- otherwise it over-loads the stuff in the web-page already.
    # -- This can mess up other loaded js such as tablesoter
    my $output = $results->output;

    # erb 5/19/2014 - older version of dotplot would generate javascript inline
    # Remove jquery and xhairs from generated html
    $output =~ s/<script src="\/CoGe\/js\/jquery-1.3.2.js"><\/script>//g;
    $output =~ s/<script src="\/CoGe\/js\/xhairs.js"><\/script>//g;

    #FIXME Return the json representation of the data instead of html
    return encode_json({ html => $output });
}

sub _filename_to_link {
    my %opts = (
        styles   => "link",
        warn     => 1,
        required => 0,
        @_,
    );
    my $file = $opts{file};
    my $url = $opts{url};
    return unless ($file or $url);

    my $link;
    if ( -r $file or $url) {
        if (!$url) {
            $url = $opts{file};
            $url =~ s/$DIR/$URL/;
        }

        $link =
            q{<span class="}
          . $opts{styles} . q{"}
          . q{onclick="window.open('}
          . $url . q{')">}
          . $opts{msg}
          . "</span><br>";

    }
    elsif ( $opts{required} ) {
        $link =
          q{<span class="alert">} . $opts{msg} . q{ (missing)} . q{</span};
    }
    else {

     #        $link = q{<span style="color:dimgray;">} . $opts{msg} . q{</span};
    }
    return $link;
}

################################################################################
# SynMap workflow
################################################################################

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

    if ($cogeweb and $write_log) {
        CoGe::Accessory::Web::write_log( "Generated organism name:",
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( " " x (2) . $title,
            $cogeweb->logfile );
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
    }
    return ( $org_name, $title );
}

sub print_debug {
    my %args = @_;

    if ( defined( $args{enabled} ) && defined( $args{msg} ) && $args{enabled} )
    {
        say STDERR "DEBUG: $args{msg}";
    }
}

# FIXME: Currently this feature is disabled.
# @by Evan Briones
# @on 3/1/2013
sub run_tandem_finder {
    my %opts    = @_;
    my $infile  = $opts{infile};    #dag file produced by dat_tools.py
    my $outfile = $opts{outfile};
    while ( -e "$outfile.running" ) {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    unless ( -r $infile && -s $infile ) {
        CoGe::Accessory::Web::write_log( "", $cogeweb->logfile );
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
    my $cmd = "$PYTHON $TANDEM_FINDER -i '$infile' > '$outfile'";
    system "/usr/bin/touch '$outfile.running'"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "run_tandem_filter: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm '$outfile.running'"
      if -r "$outfile.running";    #remove track file
    return 1 if -r $outfile;
}

#FIXME: Currently this feature is disabled
# @by Evan Briones
# @on 3/1/2013
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
    my $cmd = "$PYTHON $EVAL_ADJUST -c $cvalue '$infile' > '$outfile'";

#There is a parameter that can be passed into this to filter repetitive sequences more or less stringently:
# -c   2 gets rid of more stuff; 10 gets rid of less stuff; default is 4
#consider making this a parameter than can be adjusted from SynMap -- will need to actually play with this value to see how it works
#if implemented, this will require re-naming all the files to account for this parameter
#and updating the auto-SynMap link generator for redoing an analysis

    system "/usr/bin/touch '$outfile.running'"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log(
        "run_adjust_dagchainer_evals: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm '$outfile.running'" if -r "$outfile.running";
    ;      #remove track file
    return 1 if -r $outfile;

}

#FIXME: Currently this feature is disabled
# @by Evan Briones
# @on 6/18/2013
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
      "$PYTHON $FIND_NEARBY --diags='$infile' --all='$dag_all_file' > '$outfile'";
    system "/usr/bin/touch '$outfile.running'"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log( "run find_nearby: running\n\t$cmd",
        $cogeweb->logfile );
    `$cmd`;
    system "/bin/rm '$outfile.running'" if -r "$outfile.running";
    ;      #remove track file
    return 1 if -r $outfile;
}

sub add_reverse_match {
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
        my @line = split /\s+/;
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
    my %params;

    #print STDERR Dumper \%opts;
    $box_diags = $box_diags eq "true" ? 1 : 0;

    # base=8_8.CDS-CDS.blastn.dag_geneorder_D60_g30_A5;

    $params{flip} = $flip       if $flip;
    $params{regen} = $regen     if $regen;
    $params{width} = $width     if $width;
    $params{ksdb} = $ksdb       if $ksdb;
    $params{kstype} = $kstype   if $kstype;
    $params{log} = 1            if $kstype;
    $params{min} = $min         if $min;
    $params{max} = $max         if $max;
    $params{am} = $metric       if defined $metric;
    $params{ar} = $relationship if defined $relationship;
    $params{ct} = $color_type   if $color_type;
    $params{bd} = $box_diags    if $box_diags;
    $params{cs} = $color_scheme if defined $color_scheme;
    $params{fid1} = $fid1       if defined $fid1 && $fid1 =~ /^\d+$/;
    $params{fid2} = $fid2       if defined $fid2 && $fid2 =~ /^\d+$/;

    $url = url_for("run_dotplot.pl", %params) . "&" . $url;
    my $ua = LWP::UserAgent->new;
    $ua->timeout(10);
    my $response = $ua->get($url);

    unless ($response->is_success) {
        return "Unable to get image for dotplot: $url";
    }

    my $content = $response->decoded_content;

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

sub generate_assembly {
    my %opts = @_;
    my $gid1 = $opts{gid1};
    my $gid2 = $opts{gid2};
    my $flip = $opts{flip}; #reverse the order from default

    unless ($gid1 and $gid2) {
        return encode_json({
            success => JSON::false,
            error => "Wrong number of genomes specified"
         });
    }

    # Swap order if gid2 < gid1
    ($gid1, $gid2) = ($gid2, $gid1) if $gid2 lt $gid1;

    my $genome1 = $coge->resultset("Genome")->find($gid1);
    my $genome2 = $coge->resultset("Genome")->find($gid2);

    unless ($USER->has_access_to_genome($genome1) and
            $USER->has_access_to_genome($genome2)) {

        return encode_json({
            success => JSON::false,
            error => "User does not have access to this dataset"
        });
    }

    unless ($opts{input} and -r $opts{input}) {
        return encode_json({
            success => JSON::false,
            error => "The syntenic path assembly could not be found"
        });
    }

    my $filename = qq($gid1-$gid2-) . md5_hex($opts{input}) .".".$flip. ".tar.gz";
    my $output = catfile($DIAGSDIR, $gid1, $gid2, "assembly", $filename);

    # Submit workflow
    my $submission = generate_pseudo_assembly($JEX, $config, $opts{input}, $output, $flip);
    $output =~ s/$DIR/$URL/;

    # Fixup success to return true or false
    return encode_json({
        id => $submission->{id},
        output => $output,
        success => $submission->{success} ? JSON::true : JSON::false,
    });
}
