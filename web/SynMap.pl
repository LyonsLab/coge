#!/usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0);

use CoGeX;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Workflow;
use CoGe::Accessory::Web;
use CoGe::Algos::KsCalc;

use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use DBIxProfiler;
use Data::Dumper;
use HTML::Template;
use LWP::Simple;
use Parallel::ForkManager;
use GD;
use Digest::MD5 qw(md5_base64);
use File::Path;
use Mail::Mailer;
use Benchmark;
use DBI;
use POSIX;
use Sort::Versions;

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
    $YERBA,         $GENE_ORDER,     $PAGE_TITLE,  $KSCALC,
    $GEN_FASTA,     $RUN_ALIGNMENT,  $RUN_COVERAGE);

$DEBUG     = 1;
$YERBA     = CoGe::Accessory::Jex->new(host => "localhost", port => 5151);
$P         = CoGe::Accessory::Web::get_defaults($ENV{HOME} . 'coge.conf');
$ENV{PATH} = join ":",
  (
    $P->{COGEDIR}, $P->{BINDIR}, $P->{BINDIR} . "SynMap",
    "/usr/bin", "/usr/local/bin");
$ENV{BLASTDB}    = $P->{BLASTDB};
$ENV{BLASTMAT}   = $P->{BLASTMATRIX};
$ENV{PYTHONPATH} = "/opt/apache/CoGe/bin/dagchainer_bp";

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
    $P->{MULTI_LAST}
  . " -a $MAX_PROC --path="
  . $P->{LAST_PATH}
  . " --dbpath="
  . $P->{LASTDB};
$GZIP          = $P->{GZIP};
$GUNZIP        = $P->{GUNZIP};
$KSCALC        = $P->{KSCALC};
$GEN_FASTA     = $P->{GEN_FASTA};
$RUN_ALIGNMENT = $P->{RUN_ALIGNMENT};
$RUN_COVERAGE  = $P->{RUN_COVERAGE};

#in the web form, each sequence search algorithm has a unique number.  This table identifies those and adds appropriate options
$ALGO_LOOKUP = {
    0 => {
        algo => $BLASTN . " -task megablast",    #megablast
        opt             => "MEGA_SELECT",  #select option for html template file
        filename        => "megablast",
        displayname     => "MegaBlast",
        html_select_val => 0,
        formatdb        => 1,},
    1 => {
        algo     => $BLASTN . " -task dc-megablast",   #discontinuous megablast,
        opt      => "DCMEGA_SELECT",
        filename => "dcmegablast",
        displayname     => "Discontinuous MegaBlast",
        html_select_val => 1,
        formatdb        => 1,},
    2 => {
        algo            => $BLASTN . " -task blastn",    #blastn
        opt             => "BLASTN_SELECT",
        filename        => "blastn",
        displayname     => "BlastN",
        html_select_val => 2,
        formatdb        => 1,},
    3 => {
        algo            => $TBLASTX,                     #tblastx
        opt             => "TBLASTX_SELECT",
        filename        => "tblastx",
        displayname     => "TBlastX",
        html_select_val => 3,
        formatdb        => 1,},
    4 => {
        algo            => $LASTZ,                       #lastz
        opt             => "LASTZ_SELECT",
        filename        => "lastz",
        displayname     => "(B)lastZ",
        html_select_val => 4,},
    5 => {
        algo            => $BLASTP . " -task blastp",    #blastn
        opt             => "BLASTP_SELECT",
        filename        => "blastp",
        displayname     => "BlastP",
        html_select_val => 5,
        formatdb        => 1,},
    6 => {
        algo            => $LAST,                        #last
        opt             => "LAST_SELECT",
        filename        => "last",
        displayname     => "Last",
        html_select_val => 6,},};

$DATADIR  = $P->{DATADIR};
$DIAGSDIR = $P->{DIAGSDIR};
$FASTADIR = $P->{FASTADIR};

mkpath($TEMPDIR,     1, 0777);
mkpath($FASTADIR,    1, 0777);
mkpath($DIAGSDIR,    1, 0777);                           # mdb added 7/9/12
mkpath($P->{LASTDB}, 1, 0777);                           # mdb added 7/9/12
$BLASTDBDIR = $P->{BLASTDB};

$PYTHON        = $P->{PYTHON};                           #this was for python2.5
$PYTHON26      = $P->{PYTHON};
$DAG_TOOL      = $P->{DAG_TOOL};
$BLAST2BED     = $P->{BLAST2BED};
$GENE_ORDER    = $DIR . "/bin/SynMap/gene_order.py";
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

$DOTPLOT     = $P->{DOTPLOT} . " -cf " . $ENV{HOME} . 'coge.conf';
$SVG_DOTPLOT = $P->{SVG_DOTPLOT};

#$CONVERT_TO_GENE_ORDER = $DIR."/bin/SynMap/convert_to_gene_order.pl";
#$NWALIGN = $DIR."/bin/nwalign-0.3.0/bin/nwalign";
$NWALIGN = $P->{NWALIGN};
$|       = 1;                           # turn off buffering
$DATE    = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub {($_[5] + 1900, $_[4] + 1, $_[3]), $_[2], $_[1], $_[0]}
      ->(localtime));
$FORM       = new CGI;
$PAGE_TITLE = "SynMap";
$PAGE_NAME  = "$PAGE_TITLE.pl";

my %ajax = CoGe::Accessory::Web::ajax_func();

#$ajax{read_log}=\&read_log_test;
$DBNAME = $P->{DBNAME};
$DBHOST = $P->{DBHOST};
$DBPORT = $P->{DBPORT};
$DBUSER = $P->{DBUSER};
$DBPASS = $P->{DBPASS};
$connstr =
  "dbi:mysql:dbname=" . $DBNAME . ";host=" . $DBHOST . ";port=" . $DBPORT;
$coge = CoGeX->connect($connstr, $DBUSER, $DBPASS);

$COOKIE_NAME = $P->{COOKIE_NAME};

my ($cas_ticket) = $FORM->param('ticket');
$USER = undef;
($USER) = CoGe::Accessory::Web->login_cas(
    cookie_name => $COOKIE_NAME,
    ticket      => $cas_ticket,
    coge        => $coge,
    this_url    => $FORM->url()) if ($cas_ticket);
($USER) = CoGe::Accessory::LogUser->get_user(
    cookie_name => $COOKIE_NAME,
    coge        => $coge) unless $USER;

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
    %ajax,);

my $pj = new CGI::Ajax(%FUNCTIONS);
if ($FORM->param('jquery_ajax')) {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};

    #print STDERR Dumper \%args;
    if ($fname and defined $FUNCTIONS{$fname}) {
        if ($args{args}) {
            my @args_list = split(/,/, $args{args});
            print $FORM->header, $FUNCTIONS{$fname}->(@args_list);
        } else {
            print $FORM->header, $FUNCTIONS{$fname}->(%args);
        }
    }
} else {
    $pj->js_encode_function('escape');
    print $pj->build_html($FORM, \&gen_html);

    #   print $FORM->header; print gen_html();
}

################################################################################
# Web functions
################################################################################

sub read_log_test
{
    my %args    = @_;
    my $logfile = $args{logfile};
    my $prog    = $args{prog};
    return unless $logfile;
    $logfile .= ".log" unless $logfile =~ /log$/;
    $logfile = $TEMPDIR . "/$logfile" unless $logfile =~ /^$TEMPDIR/;
    return unless -r $logfile;
    my $str;
    open(IN, $logfile);

    while (<IN>) {
        $str .= $_;
    }
    close IN;
    return $str;
}

sub gen_html
{
    my $html;
    my ($body) = gen_body();
    my $template =
      HTML::Template->new(filename => $P->{TMPLDIR} . 'generic_page.tmpl');
    $template->param(PAGE_TITLE => 'SynMap');
    $template->param(TITLE      => 'Whole Genome Synteny');
    $template->param(HEAD       => qq{});
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER => $name);

    $template->param(LOGON => 1) unless $USER->user_name eq "public";
    $template->param(DATE => $DATE);

    #$template->param(ADJUST_BOX=>1);
    $template->param(LOGO_PNG => "SynMap-logo.png");
    $template->param(BODY     => $body);
    $template->param(HELP     => "/wiki/index.php?title=SynMap");
    $html .= $template->output;
    return $html;
}

sub gen_body
{
    my $form = shift || $FORM;
    my $template =
      HTML::Template->new(filename => $P->{TMPLDIR} . 'SynMap.tmpl');

    $template->param(MAIN => 1);

    #$template->param( EMAIL       => $USER->email )  if $USER->email;

    my $master_width = $FORM->param('w') || 0;
    $template->param(MWIDTH => $master_width);

    #set search algorithm on web-page
    if (defined($FORM->param('b'))) {
        $template->param($ALGO_LOOKUP->{$FORM->param('b')}{opt} => "selected");
    } else {
        $template->param($ALGO_LOOKUP->{6}{opt} => "selected");
    }
    my ($D, $A, $Dm, $gm, $dt, $dupdist, $cscore);
    $D       = $FORM->param('D');
    $A       = $FORM->param('A');
    $Dm      = $FORM->param('Dm');
    $gm      = $FORM->param('gm');
    $gm      = 40 unless defined $gm;
    $dt      = $FORM->param('dt');
    $cscore  = $FORM->param('csco');
    $dupdist = $FORM->param('tdd');

#   $cvalue = $FORM->param('c');       #different c value than the one for cytology.  But if you get that, you probably shouldn't be reading this code

    my $display_dagchainer_settings;
    if ($D && $A && $dt) {
        my $type;
        if ($dt =~ /gene/i) {
            $type = " genes";
            $template->param('DAG_GENE_SELECT' => 'checked');
        } else {
            $type = " bp";
            $template->param('DAG_DISTANCE_SELECT' => 'checked');
        }
        $display_dagchainer_settings =
          qq{display_dagchainer_settings([$D,$A, '$gm', $Dm],'$type');};
    } else {
        $template->param('DAG_GENE_SELECT' => 'checked');
        $display_dagchainer_settings = qq{display_dagchainer_settings();};
    }

    #   $cvalue = 4 unless defined $cvalue;
    #   $template->param( 'CVALUE'                      => $cvalue );
    $dupdist = 10 unless defined $dupdist;
    $template->param('DUPDIST' => $dupdist);
    $template->param('CSCORE' => $cscore) if defined $cscore;
    $template->param(
        'DISPLAY_DAGCHAINER_SETTINGS' => $display_dagchainer_settings);
    $template->param('MIN_CHR_SIZE' => $FORM->param('mcs'))
      if $FORM->param('mcs');

    #will the program automatically run?
    my $autogo = $FORM->param('autogo');
    $autogo = 0 unless defined $autogo;
    $template->param(AUTOGO => $autogo);

    #populate organism menus
    for (my $i = 1; $i <= 2; $i++) {
        my $dsgid = $form->param('dsgid' . $i) || 0;
        my $feattype_param = $FORM->param('ft' . $i) if $FORM->param('ft' . $i);
        my $name = $FORM->param('name' . $i) if $FORM->param('name' . $i);
        my $org_menu = gen_org_menu(
            dsgid          => $dsgid,
            num            => $i,
            feattype_param => $feattype_param,
            name           => $name);
        $template->param("ORG_MENU" . $i => $org_menu);
    }

    #set ks for coloring syntenic dots
    if ($FORM->param('ks')) {
        if ($FORM->param('ks') eq 1) {
            $template->param(KS1 => "selected");
        } elsif ($FORM->param('ks') eq 2) {
            $template->param(KS2 => "selected");
        } elsif ($FORM->param('ks') eq 3) {
            $template->param(KS3 => "selected");
        }
    } else {
        $template->param(KS0 => "selected");
    }

    #set color_scheme
    my $cs = 1;
    $cs = $FORM->param('cs') if defined $FORM->param('cs');
    $template->param("CS" . $cs => "selected");

    #set codeml min and max
    my $codeml_min;
    $codeml_min = $FORM->param('cmin') if defined $FORM->param('cmin');
    my $codeml_max;
    $codeml_max = $FORM->param('cmax') if defined $FORM->param('cmax');
    $template->param('CODEML_MIN' => $codeml_min) if defined $codeml_min;
    $template->param('CODEML_MAX' => $codeml_max) if defined $codeml_max;
    my $logks;
    $logks = $FORM->param('logks') if defined $FORM->param('logks');
    $logks = 1 unless defined $logks;    #turn on by default if not specified
    $template->param('LOGKS' => "checked") if defined $logks && $logks;

    #show non syntenic dots:  on by default
    my $snsd = 0;
    $snsd = $FORM->param('snsd') if (defined $FORM->param('snsd'));
    $template->param('SHOW_NON_SYN_DOTS' => 'checked') if $snsd;

    #are the axes flipped?
    my $flip = 0;
    $flip = $FORM->param('flip') if (defined $FORM->param('flip'));
    $template->param('FLIP' => 'checked') if $flip;

    #are the chromosomes labeled?
    my $clabel = 1;
    $clabel = $FORM->param('cl') if (defined $FORM->param('cl'));
    $template->param('CHR_LABEL' => 'checked') if $clabel;

    #are the chromosomes labeled?
    my $skip_rand = 1;
    $skip_rand = $FORM->param('sr') if (defined $FORM->param('sr'));
    $template->param('SKIP_RAND' => 'checked') if $skip_rand;

    #what is the sort for chromosome display?
    my $chr_sort_order = "S";
    $chr_sort_order = $FORM->param('cso') if (defined $FORM->param('cso'));
    if ($chr_sort_order =~ /N/i) {
        $template->param('CHR_SORT_NAME' => 'selected');
    } elsif ($chr_sort_order =~ /S/i) {
        $template->param('CHR_SORT_SIZE' => 'selected');
    }

    #set axis metric for dotplot
    if ($FORM->param('ct')) {
        if ($FORM->param('ct') eq "inv") {
            $template->param('COLOR_TYPE_INV' => 'selected');
        } elsif ($FORM->param('ct') eq "diag") {
            $template->param('COLOR_TYPE_DIAG' => 'selected');
        }
    } else {
        $template->param('COLOR_TYPE_NONE' => 'selected');
    }
    if ($FORM->param('am') && $FORM->param('am') =~ /g/i) {
        $template->param('AXIS_METRIC_GENE' => 'selected');
    } else {
        $template->param('AXIS_METRIC_NT' => 'selected');
    }

    #axis relationship:  will dotplot be forced into a square?
    if ($FORM->param('ar') && $FORM->param('ar') =~ /s/i) {
        $template->param('AXIS_RELATIONSHIP_S' => 'selected');
    } else {
        $template->param('AXIS_RELATIONSHIP_R' => 'selected');
    }

    #merge diags algorithm
    if ($FORM->param('ma')) {
        $template->param(QUOTA_MERGE_SELECT => 'selected')
          if $FORM->param('ma') eq "1";
        $template->param(DAG_MERGE_SELECT => 'selected')
          if $FORM->param('ma') eq "2";
    }
    if ($FORM->param('da')) {
        if ($FORM->param('da') eq "1") {
            $template->param(QUOTA_ALIGN_SELECT => 'selected');
        }
    }
    my $depth_org_1_ratio = 1;
    $depth_org_1_ratio = $FORM->param('do1') if $FORM->param('do1');
    $template->param(DEPTH_ORG_1_RATIO => $depth_org_1_ratio);
    my $depth_org_2_ratio = 1;
    $depth_org_2_ratio = $FORM->param('do2') if $FORM->param('do2');
    $template->param(DEPTH_ORG_2_RATIO => $depth_org_2_ratio);
    my $depth_overlap = 40;
    $depth_overlap = $FORM->param('do') if $FORM->param('do');
    $template->param(DEPTH_OVERLAP => $depth_overlap);

    $template->param('BOX_DIAGS' => "checked") if $FORM->param('bd');
    my $spa = $FORM->param('sp') if $FORM->param('sp');
    $template->param('SYNTENIC_PATH' => "checked") if $spa;
    $template->param('SHOW_NON_SYN' => "checked") if $spa && $spa =~ /2/;
    $template->param('SPA_FEW_SELECT'  => "selected") if $spa && $spa > 0;
    $template->param('SPA_MORE_SELECT' => "selected") if $spa && $spa < 0;
    my $file = $form->param('file');
    if ($file) {
        my $results = read_file($file);
        $template->param(RESULTS => $results);
    }

#place to store fids that are passed into SynMap to highlight that pair in the dotplot (if present)
    my $fid1 = 0;
    $fid1 = $FORM->param('fid1') if $FORM->param('fid1');

    $template->param('FID1' => $fid1);
    my $fid2 = 0;
    $fid2 = $FORM->param('fid2') if $FORM->param('fid2');
    $template->param('FID2'      => $fid2);
    $template->param('PAGE_NAME' => $PAGE_NAME);
    $template->param('TEMPDIR'   => $TEMPDIR);
    return $template->output;
}

sub gen_org_menu
{
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
      HTML::Template->new(filename => $P->{TMPLDIR} . 'SynMap.tmpl');
    $menu_template->param(ORG_MENU   => 1);
    $menu_template->param(NUM        => $num);
    $menu_template->param('ORG_NAME' => $name);
    $menu_template->param('ORG_DESC' => $desc);
    $menu_template->param(
        'ORG_LIST' => get_orgs(name => $name, i => $num, oid => $oid));
    my ($dsg_menu) = gen_dsg_menu(oid => $oid, dsgid => $dsgid, num => $num);
    $menu_template->param(DSG_MENU => $dsg_menu);

    if ($dsgid) {
        my ($dsg_info, $feattype_menu, $message) = get_genome_info(
            dsgid    => $dsgid,
            org_num  => $num,
            feattype => $feattype_param);
        $menu_template->param(DSG_INFO       => $dsg_info);
        $menu_template->param(FEATTYPE_MENU  => $feattype_menu);
        $menu_template->param(GENOME_MESSAGE => $message);
    }
    return $menu_template->output;
}

sub gen_dsg_menu
{
    my $t1    = new Benchmark;
    my %opts  = @_;
    my $oid   = $opts{oid};
    my $num   = $opts{num};
    my $dsgid = $opts{dsgid};
    my @dsg_menu;
    my $message;
    my $org_name;

    #   print STDERR"here!\n";

    foreach my $dsg (
        $coge->resultset('Genome')->search(
            {organism_id => $oid}, {prefetch => ['genomic_sequence_type']})
      ) {
        my $name;
        my $has_cds = 0;
        if ($dsg->restricted && !$USER->has_access_to_genome($dsg)) {
            next unless $dsgid && $dsg->id == $dsgid;
            $name = "Restricted";
        } else {
            $name .= $dsg->name . ": " if $dsg->name;
            $name .=
              $dsg->type->name . " (v" . $dsg->version . ",id" . $dsg->id . ")";
            $org_name = $dsg->organism->name unless $org_name;
            foreach my $ft (
                $coge->resultset('FeatureType')->search({
                        genome_id            => $dsg->id,
                        'me.feature_type_id' => 3},
                    {
                        join => {features => {dataset => 'dataset_connectors'}},
                        rows => 1,})
              ) {
                $has_cds = 1;
            }
        }
        push @dsg_menu, [$dsg->id, $name, $dsg, $has_cds];

    }

    my $dsg_menu = qq{
   <select id=dsgid$num onChange="\$('#dsg_info$num').html('<div class=dna_small class=loading class=small>loading. . .</div>'); get_genome_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],[handle_dsg_info])">
};
    foreach (
        sort {
                 versioncmp($b->[2]->version, $a->[2]->version)
              || $a->[2]->type->id <=> $b->[2]->type->id
              || $b->[3] cmp $a->[3]
        } @dsg_menu
      ) {
        my ($numt, $name) = @$_;
        my $selected = " selected" if $dsgid && $numt == $dsgid;
        $selected = " " unless $selected;
        $numt = 0 if $name eq "Restricted";
        $dsg_menu .= qq{
   <OPTION VALUE=$numt $selected>$name</option>
};
    }
    $dsg_menu .= "</select>";
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2, $t1));

    #    print STDERR qq{
    #-----------------
    #sub gen_dsg_menu runtime:  $time
    #-----------------
    #};
    return ($dsg_menu, $message);

}

sub read_file
{
    my $file = shift;

    my $html;
    open(IN, $TEMPDIR . $file) || die "can't open $file for reading: $!";
    while (<IN>) {
        $html .= $_;
    }
    close IN;
    return $html;
}

sub get_orgs
{
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
    } elsif ($name) {
        @db =
          $coge->resultset("Organism")
          ->search({name => {like => "%" . $name . "%"}});
    } elsif ($desc) {
        @db =
          $coge->resultset("Organism")
          ->search({description => {like => "%" . $desc . "%"}});
    } else {
        $org_count = $coge->resultset("Organism")->count;
    }

    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db) {
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
    unless (@opts && ($name || $desc)) {
        $html .= qq{<input type = hidden name="org_id$i" id="org_id$i">};
        return $html;
    }

    $html .=
      qq{<SELECT id="org_id$i" SIZE="5" MULTIPLE onChange="get_genome_info_chain($i)" >\n};
    $html .= join("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/ unless $oid;
    return $html;
}

sub get_genome_info
{
    my $t1       = new Benchmark;
    my %opts     = @_;
    my $dsgid    = $opts{dsgid};
    my $org_num  = $opts{org_num};
    my $feattype = $opts{feattype};
    $feattype = 1 unless defined $feattype;
    return " ", " ", " " unless $dsgid;
    my $html_dsg_info;

#my ($dsg) = $coge->resultset("Genome")->find({genome_id=>$dsgid},{join=>['organism','genomic_sequences'],prefetch=>['organism','genomic_sequences']});
    my ($dsg) = $coge->resultset("Genome")->find({genome_id => $dsgid},
        {join => 'organism', prefetch => 'organism'});
    return " ", " ", " " unless $dsg;
    my $org     = $dsg->organism;
    my $orgname = $org->name;
    $orgname =
        "<a href=\"OrganismView.pl?oid="
      . $org->id
      . "\" target=_new>$orgname</a>";
    my $org_desc;
    if ($org->description) {
        $org_desc = join(
            "; ",
            map {
                    qq{<span class=link onclick="\$('#org_desc}
                  . qq{$org_num').val('$_').focus();search_bar('org_desc$org_num'); timing('org_desc$org_num')">$_</span>}
              } split /\s*;\s*/,
            $org->description);
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
        $chr_count,  $plasmid,    $contig,    $scaffold) = get_dsg_info($dsg);
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
    if ($percent_gc > 0) {
        $html_dsg_info .=
          "<tr><td>DNA content: <td>GC: $percent_gc%, AT: $percent_at%, N: $percent_n%, X: $percent_x%";
    } else {
        $html_dsg_info .=
          qq{<tr><td>DNA content: <td id=gc_content$org_num class='link' onclick="get_gc($dsgid, 'gc_content$org_num')">Click to retrieve};
    }
    $html_dsg_info .= "<tr><td>Total length: <td>" . commify($chr_length);
    $html_dsg_info .= "<tr><td>Contains plasmid" if $plasmid;
    $html_dsg_info .= "<tr><td>Contains contigs" if $contig;
    $html_dsg_info .= "<tr><td>Contains scaffolds" if $scaffold;
    $html_dsg_info .= "</table>";
    if ($dsg->restricted && !$USER->has_access_to_genome($dsg)) {
        $html_dsg_info = "Restricted";
    }
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2, $t1));

    #    print STDERR qq{
    #-----------------
    #sub get_genome_info runtime:  $time
    #-----------------
    #};

    my $message;

    #create feature type menu
    my $has_cds;

    foreach my $ft (
        $coge->resultset('FeatureType')->search({
                genome_id            => $dsg->id,
                'me.feature_type_id' => 3},
            {
                join => {features => {dataset => 'dataset_connectors'}},
                rows => 1,})
      ) {
        $has_cds = 1;
    }

    my ($cds_selected, $genomic_selected) = (" ", " ");
    $cds_selected     = "selected" if $feattype eq 1 || $feattype eq "CDS";
    $genomic_selected = "selected" if $feattype eq 2 || $feattype eq "genomic";

    my $feattype_menu =
      qq{<select id="feat_type$org_num" name ="feat_type$org_num">#};
    $feattype_menu .= qq{<OPTION VALUE=1 $cds_selected>CDS</option>}
      if $has_cds;
    $feattype_menu .= qq{<OPTION VALUE=2 $genomic_selected>genomic</option>};
    $feattype_menu .= "</select>";
    $message = "<span class='small alert'>No Coding Sequence in Genome</span>"
      unless $has_cds;

    return $html_dsg_info, $feattype_menu, $message, $chr_length, $org_num,
      $dsg->organism->name, $dsg->genomic_sequence_type_id;
}

sub get_previous_analyses
{

    #FIXME:  THis whole sub needs updating or removal!  Lyons 6/12/13
    my %opts = @_;
    my $oid1 = $opts{oid1};
    my $oid2 = $opts{oid2};
    return unless $oid1 && $oid2;
    my ($org1) = $coge->resultset('Organism')->find($oid1);
    my ($org2) = $coge->resultset('Organism')->find($oid2);
    return
      if ($USER->user_name =~ /public/i
        && ($org1->restricted || $org2->restricted));
    my ($org_name1) = $org1->name;
    my ($org_name2) = $org2->name;
    ($oid1, $org_name1, $oid2, $org_name2) =
      ($oid2, $org_name2, $oid1, $org_name1)
      if ($org_name2 lt $org_name1);

    my $tmp1 = $org_name1;
    my $tmp2 = $org_name2;
    foreach my $tmp ($tmp1, $tmp2) {
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
    if (-d $dir) {
        opendir(DIR, $dir);
        while (my $file = readdir(DIR)) {
            $sqlite = 1 if $file =~ /sqlite$/;
            next unless $file =~ /\.aligncoords/;    #$/\.merge$/;
            my ($D, $g, $A) = $file =~ /D(\d+)_g(\d+)_A(\d+)/;
            my ($Dm) = $file =~ /Dm(\d+)/;
            my ($gm) = $file =~ /gm(\d+)/;
            my ($ma) = $file =~ /ma(\d+)/;
            $Dm = " " unless defined $Dm;
            $gm = " " unless defined $gm;
            $ma = 0   unless $ma;
            my $merge_algo;
            $merge_algo = "DAGChainer" if $ma && $ma == 2;

            if ($ma && $ma == 1) {
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
            next unless ($D && $g && $A);

            my ($blast) = $file =~
              /^[^\.]+\.[^\.]+\.([^\.]+)/;    #/blastn/ ? "BlastN" : "TBlastX";
            my $select_val;
            foreach my $item (values %$ALGO_LOOKUP) {
                if ($item->{filename} eq $blast) {
                    $blast      = $item->{displayname};
                    $select_val = $item->{html_select_val};
                }
            }
            my ($dsgid1, $dsgid2, $type1, $type2) =
              $file =~ /^(\d+)_(\d+)\.(\w+)-(\w+)/;
            $type1 = "CDS" if $type1 eq "protein";
            $type2 = "CDS" if $type2 eq "protein";

            #           print STDERR $file,"\n";
            #           my ($repeat_filter) = $file =~ /_c(\d+)/;
            next unless ($dsgid1 && $dsgid2 && $type1 && $type2);
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
                select_val => $select_val);
            my $geneorder = $file =~ /\.go/;
            my $dsg1 = $coge->resultset('Genome')->find($dsgid1);
            next unless $dsg1;
            next if ($dsg1->restricted && !$USER->has_access_to_genome($dsg1));
            my ($ds1) = $dsg1->datasets;
            my $dsg2 = $coge->resultset('Genome')->find($dsgid2);
            next unless $dsg2;
            next if ($dsg2->restricted && !$USER->has_access_to_genome($dsg2));
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
      . join("<TH>",
        qw(Org1 Genome1 Ver1 Genome%20Type1 Sequence%20Type1 Org2 Genome2 Ver2 Genome%20Type2 Sequence%20type2 Algo Dist%20Type Dup%20Dist Ave%20Dist(g) Max%20Dist(D) Min%20Pairs(A))
      ) . "</THEAD><TBODY>\n";
    my %seen;

    foreach my $item (
        sort {$b->{dsgid1} <=> $a->{dsgid1} || $b->{dsgid2} <=> $a->{dsgid2}}
        @items) {
        my $val = join("_",
            $item->{g},          $item->{D},       $item->{A},
            $oid1,               $item->{dsgid1},  $item->{type1},
            $oid2,               $item->{dsgid2},  $item->{type2},
            $item->{select_val}, $item->{dagtype}, $item->{tdd});
        next if $seen{$val};
        $seen{$val} = 1;
        $prev_table .=
          qq{<TR class=feat onclick="update_params('$val')" align=center><td>};
        my $ver1 = $item->{dsg1}->version;
        $ver1 = "0" . $ver1 if $ver1 =~ /^\./;
        my $ver2 = $item->{dsg2}->version;
        $ver2 = "0" . $ver2 if $ver2 =~ /^\./;
        $prev_table .= join("<td>",
            $item->{dsg1}->organism->name, $item->{genome1},
            $ver1,                         $item->{dsg1}->type->name,
            $item->{type_name1},           $item->{dsg2}->organism->name,
            $item->{genome2},              $ver2,
            $item->{dsg2}->type->name,     $item->{type_name2},
            $item->{blast},                $item->{dagtype},
            $item->{tdd},                  $item->{g},
            $item->{D},                    $item->{A})
          . "\n";
    }
    $prev_table .= qq{</TBODY></table>};
    $html .= $prev_table;
    $html .=
      "<br><span class=small>Synonymous substitution rates previously calculated</span>"
      if $sqlite;
    return "$html";
}

sub get_dsg_info
{
    my $dsg       = shift;
    my $length    = 0;
    my $chr_count = 0;

    $length = $dsg->length;    #$rs->first->get_column('total_length');
    $chr_count = $dsg->genomic_sequences->count();
    my ($gc, $at, $n, $x) = (0, 0, 0, 0);
    if ($chr_count < 100 && $length < 50000000) {
        ($gc, $at, $n, $x) = get_dsg_gc(dsg => $dsg);
    }
    my ($plasmid, $contig, $scaffold) = get_chr_types(dsg => $dsg);
    return $gc, $at, $n, $x, $length, $chr_count, $plasmid, $contig, $scaffold;
}

sub get_dsg_gc
{
    my %opts  = @_;
    my $dsg   = $opts{dsg};
    my $dsgid = $opts{dsgid};
    my $text  = $opts{text};
    $dsg = $coge->resultset('Genome')->find($dsgid) if $dsgid;
    my ($gc, $at, $n, $x) = $dsg->percent_gc;
    $gc *= 100;
    $at *= 100;
    $n  *= 100;
    $x  *= 100;

    if ($text) {
        return "GC: $gc%, AT: $at%, N: $n%, X: $x%";
    } else {
        return ($gc, $at, $n, $x);
    }
}

sub get_chr_types
{
    my %opts  = @_;
    my $dsg   = $opts{dsg};
    my $dsgid = $opts{dsgid};
    $dsg = $coge->resultset('Genome')->find($dsgid) if $dsgid;
    my $plasmid  = 0;
    my $contig   = 0;
    my $scaffold = 0;
    my @gs       = $dsg->genomic_sequences;
    if (@gs > 100) {
        return (0, 1, 0);
    }
    foreach my $chr (map {$_->chromosome} @gs) {
        $plasmid  = 1 if !$plasmid  && $chr =~ /plasmid/i;
        $contig   = 1 if !$contig   && $chr =~ /contig/i;
        $scaffold = 1 if !$scaffold && $chr =~ /scaffold/i;
    }
    return ($plasmid, $contig, $scaffold);
}

sub get_pair_info
{
    my @anno;
    foreach my $fid (@_) {
        unless ($fid =~ /^\d+$/) {
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
      . join("\n", (map {"<tr><td>" . $_ . "</td></tr>"} @anno))
      . "</table>";
    my $URL = $P->{URL};
    $output =~ s/window\.open\('(.*?)'\)/window.open('$URL$1')/g;
    return $output;
}

sub get_query_link
{
    my %url_options = @_;

    #  print STDERR Dumper \%url_options;

    my $dagchainer_D = $url_options{D};

 #    my $dagchainer_g = $url_options{g}; #depreciated -- will be a factor of -D
    my $dagchainer_A = $url_options{A};
    my $Dm           = $url_options{Dm};
    my $gm           = $url_options{gm};
    ($Dm) = $Dm =~ /(\d+)/;
    ($gm) = $gm =~ /(\d+)/;

#   my $repeat_filter_cvalue = $url_options{c}; #parameter to be passed to run_adjust_dagchainer_evals
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
        if    ($ks_type eq "ks")    {$num = 1;}
        elsif ($ks_type eq "kn")    {$num = 2;}
        elsif ($ks_type eq "kn_ks") {$num = 3;}
        $synmap_link .= ";ks=$num";
    }
    $synmap_link .= ";am=g" if $axis_metric       && $axis_metric       =~ /g/i;
    $synmap_link .= ";ar=s" if $axis_relationship && $axis_relationship =~ /s/i;
    $synmap_link .= ";ct=$color_type" if $color_type;

    my ($org_name1, $titleA) =
      gen_org_name(dsgid => $dsgid1, feat_type => $feat_type1, write_log => 0);
    my ($org_name2, $titleB) =
      gen_org_name(dsgid => $dsgid2, feat_type => $feat_type2, write_log => 0);
    my $log_msg =
        "<span class=link onclick=window.open('OrganismView.pl?dsgid="
      . $dsgid1
      . "')>$org_name1</span> v. <span class=link"
      . "onclick=window.open('OrganismView.pl?dsgid=$dsgid2')>$org_name2</span>";

    $log_msg .= " Ks" if $ks_type;
    my $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR);

    return CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $synmap_link,
        log_msg => $log_msg);
}

sub generate_basefile
{
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(tempdir => $TEMPDIR);
    return $cogeweb->basefilename;
}

sub go
{
    my %opts = @_;
    foreach my $k (keys %opts) {
        $opts{$k} =~ s/^\s+//;
        $opts{$k} =~ s/\s+$//;
    }
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};

    unless ($dsgid1 && $dsgid2) {
        return "<span class=alert>You must select two genomes.</span>";
    }

    my ($dsg1) = $coge->resultset('Genome')->find($dsgid1);
    my ($dsg2) = $coge->resultset('Genome')->find($dsgid2);

    unless ($dsg1 && $dsg2) {
        return
          "<span class=alert>Problem generating dataset group objects for ids:  $dsgid1, $dsgid2.</span>";
    }

    unless ($dsg1 && $dsg2) {
        return
          "<span class=alert>Problem generating one of the genome objects for id1: $dsgid1 or id2: $dsgid2</span>";
    }
    my ($dir1, $dir2) = sort ($dsgid1, $dsgid2);
    my $workflow_id = "$dir1-$dir2";

    ############################################################################
    # Initialize Jobs
    ############################################################################

    my $tiny_link = $opts{tiny_link};
    say STDERR "tiny_link is required for logging." unless defined($tiny_link);

    my $job = CoGe::Accessory::Web::get_job(
        tiny_link => $tiny_link,
        title     => $PAGE_TITLE,
        user_id   => $USER->id,
        db_object => $coge);

    my ($tiny_id) = $tiny_link =~ /\/(\w+)$/;
    $workflow_id .= "-$tiny_id";

    my $basename = $opts{basename};
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basename,
        tempdir  => $TEMPDIR);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Link to Regenerate Analysis",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("$tiny_link", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

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
    my $dupdist = defined($opts{tdd}) ? $opts{tdd} : 10;

    # dotplot options
    my $regen_images      = ($opts{regen_images} eq "true") ? 1 : 0;
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

    my $email = 0 if check_address_validity($opts{email}) eq 'invalid';

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
    my ($org_name1, $title1) = gen_org_name(
        dsgid     => $dsgid1,
        feat_type => $feat_type1,
        write_log => 1);

    my ($org_name2, $title2) = gen_org_name(
        dsgid     => $dsgid2,
        feat_type => $feat_type2,
        write_log => 1);

    ############################################################################
    # Generate Fasta files
    ############################################################################
    my $t0 = new Benchmark;
    my ($fasta1, $fasta2);
    my $workflow = undef;
    my $status   = undef;

    my $config = $ENV{HOME} . "coge.conf";

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
    $workflow = $YERBA->create_workflow(
        name    => "synmap-$workflow_id",
        logfile => $cogeweb->logfile);

    if ($feat_type1 eq "genomic") {
        my $genome = $coge->resultset('Genome')->find($dsgid1);
        $fasta1 = $genome->file_path;
    } else {
        my @fasta1args = ();
        $fasta1 = $FASTADIR . "/$dsgid1-$feat_type1.fasta";
        push @fasta1args, ["--config",       $config,     0];
        push @fasta1args, ["--genome_id",    $dsgid1,     1];
        push @fasta1args, ["--feature_type", $feat_type1, 1];
        push @fasta1args, ["--fasta",        $fasta1,     1];

        $workflow->add_job(
            cmd     => $GEN_FASTA,
            script  => undef,
            args    => \@fasta1args,
            inputs  => undef,
            outputs => [$fasta1]);
    }

    if ($feat_type2 eq "genomic") {
        my $genome = $coge->resultset('Genome')->find($dsgid2);
        $fasta2 = $genome->file_path;
    } else {
        $fasta2 = $FASTADIR . "/$dsgid2-$feat_type2.fasta";
        my @fasta2args = ();
        push @fasta2args, ["--config",       $config,     0];
        push @fasta2args, ["--genome_id",    $dsgid2,     1];
        push @fasta2args, ["--feature_type", $feat_type2, 1];
        push @fasta2args, ["--fasta",        $fasta2,     1];

        $workflow->add_job(
            cmd     => $GEN_FASTA,
            script  => undef,
            args    => \@fasta2args,
            inputs  => undef,
            outputs => [$fasta2]);
    }

    # Sort by genome id
    (
        $dsgid1,     $dsg1,              $org_name1,  $fasta1,
        $feat_type1, $depth_org_1_ratio, $dsgid2,     $dsg2,
        $org_name2,  $fasta2,            $feat_type2, $depth_org_2_ratio)
      = (
        $dsgid2,     $dsg2,              $org_name2,  $fasta2,
        $feat_type2, $depth_org_2_ratio, $dsgid1,     $dsg1,
        $org_name1,  $fasta1,            $feat_type1, $depth_org_1_ratio
      ) if ($org_name2 lt $org_name1);

    ############################################################################
    # Generate blastdb files
    ############################################################################
    my ($blastdb, @blastdb_files);

    if ($ALGO_LOOKUP->{$blast}{formatdb}) {
        my $write_log = 0;
        my $basename  = "$BLASTDBDIR/$dsgid2-$feat_type2";

        my @blastdbargs = ();
        push @blastdbargs, ['-p', $feat_type2 eq "protein" ? "T" : "F", 1];
        push @blastdbargs, ['-i', $fasta2, 1];
        push @blastdbargs, ['-t', '"' . $org_name2 . '"', 1];
        push @blastdbargs, ['-n', $basename, 1];

        $blastdb = $basename;
        $basename .= $feat_type2 eq "protein" ? ".p" : ".n";

        push @blastdb_files, "$basename" . "sq";
        push @blastdb_files, "$basename" . "in";
        push @blastdb_files, "$basename" . "hr";

        $workflow->add_job(
            cmd     => $FORMATDB,
            script  => undef,
            args    => \@blastdbargs,
            inputs  => [$fasta2],
            outputs => \@blastdb_files);

        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("Generating BlastDB file",
            $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("blastdb file: $blastdb",
            $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("#" x (20) . "\n",
            $cogeweb->logfile);
        CoGe::Accessory::Web::write_log(
            "creating blastdb for *"
                . $org_name2
                . "* ($blastdb): $FORMATDB",
            $cogeweb->logfile);
    } else {
        $blastdb = $fasta2;
        push @blastdb_files, $blastdb;
    }

    my $html;

    my ($orgkey1, $orgkey2) = ($title1, $title2);
    my %org_dirs = (
        $orgkey1 . "_"
          . $orgkey2 => {
            fasta    => $fasta1,
            db       => $blastdb,
            basename => $dsgid1 . "_"
              . $dsgid2
              . ".$feat_type1-$feat_type2."
              . $ALGO_LOOKUP->{$blast}{filename},
            dir => "$DIAGSDIR/$dir1/$dir2",},);

    foreach my $org_dir (keys %org_dirs) {
        my $outfile = $org_dirs{$org_dir}{dir};
        mkpath($outfile, 0, 0777) unless -d $outfile;
        warn "didn't create path $outfile: $!" unless -d $outfile;
        $outfile .= "/" . $org_dirs{$org_dir}{basename};
        $org_dirs{$org_dir}{blastfile} = $outfile;    #.".blast";
    }

    ############################################################################
    # Run Blast
    ############################################################################

    #check blast runs for problems;  Not forked in order to keep variables
    my $problem       = 0;
    my $raw_blastfile = $org_dirs{$orgkey1 . "_" . $orgkey2}{blastfile};

    foreach my $key (keys %org_dirs) {
        my $cmd =
          $ALGO_LOOKUP->{$blast}
          {algo};    #$prog =~ /tblastx/i ? $TBLASTX : $BLASTN;

        my $fasta   = $org_dirs{$key}{fasta};
        my $db      = $org_dirs{$key}{db};
        my $outfile = $org_dirs{$key}{blastfile};
        my @blastargs;

        if ($cmd =~ /lastz/i) {
            push @blastargs, ["-i", $fasta,   1];
            push @blastargs, ["-d", $db,      1];
            push @blastargs, ["-o", $outfile, 1];
        } elsif ($cmd =~ /last_wrapper/i) {
            push @blastargs, ["",   $db,      1];
            push @blastargs, ["",   $fasta,   1];
            push @blastargs, ["-o", $outfile, 1];
        } else {
            push @blastargs, ["-out",   $outfile, 1];
            push @blastargs, ["-query", $fasta,   1];
            push @blastargs, ["-db",    $db,      1];
        }

        (undef, $cmd) = CoGe::Accessory::Web::check_taint($cmd);
        push @blastdb_files, $fasta;

        $workflow->add_job(
            cmd     => $cmd,
            script  => undef,
            args    => \@blastargs,
            inputs  => \@blastdb_files,
            outputs => [$outfile]);
    }

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Running genome comparison",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log(
        "Running " . $ALGO_LOOKUP->{$blast}{displayname},
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Completed blast run(s)",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    my $t1 = new Benchmark;
    my $blast_time = timestr(timediff($t1, $t0));

    ###########################################################################
    # Converting blast to bed and finding local duplications
    ###########################################################################
    #
    # NOTES: The blast2bed program will take the input rawblast file and
    # filter it and creating a new rawblast and moving the old to
    # rawblast.orig
    #
    my $t2          = new Benchmark;
    my $query_bed   = $raw_blastfile . ".q.bed";
    my $subject_bed = $raw_blastfile . ".s.bed";

    my @blastargs = ();
    push @blastargs, ['-infile',   $raw_blastfile, 1];
    push @blastargs, ['-outfile1', $query_bed,     1];
    push @blastargs, ['-outfile2', $subject_bed,   1];

    my @bedoutputs = ();
    push @bedoutputs, $query_bed;
    push @bedoutputs, $subject_bed;
    push @bedoutputs, $raw_blastfile if ($raw_blastfile =~ /genomic/);
    push @bedoutputs, "$raw_blastfile.orig" if ($raw_blastfile =~ /genomic/);

    $workflow->add_job(
        cmd     => $BLAST2BED,
        script  => undef,
        args    => \@blastargs,
        inputs  => [$raw_blastfile],
        outputs => \@bedoutputs);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Creating .bed files", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("running: $BLAST2BED", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    ###########################################################################
    # Converting blast to raw and finding local duplications
    ###########################################################################
    my $filtered_blastfile = $raw_blastfile;
    $filtered_blastfile .= ".tdd$dupdist";
    $filtered_blastfile .= ".cs$cscore" if $cscore;
    $filtered_blastfile .= ".filtered";

    my @rawargs = ();
    push @rawargs, ["", $raw_blastfile, 1];
    push @rawargs, ["--localdups", "", 1];
    push @rawargs, ["--qbed",        $query_bed,          1];
    push @rawargs, ["--sbed",        $subject_bed,        1];
    push @rawargs, ["--tandem_Nmax", $dupdist,            1];
    push @rawargs, ["--cscore",      $cscore,             1] if $cscore;
    push @rawargs, [">",             $filtered_blastfile, 1];

    my @rawoutputs = ();
    push @rawoutputs, $filtered_blastfile;
    push @rawoutputs, "$raw_blastfile.q.localdups";
    push @rawoutputs, "$raw_blastfile.s.localdups";

    $workflow->add_job(
        cmd     => $BLAST2RAW,
        script  => undef,
        args    => \@rawargs,
        inputs  => [$raw_blastfile, $query_bed, $subject_bed],
        outputs => \@rawoutputs);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Filtering results of tandem duplicates",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("finding and removing local duplications",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("running: $BLAST2RAW", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

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
    my $local_dup_time = timestr(timediff($t2, $t1));

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
    push @dagtoolargs, ['-q', $query,              1];
    push @dagtoolargs, ['-s', $subject,            1];
    push @dagtoolargs, ['-b', $filtered_blastfile, 1];
    push @dagtoolargs, ['-c', "", 1];

    push @dagtoolargs, ['--query_dups', $query_dup_file, 1] if $query_dup_file;
    push @dagtoolargs, ['--subject_dups', $subject_dup_file, 1]
      if $subject_dup_file;
    push @dagtoolargs, ['>', $dag_file12_all, 1];

    $workflow->add_job(
        cmd     => $DAG_TOOL,
        script  => undef,
        args    => \@dagtoolargs,
        inputs  => [$filtered_blastfile],
        outputs => [$dag_file12_all]);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log(
        "Converting blast file to dagchainer input file",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("run dag_tools:\nrunning: $DAG_TOOL",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    my $t2_5 = new Benchmark;
    my $dag_tool_time = timestr(timediff($t2_5, $t2));

    ############################################################################
    # Convert to gene order
    ############################################################################
    my $dag_file12_all_geneorder = "$dag_file12_all.go";
    my $all_file;

    if ($dagchainer_type eq "geneorder") {
        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log(
            "Converting dagchainer input into gene order coordinates",
            $cogeweb->logfile);

        my @geneorderargs = ();
        push @geneorderargs, ["",           $dag_file12_all,           1];
        push @geneorderargs, ["",           $dag_file12_all_geneorder, 1];
        push @geneorderargs, ["--gid1",     $dsgid1,                   1];
        push @geneorderargs, ["--gid2",     $dsgid2,                   1];
        push @geneorderargs, ["--feature1", $feat_type1,               1];
        push @geneorderargs, ["--feature2", $feat_type2,               1];

        $workflow->add_job(
            cmd     => $GENE_ORDER,
            script  => undef,
            args    => \@geneorderargs,
            inputs  => [$dag_file12_all],
            outputs => [$dag_file12_all_geneorder]);

        my $msg = "running coversion to gene order for $dag_file12_all";
        CoGe::Accessory::Web::write_log("run dagchainer: $GENE_ORDER",
            $cogeweb->logfile);
        CoGe::Accessory::Web::write_log($msg, $cogeweb->logfile);
        $msg = "Completed conversion of gene order to file";
        $msg .= " $dag_file12_all_geneorder";
        CoGe::Accessory::Web::write_log($msg, $cogeweb->logfile);

        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

        $all_file = $dag_file12_all_geneorder;
        $dag_file12 .= ".go";
    } else {
        $all_file = $dag_file12_all;
    }

    my $t3 = new Benchmark;
    my $convert_to_gene_order_time = timestr(timediff($t3, $t2_5));

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

    my $t3_5 = new Benchmark;

    #FIXME: This step has a race condition with work_queue

    # This step will fail if the dag_file_all is larger than the system memory
    # limit. If this file does not exist, let's send a warning to the log file
    # and continue with the analysis using the dag_all file
    unless ((-r $dag_file12 && -s $dag_file12)
        || (-r $dag_file12 . ".gz" && -s $dag_file12 . ".gz")) {
        $dag_file12 = $all_file;
        CoGe::Accessory::Web::write_log(
            "WARNING: sub run_adjust_dagchainer_evals failed. "
              . "Perhaps due to Out of Memory error. "
              . "Proceeding without this step!",
            $cogeweb->logfile);
    }
    my $run_adjust_eval_time = timestr(timediff($t3_5, $t3));

    ############################################################################
    # Run dagchainer
    ############################################################################
    my ($dagchainer_file, $merged_dagchainer_file);

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Running DagChainer", $cogeweb->logfile);

    #this is for using dagchainer's merge function
    my $dag_merge_enabled = ($merge_algo == 2) ? 1 : 0;
    my $self_comparision = ($dsgid1 eq $dsgid2) ? 1 : 0;

    #length of a gap (average distance expected between two syntenic genes)
    my $gap = defined($opts{g}) ? $opts{g} : floor($dagchainer_D / 2);

    $dagchainer_file = $dag_file12;
    $dagchainer_file .= "_D$dagchainer_D" if $dagchainer_D;
    $dagchainer_file .= "_g$gap"          if $gap;
    $dagchainer_file .= "_A$dagchainer_A" if $dagchainer_A;
    $dagchainer_file .= "_Dm$Dm"          if $dag_merge_enabled;
    $dagchainer_file .= "_gm$gm"          if $dag_merge_enabled;
    $dagchainer_file .= ".aligncoords";
    $dagchainer_file .= ".ma2.dag"        if $dag_merge_enabled;

    my @dagargs = ();
    push @dagargs, ["-E", "0.05", 1];
    push @dagargs, ["-i",   $dag_file12,   1];
    push @dagargs, ["-D",   $dagchainer_D, 1] if $dagchainer_D;
    push @dagargs, ["-g",   $gap,          1] if $gap;
    push @dagargs, ["-A",   $dagchainer_A, 1] if $dagchainer_A;
    push @dagargs, ["--Dm", $Dm,           1] if $dag_merge_enabled;
    push @dagargs, ["--m",  $gm,           1] if $dag_merge_enabled;
    push @dagargs, ["--new_behavior", "", 1] if $self_comparision;

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
        print_debug(msg => "DAG-Merge is disabled.", enabled => $DEBUG);
        $merged_dagchainer_file = "$dagchainer_file.merged";
        push @dagargs, ["--merge", $merged_dagchainer_file, 1];

        $workflow->add_job(
            cmd     => $RUN_DAGCHAINER,
            script  => undef,
            args    => \@dagargs,
            inputs  => [$dag_file12],
            outputs => [$merged_dagchainer_file]);
        $post_dagchainer_file = $merged_dagchainer_file;

    } else {
        push @dagargs, [">", $dagchainer_file, 1];

        $workflow->add_job(
            cmd     => $RUN_DAGCHAINER,
            script  => undef,
            args    => \@dagargs,
            inputs  => [$dag_file12],
            outputs => [$dagchainer_file]);

        $post_dagchainer_file = $dagchainer_file;
    }

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Completed dagchainer run",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("run dagchainer\nrunning: $RUN_DAGCHAINER",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    my $t4 = new Benchmark;
    my $run_dagchainer_time = timestr(timediff($t4, $t3_5));

    ############################################################################
    # Run quota align merge
    ############################################################################
    my (
        $find_nearby_time,    $gen_ks_db_time, $dotplot_time,
        $add_gevo_links_time, $final_results_files);

    #id 1 is to specify quota align as a merge algo
    if ($merge_algo == 1) {
        $merged_dagchainer_file = "$dagchainer_file.Dm$Dm.ma1";

        my @mergeargs = ();
        push @mergeargs, ['--config',       $config,                 0];
        push @mergeargs, ['--infile',       $dagchainer_file,        1];
        push @mergeargs, ['--outfile',      $merged_dagchainer_file, 1];
        push @mergeargs, ['--max_distance', $Dm,                     1];

        $workflow->add_job(
            cmd     => $RUN_ALIGNMENT,
            script  => undef,
            args    => \@mergeargs,
            inputs  => [$dagchainer_file],
            outputs => [$merged_dagchainer_file]);

        $post_dagchainer_file = $merged_dagchainer_file;
    }
    my $post_dagchainer_file_w_nearby = $post_dagchainer_file;
    $post_dagchainer_file_w_nearby =~ s/aligncoords/all\.aligncoords/;

    print_debug(
        msg     => "Post dag chainer file: $post_dagchainer_file",
        enabled => $DEBUG);

    #add pairs that were skipped by dagchainer
    $post_dagchainer_file_w_nearby = $post_dagchainer_file;

    #TODO: Run find nearby is currently disabled
    #run_find_nearby(
    #    infile       => $post_dagchainer_file,
    #    dag_all_file => $all_file,
    #    outfile      => $post_dagchainer_file_w_nearby
    #);    #program is not working correctly.

    my $t5 = new Benchmark;
    $find_nearby_time = timestr(timediff($t5, $t4));

    ############################################################################
    # Run quota align coverage
    ############################################################################
    my ($quota_align_coverage, $grimm_stuff, $final_dagchainer_file);

    if ($depth_algo == 1)    #id 1 is to specify quota align
    {
        $quota_align_coverage = $post_dagchainer_file_w_nearby;
        $quota_align_coverage .= ".qac" . $depth_org_1_ratio . ".";
        $quota_align_coverage .= $depth_org_2_ratio . "." . $depth_overlap;

        print_debug(msg => $post_dagchainer_file_w_nearby, enabled => $DEBUG);

        my @depthargs = ();
        push @depthargs, ['--config',  $config,                        0];
        push @depthargs, ['--infile',  $post_dagchainer_file_w_nearby, 1];
        push @depthargs, ['--outfile', $quota_align_coverage,          1];
        push @depthargs, ['--depth_ratio_org1', $depth_org_1_ratio, 1];
        push @depthargs, ['--depth_ratio_org2', $depth_org_2_ratio, 1];
        push @depthargs, ['--depth_overlap',    $depth_overlap,     1];

        $workflow->add_job(
            cmd     => $RUN_COVERAGE,
            script  => undef,
            args    => \@depthargs,
            inputs  => [$post_dagchainer_file_w_nearby],
            outputs => [$quota_align_coverage]);

        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("Running Quota Align",
            $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("", $cogeweb->logfile);
        $final_dagchainer_file = $quota_align_coverage;
    } else {
        $final_dagchainer_file = $post_dagchainer_file_w_nearby;
    }

    print_debug("Final dag chainer file: $final_dagchainer_file",
        enabled => $DEBUG);

    if ($dagchainer_type eq "geneorder") {
        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log(
            "Converting gene order coordinates back to genomic coordinates",
            $cogeweb->logfile);

        my @positionargs = ();
        push @positionargs, ['', $final_dagchainer_file, 1];
        push @positionargs, ['', "$final_dagchainer_file.gcoords", 1];
        push @positionargs, ["--positional", '', 1];

        $workflow->add_job(
            cmd     => $GENE_ORDER,
            script  => undef,
            args    => \@positionargs,
            inputs  => [$final_dagchainer_file],
            outputs => ["$final_dagchainer_file.gcoords"]);

        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

        $final_dagchainer_file = $final_dagchainer_file . ".gcoords";
    }

    #generate dotplot images
    my ($org1_length, $org2_length, $chr1_count, $chr2_count) = (0) x 4;

    foreach my $gs ($dsg1->genomic_sequences) {
        $chr1_count++;
        $org1_length += $gs->sequence_length;
    }

    foreach my $gs ($dsg2->genomic_sequences) {
        $chr2_count++;
        $org2_length += $gs->sequence_length;
    }

    unless ($width) {
        my $max_chr = $chr1_count;
        $max_chr = $chr2_count if $chr2_count > $chr1_count;
        $width   = int($max_chr * 100);
        $width   = 1000 if $width > 1000;
        $width   = 500 if $width < 500;
    }

    ############################################################################
    # Create html output directory
    ############################################################################
    my ($qlead, $slead) = ("a", "b");
    my $out = $org_dirs{$orgkey1 . "_" . $orgkey2}{dir} . "/html/";
    mkpath($out, 0, 0777) unless -d $out;
    $out .= "master_";
    my ($base) = $final_dagchainer_file =~ /([^\/]*$)/;
    $out .= $base;
    $out .= "_ct$color_type" if defined $color_type;
    $out .= ".w$width";

    ############################################################################
    # KS Calculations (Slow and needs to be optimized)
    ############################################################################
    my ($ks_db, $ks_blocks_file, $svg_file, $warn);

    if ($ks_type) {
        my $check_ks = $final_dagchainer_file =~ /^(.*?CDS-CDS)/;
        $check_ks = $final_dagchainer_file =~ /^(.*?protein-protein)/ unless $check_ks;

        if ($check_ks) {
            $ks_db          = "$final_dagchainer_file.sqlite";
            $ks_blocks_file = "$final_dagchainer_file.ks";

            my @ksargs = ();
            push @ksargs, ['--config',    $config,                0];
            push @ksargs, ['--infile',    $final_dagchainer_file, 1];
            push @ksargs, ['--dbfile',    $ks_db,                 1];
            push @ksargs, ['--blockfile', $ks_blocks_file,        1];

            $workflow->add_job(
                cmd     => $KSCALC,
                script  => undef,
                args    => \@ksargs,
                inputs  => [$final_dagchainer_file],
                outputs => [$ks_blocks_file, $ks_db]);

            ####################################################################
            # Generate svg dotplot
            ####################################################################

            my @svgargs = ();
            push @svgargs, ['--dag_file', $ks_blocks_file, 1];
            push @svgargs, ['--flip', "", 1] if $flip;
            push @svgargs, ['--xhead', '"' . $org_name1 . '"', 1] if $org_name1;
            push @svgargs, ['--yhead', '"' . $org_name2 . '"', 1] if $org_name2;
            push @svgargs, ['--output', $ks_blocks_file, 1];

            $svg_file = $ks_blocks_file . ".svg";
            $workflow->add_job(
                cmd     => $SVG_DOTPLOT,
                script  => undef,
                args    => \@svgargs,
                inputs  => [$ks_blocks_file],
                outputs => [$svg_file]);

            CoGe::Accessory::Web::write_log("generate svg dotplot: $SVG_DOTPLOT",
                $cogeweb->logfile);
            CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
            CoGe::Accessory::Web::write_log("", $cogeweb->logfile);
        } else {
            $warn = "Unable to calculate Ks or Kn values due to at least"
            ." one genome lacking CDS features.";
            $ks_type = undef;
        }
    }

    my $t6 = new Benchmark;
    $gen_ks_db_time = timestr(timediff($t6, $t5));

    ############################################################################
    # Generate dot plot
    ############################################################################
    my @plotargs   = ();
    my @plotinputs = ();

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Generating dotplot", $cogeweb->logfile);

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

        push @plotargs, ['-ksdb', $ks_db,   1];
        push @plotargs, ['-kst',  $ks_type, 1];
        push @plotargs, ['-log',  $logks,   1];
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
        push @plotargs, ['-d', $dag_file12_all, 0];
        push @plotinputs, $dag_file12_all;
    } else {
        $dotfile .= ".nsd";    #no syntenic dots, yes, nomicalture is confusing.
    }

    push @plotargs, ['-a', $final_dagchainer_file, 1];
    push @plotargs, ['-b', $dotfile, 1];

    my $jsoption = "";
    $jsoption .= qq{'javascript:synteny_zoom(};
    $jsoption .= qq{"$dsgid1",};
    $jsoption .= qq{"$dsgid2",};
    $jsoption .= qq{"$basename",};
    $jsoption .= $flip ? qq{"YCHR","XCHR"} : qq{"XCHR","YCHR"};
    $jsoption .= qq{,"$ks_db"} if $ks_db;
    $jsoption .= qq{)'};

    push @plotargs, ['-l',    $jsoption, 0];
    push @plotargs, ['-dsg1', $dsgid1,   1];
    push @plotargs, ['-dsg2', $dsgid2,   1];
    push @plotargs, ['-w',    $width,    1];
    push @plotargs, ['-lt', 2, 1];
    push @plotargs, ['-assemble', $assemble,    1] if $assemble;
    push @plotargs, ['-am',       $axis_metric, 1] if $axis_metric;
    push @plotargs, ['-fb', '', 1]
      if $axis_relationship && $axis_relationship =~ /s/;
    push @plotargs, ['-mcs', $min_chr_size, 1] if $min_chr_size;
    push @plotargs, ['-cdt', $color_type,   1] if $color_type;
    push @plotargs, ['-bd', 1, 1] if $box_diags;
    push @plotargs, ['-fid1', $fid1, 1] if $fid1;
    push @plotargs, ['-fid2', $fid2, 1] if $fid2;
    push @plotargs, ['-f',      1, 1] if $flip;
    push @plotargs, ['-labels', 0, 1] if $clabel eq 0;
    push @plotargs, ['-sr',     1, 1] if $skip_rand;
    push @plotargs, ['-color_scheme', $color_scheme, 1]
      if defined $color_scheme;
    push @plotargs, ['-chr_sort_order', $chr_sort_order, 1]
      if defined $chr_sort_order;
    push @plotargs, ['-min', $codeml_min, 1] if defined $codeml_min;
    push @plotargs, ['-max', $codeml_max, 1] if defined $codeml_max;

    push @plotinputs, $final_dagchainer_file;

    my $hist = $dotfile . ".hist.png";

    # this would be generated by the DOTPLOT program is Syntenic path assembly
    # was requested
    my $spa_file = $dotfile . ".spa_info.txt";
    my @plotoutputs =
      ("$dotfile.html", "$dotfile.png", "$dotfile.x.png", "$dotfile.y.png");
    push @plotoutputs, $hist     if $ks_db;
    push @plotoutputs, $spa_file if $assemble;

    $workflow->add_job(
        cmd       => $DOTPLOT,
        script    => undef,
        args      => \@plotargs,
        inputs    => \@plotinputs,
        outputs   => \@plotoutputs,
        overwrite => $regen_images);

    CoGe::Accessory::Web::write_log("generate dotplot: running\n\t$DOTPLOT",
        $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    $status = $YERBA->submit_workflow($workflow);
    $YERBA->wait_for_completion($workflow->name);
    $out = $dotfile;

    my $t7 = new Benchmark;
    $dotplot_time = timestr(timediff($t7, $t6));

    ############################################################################
    # Generate html
    ############################################################################

    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Adding GEvo links to final output files",
        $cogeweb->logfile);
    add_GEvo_links(
        infile => $final_dagchainer_file,
        dsgid1 => $dsgid1,
        dsgid2 => $dsgid2);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);
    my $t8 = new Benchmark;
    $add_gevo_links_time = timestr(timediff($t8, $t7));

    if (-r "$out.html") {
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
        open(IN, "$out.html") || warn "problem opening $out.html for reading\n";
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

        if (-r $y_lab) {
            $y_lab =~ s/$DIR/$URL/;
            $html .= qq{ <div><img src="$y_lab"></div><td> };
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
        if (-r $x_lab) {
            $x_lab =~ s/$DIR/$URL/;
            $html .= qq{ <br><img src="$x_lab"> };
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
        $html .= qq{<TR valign=top><td>Homolog search<td>Diagonals<td>Results};
        $html .= qq{<tr valign=top><td>};

  #       $html .= qq{<span class='small link' onclick=window.open('')></span>};
        $fasta1 =~ s/$DIR/$URL/;
        $html .=
          qq{<span class='small link' onclick=window.open('$fasta1')>Fasta file for $org_name1: $feat_type1</span><br>};
        $fasta2 =~ s/$DIR/$URL/;
        $html .=
          qq{<span class='small link' onclick=window.open('$fasta2')>Fasta file for $org_name2: $feat_type2</span><br>};
        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("Processing Tandem Duplicate File",
            $cogeweb->logfile);
        my $org1_localdups = process_local_dups_file(
            infile  => $raw_blastfile . ".q.localdups",
            outfile => $raw_blastfile . ".q.tandems");
        my $org2_localdups = process_local_dups_file(
            infile  => $raw_blastfile . ".s.localdups",
            outfile => $raw_blastfile . ".s.tandems");
        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

        my $final_dagchainer_file_condensed =
          $final_dagchainer_file . ".condensed";
        my $qa_file = $merged_dagchainer_file;
        $qa_file =~ s/\.ma\d$/\.qa/;
        my $qa_merged_file  = $qa_file . ".merged";
        my $qa_coverage_tmp = $quota_align_coverage . ".tmp"
          if $quota_align_coverage;
        my $qa_coverage_qa = $quota_align_coverage . ".qa"
          if $quota_align_coverage;

        ########################################################################
        # Compress Results
        ########################################################################

        my $file_list = [
            \$raw_blastfile,          \$filtered_blastfile,
            \$query_bed,              \$subject_bed,
            \$org1_localdups,         \$org2_localdups,
            \$dag_file12_all,         \$dag_file12_all_geneorder,
            \$dag_file12,             \$dagchainer_file,
            \$final_dagchainer_file,  \$final_dagchainer_file_condensed,
            \$merged_dagchainer_file, \$qa_file,
            \$qa_merged_file,         \$ks_blocks_file,
            \$quota_align_coverage,   \$qa_coverage_qa];    #, \$spa_file];
        my $pm = new Parallel::ForkManager($MAX_PROC);

        foreach my $item (@$file_list) {
            $pm->start and next;

            #$$item = CoGe::Accessory::Web::gzip($$item);
            $pm->finish;
        }
        $pm->wait_all_children();

        foreach my $item (@$file_list) {

            #$$item = CoGe::Accessory::Web::gzip($$item);
        }

        $raw_blastfile      =~ s/$DIR/$URL/;
        $filtered_blastfile =~ s/$DIR/$URL/;
        $org1_localdups     =~ s/$DIR/$URL/;
        $org2_localdups     =~ s/$DIR/$URL/;
        $dag_file12_all     =~ s/$DIR/$URL/;
        $dag_file12         =~ s/$DIR/$URL/;
        $dagchainer_file    =~ s/$DIR/$URL/;

        $html .=
          "<span class='link small' onclick=window.open('$raw_blastfile')>Unfiltered $algo_name results</span><br>";
        $html .=
          "<span class='link small' onclick=window.open('$filtered_blastfile')>Filtered $algo_name results (no tandem duplicates)</span><br>";
        $html .=
          "<span class='link small' onclick=window.open('$org1_localdups')>Tandem Duplicates for $org_name1</span><br>";
        $html .=
          "<span class='link small' onclick=window.open('$org2_localdups')>Tandem Duplicates for $org_name2</span><br>";
        $html .= "<td>";
        $html .=
          qq{<span class='small link' onclick=window.open('$dag_file12_all')>DAGChainer Initial Input file</span><br>};

        if ($dag_file12_all_geneorder && -r $dag_file12_all_geneorder) {
            $dag_file12_all_geneorder =~ s/$DIR/$URL/;
            $html .=
              qq{<span class='small link' onclick=window.open('$dag_file12_all_geneorder')>DAGChainer Input file converted to gene order</span><br>};
        }
        $html .=
          qq{<span class='small link' onclick=window.open('$dag_file12')>DAGChainer Input file post repetivive matches filter</span><br>};
        $html .= "<td>";
        $html .=
          qq{<span class='small link' onclick=window.open('$dagchainer_file')>DAGChainer output</span>};

        if (-r $merged_dagchainer_file) {
            $merged_dagchainer_file =~ s/$DIR/$URL/;
            $html .=
              qq{<br><span class='small link' onclick=window.open('$merged_dagchainer_file')>Merged DAGChainer output</span>};
        }

#       $final_dagchainer_file =~ s/$DIR/$URL/;
#       if ($final_dagchainer_file=~/gcoords/)
#         {
#       $post_dagchainer_file_w_nearby =~ s/$DIR/$URL/;
#       $html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results with nearby genes added</span><br>};
#       $html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results converted back to genomic coordinates</span>};
#         }
#       else
#         {
#       $post_dagchainer_file_w_nearby =~ s/$DIR/$URL/;
#       $html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results with nearby genes added</span>};
#         }
        if ($quota_align_coverage && -r $quota_align_coverage) {
            $quota_align_coverage =~ s/$DIR/$URL/;
            $html .=
              qq{<br><span class='small link' onclick=window.open('$quota_align_coverage')>Quota Alignment Results</span>};
        }

        if ($final_dagchainer_file =~ /gcoords/) {

            #my $tmp= $final_dagchainer_file;
            $final_dagchainer_file =~ s/$DIR/$URL/;

            #$tmp =~ s/\.gcoords//;

#$html .= qq{<br><span class='small link' onclick=window.open('$tmp')>DAGChainer output in gene coordinates</span>};
            $html .=
              qq{<br><span class='small link' onclick=window.open('$final_dagchainer_file')>DAGChainer output in genomic coordinates</span>};

#$html .= qq{<span class='small link' onclick=window.open('$post_dagchainer_file_w_nearby')>Results converted back to genomic coordinates</span>};
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
        if (-r $final_dagchainer_file_condensed) {
            $final_dagchainer_file_condensed =~ s/$DIR/$URL/;
            $html .=
              "<br><span class='small link' onclick=window.open('$final_dagchainer_file_condensed')>Condensed syntelog file with GEvo links</span>";
        }

        if ($svg_file && -r $svg_file) {
            $svg_file =~ s/$DIR/$URL/;
            $html .=
              "<br><span class='small link' onclick=window.open('$svg_file')>SVG Version of Syntenic Dotplot</span>";
        }

        if ($spa_file && -r $spa_file) {
            $spa_file =~ s/$DIR/$URL/;
            $html .=
              "<br><span class='small link' onclick=window.open('$spa_file')>Syntenic Path Assembly mapping</span>";
        }

        $html .= "<tr><td>";
        my $conffile = $ENV{HOME} . 'coge.conf';
        $dagchainer_file =~ s/^$URL/$DIR/;
        $html .= "<br>"
          . qq{<span class="small link" id="" onClick="window.open('bin/SynMap/order_contigs_to_chromosome.pl?f=$dagchainer_file&cf=$conffile;l=$tiny_link');" >Generate Pseudo-Assembled Genomic Sequence</span>}
          if $assemble;
        $html .= qq{</table>};

        ########################################################################
        # Regenerate Analysis Link - HTML
        ########################################################################

        $html .= "<a href='$tiny_link' class='ui-button ui-corner-all'";
        $html .=
          " style='color: #000000' target=_new_synmap>Regenerate this analysis: $tiny_link</a>";

        if ($ks_type) {
            $html .=
              qq{<span  class='ui-button ui-corner-all' onclick="window.open('SynSub.pl?dsgid1=$dsgid1;dsgid2=$dsgid2')">Generate Substitution Matrix of Syntelogs</span>};
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
        }
        $html .= "<br>";
    }
    ##print out all the datafiles created
    $html .= "<br>";

    if ($problem) {
        $html .=
          qq{<span class=alert>There was a problem running your analysis.  Please check the log file for details.</span><br>};
    }

    if ($warn) {
        $html .=
          qq{<span class=alert>Warning: $warn</span><br>};
    }

    ############################################################################
    # Email results, output benchmark and return results
    ############################################################################

    email_results(
        email    => $email,
        html     => $html,
        org1     => $org_name1,
        org2     => $org_name2,
        jobtitle => $job_title,
        link     => $tiny_link) if $email;

    my $benchmarks = qq{
            $org_name1 versus $org_name2
            Benchmarks:
            $algo_name:               $blast_time
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
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("Benchmark", $cogeweb->logfile);
    CoGe::Accessory::Web::write_log($benchmarks, $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
    CoGe::Accessory::Web::write_log("", $cogeweb->logfile);

    # Need to remove this from the output from dotplot
    # -- otherwise it over-loads the stuff in the web-page already.
    # -- This can mess up other loaded js such as tablesoter

    # Update status to completed.
    $job->update({status => 2}) if defined($job);

    $html =~ s/<script src="\/CoGe\/js\/jquery-1.3.2.js"><\/script>//g;
    return $html;
}

################################################################################
# SynMap workflow
################################################################################

sub gen_org_name
{
    my %opts      = @_;
    my $dsgid     = $opts{dsgid};
    my $feat_type = $opts{feat_type} || 1;
    my $write_log = $opts{write_log} || 0;
    my ($dsg) = $coge->resultset('Genome')->search({genome_id => $dsgid},
        {join => 'organism', prefetch => 'organism'});

    my $org_name = $dsg->organism->name;
    my $title =
        $org_name . " (v"
      . $dsg->version
      . ", dsgid"
      . $dsgid . ") "
      . $feat_type;
    $title =~ s/(`|')//g;

    if ($write_log) {
        CoGe::Accessory::Web::write_log("#" x (20), $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("Generating organism name:",
            $cogeweb->logfile);
        CoGe::Accessory::Web::write_log($title, $cogeweb->logfile);
        CoGe::Accessory::Web::write_log("#" x (20) . "\n", $cogeweb->logfile);
    }
    return ($org_name, $title);
}

sub print_debug
{
    my %args = @_;

    if (defined($args{enabled}) && defined($args{msg}) && $args{enabled}) {
        say STDERR "DEBUG: $args{msg}";
    }
}

sub process_local_dups_file
{
    my %opts    = @_;
    my $infile  = $opts{infile};
    my $outfile = $opts{outfile};
    if (   (-r $outfile && -s $outfile)
        || (-r $outfile . ".gz" && -s $outfile . ".gz")) {
        CoGe::Accessory::Web::write_log(
            "Processed tandem duplicate file found: $outfile",
            $cogeweb->logfile);
        return $outfile;
    }
    CoGe::Accessory::Web::gunzip($infile . ".gz") if -r $infile . ".gz";
    return unless -r $infile;
    CoGe::Accessory::Web::write_log(
        "Adding coge links to tandem duplication file.  Infile $infile : Outfile $outfile",
        $cogeweb->logfile);
    $/ = "\n";
    open(IN,  $infile);
    open(OUT, ">$outfile");
    print OUT "#",
      join("\t",
        "FeatList_link", "GEvo_link", "FastaView_link",
        "chr||start||stop||name||strand||type||database_id||gene_order"),
      "\n";

    while (<IN>) {
        chomp;
        next unless $_;
        my @line = split /\t/;
        my %fids;
        foreach (@line) {
            my @item = split /\|\|/;
            next unless $item[6];
            $fids{$item[6]} = 1;
        }
        next unless keys %fids;
        my $featlist = $BASE_URL . "FeatList.pl?";
        map {$featlist .= "fid=$_;"} keys %fids;
        my $fastaview = $BASE_URL . "FastaView.pl?";
        map {$fastaview .= "fid=$_;"} keys %fids;
        my $gevo  = $BASE_URL . "GEvo.pl?";
        my $count = 1;

        foreach my $id (keys %fids) {
            $gevo .= "fid$count=$id;";
            $count++;
        }
        $gevo .= "num_seqs=" . scalar keys %fids;
        print OUT join("\t", $featlist, $gevo, $fastaview, @line), "\n";
    }
    close OUT;
    close IN;

    #`/bin/rm $infile`;
    #$infile =~ s/localdups/nolocaldups.bed/;
    #`/bin/rm $infile`;
    return $outfile;
}

# FIXME: Currently this feature is not called.
# @by Evan Briones
# @on 3/1/2013
#sub run_tandem_finder
#{
#   my %opts    = @_;
#   my $infile  = $opts{infile};    #dag file produced by dat_tools.py
#   my $outfile = $opts{outfile};
#   while ( -e "$outfile.running" )
#   {
#       print STDERR "detecting $outfile.running.  Waiting. . .\n";
#       sleep 60;
#   }
#   unless ( -r $infile && -s $infile )
#   {
#       CoGe::Accessory::Web::write_log( "WARNING:   Cannot run tandem finder! DAGChainer input file ($infile) contains no data!", $cogeweb->logfile );
#       return 0;
#   }
#   if ( -r $outfile )
#   {
#       CoGe::Accessory::Web::write_log( "run_tandem_filter: file $outfile already exists", $cogeweb->logfile );
#       return 1;
#   }
#   my $cmd = "$PYTHON $TANDEM_FINDER -i $infile > $outfile";
#   system "/usr/bin/touch $outfile.running";    #track that a blast anlaysis is running for this
#   CoGe::Accessory::Web::write_log( "run_tandem_filter: running\n\t$cmd", $cogeweb->logfile );
#   `$cmd`;
#   system "/bin/rm $outfile.running" if -r "$outfile.running";    #remove track file
#   return 1                          if -r $outfile;
#}

#FIXME: Currently this feature is disabled
# @by Evan Briones
# @on 3/1/2013
#sub run_adjust_dagchainer_evals
#{
#   my %opts    = @_;
#   my $infile  = $opts{infile};
#   my $outfile = $opts{outfile};
#   my $cvalue  = $opts{cvalue};
#   $cvalue = 4 unless defined $cvalue;
#   while ( -e "$outfile.running" )
#   {
#       print STDERR "detecting $outfile.running.  Waiting. . .\n";
#       sleep 60;
#   }
#   if ( -r $outfile || -r $outfile . ".gz" )
#   {
#       CoGe::Accessory::Web::write_log( "run_adjust_dagchainer_evals: file $outfile already exists", $cogeweb->logfile );
#       return 1;
#   }
#   CoGe::Accessory::Web::gunzip( $infile . ".gz" ) if -r $infile . ".gz";
#   unless ( -r $infile && -s $infile )
#   {
#       CoGe::Accessory::Web::write_log( "WARNING:   Cannot adjust dagchainer evals! DAGChainer input file ($infile) contains no data!", $cogeweb->logfile );
#       return 0;
#   }
#   my $cmd = "$PYTHON $EVAL_ADJUST -c $cvalue $infile > $outfile";
#
#   #There is a parameter that can be passed into this to filter repetitive sequences more or less stringently:
#   # -c   2 gets rid of more stuff; 10 gets rid of less stuff; default is 4
#   #consider making this a parameter than can be adjusted from SynMap -- will need to actually play with this value to see how it works
#   #if implemented, this will require re-naming all the files to account for this parameter
#   #and updating the auto-SynMap link generator for redoing an analysis
#
#   system "/usr/bin/touch $outfile.running";    #track that a blast anlaysis is running for this
#   CoGe::Accessory::Web::write_log( "run_adjust_dagchainer_evals: running\n\t$cmd", $cogeweb->logfile );
#   `$cmd`;
#   system "/bin/rm $outfile.running" if -r "$outfile.running";
#   ;                                            #remove track file
#   return 1 if -r $outfile;
#
#}

#FIXME: Currently this feature is disabled
# @by Evan Briones
# @on 6/18/2013
sub run_find_nearby
{
    my %opts         = @_;
    my $infile       = $opts{infile};
    my $dag_all_file = $opts{dag_all_file};
    my $outfile      = $opts{outfile};
    while (-e "$outfile.running") {
        print STDERR "detecting $outfile.running.  Waiting. . .\n";
        sleep 60;
    }
    if (-r $outfile) {
        CoGe::Accessory::Web::write_log(
            "run find_nearby: file $outfile already exists",
            $cogeweb->logfile);
        return 1;
    }
    my $cmd =
      "$PYTHON $FIND_NEARBY --diags=$infile --all=$dag_all_file > $outfile";
    system "/usr/bin/touch $outfile.running"
      ;    #track that a blast anlaysis is running for this
    CoGe::Accessory::Web::write_log("run find_nearby: running\n\t$cmd",
        $cogeweb->logfile);
    `$cmd`;
    system "/bin/rm $outfile.running" if -r "$outfile.running";
    ;      #remove track file
    return 1 if -r $outfile;
}

sub add_GEvo_links
{
    my %opts   = @_;
    my $infile = $opts{infile};
    my $dsgid1 = $opts{dsgid1};
    my $dsgid2 = $opts{dsgid2};
    $/ = "\n";
    return if (-r "$infile.condensed" || -r "$infile.condensed.gz"); #check this
    CoGe::Accessory::Web::gunzip($infile . ".gz") if -r $infile . ".gz";
    open(IN,  $infile);
    open(OUT, ">$infile.tmp");
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
        my ($fid1, $fid2);

        if ($feat1[6]) {
            $fid1 = $feat1[6];
            $link .= "fid1=" . $fid1;
        } else {
            my ($xmin) = sort ($feat1[1], $feat1[2]);
            my $x = sprintf("%.0f", $xmin + abs($feat1[1] - $feat1[2]) / 2);
            $link .= "chr1=" . $feat1[0] . ";x1=" . $x;
        }
        if ($feat2[6]) {
            $fid2 = $feat2[6];
            $link .= ";fid2=" . $fid2;
        } else {
            my ($xmin) = sort ($feat2[1], $feat2[2]);
            my $x = sprintf("%.0f", $xmin + abs($feat2[1] - $feat2[2]) / 2);
            $link .= ";chr2=" . $feat2[0] . ";x2=" . $x;
        }
        $link .= ";dsgid1=" . $dsgid1;
        $link .= ";dsgid2=" . $dsgid2;

        if ($fid1 && $fid2) {
            $condensed{$fid1 . "_" . $dsgid1}{$fid2 . "_" . $dsgid2} = 1;
            $condensed{$fid2 . "_" . $dsgid2}{$fid1 . "_" . $dsgid1} = 1;
            $names{$fid1} = $feat1[3];
            $names{$fid2} = $feat2[3];
        }

#   accn1=".$feat1[3]."&fid1=".$feat1[6]."&accn2=".$feat2[3]."&fid2=".$feat2[6] if $feat1[3] && $feat1[6] && $feat2[3] && $feat2[6];
        print OUT $_;
        print OUT "\t", $link;
        print OUT "\n";
    }
    close IN;
    close OUT;
    if ($previously_generated) {
        `/bin/rm $infile.tmp` if -r "$infile.tmp";
    } else {
        my $cmd = "/bin/mv $infile.tmp $infile";
        `$cmd`;
    }
    if (keys %condensed
        && !(-r "$infile.condensed" || -r "$infile.condensed.gz")) {

        #take into account transitivity
        foreach my $id1 (keys %condensed) {
            foreach my $id2 (keys %{$condensed{$id1}}) {
                foreach my $id3 (keys %{$condensed{$id2}}) {
                    next if $id1 eq $id2;
                    $condensed{$id1}{$id3} = 1;
                    $condensed{$id3}{$id1} = 1;
                }
            }
        }

        open(OUT, ">$infile.condensed");
        print OUT
          join("\t",
            qw(COUNT GEVO MASKED_GEVO FASTA_LINK GENE_LIST GENE_NAMES)), "\n";
        my %seen;
        foreach my $id1 (
            sort {
                scalar(keys %{$condensed{$b}}) <=>
                  scalar(keys %{$condensed{$a}})
            } keys %condensed
          ) {
            my ($fid1, $dsgid1) = split /_/, $id1;
            next if $seen{$fid1};
            $seen{$fid1} = 1;
            my @names     = $names{$fid1};
            my $gevo_link = $BASE_URL . "GEvo.pl?fid1=$fid1;dsgid1=$dsgid1";
            my $fids      = "fid=$fid1";
            my $count     = 2;

            foreach my $id2 (sort keys %{$condensed{$id1}}) {
                my ($fid2, $dsgid2) = split /_/, $id2, 2;
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
            for my $i (1 .. $count) {
                $gevo_link2 .= ";mask$i=non-cds";
            }
            $gevo_link2 .= ";pad_gs=200000";
            $gevo_link2 .= ";autogo=1";
            my $fasta_link    = $BASE_URL . "FastaView.pl?$fids";
            my $featlist_link = $BASE_URL . "FeatList.pl?$fids";
            print OUT join("\t",
                $count, $gevo_link, $gevo_link2, $fasta_link, $featlist_link,
                @names),
              "\n";
        }
        close OUT;
    }
}

sub add_reverse_match
{
    my %opts   = @_;
    my $infile = $opts{infile};
    $/ = "\n";
    open(IN, $infile);
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
            $stuff .= join(" ", @line) . "\n";
            next;
        }
        my $chr1 = $line[0];
        my $chr2 = $line[4];
        $chr1 =~ s/^a//;
        $chr2 =~ s/^b//;
        next if $chr1 eq $chr2;
        my @tmp1 = @line[1 .. 3];
        my @tmp2 = @line[5 .. 7];
        @line[1 .. 3] = @tmp2;
        @line[5 .. 7] = @tmp1;
        $line[0]      = "a" . $chr2;
        $line[4]      = "b" . $chr1;
        $stuff .= join("\t", @line) . "\n";
    }
    return if $skip;
    close IN;
    open(OUT, ">>$infile");
    print OUT $stuff;
    close OUT;

}

sub check_address_validity
{
    my $address = shift;
    return 'valid' unless $address;
    my $validity =
      $address =~
      /^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/
      ? 'valid'
      : 'invalid';
    return $validity;
}

sub commify
{
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub get_dotplot
{
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

    $url = $P->{SERVER} . "run_dotplot.pl?" . $url;
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
    my ($w, $h) = $img->getBounds();
    $w += 600;
    $h += 250;

    if ($loc) {
        return ($url, $loc, $w, $h);
    }
    my $html =
      qq{<iframe src=$url frameborder=0 width=$w height=$h scrolling=no></iframe>};
    return $html;
}
