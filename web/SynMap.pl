#! /usr/bin/perl -w
use v5.10;
use strict;
no warnings 'redefine';
umask(0); # what is this for? (mdb 10/19/16)

use CoGeX;
use CoGeX::Result::Genome qw(ERROR LOADING);
use CoGe::Accessory::Web qw(url_for api_url_for get_command_path);
use CoGe::Accessory::Utils qw( commify html_escape );
use CoGe::Builder::Tools::SynMap;
use CoGe::Core::Genome qw(genomecmp genomecmp2);
use CoGe::Core::Favorites;
use CoGeDBI qw(get_feature_counts);
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
#use DBIxProfiler;
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
use HTML::Template;
use JSON::XS;
use LWP::UserAgent;
#use Parallel::ForkManager;
use GD;
use File::Path;
use File::Spec::Functions;
use Mail::Mailer;
use Benchmark;
use DBI;
use POSIX;
#use Sort::Versions;

our (
	$config,        $DIR,
	$URL,           $SERVER,        $USER,
	$FORM,          $coge,          $cogeweb,
	$PAGE_NAME,     $PAGE_TITLE,
	$FASTADIR,      $BLASTDBDIR,    $DIAGSDIR,
	$PYTHON,
	$TANDEM_FINDER, $RUN_DAGCHAINER,
	$EVAL_ADJUST,   $FIND_NEARBY,
	$NWALIGN,       $QUOTA_ALIGN,
	$CLUSTER_UTILS, $BASE_URL,
	$BLAST2BED,     $SYNTENY_SCORE, $TEMPDIR,
	$TEMPURL,       $ALGO_LOOKUP,
#	$GZIP,$GUNZIP,
	%FUNCTIONS,     $SCRIPTDIR,    $LINK
);

$| = 1;    # turn off buffering

$FORM       = new CGI;
$PAGE_TITLE = "SynMap";
$PAGE_NAME  = "$PAGE_TITLE.pl";

( $coge, $USER, $config, $LINK ) = CoGe::Accessory::Web->init(
	cgi        => $FORM,
	page_title => $PAGE_TITLE,
);

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
$TEMPDIR  = catdir($config->{TEMPDIR}, 'SynMap');
$TEMPURL  = $config->{TEMPURL} . "SynMap";
#$BLAST    = $config->{BLAST} . " -a " . $MAX_PROC . " -K 80 -m 8 -e 0.0001";

# mdb removed 9/20/13 issue 213
#  . " --dbpath="
#  . $config->{LASTDB};

#$GZIP          = $config->{GZIP};
#$GUNZIP        = $config->{GUNZIP};

$ALGO_LOOKUP = algo_lookup();

$DIAGSDIR = $config->{DIAGSDIR};
$FASTADIR = $config->{FASTADIR};

mkpath( $FASTADIR,         0, 0777 );
mkpath( $DIAGSDIR,         0, 0777 );    # mdb added 7/9/12
mkpath( $config->{LASTDB}, 0, 0777 );    # mdb added 7/9/12
$BLASTDBDIR = $config->{BLASTDB};

$PYTHON        = $config->{PYTHON};                    #this was for python2.5
$TANDEM_FINDER = $config->{TANDEM_FINDER} . " -d 5 -s -r"; #-d option is the distance (in genes) between dups -- not sure if the -s and -r options are needed -- they create dups files based on the input file name

$EVAL_ADJUST    = $config->{EVALUE_ADJUST};

$FIND_NEARBY = $config->{FIND_NEARBY} . " -d 20"; #the parameter here is for nucleotide distances -- will need to make dynamic when gene order is selected -- 5 perhaps?

#programs to run Haibao Tang's quota_align program for merging diagonals and mapping coverage
$QUOTA_ALIGN = 'nice ' . $config->{QUOTA_ALIGN};    #the program
$CLUSTER_UTILS = $config->{CLUSTER_UTILS};    #convert dag output to quota_align input
$SYNTENY_SCORE = $config->{SYNTENY_SCORE};

#$CONVERT_TO_GENE_ORDER = $DIR."/bin/SynMap/convert_to_gene_order.pl";
$NWALIGN = get_command_path('NWALIGN');

my %ajax = CoGe::Accessory::Web::ajax_func();

%FUNCTIONS = (
	get_orgs               => \&get_orgs,
	get_genome_info        => \&get_genome_info,
#	get_previous_analyses  => \&get_previous_analyses,
	get_pair_info          => \&get_pair_info,
	check_address_validity => \&check_address_validity,
#	generate_basefile      => \&generate_basefile,
	get_dotplot            => \&get_dotplot,
	gen_dsg_menu           => \&gen_dsg_menu,
	get_dsg_gc             => \&get_dsg_gc,
	#read_log               => \&CoGe::Accessory::Web::read_log,
	get_results       => \&get_results,
	generate_assembly => \&generate_assembly,
	%ajax,
);

my $pj = new CGI::Ajax(%FUNCTIONS);
if ( $FORM->param('jquery_ajax') ) {
	my %args  = $FORM->Vars;
	my $fname = $args{'fname'};

	if ( $fname and defined $FUNCTIONS{$fname} ) {
	    my $header = $FORM->header;
		if ( $args{args} ) {
			my @args_list = split( /,/, $args{args} );
			print $header, $FUNCTIONS{$fname}->(@args_list);
		}
		else {
			print $header, $FUNCTIONS{$fname}->(%args);
		}
	}
}
else {
	$pj->js_encode_function('escape');
	print $pj->build_html( $FORM, \&gen_html );
}

################################################################################
# Web functions
################################################################################

#sub read_log_test {
#	my %args    = @_;
#	my $logfile = $args{logfile};
#	my $prog    = $args{prog};
#	return unless $logfile;
#	$logfile .= ".log" unless $logfile =~ /log$/;
#	$logfile = $TEMPDIR . "/$logfile" unless $logfile =~ /^$TEMPDIR/;
#	return unless -r $logfile;
#	my $str;
#	open( IN, $logfile );
#
#	while (<IN>) {
#		$str .= $_;
#	}
#	close IN;
#	return $str;
#}

sub gen_html {
	my $template = HTML::Template->new( filename => $config->{TMPLDIR} . 'generic_page.tmpl' );
	$template->param(
		PAGE_TITLE => 'SynMap',
		TITLE      => 'SynMap: Whole Genome Synteny Analysis',
		PAGE_LINK  => $LINK,
		SUPPORT_EMAIL => $config->{SUPPORT_EMAIL},
		HEAD       => qq{},
		USER       => $USER->display_name || '',
		NO_DOCTYPE => 1
	);

	$template->param( LOGON => 1 ) unless $USER->user_name eq "public";

    my ($body) = gen_body();
	$template->param( BODY => $body );
	$template->param(
		HOME       => $config->{SERVER},
		HELP       => 'SynMap',
		WIKI_URL   => $config->{WIKI_URL} || '',
		ADMIN_ONLY => $USER->is_admin,
		CAS_URL    => $config->{CAS_URL} || '',
		COOKIE_NAME => $config->{COOKIE_NAME} || ''
	);
	return $template->output;
}

sub gen_body {
	my $template = HTML::Template->new( filename => $config->{TMPLDIR} . 'SynMap.tmpl' );

	$template->param(
		MAIN => 1,
		API_BASE_URL  => $config->{SERVER} . 'api/v1/', #TODO move into config file or module
        MWIDTH => $FORM->param('w') || 0,
        SUPPORT_EMAIL => $config->{SUPPORT_EMAIL},
        USER_NAME => $USER->user_name
    );

	#set search algorithm on web-page
	my $b = $FORM->param('b');
	if ( defined($b) ) {
		$template->param( $ALGO_LOOKUP->{$b}{opt} => "selected" );
	}
	else {
		$template->param( $ALGO_LOOKUP->{6}{opt} => "selected" );
	}
	my ( $D, $A, $Dm, $gm, $dt, $vis, $dupdist, $cscore ); #AKB Added '$vis' 2016-10-18
	$D  = $FORM->param('D');
	$A  = $FORM->param('A');
	$Dm = $FORM->param('Dm');
	$gm = $FORM->param('gm');
	$gm //= 40;    #/
	$dt     = $FORM->param('dt');
	$vis = $FORM->param('vis');  #AKB Added 2016-10-18
	$cscore = $FORM->param('csco');
	$cscore //= 0;    #/
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
	$template->param( 'DISPLAY_DAGCHAINER_SETTINGS' => $display_dagchainer_settings );
	my $mcs = $FORM->param('mcs');
	$template->param( 'MIN_CHR_SIZE' => $mcs ) if $mcs;

	# Set Visualizer Option # AKB Added 2016-10-18
	if ($vis && $vis =~ /legacy/i ) {
		$template->param( 'LEGACY_SELECT' => 'checked' );
	}
    else {
		$template->param( 'SYNMAP2_SELECT' => 'checked' );
	}

	#will the program automatically run?
	my $autogo = $FORM->param('autogo');
	$autogo = 0 unless defined $autogo;
	$template->param( AUTOGO => $autogo );

    #if the page is loading with genomes, there will be a check for whether the genome is rest
    #populate organism menus
	my $error = 0;

	for ( my $i = 1 ; $i <= 2 ; $i++ ) {
		my $dsgid = 0;
		$dsgid = $FORM->param( 'dsgid' . $i ) if $FORM->param( 'dsgid' . $i );    #old method for specifying genome
		$dsgid = $FORM->param( 'gid' . $i ) if $FORM->param( 'gid' . $i );      #new method for specifying genome
		my $feattype_param = $FORM->param( 'ft' . $i ) if $FORM->param( 'ft' . $i );
		my $name = $FORM->param( 'name' . $i ) if $FORM->param( 'name' . $i );
		my $org_menu = gen_org_menu(
			dsgid          => $dsgid,
			num            => $i,
			feattype_param => $feattype_param,
			name           => $name
		);
		$template->param( "ORG_MENU" . $i => $org_menu );

		my ($dsg) = $coge->resultset('Genome')->find($dsgid);
		if ( $dsgid > 0 and !$USER->has_access_to_genome($dsg) ) {
			$error = 1;
		}
	}

	if ($error) {
		$template->param(
			"error" => 'The genome was not found or is restricted.' );
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
	
	my $fb = $FORM->param('fb');
	$template->param( FRAC_BIAS => "checked" ) if (defined $fb && $fb);
	my $fb_window_size = 100;
	$fb_window_size = $FORM->param('fb_ws') if $FORM->param('fb_ws');
	$template->param( FB_WINDOW_SIZE => $fb_window_size );
	if ($FORM->param('fb_tg')) {
		$template->param( FB_TARGET_GENES => "checked" );
	}
	else {
		$template->param ( FB_ALL_GENES => "checked" );
	}
	my $fb_numquerychr = 25;
	$fb_numquerychr = $FORM->param('fb_nqc') if $FORM->param('fb_nqc');
	$template->param( FB_NUMQUERYCHR => $fb_numquerychr );
	my $fb_numtargetchr = 25;
	$fb_numtargetchr = $FORM->param('fb_ntc') if $FORM->param('fb_ntc');
	$template->param( FB_NUMTARGETCHR => $fb_numtargetchr );
	my $fb_rru = 1;
	$fb_rru = $FORM->param('fb_rru') if defined $FORM->param('fb_rru');
	if ($fb_rru) {
		$template->param( FB_REMOVE_RANDOM_UNKNOWN => "checked" );
	}

	$template->param( 'BOX_DIAGS' => "checked" ) if $FORM->param('bd');
	my $spa = $FORM->param('sp') if $FORM->param('sp');
	$template->param( 'SYNTENIC_PATH' => "checked" ) if $spa;
	$template->param( 'SHOW_NON_SYN' => "checked" ) if $spa && $spa =~ /2/;
	$template->param( 'SPA_FEW_SELECT'  => "selected" ) if $spa && $spa > 0;
	$template->param( 'SPA_MORE_SELECT' => "selected" ) if $spa && $spa < 0;

	my $file = $FORM->param('file');
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
    # $template->param(
	# 	PAGE_TITLE         => $PAGE_TITLE,
	# 	SPLASH_COOKIE_NAME => $PAGE_TITLE . '_splash_disabled',
    #     SPLASH_CONTENTS    => 'This page allows you to compare synteny between two genomes.'
	# 	HELP_URL           => 'https://genomevolution.org/wiki/index.php/SynMap'
	# );
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

	if ( $dsg and $USER->has_access_to_genome($dsg) ) {
		my $org = $dsg->organism;
		$oid = $org->id;

		my ( $dsg_info, $feattype_menu, $message, $total_length, undef, undef, $seq_type ) = get_genome_info(
			dsgid    => $dsgid,
			org_num  => $num,
			feattype => $feattype_param
		);

		$template->param(
			DSG_INFO       => $dsg_info,
			FEATTYPE_MENU  => $feattype_menu,
			GENOME_MESSAGE => $message,
			 # mdb added 11/3/16 COGE-768 -- set values needed for block check because handle_dsg_info() doesn't get called when gid's passed in URL
			ORG_LENGTH     => $total_length,
			SEQ_TYPE       => $seq_type
		);
	}
	else {
		$oid   = 0;
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

    my $favorites = CoGe::Core::Favorites->new(user => $USER);
    
	foreach my $dsg ( sort { genomecmp2($a, $b, $favorites) }
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
		elsif ( $dsg->deleted ) {
			if ( $dsgid && $dsgid == $dsg->id ) {
				$name =
				    "DELETED: "
				  . $dsg->type->name . " (v"
				  . $dsg->version . ",id"
				  . $dsg->id . ")";
			}
			else {
				next;
			}
		}
		else {
		    $name .= "&#11088; " if ($favorites->is_favorite($dsg));
			$name .= "&#x2705; " if $dsg->certified;
		    $name .= "&#x1f512; " if $dsg->restricted;
			$name .= $dsg->name . ": " if $dsg->name;
			$name .= $dsg->type->name . " (v" . $dsg->version . ",id" . $dsg->id . ")";
			$org_name = $dsg->organism->name unless $org_name;
			foreach my $ft (
				$coge->resultset('FeatureType')->search(
					{	genome_id            => $dsg->id,
						'me.feature_type_id' => 3
					},
					{	join =>
						  { features => { dataset => 'dataset_connectors' } },
						rows => 1,
					}
				))
			{
				$has_cds = 1;
			}
		}

		push @dsg_menu, [ $dsg->id, $name, $dsg, $has_cds ];
	}

	return ( qq{<span id="dsgid$num" class="hidden"></span>}, '' )
	  unless (@dsg_menu);

    #my $dsg_menu = qq{<select id="dsgid$num" onChange="\$('#dsg_info$num').html('<div class=dna_small class="loading" class="small">loading. . .</div>'); get_genome_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],[handle_dsg_info])">};
	my $dsg_menu =
	    qq{<div class="coge-padded-top inline bottom">}
	  . qq{<span class="small text">Genomes: </span>}
	  . qq{<select id="dsgid$num" style="max-width:400px;height:2em;" onChange="get_genome_info(['args__dsgid','dsgid$num','args__org_num','args__$num'],[handle_dsg_info])">};

	foreach (@dsg_menu) {
		my ( $numt, $name ) = @$_;
		my $selected = " selected" if $dsgid && $numt == $dsgid;
		$selected = " " unless $selected;
		$numt = 0 if $name eq "Restricted";
		$dsg_menu .= qq{<OPTION VALUE=$numt $selected>$name</option>};
	}
	$dsg_menu .= "</select>";
	$dsg_menu .= "</div>";

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
	my %opts   = @_;
	my $search = $opts{search};
	my $oid    = $opts{oid};
	my $i      = $opts{i};

	#get rid of trailing white-space
	$search =~ s/^\s+//g if $search;
	$search =~ s/\s+$//g if $search;
	$search = ""
	  if $search && $search =~ /Search/;    #need to clear to get full org count

	my @organisms;
	my $org_count;

	# Create terms for search
	my @terms = split /\s+/, $search if defined $search;

	if ( scalar @terms or $oid ) {
		my @constraints = map {
			-or => [
				{ name        => { like => qq{%$_%} } },
				{ description => { like => qq{%$_%} } }
			  ]
		} @terms;

		@organisms = $coge->resultset("Organism")->search(
			{
				-or => [
					-and => \@constraints,
					{ organism_id => $oid },
				]
			}
		);
	}
	else {
		$org_count = $coge->resultset("Organism")->count;
	}

	my @opts;
	foreach my $item ( sort { uc( $a->name ) cmp uc( $b->name ) } @organisms ) {
		my $option = "<OPTION value=\"" . $item->id . "\"";
		$option .= " selected" if $oid && $oid == $item->id;
		$option .= ">" . $item->name . " (id" . $item->id . ")</OPTION>";
		push @opts, $option;
	}

	unless ( @opts && @organisms ) {
		return qq{<span name="org_id$i" id="org_id$i"></span>};
	}

	$org_count = scalar @opts unless $org_count;
	my $html =
	    qq{<span class="small info">Organisms: (}
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
	return ( "<div class='small note indent'>No matching results found</div>",
		" ", " ", '', $org_num, '', '' )
	  unless ( $dsgid && $dsgid =~ /\d+/ );

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
	$html_dsg_info .= qq{<tr><td>Description:</td><td class="link" onclick=window.open('GenomeInfo.pl?gid=$dsgid')>}
	  . $dsg->info
	  . qq{</td></tr>};

    #$html_dsg_info .= qq{<tr><td>Organism:</td><td>$orgname</td></tr>}; # redundant
	$html_dsg_info .= qq{<tr><td>Taxonomy:</td><td>$org_desc</td></tr>};
	$html_dsg_info .= "<tr><td>Name: <td>" . $dsg->name if $dsg->name;
	$html_dsg_info .= "<tr><td>Description: <td>" . $dsg->description if $dsg->description;
	$html_dsg_info .=
	    "<tr><td>Source:  <td><a href=" 
	  . $link
	  . " target=_new>"
	  . $ds->data_source->name . "</a>";
	$html_dsg_info .= 
	    '&nbsp;&nbsp;'
	  . HTML::Template->new( filename => $config->{TMPLDIR} . 'widgets/Certified.tmpl' )->output if $dsg->certified;

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
	$html_dsg_info .= qq{<tr class="alert"><td>Restricted:</td><td>Yes}
	  if $dsg->restricted;
	$html_dsg_info .= "</table>";
	if ( $dsg->restricted && !$USER->has_access_to_genome($dsg) ) {
		$html_dsg_info = "Restricted";
	}
	if ( $dsg->deleted ) {
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

	my $feattype_menu = qq{<select id="feat_type$org_num" name="feat_type$org_num" style="height:2em;">#};
	$feattype_menu .= qq{<OPTION VALUE=1 $cds_selected>CDS</option>} if $has_cds;
	$feattype_menu .= qq{<OPTION VALUE=2 $genomic_selected>genomic</option>};
	$feattype_menu .= "</select>";
	$message = "<span class='small alert'>No Coding Sequence in Genome</span>" unless $has_cds;
	$message = "<span class='small alert'>Genome is still being loaded</span>" if ($dsg->is_loading());
	$message = "<span class='small alert'>Genome is in invalid state</span>"   if ($dsg->is_error());

	return $html_dsg_info, $feattype_menu, $message, $chr_length, $org_num,
	  $dsg->organism->name, $dsg->genomic_sequence_type_id;
}

#sub get_previous_analyses {
#
#	#FIXME:  THis whole sub needs updating or removal!  Lyons 6/12/13
#	my %opts = @_;
#	my $oid1 = $opts{oid1};
#	my $oid2 = $opts{oid2};
#	return unless $oid1 && $oid2;
#	my ($org1) = $coge->resultset('Organism')->find($oid1);
#	my ($org2) = $coge->resultset('Organism')->find($oid2);
#	return
#	  if ( $USER->user_name =~ /public/i
#		&& ( $org1->restricted || $org2->restricted ) );
#	my ($org_name1) = $org1->name;
#	my ($org_name2) = $org2->name;
#	( $oid1, $org_name1, $oid2, $org_name2 ) =
#	  ( $oid2, $org_name2, $oid1, $org_name1 )
#	  if ( $org_name2 lt $org_name1 );
#
#	my $tmp1 = $org_name1;
#	my $tmp2 = $org_name2;
#	foreach my $tmp ( $tmp1, $tmp2 ) {
#		$tmp =~ s/\///g;
#		$tmp =~ s/\s+/_/g;
#		$tmp =~ s/\(//g;
#		$tmp =~ s/\)//g;
#		$tmp =~ s/://g;
#		$tmp =~ s/;//g;
#		$tmp =~ s/#/_/g;
#		$tmp =~ s/'//g;
#		$tmp =~ s/"//g;
#	}
#
#	my $dir = $tmp1 . "/" . $tmp2;
#	$dir = "$DIAGSDIR/" . $dir;
#	my $sqlite = 0;
#	my @items;
#	if ( -d $dir ) {
#		opendir( DIR, $dir );
#		while ( my $file = readdir(DIR) ) {
#			$sqlite = 1 if $file =~ /sqlite$/;
#			next unless $file =~ /\.aligncoords/;    #$/\.merge$/;
#			my ( $D, $g, $A ) = $file =~ /D(\d+)_g(\d+)_A(\d+)/;
#			my ($Dm) = $file =~ /Dm(\d+)/;
#			my ($gm) = $file =~ /gm(\d+)/;
#			my ($ma) = $file =~ /ma(\d+)/;
#			$Dm = " " unless defined $Dm;
#			$gm = " " unless defined $gm;
#			$ma = 0   unless $ma;
#			my $merge_algo;
#			$merge_algo = "DAGChainer" if $ma && $ma == 2;
#
#			if ( $ma && $ma == 1 ) {
#				$merge_algo = "Quota Align";
#				$gm         = " ";
#			}
#			unless ($ma) {
#				$merge_algo = "--none--";
#				$gm         = " ";
#				$Dm         = " ";
#			}
#
#			#       $Dm = 0 unless $Dm;
#			#       $gm = 0 unless $gm;
#			next unless ( $D && $g && $A );
#
#			my ($blast) = $file =~
#			  /^[^\.]+\.[^\.]+\.([^\.]+)/;    #/blastn/ ? "BlastN" : "TBlastX";
#			my $select_val;
#			foreach my $item ( values %$ALGO_LOOKUP ) {
#				if ( $item->{filename} eq $blast ) {
#					$blast      = $item->{displayname};
#					$select_val = $item->{html_select_val};
#				}
#			}
#			my ( $dsgid1, $dsgid2, $type1, $type2 ) =
#			  $file =~ /^(\d+)_(\d+)\.(\w+)-(\w+)/;
#			$type1 = "CDS" if $type1 eq "protein";
#			$type2 = "CDS" if $type2 eq "protein";
#
#			#           print STDERR $file,"\n";
#			#           my ($repeat_filter) = $file =~ /_c(\d+)/;
#			next unless ( $dsgid1 && $dsgid2 && $type1 && $type2 );
#			my ($dupdist) = $file =~ /tdd(\d+)/;
#			my %data = (
#
#		   #                                    repeat_filter => $repeat_filter,
#				tdd        => $dupdist,
#				D          => $D,
#				g          => $g,
#				A          => $A,
#				Dm         => $Dm,
#				gm         => $gm,
#				ma         => $ma,
#				merge_algo => $merge_algo,
#				blast      => $blast,
#				dsgid1     => $dsgid1,
#				dsgid2     => $dsgid2,
#				select_val => $select_val
#			);
#			my $geneorder = $file =~ /\.go/;
#			my $genome1 = $coge->resultset('Genome')->find($dsgid1);
#			next unless $genome1;
#			next
#			  if ( $genome1->restricted
#				&& !$USER->has_access_to_genome($genome1) );
#			my ($ds1) = $genome1->datasets;
#			my $genome2 = $coge->resultset('Genome')->find($dsgid2);
#			next unless $genome2;
#			next
#			  if ( $genome2->restricted
#				&& !$USER->has_access_to_genome($genome2) );
#			my ($ds2) = $genome2->datasets;
#			$data{dsg1} = $genome1;
#			$data{dsg2} = $genome2;
#			$data{ds1}  = $ds1;
#			$data{ds2}  = $ds2;
#
#			$genome1 .= $genome1->name if $genome1->name;
#			$genome1 .= ": "           if $genome1;
#			$genome1 .= $ds1->data_source->name;
#			$genome2 .= $genome2->name if $genome2->name;
#			$genome2 .= ": "           if $genome2;
#			$genome2 .= $ds2->data_source->name;
#			$data{genome1}    = $genome1;
#			$data{genome2}    = $genome2;
#			$data{type_name1} = $type1;
#			$data{type_name2} = $type2;
#			$type1 = $type1 eq "CDS" ? 1 : 2;
#			$type2 = $type2 eq "CDS" ? 1 : 2;
#			$data{type1}   = $type1;
#			$data{type2}   = $type2;
#			$data{dagtype} = $geneorder ? "Ordered genes" : "Distance";
#			push @items, \%data;
#		}
#		closedir(DIR);
#	}
#	return unless @items;
#	my $size = scalar @items;
#	$size = 8 if $size > 8;
#	my $html;
#	my $prev_table = qq{<table id=prev_table class="small resultborder">};
#	$prev_table .= qq{<THEAD><TR><TH>}
#	  . join( "<TH>",
#		qw(Org1 Genome1 Ver1 Genome%20Type1 Sequence%20Type1 Org2 Genome2 Ver2 Genome%20Type2 Sequence%20type2 Algo Dist%20Type Dup%20Dist Ave%20Dist(g) Max%20Dist(D) Min%20Pairs(A))
#	  ) . "</THEAD><TBODY>\n";
#	my %seen;
#
#	foreach my $item (
#		sort { $b->{dsgid1} <=> $a->{dsgid1} || $b->{dsgid2} <=> $a->{dsgid2} }
#		@items )
#	{
#		my $val = join( "_",
#			$item->{g},          $item->{D},       $item->{A},
#			$oid1,               $item->{dsgid1},  $item->{type1},
#			$oid2,               $item->{dsgid2},  $item->{type2},
#			$item->{select_val}, $item->{dagtype}, $item->{tdd} );
#		next if $seen{$val};
#		$seen{$val} = 1;
#		$prev_table .=
#		  qq{<TR class=feat onclick="update_params('$val')" align=center><td>};
#		my $ver1 = $item->{dsg1}->version;
#		$ver1 = "0" . $ver1 if $ver1 =~ /^\./;
#		my $ver2 = $item->{dsg2}->version;
#		$ver2 = "0" . $ver2 if $ver2 =~ /^\./;
#		$prev_table .= join( "<td>",
#			$item->{dsg1}->organism->name, $item->{genome1},
#			$ver1,                         $item->{dsg1}->type->name,
#			$item->{type_name1},           $item->{dsg2}->organism->name,
#			$item->{genome2},              $ver2,
#			$item->{dsg2}->type->name,     $item->{type_name2},
#			$item->{blast},                $item->{dagtype},
#			$item->{tdd},                  $item->{g},
#			$item->{D},                    $item->{A} )
#		  . "\n";
#	}
#	$prev_table .= qq{</TBODY></table>};
#	$html .= $prev_table;
#	$html .=
#"<br><span class=small>Synonymous substitution rates previously calculated</span>"
#	  if $sqlite;
#	return "$html";
#}

sub get_dsg_info {
	my $dsg       = shift;
	my $length    = 0;
	my $chr_count = 0;

	$length    = $dsg->length;          #$rs->first->get_column('total_length');
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
	my $plasmid     = 0;
	my $contig      = 0;
	my $scaffold    = 0;
	my @chromosomes = $dsg->chromosomes;
	if ( @chromosomes > 100 ) {
		return ( 0, 1, 0 );
	}
	foreach my $chr (@chromosomes) {
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

#sub generate_basefile {
#	$cogeweb = CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
#	return $cogeweb->basefilename;
#}

sub get_results {
	my %opts = @_;
	foreach my $k ( keys %opts ) {
		$opts{$k} =~ s/^\s+//;
		$opts{$k} =~ s/\s+$//;
	}
	my $dsgid1 = $opts{dsgid1};
	my $dsgid2 = $opts{dsgid2};

	return encode_json( { error => "You must select two genomes." } ) unless ( $dsgid1 && $dsgid2 );

	my ($genome1) = $coge->resultset('Genome')->find($dsgid1);
	my ($genome2) = $coge->resultset('Genome')->find($dsgid2);

	return encode_json({
		error => "Problem generating dataset group objects for ids:  $dsgid1, $dsgid2."
	}) unless ( $genome1 && $genome2 );

	############################################################################
	# Initialize Job info
	############################################################################
	my $tiny_link = get_query_link($config, $coge, @_);

	say STDERR "tiny_link is required for logging." unless defined($tiny_link);

	#    my ($tiny_id) = $tiny_link =~ /\/(\w+)$/;
	#    my $workflow_name .= "-$tiny_id";

	my $basename = $opts{basename};
	$cogeweb = CoGe::Accessory::Web::initialize_basefile(
		basename => $basename,
		tempdir  => $TEMPDIR
	);
	
	my $result_path = get_result_path($DIAGSDIR, $dsgid1, $dsgid2);
	my $log_path = get_log_file_path($result_path, $tiny_link);
	$cogeweb->logfile($log_path);

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
	my $regen             = $opts{regen_images};
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
		db		  => $coge,
		genome_id     => $dsgid1,
		feat_type => $feat_type1,
	);

	my ( $org_name2, $title2 ) = gen_org_name(
		db		  => $coge,
		genome_id     => $dsgid2,
		feat_type => $feat_type2,
	);

	############################################################################
	# Generate Fasta files
	############################################################################
	my ( $fasta1, $fasta2 );
#	my $workflow = undef;
#	my $status   = undef;

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
		$dsgid1,     $genome1,           $org_name1,  $fasta1,
		$feat_type1, $depth_org_1_ratio, $dsgid2,     $genome2,
		$org_name2,  $fasta2,            $feat_type2, $depth_org_2_ratio
	  )
	  = (
		$dsgid2,     $genome2,           $org_name2,  $fasta2,
		$feat_type2, $depth_org_2_ratio, $dsgid1,     $genome1,
		$org_name1,  $fasta1,            $feat_type1, $depth_org_1_ratio
	  ) if ( $dsgid2 lt $dsgid1 );

	############################################################################
	# Generate blastdb files
	############################################################################
	my ( $blastdb, @blastdb_files );

	if ( $ALGO_LOOKUP->{$blast}{formatdb} ) {
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
			dir => $result_path
		  }
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
	$raw_blastfile .= ".new" if $raw_blastfile =~ /genomic/;

	###########################################################################
	# Converting blast to bed and finding local duplications
	###########################################################################
	#
	# NOTES: The blast2bed program will take the input rawblast file and
	# filter it and creating a new rawblast and moving the old to
	# rawblast.orig
	#
#	my $query_bed   = $raw_blastfile . ".q.bed";
#	my $subject_bed = $raw_blastfile . ".s.bed";

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
#	my $query_dup_file   = $opts{query_dup_files};
#	my $subject_dup_file = $opts{subject_dup_files};
#	my $query            = "a" . $dsgid1;
#	my $subject          = "b" . $dsgid2;

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
		my $max_chr    = $chr1_count;
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
			$svg_file       = $ks_blocks_file . ".svg";
		}
		else {
			$warn = "Unable to calculate Ks or Kn values due to at least"
			  . " one genome lacking CDS features.";
			$ks_type = undef;
		}
	}
	else {
		$svg_file = "$final_dagchainer_file.svg";
	}

	############################################################################
	# Generate dot plot
	############################################################################
	($basename) = $final_dagchainer_file =~ /([^\/]*aligncoords.*)/;    #.all.aligncoords/;
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
	my $spa_file       = $dotfile . ".spa_info.txt";
	my $json_file      = "$json_basename.json";
	my $all_json_file  = "$json_basename.all.json";
	my $hist_json_file = "$json_basename.datasets.json";

	$out = $dotfile;

	############################################################################
	# Generate html
	############################################################################
	my $results = HTML::Template->new( filename => $config->{TMPLDIR} . 'partials/synmap_results.tmpl' );

	my ( $x_label, $y_label );

	if ($clabel) {
		$y_label = "$out.y.png";
		$x_label = "$out.x.png";
		$warn .= qq{Unable to display the y-axis.} unless -r $y_label;
		$warn .= qq{Unable to display the x-axis.} unless -r $x_label;
	}

	my $final_dagchainer_file_gevolinks = $final_dagchainer_file . ".gevolinks"; # mdb added 12/2/16 COGE-794
	my $final_dagchainer_file_condensed = $final_dagchainer_file_gevolinks . ".condensed";

	my $problem;
	if ( -r "$out.html" &&
	     -r $final_dagchainer_file &&           # mdb added 12/2/16 COGE-794
		 -r $final_dagchainer_file_gevolinks && # mdb added 12/2/16 COGE-794
		 -r $final_dagchainer_file_condensed )  # mdb added 12/2/16 COGE-794
	{
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
# not sure why this test is done
#		if ($DIR =~ /$URL/) {
    		$out_url =~ s/$DIR/$URL/;
    		$y_label =~ s/$DIR/$URL/;
    		$x_label =~ s/$DIR/$URL/;
#		}
#		else {
#		    print STDERR "get_results: ERROR !!!!, cannot perform URL substitution: out_url=$out_url DIR=$DIR URL=$URL\n";
#		}

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
			$results->param( yorg_name => html_escape($org_name1) );
			$results->param( xorg_name => html_escape($org_name2) );
		}
		else {
			$results->param( yorg_name => html_escape($org_name2) );
			$results->param( xorg_name => html_escape($org_name1) );
		}

		$results->param( dotplot   => $tmp );
		$results->param( algorithm => $algo_name );

		if ( $hist and $ks_type ) {
			if ( -r $hist and -s $hist ) {
				$results->param( histogram => $out_url . '.hist.png' );
				$results->param( ks_type   => $ks_type );
			}
			else {
				$warn =
				  qq{The histogram was not generated no ks or kn data found.};
				warn "problem reading $hist or is empty";
			}
		}

		# dotplot
		if ($ks_type) {
			my $ks_blocks_file_url = $ks_blocks_file;
			$ks_blocks_file_url =~ s/$DIR/$URL/;
			$results->param( file_url => $ks_blocks_file_url );
		}
		else {
			my $final_dagchainer_url = $final_dagchainer_file;
			$final_dagchainer_url =~ s/$DIR/$URL/;
			$results->param( file_url => $final_dagchainer_url );
		}
		$results->param( dsgid1 => $dsgid1 );
		$results->param( dsgid2 => $dsgid2 );
	    my $chromosomes1 = $genome1->chromosomes_all;
	    my $feature_counts = get_feature_counts($coge->storage->dbh, $genome1->id);
	    foreach (@$chromosomes1) {
	        $_->{gene_count} = $feature_counts->{$_->{name}}{1}{count} ? int($feature_counts->{$_->{name}}{1}{count}) : 0;
	    }
	    $results->param( chromosomes1 => encode_json($chromosomes1) );
	    my $chromosomes2 = $genome2->chromosomes_all;
		$feature_counts = get_feature_counts($coge->storage->dbh, $genome2->id);
	    foreach (@$chromosomes2) {
	        $_->{gene_count} = $feature_counts->{$_->{name}}{1}{count} ? int($feature_counts->{$_->{name}}{1}{count}) : 0;
	    }
	    $results->param( chromosomes2 => encode_json($chromosomes2) );

		# fractionation bias
		# my $gff_sort_output_file;
		my $synmap_dictionary_output_file;
		my $fract_bias_raw_output_file;
		my $fract_bias_results_file;
		if ( $opts{frac_bias} =~ /true/i ) {
			my $output_url = $result_path;
			$output_url =~ s/$DIR/$URL/;
			my $query_id;
			my $target_id;
			if ( $depth_org_1_ratio < $depth_org_2_ratio ) {
				$query_id      = $dsgid2;
				$target_id     = $dsgid1;
				$results -> param( target_genome => html_escape($org_name1))
			}
			else {
				$query_id      = $dsgid1;
				$target_id     = $dsgid2;
				$results -> param( target_genome => html_escape($org_name2))
			}
			my $all_genes = ($opts{'fb_target_genes'} eq 'true') ? 'False' : 'True';
			my $rru = $opts{'fb_remove_random_unknown'} ? 'True' : 'False';
			my $syn_depth = $depth_org_1_ratio . 'to' . $depth_org_2_ratio;
			my $fb_prefix = substr($final_dagchainer_file, length($result_path)) . '_tc' . $opts{fb_numtargetchr} . '_qc' . $opts{fb_numquerychr} . '_sd' . $syn_depth . '_ag' . $all_genes . '_rr' . $rru . '_ws' . $opts{fb_window_size};
			my $fb_json_file = $fb_prefix . '.fractbias-fig.json';
			if (! -r catfile($result_path, $fb_json_file)) {
				return encode_json( { error => "The fractionation bias data could not be found." } );
			}
			$results->param( frac_bias => catfile($output_url, $fb_json_file) );
			# $gff_sort_output_file = _filename_to_link(
			# 	file => catfile($result_path, 'gff_sort.txt'),
			# 	msg  => qq{GFF Sort output file},
			# 	required => 1
			# );
			$synmap_dictionary_output_file = _filename_to_link(
				file => catfile($result_path, $fb_prefix . '.fractbias-synmap-data.json'),
				msg  => qq{SynMap dictionary output file},
				required => 1
			);
			$fract_bias_raw_output_file = _filename_to_link(
				file => catfile($result_path, $fb_prefix . '.fractbias-genes.csv'),
				msg  => qq{Fractionation Bias synteny report},
				required => 1
			);
			$fract_bias_results_file = _filename_to_link(
				file => catfile($result_path, $fb_prefix . '.fractbias-results.csv'),
				msg  => qq{Fractionation Bias sliding window results},
				required => 1
			);
		}

		my $qa_file = $merged_dagchainer_file;
		$qa_file =~ s/\.ma\d$/\.qa/ if $qa_file;
#		my $qa_merged_file = $qa_file . ".merged" if $qa_file;
#		my $qa_coverage_tmp = $quota_align_coverage . ".tmp" if $quota_align_coverage;
#		my $qa_coverage_qa = $quota_align_coverage . ".qa" if $quota_align_coverage;

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

		my $sequence_url1 = api_url_for("genomes/$dsgid1/sequence"); #"api/v1/legacy/sequence; # mdb changed 2/12/16 for hypnotoad
        my $sequence_url2 = api_url_for("genomes/$dsgid2/sequence");

		my $fasta1_url = _filename_to_link(
			url =>
			  ( $feat_type1 eq "genomic" ? $sequence_url1 : undef ),
			file => (
				$feat_type1 eq "genomic"
				? undef
				: $FASTADIR . "/$dsgid1-$feat_type1.fasta"
			),
			msg => qq{Fasta file for $org_name1: $feat_type1}
		);
		my $fasta2_url = _filename_to_link(
			url =>
			  ( $feat_type2 eq "genomic" ? $sequence_url2 : undef ),
			file => (
				$feat_type2 eq "genomic"
				? undef
				: $FASTADIR . "/$dsgid2-$feat_type2.fasta"
			),
			msg => qq{Fasta file for $org_name2: $feat_type2}
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
		if ( $merged_dagchainer_file and -r $merged_dagchainer_file ) {
			$dagchainer_url .= $merged_dag_url;
		}

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
			file => $final_dagchainer_file_gevolinks,
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

		if ( $spa_url and $assemble ) {
			#added by EHL: 6/17/15
			#fixing problem where
			my $genome1_chr_count = $genome1->chromosome_count();
			my $genome2_chr_count = $genome2->chromosome_count();
			my $flip              = 0;
			if ( $genome1_chr_count <= $genome2_chr_count && $assemble < 0 ) {
				$flip = 1;
			}
			$spa_result =
			    $spa_url
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
				result   => $spa_result,
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
		if ( $opts{frac_bias} =~ /true/i ) {
			push @$rows,
			  {
				general  => undef,
				homolog  => undef,
				diagonal => undef,
				result   => $synmap_dictionary_output_file
			  };
			push @$rows,
			  {
				general  => undef,
				homolog  => undef, #$gff_sort_output_file,
				diagonal => undef,
				result   => $fract_bias_raw_output_file
			  };
			push @$rows,
			  {
				general  => undef,
				homolog  => undef,
				diagonal => undef,
				result   => $fract_bias_results_file
			  };
		}
		$results->param( files => $rows );

		########################################################################
		# SynMap3D Link
		########################################################################
		my $syn3d = $BASE_URL . "SynMap3D.pl";
		my $threedlink = $syn3d . "?x_gid=" . $dsgid1 . ";y_gid=" . $dsgid2;
		#print STDERR $threedlink . "\n";
		$results->param( syn3dlink => $threedlink) ;


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
            <span class="coge-button" id = "grimm_link" onclick="post_to_grimm('$seq1','$seq2')" > Rearrangement Analysis</span> <a class="small" href=http://grimm.ucsd.edu/GRIMM/index.html target=_new>(Powered by GRIMM!)</a>
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
		return encode_json(
			{ error => "The output $out.html could not be found." } );
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
	$log            =~ s/$DIR/$URL/;
	$json_file      =~ s/$DIR/$URL/;
	$all_json_file  =~ s/$DIR/$URL/;
	$hist_json_file =~ s/$DIR/$URL/;

	$results->param( error   => $problem ) if $problem;
	$results->param( warning => $warn )    if $warn;
	$results->param( log     => $log );

	#$results->param( json     => $json_file );
	#$results->param( allpairs => $all_json_file );
	#$results->param( hist     => $hist_json_file );

	##print out all the datafiles created
	$html .= "<br>";
	$html .= qq{<span id="clear" style="font-size: 0.8em" class="coge-button" onClick="\$('#results').hide(); \$(this).hide(); \$('#intro').fadeIn();" >Hide Results</span>};
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
	return encode_json( { html => $output } );
}

sub _filename_to_link {
	my %opts = (
		styles   => "link",
		required => 0,
		@_,
	);
	my $file = $opts{file};
	my $url  = $opts{url};
	return unless ( $file or $url );

	my $link;
	if ( -r $file or $url ) {
		if ( !$url ) {
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
		$link = q{<span class="alert">} . $opts{msg} . q{ (missing)} . q{</span};
		warn "missing file: $file";
		my ($package, $filename, $line) = caller;
		warn 'called from line: ' . $line;
	}
	else {

	 #        $link = q{<span style="color:dimgray;">} . $opts{msg} . q{</span};
	}
	return $link;
}

################################################################################
# SynMap workflow
################################################################################

#sub add_reverse_match {
#	my %opts   = @_;
#	my $infile = $opts{infile};
#	$/ = "\n";
#	open( IN, $infile );
#	my $stuff;
#	my $skip = 0;
#	while (<IN>) {
#		chomp;
#		s/^\s+//;
#		$skip = 1
#		  if /GEvo\.pl/
#		; #GEvo links have been added, this file was generated on a previous run.  Skip!
#		last if ($skip);
#		next unless $_;
#		my @line = split /\s+/;
#		if (/^#/) {
#			my $chr1 = $line[2];
#			my $chr2 = $line[4];
#			$chr1 =~ s/^a//;
#			$chr2 =~ s/^b//;
#			next if $chr1 eq $chr2;
#			$line[2] = "b" . $chr2;
#			$line[4] = "a" . $chr1;
#			$stuff .= join( " ", @line ) . "\n";
#			next;
#		}
#		my $chr1 = $line[0];
#		my $chr2 = $line[4];
#		$chr1 =~ s/^a//;
#		$chr2 =~ s/^b//;
#		next if $chr1 eq $chr2;
#		my @tmp1 = @line[ 1 .. 3 ];
#		my @tmp2 = @line[ 5 .. 7 ];
#		@line[ 1 .. 3 ] = @tmp2;
#		@line[ 5 .. 7 ] = @tmp1;
#		$line[0]        = "a" . $chr2;
#		$line[4]        = "b" . $chr1;
#		$stuff .= join( "\t", @line ) . "\n";
#	}
#	return if $skip;
#	close IN;
#	open( OUT, ">>$infile" );
#	print OUT $stuff;
#	close OUT;
#
#}

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

	$params{flip}   = $flip         if $flip;
	$params{regen}  = $regen        if $regen;
	$params{width}  = $width        if $width;
	$params{ksdb}   = $ksdb         if $ksdb;
	$params{kstype} = $kstype       if $kstype;
	$params{log}    = 1             if $kstype;
	$params{min}    = $min          if $min;
	$params{max}    = $max          if $max;
	$params{am}     = $metric       if defined $metric;
	$params{ar}     = $relationship if defined $relationship;
	$params{ct}     = $color_type   if $color_type;
	$params{bd}     = $box_diags    if $box_diags;
	$params{cs}     = $color_scheme if defined $color_scheme;
	$params{fid1}   = $fid1         if defined $fid1 && $fid1 =~ /^\d+$/;
	$params{fid2}   = $fid2         if defined $fid2 && $fid2 =~ /^\d+$/;

	$url = url_for( "run_dotplot.pl", %params ) . "&" . $url;
	my $ua = LWP::UserAgent->new;
	$ua->timeout(10);
	my $response = $ua->get($url);
	unless ( $response->is_success ) {
		return "Unable to get image for dotplot (failed): $url";
	}

	my $content = $response->decoded_content;
	unless ( $content ) {
        return "Unable to get image for dotplot (no content): $url";
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
	
	my $html = qq{<iframe src=$url frameborder=0 width=$w height=$h scrolling=no></iframe>};
	return $html;
}

sub generate_assembly {
	my %opts = @_;
	my $gid1 = $opts{gid1};
	my $gid2 = $opts{gid2};
	my $flip = $opts{flip};    #reverse the order from default

	unless ( $gid1 and $gid2 ) {
		return encode_json(
			{
				success => JSON::false,
				error   => "Wrong number of genomes specified"
			}
		);
	}

	# Swap order if gid2 < gid1
	( $gid1, $gid2 ) = ( $gid2, $gid1 ) if $gid2 lt $gid1;

	my $genome1 = $coge->resultset("Genome")->find($gid1);
	my $genome2 = $coge->resultset("Genome")->find($gid2);

	unless ($USER->has_access_to_genome($genome1)
		and $USER->has_access_to_genome($genome2) )
	{

		return encode_json(
			{
				success => JSON::false,
				error   => "User does not have access to this dataset"
			}
		);
	}

	unless ( $opts{input} and -r $opts{input} ) {
		return encode_json(
			{
				success => JSON::false,
				error   => "The syntenic path assembly could not be found"
			}
		);
	}

	my $filename =
	  qq($gid1-$gid2-) . md5_hex( $opts{input} ) . "." . $flip . ".tar.gz";
	my $output = catfile( $DIAGSDIR, $gid1, $gid2, "assembly", $filename );

	# Submit workflow
	my $submission =
	  generate_pseudo_assembly( $config, $opts{input}, $output, $flip );
	$output =~ s/$DIR/$URL/;

	# Fixup success to return true or false
	return encode_json(
		{
			id      => $submission->{id},
			output  => $output,
			success => $submission->{success} ? JSON::true : JSON::false,
		}
	);
}
