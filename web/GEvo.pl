#!/usr/bin/perl -w

use strict;
use warnings;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Cookie;
use CGI::Ajax;
use URI::Escape;
use HTML::Template;
use Digest::MD5 qw(md5_base64);
use Data::Dumper;
use File::Basename;
use File::Temp;
use File::Path;

use CoGeX;
use CoGeX::Result::Feature;
use CoGe::Accessory::GenBank;
use CoGe::Accessory::Web;
use CoGe::Accessory::Utils qw( commify );
use CoGe::Accessory::bl2seq_report;
use CoGe::Accessory::blastz_report;
use CoGe::Accessory::lagan_report;
use CoGe::Accessory::chaos_report;
use CoGe::Accessory::GenomeThreader_report;
use CoGe::Accessory::dialign_report;
use CoGe::Accessory::dialign_report::anchors;
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Graphics::Feature::GAGA;
use CoGe::Graphics::Feature::Exon_motifs;
use CoGe::Graphics::Feature::AminoAcid;
use CoGe::Graphics::Feature::Domain;
use CoGe::Graphics::Feature::HSP;
use CoGe::Graphics::Feature::Block;
use CoGe::Graphics::Feature::Outline;

#use CoGe::Graphics::Feature::Line;
use Sort::Versions;
use DBIxProfiler;
use DBI;
use Parallel::ForkManager;
use Statistics::Basic;
use Benchmark qw(:all);

#use Mail::Mailer;
use Mail::Mailer;    # ("mail);
use Digest::MD5 qw(md5_hex);
use GD;
use LWP::Simple;     #needed for the merge URL function in GEvo
no warnings 'redefine';

# for security purposes

delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
use vars qw($P $PAGE_TITLE $PAGE_NAME $LINK $DATE
  $DEBUG $BL2SEQ $BLASTZ $LAGAN $CHAOS $DIALIGN $GENOMETHREADER
  $TEMPDIR $TEMPURL $USER $FORM $cogeweb $BENCHMARK $coge
  $NUM_SEQS $MAX_SEQS $MAX_PROC %FUNCTION);

$FORM                 = new CGI;
$CGI::POST_MAX        = 60 * 1024 * 1024;    # 24MB
$CGI::DISABLE_UPLOADS = 0;

$PAGE_TITLE = 'GEvo';
$PAGE_NAME  = "$PAGE_TITLE.pl";

( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

$ENV{PATH} = $P->{COGEDIR};

#print Dumper $P;
$BL2SEQ     = $P->{BL2SEQ};
$BLASTZ     = get_command_path('LASTZ');
$BLASTZ     .= " --ambiguous=iupac";
$LAGAN          = $P->{LAGAN};
$CHAOS          = $P->{CHAOS};
$GENOMETHREADER = $P->{GENOMETHREADER};
$DIALIGN        = $P->{DIALIGN};
$TEMPDIR        = $P->{TEMPDIR} . "GEvo";
$TEMPURL        = $P->{TEMPURL} . "GEvo";
$MAX_PROC       = $P->{MAX_PROC};

$DATE = sprintf(
    "%04d-%02d-%02d %02d:%02d:%02d",
    sub { ( $_[5] + 1900, $_[4] + 1, $_[3] ), $_[2], $_[1], $_[0] }
      ->(localtime)
);

#for chaos
$ENV{'LAGAN_DIR'} = $P->{LAGANDIR};

#for dialign
$ENV{'DIALIGN2_DIR'} = $P->{DIALIGN2_DIR};

# set this to 1 to print verbose messages to logs
$DEBUG     = 0;
$BENCHMARK = 1;
$NUM_SEQS  = 2;         #SHABARI EDIT
$MAX_SEQS  = 30;
$|         = 1;         # turn off buffering

my %ajax = CoGe::Accessory::Web::ajax_func();

#$ajax{dataset_search} = \&dataset_search; #override this method from Accessory::Web for restricted organisms
#$ajax{feat_search} = \&feat_search;
#print STDERR join ("\n", keys %ajax),"\n";

%FUNCTION = (
    %ajax,
    run                    => \&run,
    loading                => \&loading,
    merge_previous         => \&merge_previous,
    add_seq                => \&add_seq,
    get_file               => \&get_file,
    gen_go_run             => \&gen_go_run,
    gen_hsp_colors         => \&gen_hsp_colors,
    save_settings_gevo     => \&save_settings_gevo,
    reset_settings_gevo    => \&reset_settings_gevo,
    check_address_validity => \&check_address_validity,
    genome_search          => \&genome_search,
    get_tiny_url           => \&get_tiny_url,
    add_to_user_history    => \&add_to_user_history,
    get_image_info         => \&get_image_info,
    dataset_search         => \&dataset_search,
    feat_search            => \&feat_search,
);
my $pj = new CGI::Ajax(%FUNCTION);
$pj->JSDEBUG(0);
$pj->DEBUG(0);

#$pj->js_encode_function('escape');
if ( $FORM->param('jquery_ajax') ) {
    dispatch();
}
else {
    print $pj->build_html( $FORM, \&gen_html );
}

#print $pj->build_html($FORM, \&gen_html);

#print $FORM->header;gen_html();

sub dispatch {
    my %args  = $FORM->Vars;
    my $fname = $args{'fname'};
    if ($fname) {

        #my %args = $cgi->Vars;
        #print STDERR Dumper \%args;
        if ( $args{args} ) {
            my @args_list = split( /,/, $args{args} );
            print $FORM->header, $FUNCTION{$fname}->(@args_list);
        }
        else {
            print $FORM->header, $FUNCTION{$fname}->(%args);
        }
    }
}

sub loading {
    my $message = shift || "Generating results. . .";
    return qq{<div class="dna"><div id="loading">$message</div></div>};
}

sub gen_html {
    my $html;    # =  "Content-Type: text/html\n\n";
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( TITLE      => 'GEvo: Genome Evolution Analysis',
    				  PAGE_TITLE => 'GEvo',
    				  PAGE_LINK  => $LINK,
    				  SUPPORT_EMAIL => $P->{SUPPORT_EMAIL},
    				  HOME       => $P->{SERVER},
                      HELP       => 'GEvo',
                      WIKI_URL   => $P->{WIKI_URL} || '' );
    $template->param( USER   => $USER->display_name || '' );
    $template->param( LOGON  => 1 ) unless $USER->user_name eq "public";
    $template->param( NO_BOX => 1 );
    $template->param( BODY   => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin,
                      CAS_URL    => $P->{CAS_URL} || '',
                      COOKIE_NAME => $P->{COOKIE_NAME} || '' );
    my $prebox = HTML::Template->new( filename => $P->{TMPLDIR} . 'GEvo.tmpl' );
    $prebox->param( RESULTS_DIV => 1 );
    $template->param( PREBOX     => $prebox->output );
    $html .= $template->output;
    return $html;
}

sub gen_body {
    my $form  = $FORM;
    my $prefs = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    my $num_seqs = $form->param('num_seqs');#= get_opt(params=>$prefs, form=>$form, param=>'num_seqs');
    $num_seqs = $NUM_SEQS unless defined $num_seqs;

    #$MAX_SEQS = 20 if $form->param('override');
    my $message;
    if ( !( $num_seqs =~ /^\d+$/ ) ) {
        $message .= "Problem with requested number of sequences: '$num_seqs'.  Defaulting to $NUM_SEQS input sequences.";
        $num_seqs = $NUM_SEQS;
    }
    elsif ( $num_seqs < 2 ) {
        $message .= "Minimum number of sequences to compare is two.";
        $num_seqs = 2;
    }
    elsif ( $num_seqs > $MAX_SEQS ) {
        $message .= "Maximum number of sequence is $MAX_SEQS.";
        $num_seqs = $MAX_SEQS;
    }
    my @seq_nums;
    my @seq_sub;

    my $index = 1;
    for ( my $i = 1 ; $i <= $num_seqs ; $i++ ) {
        my (
            $draccn, $pos,   $chr,  $fid, $dsid,
            $dsgid,  $gstid, $mask, $display_order
        );

        #order by which genomi regions are displayed, top to bottom
        $display_order = $index unless $display_order;
        $draccn = $form->param( "accn" . $i ) if $form->param( "accn" . $i );
        $pos    = $form->param( "x" . $i )    if $form->param( "x" . $i );
        $chr = $form->param( "chr" . $i )     if defined $form->param( "chr" . $i );
        $fid = $form->param( "fid" . $i )     if $form->param( "fid" . $i );
        $dsid  = $form->param( 'dsid' . $i )  if $form->param( 'dsid' . $i );
        $dsgid = $form->param( 'dsgid' . $i ) if $form->param( 'dsgid' . $i );
        $dsgid = $form->param( 'gid' . $i )   if $form->param( 'gid' . $i );
        $gstid = $form->param( 'gstid' . $i ) if $form->param( 'gstid' . $i );
        $mask = $form->param( 'mask' . $i )   if ( $form->param( 'mask' . $i ) );
        $mask = undef if $mask && $mask eq "--None--";
        $mask = "non-cds"
          if $form->param( 'maskncs' . $i );    #backwards compatibility

        if ( $fid && $fid =~ /_/ ) {
            ( $fid, $gstid ) = split( /_/, $fid );
        }

        if ( $fid && !$draccn ) {
            my ($feat) = $coge->resultset('Feature')->find($fid);
            ($draccn) = $feat->names if $feat;

            # Check for available genomes
            next unless grep {
                #$_ && !$_->deleted && $USER->has_access_to_genome($_) # mdb removed 7/16/15
                $USER->has_access_to_genome($_) # mdb added 7/16/15
            } $feat->dataset->genomes;

            unless ($draccn)                    #no name!  This is a problem
            {
                $pos  = $feat->start;
                $chr  = $feat->chromosome;
                $dsid = $feat->dataset_id;
            }
        }

        # Check if all genomes have been deleted
        if ($dsid) {
            my $dataset = $coge->resultset("Dataset")->find($dsid);

            if ($dataset) {
                next unless grep {
                    #$_ && !$_->deleted && $USER->has_access_to_genome($_) # mdb removed 7/16/15
                    $USER->has_access_to_genome($_) # mdb added 7/16/15 
                } $dataset->genomes;
            }
        }

        ## Check if the genome was deleted
        if ($dsgid) {
            my $genome = $coge->resultset('Genome')->find($dsgid);
            #next unless $genome && !$genome->deleted && $USER->has_access_to_genome($genome); # mdb removed 7/16/15
            next unless $genome && $USER->has_access_to_genome($genome); # mdb added 7/16/15
        }

        if ($draccn) {
            my $rs = $coge->resultset('Dataset')->search(
                { 'feature_names.name' => uri_unescape($draccn), },
                {
                    'join' => { 'features' => 'feature_names', },
                }
            );

            if ($rs) {
                while (my $ds = $rs->next()) {
                    if ($ds) {
                        next unless grep {
                            #$_ && !$_->deleted && $USER->has_access_to_genome($_) # mdb removed 7/16/15
                            $USER->has_access_to_genome($_) # mdb added 7/16/15
                        } $ds->genomes;
                    }
                }
            }
        }
        
        my $drup = $form->param( 'dr' . $i . 'up' ) if defined $form->param( 'dr' . $i . 'up' );
        my $drdown = $form->param( 'dr' . $i . 'down' ) if defined $form->param( 'dr' . $i . 'down' );
        $drup = $form->param( 'drup' . $i ) if defined $form->param( 'drup' . $i );
        $drdown = $form->param( 'drdown' . $i ) if defined $form->param( 'drdown' . $i );
        $drup   = 10000 unless defined $drup;
        $drdown = 10000 unless defined $drdown;

        #$drup += $pad_gs if $pad_gs;
        #$drdown += $pad_gs if $pad_gs;
        my $gbaccn = $form->param( "gbaccn" . $i )
          if $form->param( "gbaccn" . $i );
        my $gbstart = $form->param( "gbstart" . $i )
          if $form->param( "gbstart" . $i );
        $gbstart = 1 unless defined $gbstart;
        my $gblength = $form->param( "gblength" . $i )
          if defined $form->param( "gblength" . $i );
	    $gblength = 0 unless defined $gblength;
        my $revy = "checked" if $form->param( 'rev' . $i );
        my $revn = "checked" unless $revy;

        $form->param( 'ref' . $i, 1 ) unless defined $form->param( 'ref' . $i ); #when not specified, set to 1
        $form->param( 'ref' . $i, 0 ) if $form->param( 'nref' . $i );    #backward compatiblity

        my $refy = "checked" if $form->param( 'ref' . $i );
        my $refn = "checked" unless $refy;
        my $org_title;
        my $dsg_menu = qq{<input type="hidden" id="dsgid$index"};
        $dsg_menu .= qq{value ="$dsgid"} if $dsgid;
        $dsg_menu .= qq{>};

        ( $org_title, $dsid, $gstid, $dsgid, $dsg_menu ) = get_org_info(
            dsid    => $dsid,
            chr     => $chr,
            gstid   => $gstid,
            dsgid   => $dsgid,
            seq_num => $index,
        ) if $pos;

        # mdb added 11/20/13 issue 254
        my ($start, $stop);
        if ($fid) {
        	my ($feat) = $coge->resultset('Feature')->find($fid);
        	if ($feat) {
        		$start = $feat->start - $drup;
        		$stop = $feat->stop + $drdown;
        	}
        }
        elsif (defined $gbstart) {
        	$start = $gbstart - $drup;
        	$stop = $gbstart + $gblength-1 + $drdown;
        }

        push @seq_nums, { SEQ_NUM => $index, };
        my %opts = (
            SEQ_NUM   => $display_order,
            REV_YES   => $revy,
            REV_NO    => $revn,
            REF_YES   => $refy,
            REF_NO    => $refn,
            DRUP      => $drup,
            DRDOWN    => $drdown,
            DRACCN    => $draccn,
            DSID      => $dsid,
            DSGID     => $dsgid,
            DSG_STUFF => $dsg_menu,
            GSTID     => $gstid,
            GBACCN    => $gbaccn,
            GBSTART   => $gbstart,
            GBLENGTH  => $gblength,
            POS       => $pos,
            ORGINFO   => $org_title,
            CHR       => $chr,
            FEATID    => $fid,
            START     => $start, # mdb added 11/20/13 issue 254
            STOP      => $stop, # mdb added 11/20/13 issue 254
        );

        if ($mask) {
            $opts{MASK_CDS}    = "selected" if $mask eq "cds";
            $opts{MASK_RNA}    = "selected" if $mask eq "rna";
            $opts{MASK_NCDS}   = "selected" if $mask eq "non-cds";
            $opts{MASK_NGENIC} = "selected" if $mask eq "non-genic";
        }
        $opts{COGEPOS} = qq{<option value="cogepos$i" selected="selected">CoGe Database Position</option>}
          if $pos;
        push @seq_sub, {%opts};

        $index++;
    }
    @seq_sub =
      sort { $a->{SEQ_NUM} <=> $b->{SEQ_NUM} } @seq_sub;  #sort based on seq_num
                                                          #page preferences
    my $pad_gs = $form->param("pad_gs") if $form->param("pad_gs");
    $pad_gs = 0 unless $pad_gs;
    my $apply_all =
      get_opt( params => $prefs, form => $form, param => 'apply_all' );
    my $prog = get_opt( params => $prefs, form => $form, param => 'prog' );
    $prog = "blastz" unless $prog;
    my $image_width = get_opt( params => $prefs, form => $form, param => 'iw' );
    $image_width = 1000 unless $image_width;
    my $feature_height =
      get_opt( params => $prefs, form => $form, param => 'fh' );
    $feature_height = 20 unless $feature_height;
    $feature_height = 10 if $num_seqs > 6;
    my $padding =
      get_opt( params => $prefs, form => $form, param => 'padding' );
    $padding = 2 unless defined $padding;
    $padding = 1 if $feature_height <= 10;
    my $gc_color = get_opt( params => $prefs, form => $form, param => 'gc' );
    $gc_color = 0 unless $gc_color;
    my $nt_color = get_opt( params => $prefs, form => $form, param => 'nt' );
    $nt_color = 1 unless defined $nt_color;
    my $cbc_color = get_opt( params => $prefs, form => $form, param => 'cbc' );
    $cbc_color = 0 unless $cbc_color;
    my $skip_feat_overlap_adjust =
      get_opt( params => $prefs, form => $form, param => 'skip_feat_overlap' );
    $skip_feat_overlap_adjust = 1 unless defined $skip_feat_overlap_adjust;
    my $skip_hsp_overlap_adjust =
      get_opt( params => $prefs, form => $form, param => 'skip_hsp_overlap' );
    $skip_hsp_overlap_adjust = 1 unless defined $skip_hsp_overlap_adjust;
    my $hiqual = get_opt( params => $prefs, form => $form, param => 'hiqual' );
    $hiqual = 0 unless $hiqual;
    my $color_anchor = get_opt( params => $prefs, form => $form, param => 'ca' );
    $color_anchor = 1 unless defined $color_anchor;
    my $comp_adj =
      get_opt( params => $prefs, form => $form, param => 'comp_adj' );
    $comp_adj = 0 unless $comp_adj;
    my $hsp_track =
      get_opt( params => $prefs, form => $form, param => 'hsp_track' );
    $hsp_track = 0 unless $hsp_track;
    my $hsp_top =
      get_opt( params => $prefs, form => $form, param => 'hsp_top' );
    $hsp_top = 0 unless $hsp_top;
    my $hsp_single_color =
      get_opt( params => $prefs, form => $form, param => 'hspsc' );
    $hsp_single_color = 0 unless $hsp_single_color;
    my $color_hsp =
      get_opt( params => $prefs, form => $form, param => 'color_hsp' );
    $color_hsp = 0 unless $color_hsp;
    my $color_feat =
      get_opt( params => $prefs, form => $form, param => 'colorfeat' );
    $color_feat = 0 unless $color_feat;
    my $hsp_label =
      get_opt( params => $prefs, form => $form, param => 'hsp_labels' );
    $hsp_label = undef unless defined $hsp_label;

#    my $hsp_limit = get_opt(params=>$prefs, form=>$form, param=>'hsplim');
#    $hsp_limit = 0 unless $hsp_limit;
#    my $hsp_limit_num = get_opt(params=>$prefs, form=>$form, param=>'hsplimnum');
#    $hsp_limit_num = 20 unless defined $hsp_limit_num;
    my $draw_model =
      get_opt( params => $prefs, form => $form, param => 'draw_model' );
    $draw_model = "full" unless $draw_model;
    my $hsp_overlap_limit =
      get_opt( params => $prefs, form => $form, param => 'hsp_overlap_limit' );
    $hsp_overlap_limit = 0 unless $hsp_overlap_limit;
    my $hsp_size_limit =
      get_opt( params => $prefs, form => $form, param => 'hsp_size_limit' );
    $hsp_size_limit = 0 unless $hsp_size_limit;
    my $show_cns =
      get_opt( params => $prefs, form => $form, param => 'show_cns' );
    $show_cns = 0 unless $show_cns;
    my $show_ofeat =
      get_opt( params => $prefs, form => $form, param => 'show_ofeat' );
    $show_ofeat = 0 unless $show_ofeat;
    my $feat_labels =
      get_opt( params => $prefs, form => $form, param => 'feat_labels' );
    $feat_labels = 0 unless $feat_labels;
    my $show_gene_space =
      get_opt( params => $prefs, form => $form, param => 'show_gene_space' );
    my $show_contigs =
      get_opt( params => $prefs, form => $form, param => 'show_contigs' );
    $show_gene_space = 0 unless $show_gene_space;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GEvo.tmpl' );

    # Check if a genome was specified
    my $error = (scalar @seq_sub < $num_seqs) ? 1 : 0;

    $template->param( ERROR             => $error );
    $template->param( PAD_GS            => $pad_gs );
    $template->param( APPLY_ALL         => $apply_all );
    $template->param( IMAGE_WIDTH       => $image_width );
    $template->param( FEAT_HEIGHT       => $feature_height );
    $template->param( PADDING           => $padding );
    $template->param( HSP_OVERLAP_LIMIT => $hsp_overlap_limit );
    $template->param( HSP_SIZE_LIMIT    => $hsp_size_limit );
    if   ($gc_color) { $template->param( GC_COLOR_YES => "checked" ); }
    else             { $template->param( GC_COLOR_NO  => "checked" ); }
    if   ($nt_color) { $template->param( NT_COLOR_YES => "checked" ); }
    else             { $template->param( NT_COLOR_NO  => "checked" ); }
    if   ($cbc_color) { $template->param( CBC_YES => "checked" ); }
    else              { $template->param( CBC_NO  => "checked" ); }
    if   ($show_contigs) { $template->param( SHOW_CONTIGS_YES => "checked" ); }
    else                 { $template->param( SHOW_CONTIGS_NO  => "checked" ); }

    if ($skip_feat_overlap_adjust) {
        $template->param( FEAT_OVERLAP_NO => "checked" );
    }
    else { $template->param( FEAT_OVERLAP_YES => "checked" ); }
    if ($skip_hsp_overlap_adjust) {
        $template->param( HSP_OVERLAP_NO => "checked" );
    }
    else { $template->param( HSP_OVERLAP_YES => "checked" ); }
    if   ($hiqual) { $template->param( HIQUAL_YES => "checked" ); }
    else           { $template->param( HIQUAL_NO  => "checked" ); }
    if   ($color_anchor) { $template->param( CA_YES => "checked" ); }
    else           { $template->param( CA_NO  => "checked" ); }
    if   ($comp_adj) { $template->param( COMP_ADJ_YES => "checked" ); }
    else             { $template->param( COMP_ADJ_NO  => "checked" ); }
    if   ($hsp_top) { $template->param( HSP_TOP_YES => "checked" ); }
    else            { $template->param( HSP_TOP_NO  => "checked" ); }

    if ($hsp_single_color) {
        $template->param( HSP_SINGLE_COLOR_YES => "checked" );
    }
    else { $template->param( HSP_SINGLE_COLOR_NO => "checked" ); }
    if   ($hsp_track) { $template->param( HSP_TRACK_YES => "checked" ); }
    else              { $template->param( HSP_TRACK_NO  => "checked" ); }
    if   ($color_hsp) { $template->param( COLOR_HSP_YES => "checked" ); }
    else              { $template->param( COLOR_HSP_NO  => "checked" ); }
    if   ($color_feat) { $template->param( COLOR_FEAT_YES => "checked" ); }
    else               { $template->param( COLOR_FEAT_NO  => "checked" ); }

    if   ($show_cns) { $template->param( SHOW_CNS_YES => "checked" ); }
    else             { $template->param( SHOW_CNS_NO  => "checked" ); }

    if   ($show_ofeat) { $template->param( SHOW_OFEAT_YES => "checked" ); }
    else               { $template->param( SHOW_OFEAT_NO  => "checked" ); }

    if ($show_gene_space) {
        $template->param( SHOW_GENESPACE_YES => "checked" );
    }
    else { $template->param( SHOW_GENESPACE_NO => "checked" ); }
    if ( $hsp_label && $hsp_label eq "staggered" ) {
        $template->param( HSP_LABELS_STAG => "selected" );
    }
    elsif ( $hsp_label && $hsp_label eq "linear" ) {
        $template->param( HSP_LABELS_LIN => "selected" );
    }
    else { $template->param( HSP_LABELS_NO => "selected" ); }
    if ( $feat_labels && $feat_labels eq "staggered" ) {
        $template->param( FEAT_LABELS_STAG => "selected" );
    }
    elsif ( $feat_labels && $feat_labels eq "linear" ) {
        $template->param( FEAT_LABELS_LIN => "selected" );
    }
    else { $template->param( FEAT_LABELS_NO => "selected" ); }

    #    if ($feat_labels) {$template->param(FEAT_LABELS_YES=>"checked");}
    #    else {$template->param(FEAT_LABELS_NO=>"checked");}
    #    if ($hsp_limit) {$template->param(HSP_LIMIT_YES=>"checked");}
    #    else {$template->param(HSP_LIMIT_NO=>"checked");}
    #    $template->param(HSP_LIMIT_NUM=>$hsp_limit_num);
    if ( $draw_model eq "full" ) {
        $template->param( DRAW_MODEL_FULL => "selected" );
    }
    elsif ( $draw_model eq "gene" ) {
        $template->param( DRAW_MODEL_GENE => "selected" );
    }
    elsif ( $draw_model eq "mRNA" ) {
        $template->param( DRAW_MODEL_mRNA => "selected" );
    }
    elsif ( $draw_model eq "CDS" ) {
        $template->param( DRAW_MODEL_CDS => "selected" );
    }
    elsif ( $draw_model eq "RNA" ) {
        $template->param( DRAW_MODEL_RNA => "selected" );
    }

#elsif ($draw_model eq "Gene_space") {$template->param(DRAW_MODEL_GENE_SPACE=>"selected");}
    else { $template->param( DRAW_MODEL_NO => "selected" ); }

    #set algorithm parameters
    my $bnW = 7;
    $bnW = $form->param('bnW') if defined $form->param('bnW');
    $template->param( BLAST_WORDSIZE => $bnW );
    my $bnG = 5;
    $bnG = $form->param('bnG') if defined $form->param('bnG');
    $template->param( BLAST_GAPOPEN => $bnG );
    my $bnE = 2;
    $bnE = $form->param('bnE') if defined $form->param('bnE');
    $template->param( BLAST_GAPEXT => $bnE );
    my $bnq = -2;
    $bnq = $form->param('bnq') if defined $form->param('bnq');
    $template->param( BLAST_MISMATCH => $bnq );
    my $bnr = 1;
    $bnr = $form->param('bnr') if defined $form->param('bnr');
    $template->param( BLAST_MATCH => $bnr );
    my $bne = 30;
    $bne = $form->param('bne') if defined $form->param('bne');
    $template->param( BLAST_EVAL => $bne );

    #param for whether the dust/seg filter is on or off
    my $bnF = "F";
    $bnF = $form->param('bnF') if defined $form->param('bnF');
    if ( $bnF eq "F" ) {
        $template->param( BLAST_FILTER_NO => "selected" );
    }
    else {
        $template->param( BLAST_FILTER_YES => "selected" );
    }

   #param for whether or not hsps with stop codons are shown during tblastx runs
    my $hide_stop;
    $hide_stop = $form->param('hs') if defined $form->param('hs');

    if ($hide_stop) {
        $template->param( HIDE_STOP_YES => "checked" );
    }
    else {
        $template->param( HIDE_STOP_NO => "checked" );
    }
    my $bzW = 8;
    $bzW = $form->param('bzW') if defined $form->param('bzW');
    $template->param( BLASTZ_WORDSIZE => $bzW );
    my $bzO = 400;
    $bzO = $form->param('bzO') if defined $form->param('bzO');
    $template->param( BLASTZ_GAPSTART => $bzO );
    my $bzE = 30;
    $bzE = $form->param('bzE') if defined $form->param('bzE');
    $template->param( BLASTZ_GAPEXT => $bzE );
    my $bzK = 3000;
    $bzK = $form->param('bzK') if defined $form->param('bzK');
    $template->param( BLASTZ_SCORETHRESH => $bzK );
    my $bzM = 0;
    $bzM = $form->param('bzM') if defined $form->param('bzM');
    $template->param( BLASTZ_MASKTHRESH => $bzM );
    my $bzC = "0";
    $bzC = $form->param('bzC') if defined $form->param('bzC');

    if ( $bzC == 0 ) {
        $template->param( BLASTZ_CHAIN_NO => "selected" );
    }
    elsif ( $bzC == 1 ) {
        $template->param( BLASTZ_CHAIN_OUT => "selected" );
    }
    elsif ( $bzC == 2 ) {
        $template->param( BLASTZ_CHAIN_EXT => "selected" );
    }
    elsif ( $bzC == 3 ) {
        $template->param( BLASTZ_CHAIN_HSP => "selected" );
    }

    my $box = HTML::Template->new( filename => $P->{TMPLDIR} . 'box.tmpl' );

    #generate sequence submission selector
    $template->param( SEQ_SELECT      => 1 );
    $template->param( SEQ_SELECT_LOOP => \@seq_sub );
    my $seq_submission = $template->output;
    $template->param( SEQ_SELECT => 0 );

    #generate the hsp color option
    my ( $hsp_colors, $num_colors ) =
      gen_hsp_colors( num_seqs => $num_seqs, prefs => $prefs );
    my $spike_len = 15;
    $spike_len = $form->param('spike_len') if defined $form->param('spike_len');
    $template->param( SPIKE_LEN     => $spike_len );
    $template->param( SEQ_RETRIEVAL => 1 );
    $template->param( NUM_SEQS      => scalar @seq_sub);

    # $template->param(COLOR_NUM=>$num_colors);
    $message .= "<BR/>" if $message;
    $template->param( MESSAGE   => $message );
    $template->param( SEQ_SUB   => $seq_submission );
    $template->param( HSP_COLOR => $hsp_colors );
    $template->param( GO_RUN    => gen_go_run($num_seqs) );
    #my $cmd          = "/usr/bin/svnversion " . $P->{COGEDIR} . "gobe/flash";
    #my $gobe_version = `$cmd`;
    #$gobe_version =~ s/\n//g;
    #$template->param( GOBE_VERSION => $gobe_version );

    $template->param( OPTIONS            => 1 );
    $template->param( ALIGNMENT_PROGRAMS => algorithm_list($prog) );
    $template->param( SAVE_SETTINGS      => gen_save_settings($num_seqs) )
      unless !$USER || $USER->user_name =~ /public/i;
    $template->param( 'TEMPDIR' => $TEMPDIR );

    $box->param( BODY => $template->output );

    my $html;
    $html .= $box->output;
    return $html;
}

sub run {
    my %opts              = @_;
    my $num_seqs          = $opts{num_seqs} || $NUM_SEQS;
    my $spike_len         = $opts{spike};
    my $iw                = $opts{iw};
    my $feat_h            = $opts{fh};
    my $show_gc           = $opts{gc};
    my $show_nt           = $opts{nt};
    my $show_cbc          = $opts{cbc};
    my $color_hsp         = $opts{color_hsp};
    my $comp_adj          = $opts{comp_adj};
    my $hsp_top           = $opts{hsp_top};
    my $hsp_track         = $opts{hsp_track};
    my $hsp_labels        = $opts{hsp_labels};
    my $hsp_single_color  = $opts{hsp_single_color};
    my $feat_labels       = $opts{feat_labels};
    my $draw_model        = $opts{draw_model};
    my $hsp_overlap_limit = $opts{hsp_overlap_limit};
    my $hsp_size_limit    = $opts{hsp_size_limit};
    my $hiqual            = $opts{hiqual};

    #my $hsp_limit = $opts{hsplim};
    #my $hsp_limit_num = $opts{hsplimnum};
    my $show_hsps_with_stop_codon = $opts{showallhsps};
    my $padding                   = $opts{padding};
    my ( $analysis_program, $param_string, $parser_opts, $add_gevo_link ) =
      get_algorithm_options(%opts);
    my $pad_gs                    = $opts{pad_gs} || 0;
    my $color_overlapped_features = $opts{color_overlapped_features};
    my $hsp_overlap_length        = $opts{hsp_overlap_length};
    my $email_address             = $opts{email};
    my $show_cns                  = $opts{show_cns};
    my $show_ofeat                = $opts{show_ofeat};
    my $show_gene_space           = $opts{show_gene_space};
    my $show_contigs              = $opts{show_contigs};
    my $skip_feat_overlap_search  = $opts{skip_feat_overlap};
    my $skip_hsp_overlap_search   = $opts{skip_hsp_overlap};
    my $font_size                 = $opts{font_size};
    my $color_anchor              = $opts{ca};
    $color_anchor = 1 unless defined $color_anchor; #default to 1 which is the original behavior of GEvo, EHL 5/5/15
    my $message;
    my $gen_prot_sequence =
      0;    #flag for generating fasta file of protein_sequence;
    $gen_prot_sequence = 1 if $analysis_program eq "GenomeThreader";
    my $basefilename = $opts{basefile};
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basefilename,
        tempdir  => $TEMPDIR
    );
###### working on basefile name initialization problems
    if ( !$basefilename || $basefilename eq "undefined" ) {
        $cogeweb =
          CoGe::Accessory::Web::initialize_basefile( tempdir => $TEMPDIR );
        $basefilename = $cogeweb->basefilename;
    }
######
    #print STDERR "Running GEvo:  basefile:  $basefilename\n";
    CoGe::Accessory::Web::write_log(
        "Beginning GEvo analysis.  Basefilename: $basefilename",
        $cogeweb->logfile );
    my @hsp_colors;
    for ( my $i = 1 ; $i <= num_colors($num_seqs) ; $i++ ) {
        my $rgb = $opts{"rgb$i"};
        my @tmp;
        my ( $r, $g, $b ) =
          $rgb =~ /^rgb\(\s*(\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})\)$/;
        foreach my $color ( $r, $g, $b ) {
            $color = 0 unless $color && $color =~ /^\d+$/;
            $color = 0   if $color < 0;
            $color = 255 if $color > 255;
            push @tmp, $color;
        }
        push @hsp_colors, \@tmp;
    }

    #my $stagger_label = $hsp_label && $hsp_label =~ /staggered/i ? 1 : 0;
    #my $feature_labels = !$hsp_label ? 0 : 1;
    my $form      = $FORM;
    my $gevo_link = $form->url . "?prog=$analysis_program";
    $gevo_link .= ";show_cns=1"        if $show_cns;
    $gevo_link .= ";show_ofeat=1"      if $show_ofeat;
    $gevo_link .= ";show_gene_space=1" if $show_gene_space;
    $gevo_link .= ";show_contigs=1"    if $show_contigs;
    $gevo_link .= ";iw=$iw";
    $gevo_link .= ";fh=$feat_h";
    $gevo_link .= ";padding=$padding";
    $gevo_link .= ";gc=$show_gc"       if $show_gc;
    $gevo_link .= ";color_hsp=1"       if $color_hsp;
    $gevo_link .= ";hsp_top=1"         if $hsp_top;
    $gevo_link .= ";hspsc=1"           if $hsp_single_color;
    $gevo_link .= ";hsp_track=1"       if $hsp_track;
    $gevo_link .= ";comp_adj=1"        if $comp_adj;
    $gevo_link .= ";colorfeat=1"       if $color_overlapped_features;
    $gevo_link .= ";nt=$show_nt";
    $gevo_link .= ";cbc=$show_cbc";
    $gevo_link .= ";spike_len=$spike_len";
    $gevo_link .= ";ca=$color_anchor";
    $gevo_link .= ";skip_feat_overlap=$skip_feat_overlap_search"
      if defined $skip_feat_overlap_search;
    $gevo_link .= ";skip_hsp_overlap=$skip_hsp_overlap_search"
      if defined $skip_hsp_overlap_search;
    $gevo_link .= ";hs=$show_hsps_with_stop_codon"
      if defined $show_hsps_with_stop_codon;
    $gevo_link .= $add_gevo_link;

    my @gevo_link_seqs;
    my @coge_seqs
      ; #place to store stuff for parallel creation of sequence file from genome database
    my @sets;
    my $html;
    my $t1 = new Benchmark;

    for ( my $i = 1 ; $i <= $num_seqs ; $i++ ) {
        my $display_order = $opts{"display_order$i"};
        ($display_order) = $display_order =~ /(\d+)/ if $display_order;
        $display_order = $i unless $display_order;

        my $skip_seq = $opts{"skip_seq$i"};
        next if $skip_seq;
        my $accn = $opts{"draccn$i"};
        $accn = uri_unescape($accn);    # mdb added 10/9/12 issue #11

        #print STDERR "run: accn=$accn\n";
        #trim whitespace at ends of accn
        $accn =~ s/^\s+// if $accn;
        $accn =~ s/\s+$// if $accn;
        my $featid = $opts{"featid$i"};
        my $gstid =
          $opts{"gstid$i"}; #currently not used, value passed attached to featid
        ( $featid, $gstid ) = split( /_/, $featid ) if ( $featid =~ /_/ );
        my $feat = $coge->resultset('Feature')->find($featid) if $featid;
        my $dsgid = $opts{"dsgid$i"};
        my $dsid;
        $dsid = $feat->dataset_id if $feat;
        $dsid = $opts{"dsid$i"} unless $dsid;
        my $chr;
        $chr = $feat->chromosome if $feat;
        $chr = $opts{"chr$i"} unless defined $chr;

        my $drup   = $opts{"drup$i"};
        my $drdown = $opts{"drdown$i"};
        $drup   += $pad_gs;
        $drdown += $pad_gs;
        my $pos    = $opts{"pos$i"};
        my $gbaccn = $opts{"gbaccn$i"};
        $gbaccn =~ s/^\s+//g;
        $gbaccn =~ s/\s+$//g;
        my $gbstart = $opts{"gbstart$i"};
        $gbstart = 1 unless defined $gbstart;
        my $gblength = $opts{"gblength$i"};
        $gblength =~ s/\s+//g;
        my $dirseq    = $opts{"dirseq$i"};
        my $dirstart  = $opts{"dirstart$i"};
        my $dirlength = $opts{"dirlength$i"};

        my $rev  = $opts{"rev$i"};
        my $mask = $opts{"mask$i"};
        $mask = undef if $mask eq "--None--";
        my ( $up,   $down );
        my ( $file, $obj );
        my $reference_seq = $opts{"ref_seq$i"};
        next unless $accn || $featid || $gbaccn || $dirseq || $pos;
        my %gevo_link_info;
        $gevo_link_info{accn}     = CGI::escape($accn)   if $accn;
        $gevo_link_info{x}        = CGI::escape($pos)    if $pos;
        $gevo_link_info{fid}      = CGI::escape($featid) if $featid;
        $gevo_link_info{dsid}     = CGI::escape($dsid)   if $dsid;
        $gevo_link_info{dsgid}    = CGI::escape($dsgid)  if $dsgid;
        $gevo_link_info{gstid}    = CGI::escape($gstid)  if $gstid;
        $gevo_link_info{chr}      = CGI::escape($chr)    if defined $chr;
        $gevo_link_info{drup}     = $drup                if defined $drup;
        $gevo_link_info{drdown}   = $drdown              if defined $drdown;
        $gevo_link_info{gbaccn}   = CGI::escape($gbaccn) if $gbaccn;
        $gevo_link_info{gbstart}  = $gbstart             if $gbstart;
        $gevo_link_info{gblength} = $gblength            if $gblength;
        $gevo_link_info{rev}      = 1                    if $rev;
        $gevo_link_info{ref} = $reference_seq ? 1 : 0;
        $gevo_link_info{mask} = $mask          if $mask;
        $gevo_link_info{do}   = $display_order if $display_order;
        push @gevo_link_seqs, \%gevo_link_info;

        if ( $featid || $pos ) {
            my $error;
            ( $obj, $error ) = get_obj_from_genome_db(
                accn              => $accn,
                featid            => $featid,
                pos               => $pos,
                dsid              => $dsid,
                rev               => $rev,
                up                => $drup,
                down              => $drdown,
                chr               => $chr,
                gstid             => $gstid,
                mask              => $mask,
                dsgid             => $dsgid,
                gen_prot_sequence => $gen_prot_sequence,
                color_anchor		  => $color_anchor #color anchor gene yellow flag.  Defaults to 1 in the sub.
            );
            if ($obj)
            { #going to generalize this in parallel after all sequences to be retrieved from coge's db are specified.
                push @coge_seqs, { obj => $obj, mask => $mask };
                $up   = $drup;
                $down = $drdown;
            }
            else {
                $message .=
                  "Unable to generate sequence for sequence $i.  (Skipped.)\n";
                $message .= " Error message: $error\n" if $error;
            }
        }
        elsif ($dirseq) {
            ($obj) = generate_obj_from_seq(
                seq    => $dirseq,
                num    => $i,
                rc     => $rev,
                start  => $dirstart,
                length => $dirlength
            );
            $dirlength = length( $obj->sequence ) - $dirstart + 1
              unless $dirlength;
            if ($obj) {

                #add an anchor
                my $anchor_stop =
                    $dirlength
                  ? $dirlength
                  : length( $obj->sequence ) - $dirstart + 1;
                $obj->add_feature(
                    type => "direct sequence submission",

                    #location=>(1-$dirstart*2+1)."..".(1-$dirstart*2+1),
                    location   => ( 1 - $dirstart + 1 ) . ".." . $anchor_stop,
                    strand     => 1,
                    qualifiers => {
                        type  => "anchor",
                        names => [ $obj->accn ],
                    },
                    force => 1,
                );
                $file = write_fasta(
                    obj      => $obj,
                    mask     => $mask,
                    startpos => $dirstart,
                    length   => $dirlength,
                    force    => 1,
                );
                $obj->start($dirstart);
                $obj->stop( $dirstart + length( $obj->sequence ) - 1 );
                $obj->chromosome(1);
                $up   = $dirstart;
                $down = $dirlength;
            }
            else {
                $message .= "Problem with direct sequence submission\n";
            }
        }
        elsif ($gbaccn) {
            my $got = 0;
            my $try = 0;
            my ($tmp) = $gbaccn =~ /(.*)/;
            print STDERR "GBACCN: $gbaccn\n";
            $tmp =~ s/^\s+//;
            $tmp =~ s/\s+$//;
            my $gbfile = $TEMPDIR . "/" . uc($tmp) . ".gbk";
            $obj = new CoGe::Accessory::GenBank;

            while ( !$got && $try < 5 ) {
                if ( $try > 1 && -r $gbfile ) {
                    unlink $gbfile;
                }
                my ( $res, $error ) = $obj->get_genbank_from_ncbi(
                    file       => $gbfile,
                    accn       => $tmp,
                    rev        => $rev,
                    start      => $gbstart,
                    length     => $gblength,
                    retmax     => 1,
                    complexity => 1
                );
                $message .= $error . "\n" unless $res;
                $got = 1 if $obj->accn;
                $try++;
            }
            $obj->srcfile( $TEMPDIR . "/" . uc($tmp) . ".faa" );

      #	    $obj->add_gene_models(1); #may want to make a user selectable option
            if ( $obj->accn ) {

                #add an anchor
                my $anchor_stop =
                    $gblength
                  ? $gblength
                  : length( $obj->sequence ) - $gbstart + 1;
                $obj->add_feature(
                    type       => "genbank entry",
                    location   => ( 1 - $gbstart + 1 ) . ".." . $anchor_stop,
                    strand     => 1,
                    qualifiers => {
                        type  => "anchor",
                        names => [$gbaccn],
                    },
                    force => 1
                );
                $gbstart = $obj->seq_length - $gbstart + 1 if $rev;
                $gbstart -= $gblength if $rev && $gblength;
                $file = write_fasta(
                    obj      => $obj,
                    mask     => $mask,
                    startpos => $gbstart,
                    length   => $gblength,
                    force    => 1,
                );
                $obj->start($gbstart);
                $obj->stop( $gbstart + length( $obj->sequence ) - 1 );
                $obj->chromosome("?") unless $obj->chromosome;
                $up   = $gbstart;
                $down = $gblength;
            }
            else {
                $message .= "No GenBank entry found for $gbaccn\n";
            }
        }
        next unless $obj;
        if ( $obj->sequence && $obj->start ne $obj->stop )
        { #need to check for duplicate accession names -- sometimes happens and major pain in the ass for other parts of the code
            my $accn  = $obj->accn;
            my $count = 0;
            foreach my $accn2 ( map { $_->{obj}->accn() } @sets ) {
                $accn2 =~ s/\*\*\d+\*\*$//;
                $count++ if $accn eq $accn2;
            }
            $accn .= "**" . ( $count + 1 ) . "**" if $count;
            $obj->accn($accn);
            push @sets,
              {
                obj           => $obj,
                file          => $obj->srcfile,
                accn          => $accn,
                rev           => $rev,
                up            => $up,
                down          => $down,
                spike_length  => $spike_len,
                reference_seq => $reference_seq,
                seq_num       => $display_order,
              };
        }
        else {
            #push @sets, {seq_num=>$i};
        }
    }
    @sets = sort { $a->{seq_num} <=> $b->{seq_num} } @sets;
    my $pm = new Parallel::ForkManager($MAX_PROC);
    foreach my $param (@coge_seqs) {
        $pm->start and next;
        write_fasta(%$param);
        $pm->finish;
    }
    $pm->wait_all_children;

    my $i = 1;
    foreach my $item ( sort { $a->{do} <=> $b->{do} } @gevo_link_seqs ) {
        $gevo_link .= ";accn$i=" . $item->{accn}   if $item->{accn};
        $gevo_link .= ";x$i=" . $item->{x}         if $item->{x};
        $gevo_link .= ";fid$i=" . $item->{fid}     if $item->{fid};
        $gevo_link .= ";dsid$i=" . $item->{dsid}   if $item->{dsid};
        $gevo_link .= ";dsgid$i=" . $item->{dsgid} if $item->{dsgid};
        $gevo_link .= ";gstid$i=" . $item->{gstid}
          if $item->{gstid} && !$item->{dsgid};
        $gevo_link .= ";chr$i=" . $item->{chr} if $item->{chr};
        $gevo_link .= ";dr$i" . "up=" . $item->{drup} if defined $item->{drup};
        $gevo_link .= ";dr$i" . "down=" . $item->{drdown}
          if defined $item->{drdown};
        $gevo_link .= ";gbaccn$i=" . $item->{gbaccn} if $item->{gbaccn};
        $gevo_link .= ";gbstart$i=" . $item->{gbstart}
          if $item->{gbaccn} && $item->{gbstart};
        $gevo_link .= ";gblength$i=" . $item->{gblength}
          if $item->{gbaccn} && $item->{gblength};
        $gevo_link .= ";rev$i=1"                 if $item->{rev};
        $gevo_link .= ";ref$i=" . $item->{ref};
        $gevo_link .= ";mask$i=" . $item->{mask} if $item->{mask};

        #$gevo_link .= ";do$i=".$item->{d} if $item->{do};
        $i++;
    }
    $i--;
    $gevo_link .= ";num_seqs=$i";
    $gevo_link .= ";hsp_overlap_limit=" . $hsp_overlap_limit
      if defined $hsp_overlap_limit;
    $gevo_link .= ";hsp_size_limit=" . $hsp_size_limit
      if defined $hsp_size_limit;
    unless ( @sets > 1 ) {
        $message .=
"Problem retrieving information, please check submissions.  At least 2 sequences should be specified.\n";
        return '', '', '', '', 0, '', '', $message;
    }

    # mdb added 11/12/13 - log it
    my $tiny_link = CoGe::Accessory::Web::get_tiny_link(url => $gevo_link);
    my $log = CoGe::Accessory::Web::log_history(
        db      => $coge,
        user_id => $USER->id,
        description => "Compared $num_seqs regions",
        page    => $PAGE_TITLE,
        link => $tiny_link,
    );

    my $job = CoGe::Accessory::Web::get_job(
        tiny_link => $tiny_link,
        title     => $PAGE_TITLE,
        user_id   => $USER->id,
        log_id    => $log->id,
        db_object => $coge
    );
    CoGe::Accessory::Web::schedule_job(job => $job);

    my $t2 = new Benchmark;

    # set up output page

    # run bl2seq
    my $analysis_reports;

    if ( $analysis_program eq "blastz" ) {
        $analysis_reports = run_blastz(
            sets        => \@sets,
            params      => $param_string,
            parser_opts => $parser_opts,
            comp_adj    => $comp_adj
        );
    }
    elsif ( $analysis_program eq "LAGAN" ) {
        $analysis_reports = run_lagan(
            sets        => \@sets,
            params      => $param_string,
            parser_opts => $parser_opts,
            comp_adj    => $comp_adj
        );
    }
    elsif ( $analysis_program eq "CHAOS" ) {
        $analysis_reports = run_chaos(
            sets        => \@sets,
            params      => $param_string,
            parser_opts => $parser_opts,
            comp_adj    => $comp_adj
        );
    }
    elsif ( $analysis_program eq "DiAlign_2" ) {
        ( $analysis_reports, $analysis_program, $message ) = run_dialign(
            sets        => \@sets,
            params      => $param_string,
            parser_opts => $parser_opts,
            comp_adj    => $comp_adj
        );
    }
    elsif ( $analysis_program eq "GenomeThreader" ) {
        ($analysis_reports) = run_genomethreader(
            sets        => \@sets,
            params      => $param_string,
            parser_opts => $parser_opts,
            comp_adj    => $comp_adj
        );
    }
    elsif ( $analysis_program eq "blastn" || $analysis_program eq "tblastx" ) {
        $analysis_reports = run_bl2seq(
            sets          => \@sets,
            params        => $param_string,
            parser_opts   => $parser_opts,
            blast_program => $analysis_program,
            comp_adj      => $comp_adj
        );
    }

    $analysis_reports = [] unless ref($analysis_reports) =~ /ARRAY/i;
    CoGe::Accessory::Web::write_log( $message, $cogeweb->logfile ) if $message;

    #sets => array or data for blast
    #blast_reports => array of arrays (report, accn1, accn2, parsed data)
    my $t3 = new Benchmark;
    initialize_sqlite();
    my $bitscore_cutoff = get_bit_score_cutoff(
        seq_length => $spike_len,
        match      => $opts{blast_match},
        mismatch   => $opts{blast_mmatch}
    ) if $spike_len;
    my @gfxs;

    foreach my $item (@sets) {
        my $obj      = $item->{obj};
        my $filename = $cogeweb->basefile . "_" . $item->{seq_num} . ".png";
        $filename = CoGe::Accessory::Web::check_filename_taint($filename);
        $item->{png_filename} = $filename;
        my $image = basename($filename);
        $item->{image} = $image;
    }
    my $count = 1;
    foreach my $item (@sets) {
        $pm->start and next;
        my $obj = $item->{obj};
        next unless $obj->sequence;
        if ($obj) {
            CoGe::Accessory::Web::write_log(
                "generating image ($count/"
                  . scalar @sets . ")for "
                  . $obj->accn . ": "
                  . $item->{png_filename},
                $cogeweb->logfile
            );
            $count++;
            my ($gfx) = generate_image(
                obj         => $item->{obj},
                start       => 1,
                stop        => length( $obj->sequence ),
                data        => $analysis_reports,
                iw          => $iw,
                fh          => $feat_h,
                show_gc     => $show_gc,
                show_nt     => $show_nt,
                show_cbc    => $show_cbc,
                hsp_labels  => $hsp_labels,
                feat_labels => $feat_labels,

                #			     hsp_limit=>$hsp_limit,
                #			     hsp_limit_num=>$hsp_limit_num,
                color_hsp        => $color_hsp,
                hsp_top          => $hsp_top,
                hsp_track        => $hsp_track,
                hsp_colors       => \@hsp_colors,
                hsp_single_color => $hsp_single_color
                ,    #just use the first hsp_color for all of them
                show_hsps_with_stop_codon => $show_hsps_with_stop_codon,
                hiqual                    => $hiqual,
		color_anchor		  => $color_anchor,
                padding                   => $padding,
                reverse_image             => $item->{rev},
                draw_model                => $draw_model,
                hsp_overlap_limit         => $hsp_overlap_limit,
                hsp_size_limit            => $hsp_size_limit,
                color_overlapped_features => $color_overlapped_features,
                hsp_overlap_length        => $hsp_overlap_length,
                show_cns                  => $show_cns,
                show_ofeat                => $show_ofeat,
                show_gene_space           => $show_gene_space,
                show_contigs              => $show_contigs,
                bitscore_cutoff           => $bitscore_cutoff,
                skip_feat_overlap_search  => $skip_feat_overlap_search,
                skip_hsp_overlap_search   => $skip_hsp_overlap_search,
                font_size                 => $font_size,
                analysis_program          => $analysis_program,
            );
            $gfx->generate_png( file => $item->{png_filename} );
            generate_image_db( set => $item, gfx => $gfx );
        }
        $pm->finish;
    }
    $pm->wait_all_children;

    #get combined height of all the images
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    my $query = qq{select sum(px_height) from image_info;};
    my $sth   = $dbh->prepare($query);
    $sth->execute;
    my ($frame_height) = $sth->fetchrow_array;
    $sth->finish();
    undef $sth;
    $dbh->disconnect();
    undef $dbh;
    my $t3_5 = new Benchmark;
    image_db_create_hsp_pairs();

    my $t4 = new Benchmark;

    #set up buttons for gobe

    my $gobe_buttons = qq{
<table>
<tr>
<td><span class='coge-button coge-button-sm' id="clear_lines" onclick="Gobe.clear()">Clear Connectors</span>

<td><span class='coge-button coge-button-sm drawline' id="set_lines" onclick="\$('.drawline').hide();\$('#set_wedges').show();\$('.lineopt').show();Gobe.set_connector('line')">Set connector as Lines</span>

<span style="display: none" class='coge-button coge-button-sm lineopt' id="set_wedges" onclick="\$('.drawline').show();\$('.lineopt').hide();Gobe.set_connector('wedge')">Set connector as Wedges</span>

<td><div class=lineopt style="float: left; display: none">
 <span class='coge-button coge-button-sm' id="">Line Width</span>
 <input type=textbox size=2 class="backbox line_width" id=line_width_val value=3 readonly>
 <span class='coge-button coge-button-sm' id="" onclick="update_line_width(1)">+</span>
 <span class='coge-button coge-button-sm' id="" onclick="update_line_width(-1)">-</span>
</div>
};
    #$gobe_buttons .=
#qq{<td><a href="javascript:void(0);" id="history_dialog_button" class='ui-button ui-corner-all ui-button-icon-left' onClick="save_GEvo_results()"><span class="ui-icon ui-icon-newwin"></span>Save Results</a>}
      #unless $USER->user_name eq 'public';
    $gobe_buttons .= "</table>";
    $html         .= $gobe_buttons;
    $html         .= qq{<DIV id=flashcontent></DIV>};
    $html .=
qq{<br><a href="http://genomevolution.org/wiki/index.php/Gobe" class="small" style="color: red" target=_new>Click here for help!</a>  <a href="http://get.adobe.com/flashplayer/" class="small" target=_new >No results?  Rerun by pressing "Run GEvo Analysis!" again.  Still no results? Try installing the latest version of Flash</a>.};
    $html .= $gobe_buttons;
    $html .= qq{<table class=small>};
    $html .= qq{<tr valign=top><td><span class=bold>Alignment reports</span>};

    #    my $stats_file = $cogeweb->basefile."_stats.txt";
    #    $stats_file = CoGe::Accessory::Web::check_filename_taint($stats_file);

    if ( $analysis_reports && @$analysis_reports ) {
        foreach my $item (@$analysis_reports) {
            my $report         = $item->[0];
            my $accn1          = $item->[1];
            my $accn2          = $item->[2];
            my $basereportname = basename($report);
            $basereportname = $TEMPURL . "/$basereportname\n";
            $html .=
"<div><A HREF=\"$basereportname\" target=_new>$accn1 versus $accn2</A></font></DIV>\n";
        }
    }
    else {
        $html .= "<divl>No alignment reports were generated</DIV>\n";
    }
    $html .= qq{<td class=dropmenu><td><span class=bold>Fasta files</span>};
    my $all_file = $cogeweb->basefile() . ".all.fasta";
    foreach my $item (@sets) {
        next unless $item->{file};
        my $basename = $TEMPURL . "/" . basename( $item->{file} );
        print STDERR "basename is undefined: $basename\n"
          if $basename =~ /defined/i;
        my $accn = $item->{accn};
        $html .=
          "<div><A HREF=\"$basename\" target=_new>$accn</A></font></DIV>\n";
        my $x;
        ( $x, $all_file ) = CoGe::Accessory::Web::check_taint($all_file);
        my $seq_file = $item->{file};
        ( $x, $seq_file ) = CoGe::Accessory::Web::check_taint($seq_file);
        my $cmd = "/bin/cat '$seq_file' >> '$all_file'";
        `$cmd`;
    }

    $html .=
        "<div><A HREF=\""
      . $TEMPURL . "/"
      . basename($all_file)
      . "\" target=_new>all sequences</A></font></DIV>\n";
    $html .=
qq{<td class="dropmenu"><td><span class="bold">Third party system annotation files</span>};
    foreach my $item (@sets) {
        my $anno_file = generate_annotation(%$item);
        next unless $anno_file;
        my $basename = $TEMPURL . "/" . basename($anno_file);
        my $accn     = $item->{accn};
        $html .=
qq{<div><a href = "http://genome.lbl.gov/vista/mvista/instructions.shtml">VISTA</a> Annotation: <A HREF="$basename" target=_new>$accn</A></font></DIV>\n};
    }
    my ( $synfile, $annofile ) = generate_mGSV_files( $cogeweb->sqlitefile );
    $synfile  =~ s/$TEMPDIR/$TEMPURL/;
    $annofile =~ s/$TEMPDIR/$TEMPURL/;
    $html .=
qq{<div><a href="http://cas-bioinfo.cas.unt.edu/mgsv/index.php" target=_new>mGSV</a> <a href="$annofile" target=_new>Annotation File</a></div>};
    $html .=
qq{<div><a href="http://cas-bioinfo.cas.unt.edu/mgsv/index.php" target=_new>mGSV</a> <a href="$synfile" target=_new>Synteny File</a></div>};

    $html .= qq{<td class="dropmenu"><td><span class="bold">Image Files</span>};
    foreach my $item (@sets) {
        my $png = $TEMPURL . "/" . basename( $item->{png_filename} );
        $html .=
          qq{<br><a href="$png" target=_new>} . $item->{obj}->accn . "</a>";
    }
    $html .= qq{<td class="dropmenu"><td><span class="bold">SQLite db</span>};
    my $dbname = $TEMPURL . "/" . basename( $cogeweb->sqlitefile );

    $html .= "<div><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
    $html .= qq{<td class="dropmenu"><td><span class="bold">Log File</span>};
    my $logfile = $TEMPURL . "/" . basename( $cogeweb->logfile );
    $html .= "<div><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
    $html .= qq{<td class="dropmenu"><td><span class="bold">Return to this analysis</span>};
    $html .= qq{<div id="tiny_link"></div>};
    #$html .= qq{<div><a href="GEvo_direct.pl?name=$basefilename" target=_new>Results only</a></div></td>}; # mdb removed 8/20/14 issue 467

    my (@ncbi_links);
    foreach my $item (@sets) {
        if ( $item->{'obj'}->ncbi_link ) {
            my $ncbi_link = $item->{'obj'}->ncbi_link;
            push @ncbi_links,
                qq{<div class="link" onclick="window.open('$ncbi_link')">}
              . $item->{'obj'}->accn
              . "</div>";
        }
    }
    if (@ncbi_links) {
        $html .= qq{<td>NCBI links\n};
        $html .= join( "\n", @ncbi_links );
    }

    $html .= qq{</table>};

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'box.tmpl' );
    my $results_name = "Results: $analysis_program";
    $results_name .=
      qq{ <span class="small">(spike sequence filter length: $spike_len)</span>}
      if $spike_len;
    $template->param( BODY     => $html );
    my $outhtml     = $template->output;
    my $t5          = new Benchmark;
    my $db_time     = timestr( timediff( $t2, $t1 ) );
    my $blast_time  = timestr( timediff( $t3, $t2 ) );
    my $image_time  = timestr( timediff( $t3_5, $t3 ) );
    my $db_hsp_time = timestr( timediff( $t4, $t3_5 ) );
    my $html_time   = timestr( timediff( $t5, $t4 ) );
    my $total_time  = timestr( timediff( $t5, $t1 ) );
    my $bench       = qq{
GEvo Benchmark: $DATE
  Time to get sequence                              : $db_time
  Time to run $analysis_program                     : $blast_time
  Time to generate images, maps, and sqlite database: $image_time
  Time to find and update sqlite database for HSPs  : $db_hsp_time
  Time to process html                              : $html_time
Total time                                          : $total_time

};

    #    print STDERR "\n",$gevo_link,"\n";
    #    print STDERR $bench if $BENCHMARK;
    #    print STDERR "GEvo run time: $total_time\n";
    CoGe::Accessory::Web::write_log( $bench,      $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "Finished!", $cogeweb->logfile );
    CoGe::Accessory::Web::write_log( "GEvo link: $gevo_link",
        $cogeweb->logfile );

    $job->update( {
        status => 2,
        end_time => \"current_timestamp"
    } ) if defined($job);

    #    CoGe::Accessory::Web::write_log("Tiny url: $tiny", $cogeweb->logfile);
    email_results(
        email         => $email_address,
        basefile      => $basefilename,
        full_gevo_url => $gevo_link
    ) if $email_address;
    $iw += 10;    #extra padding to make it easier to grab slider bars
    return $outhtml, $iw, $frame_height, $cogeweb->basefilename, scalar(@sets),
      $gevo_link, $basefilename, $message;

}

sub generate_mGSV_files {
    my $db = shift;
    $db = $TEMPDIR . "/" . $db unless $db =~ /^$TEMPDIR/;
    my $basefile = $db;
    $basefile =~ s/sqlite//;
    return "error, unable to find db: $db" unless -r $db;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=" . $db, "", "" );

    my $image_query = qq{select display_id, title from image_info};
    my %genome_info;
    my $count = 1;
    foreach my $item ( @{ $dbh->selectall_arrayref($image_query) } ) {

        #my ($name) = $item->[1];#=~ /^(\w{10})/;
        #$name =~ s/\s|\(|\)/_/g;
        my $name = "Org" . $count;
        $name .= " $item->[1]";
        $name =~ s/\s|\(|\)/_/g;
        $genome_info{ $item->[0] } = $name;
        $count++;
    }
    my $query =
qq{select name, type, bpmin, bpmax, image_id, pair_id, link, annotation, strand from image_data};

    my %annotations;
    my %synpairs;
    foreach my $item ( @{ $dbh->selectall_arrayref($query) } ) {

        #	print join ("\t", @$item[0..6]),"\n";
        if ( $item->[5] == -99 ) {
            my ($fid) = $item->[6] =~ /fid=(\d+)/;
            push @{ $annotations{ $item->[4] }{ $item->[1] } },
              {
                start  => $item->[2],
                stop   => $item->[3],
                fid    => $fid,
                strand => $item->[8],
                anno   => $item->[7],
              };

        }
        else {
            push @{ $synpairs{ $item->[0] } },
              {
                genome_id => $item->[4],
                start     => $item->[2],
                stop      => $item->[3],
                name      => $item->[0],
                link      => $item->[6],
                strand    => $item->[8],
                anno      => $item->[7],
              };
        }
    }
    my $synfile = $basefile . "syn_pairs.mGSV.txt";
    open( OUT, ">" . $synfile );
    print OUT "#",
      join( "\t",
        qw (org1   org1_start      org1_end        org2    org2_start      org2_end        score   evalue)
      ),
      "\n";

    foreach my $key ( sort keys %synpairs ) {
        my ( $item1, $item2 ) = @{ $synpairs{$key} };
        my $start = $item2->{start};
        my $stop  = $item2->{stop};
        ( $start, $stop ) = ( $stop, $start )
          if $item1->{strand} == -1 || $item2->{strand} == -1;
        print OUT join( "\t",
            $genome_info{ $item1->{genome_id} },
            $item1->{start}, $item1->{stop},
            $genome_info{ $item2->{genome_id} },
            $start, $stop, 10, 0 ),
          "\n";
    }
    close OUT;
    my $annofile = $basefile . "anno.mGSV.txt";
    open( OUT, ">$annofile" );
    print OUT "#",
      join( "\t",
        qw(org_id start   end     strand  feature_name    feature_value   track_name      track_shape     track_color)
      ),
      "\n";
    foreach my $genome_id ( sort keys %annotations ) {
        foreach my $type ( sort keys %{ $annotations{$genome_id} } ) {
            foreach my $item ( @{ $annotations{$genome_id}{$type} } ) {
                my $strand = $item->{strand} eq 1 ? "+" : "-";
                my ($feat) = $coge->resultset('Feature')->find( $item->{fid} );
                my ($name) = $feat->names if $feat;
                $name = "NA" unless $name;
                if ( $item->{fid} ) {
                    $name =
                        "<a href=\""
                      . $P->{SERVER}
                      . "/FeatView.pl?fid="
                      . $item->{fid} . "\">";
                    $name .= $name ? $name : $item->{fid};
                    $name .= "</a>";
                }
                print OUT join(
                    "\t",
                    $genome_info{$genome_id},
                    $item->{start},
                    $item->{stop},
                    $strand,
                    $name,
                    "",    #value
                    $type,
                    "arrow",
                    "green",
                  ),
                  "\n";
            }
        }
    }
    close OUT;
    return ( $synfile, $annofile );
}

sub generate_image {
    my %opts          = @_;
    my $gbobj         = $opts{obj};
    my $start         = $opts{start};
    my $stop          = $opts{stop};
    my $data          = $opts{data};
    my $iw            = $opts{iw} || 1600;
    my $fh            = $opts{fh} || 25;
    my $show_gc       = $opts{show_gc};
    my $show_nt       = $opts{show_nt};
    my $show_cbc      = $opts{show_cbc};
    my $reverse_image = $opts{reverse_image};
    my $hsp_labels    = $opts{hsp_labels};
    my $feat_labels   = $opts{feat_labels};

    #    my $hsp_limit = $opts{hsp_limit};
    #    my $hsp_limit_num = $opts{hsp_limit_num};
    my $hsp_top                   = $opts{hsp_top};
    my $hsp_single_color          = $opts{hsp_single_color};
    my $hsp_track                 = $opts{hsp_track};
    my $color_hsp                 = $opts{color_hsp};
    my $bitscore_cutoff           = $opts{bitscore_cutoff};
    my $hsp_colors                = $opts{hsp_colors};
    my $show_hsps_with_stop_codon = $opts{show_hsps_with_stop_codon};
    my $hiqual                    = $opts{hiqual};
    my $color_anchor		  = $opts{color_anchor};
    my $padding                   = $opts{padding} || 5;
    my $draw_model                = $opts{draw_model};
    my $hsp_overlap_limit         = $opts{hsp_overlap_limit};
    my $hsp_size_limit            = $opts{hsp_size_limit};
    my $color_overlapped_features = $opts{color_overlapped_features};
    my $hsp_overlap_length        = $opts{hsp_overlap_length};
    my $show_cns                  = $opts{show_cns};
    my $show_ofeat                = $opts{show_ofeat};
    my $show_gene_space           = $opts{show_gene_space};
    my $show_contigs              = $opts{show_contigs};
    my $graphic                   = new CoGe::Graphics;
    my $gfx                       = new CoGe::Graphics::Chromosome;
    my $skip_feat_overlap_search  = $opts{skip_feat_overlap_search};
    my $skip_hsp_overlap_search   = $opts{skip_hsp_overlap_search};
    my $font_size                 = $opts{font_size};
    my $analysis_program          = $opts{analysis_program};
    $skip_hsp_overlap_search = 1 unless defined $skip_hsp_overlap_search;
    $graphic->initialize_c(
        c                 => $gfx,
        iw                => $iw,
        start             => $start,
        stop              => $stop,
        draw_chr          => 1,
        draw_ruler        => 1,
        draw_chr_end      => 1,
        feature_height    => $fh,
        chr_length        => length( $gbobj->sequence ),
        fill_labels       => 1,
        minor_tick_labels => 1,
        feature_labels    => $feat_labels,
        draw_hi_qual      => $hiqual,
        padding           => $padding,
    );
    $gfx->overlap_adjustment(1);
    $gfx->top_padding(15);
    $gfx->skip_duplicate_features(1);
    $gfx->DEBUG(0);
    $gfx->major_tick_labels(0);
    $graphic->process_nucleotides(
        c      => $gfx,
        seq    => $gbobj->sequence,
        layers => { gc => $show_gc, nt => $show_nt },
        n_color => [ 255, 155, 0 ]
    );
    process_hsps(
        c          => $gfx,
        data       => $data,
        accn       => $gbobj->accn,
        rev        => $reverse_image,
        seq_length => length( $gbobj->sequence ),
        hsp_labels => $hsp_labels,

#			     hsp_limit=>$hsp_limit,
#			     hsp_limit_num=>$hsp_limit_num, #number of hsps for which to draw labels
        gbobj                     => $gbobj,
        color_hsp                 => $color_hsp,
        hsp_top                   => $hsp_top,
        hsp_single_color          => $hsp_single_color,
        hsp_track                 => $hsp_track,
        colors                    => $hsp_colors,
        show_hsps_with_stop_codon => $show_hsps_with_stop_codon, #tblastx option
        hsp_overlap_limit         => $hsp_overlap_limit
        , #if an hsp overlaps with other hsps this many times, don't draw (delete from graphics object)
        hsp_size_limit => $hsp_size_limit,    #minimum length of hsp to be shown
        hsp_overlap_length => $hsp_overlap_length
        , #length that an hsp and a genomic feature must overlap in order for the feature to be flagged as 'hit' by an HSP
        bitscore_cutoff           => $bitscore_cutoff,
        color_overlapped_features => $color_overlapped_features,
        skip_overlap_search       => $skip_hsp_overlap_search,
        font_size                 => $font_size,
        analysis_program          => $analysis_program,
    );

    my ($feat_counts) = process_features(
        c                         => $gfx,
        obj                       => $gbobj,
        start                     => $start,
        stop                      => $stop,
        skip_overlap_search       => $skip_feat_overlap_search,
        draw_model                => $draw_model,
        color_overlapped_features => $color_overlapped_features,
        cbc                       => $show_cbc,
        cns                       => $show_cns,
        ofeat                     => $show_ofeat,
        gene_space                => $show_gene_space,
        show_contigs              => $show_contigs,
        feat_labels               => $feat_labels,
        font_size                 => $font_size,
	color_anchor		  => $color_anchor,
    );

    return ($gfx);
}

sub initialize_sqlite {
    my %opts   = @_;
    my $dbfile = $cogeweb->sqlitefile;
    return if -r $dbfile;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfile", "", "" );
    unless (defined $dbh) {
        print STDERR "Gevo.pl ERROR connecting to sqlite file\n";
        return;
    }
    my $create = qq{
CREATE TABLE image_data
(
id INTEGER PRIMARY KEY,
name varchar,
type varchar,
xmin integer,
xmax integer,
ymin integer,
ymax integer,
bpmin integer,
bpmax integer,
image_track integer,
image_id integer,
pair_id integer,
link varchar,
annotation text,
color varchar,
strand varchar
)
};
    $dbh->do($create);
    $dbh->do('CREATE INDEX type ON image_data (type)');
    $dbh->do('CREATE INDEX image_id ON image_data (image_id)');
    $dbh->do('CREATE INDEX x_idx ON image_data (xmin,xmax)');
    $dbh->do('CREATE INDEX y_idx ON image_data (ymin,ymax)');
    $dbh->do('CREATE INDEX image_track ON image_data (image_track)');
    $dbh->do('CREATE INDEX pair_id ON image_data (pair_id)');

    #id INTEGER PRIMARY KEY AUTOINCREMENT,
    $create = qq{
CREATE TABLE image_info
(
id INTEGER,
display_id INTEGER,
iname varchar,
title varchar,
px_width integer,
bpmin integer,
bpmax integer,
dsid integer,
chromosome varchar,
reverse_complement integer,
px_height integer
)
};

#TODO: make sure to populate the bpmin and pbmax and image width!!!!!!!!!!!!!!!!!
    $dbh->do($create);
    $dbh->do('CREATE INDEX id ON image_info (id)');
    system "/bin/chmod +rw '$dbfile'";
    $dbh->disconnect();
    undef $dbh;
}

sub generate_image_db {
    my %args = @_;
    my $gfx  = $args{gfx};
    my $set  = $args{set};
    next unless $gfx;
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    $dbh->do("begin exclusive transaction") or die $dbh->errstr;
    my $image = $set->{image};
    my $accn  = $set->{accn};
    my $title;
    $title = $set->{obj}->organism() if $set->{obj}->organism();
    $title .= " " if $title;
    $title .= $accn;
    $title .=
        " (chr: "
      . $set->{obj}->chromosome . " "
      . $set->{obj}->start . "-"
      . $set->{obj}->stop . ")"
      if defined $set->{up};
    $title .= qq! Reverse Complement! if $set->{rev};
    my $width  = $gfx->image_width;
    my $height = $gfx->image_height;
    my $dsid   = $set->{obj}->dataset;
    my ($chr)  = $set->{obj}->chromosome;    # =~ /(\d+)/;
    $chr  = "NULL" unless defined $chr;
    $dsid = "NULL" unless $dsid;
    my $image_start = $set->{obj}->start;
    my $image_stop  = $set->{obj}->stop;
    my $image_id    = $set->{seq_num};
    my $rc          = $set->{rev};
    $rc = 0 unless defined $rc;
    my $statement = qq{
INSERT INTO image_info (id, display_id, iname, title, px_width, px_height, dsid, chromosome, bpmin, bpmax, reverse_complement) values ($image_id, $image_id, "$image", "$title", $width, $height, "$dsid", "$chr", $image_start, $image_stop, $rc)
};
    my $try           = 1;
    my $run_statement = $dbh->do($statement);
    print STDERR $cogeweb->sqlitefile, "\n", $statement, "\n"
      unless $run_statement;

    unless ( $run_statement || $try > 100 ) {
        sleep(1);
        $run_statement = $dbh->do($statement);
        $try++;
    }

#    print STDERR $cogeweb->sqlitefile, "\n", $statement unless $dbh->do($statement);
    foreach my $feat ( $gfx->get_feats ) {

        #		print STDERR Dumper $feat unless $feat->color;
        if ( $feat->fill ) {
            next unless $feat->{"anchor"} || $feat->type eq "anchor";
        }
        next
          unless $feat->image_coordinates
              || $feat->{"anchor"}
              || $feat->type eq "anchor";
        my $type    = $feat->type;
        my $pair_id = "-99";
        my $coords  = $feat->image_coordinates;
        $coords =~ s/\s//g;
        next if $feat->type eq "unknown";
        my $name = $feat->type =~ /HSP/i ? $feat->alt : $feat->label;
        $name .= "_" . $feat->type;    # unless $name;

        my $color = "NULL";
        if ( $feat->type =~ /HSP/ ) {
            $color = "#";
            foreach my $c ( @{ $feat->color } ) {
                $color .= sprintf( "%X", $c );
            }
        }

        #generate link
        my $link = $feat->link;
        $link = " " unless $feat->link;
        $link =~ s/'//g if $link;

        #generate image track
        my $image_track = $feat->track;
        $image_track = "-" . $image_track if $feat->strand =~ /-/;
        my ( $xmin, $ymin, $xmax, $ymax ) = ( -1, -1, -1, -1 );
        ( $xmin, $ymin, $xmax, $ymax ) = split( /,/, $coords ) if $coords;
        $xmin++;
        $xmax++;
        my $anno = $feat->description;

        #	    $anno =~ s/'|"//g;
        #	$anno =~ s/'//g if $anno;
        #	$anno =~ s/<br\/?>/&#10;/ig if $anno;
        #	$anno =~ s/\n/&#10;/g if $anno;
        #	$anno =~ s/[\[\]\(\)]/ /g if $anno;
        #	$anno = " " unless $anno;
        my $start      = $feat->start;
        my $stop       = $feat->stop;
        my $strand     = $feat->strand;
        my $length_nt  = $stop - $start + 1;
        my $length_pix = $xmax - $xmin;
        next if !$feat->{anchor} && ( $length_nt == 0 );
        $type = "anchor" if $feat->{anchor};
        my $bpmin = $set->{obj}->start + $start - 1;
        my $bpmax = $set->{obj}->start + $stop - 1;
        $anno =~ s/'//g;    #need to remove these or the sql will get messed up
        $statement = qq{
INSERT INTO image_data (name, type, xmin, xmax, ymin, ymax, bpmin,bpmax,image_id, image_track,pair_id, link, annotation, color, strand) values ("$name", "$type", $xmin, $xmax, $ymin, $ymax, $bpmin, $bpmax,$image_id, "$image_track",$pair_id, '$link', '$anno', '$color', '$strand')
};
        my $try           = 1;
        my $run_statement = $dbh->do($statement);
        print STDERR $cogeweb->sqlitefile, "\n", $statement, "\n"
          unless $run_statement;

        unless ( $run_statement || $try > 100 ) {
            sleep(1);
            $run_statement = $dbh->do($statement);
            $try++;
        }
    }
    $dbh->do("commit transaction") or die $dbh->errstr;
    $dbh->disconnect();
    undef $dbh;
}

sub image_db_create_hsp_pairs {
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    $dbh->do("begin exclusive transaction") or die $dbh->errstr;
    my $query = qq{select id,name from image_data where type = "HSP"};
    my %ids;
    foreach my $item ( @{ $dbh->selectall_arrayref($query) } ) {
        push @{ $ids{ $item->[1] } }, $item->[0];
    }
    foreach my $name ( keys %ids ) {
        if ( @{ $ids{$name} } != 2 ) {
            print STDERR "Problem with $name.  Retrieved", Dumper $ids{$name}
              if $DEBUG;
            next;
        }
        my ( $id1, $id2 ) = @{ $ids{$name} };
        my $statement = "update image_data set pair_id = $id1 where id = $id2";

        my $try           = 1;
        my $run_statement = $dbh->do($statement);
        print STDERR $cogeweb->sqlitefile, "\n", $statement, "\n"
          unless $run_statement;
        unless ( $run_statement || $try > 100 ) {
            sleep(1);
            $run_statement = $dbh->do($statement);
            $try++;
        }

        #	$dbh->do($statement);
        $statement     = "update image_data set pair_id = $id2 where id = $id1";
        $try           = 1;
        $run_statement = $dbh->do($statement);
        print STDERR $cogeweb->sqlitefile, "\n", $statement, "\n"
          unless $run_statement;
        unless ( $run_statement || $try > 100 ) {
            sleep(1);
            $run_statement = $dbh->do($statement);
            $try++;
        }

        #	$dbh->do($statement);
    }
    $dbh->do("commit transaction") or die $dbh->errstr;
    $dbh->disconnect();
    undef $dbh;
}

sub process_features {

    #process features
    my %opts                      = @_;
    my $c                         = $opts{c};
    my $obj                       = $opts{obj};
    my $start                     = $opts{start};
    my $stop                      = $opts{stop};
    my $skip_overlap_search       = $opts{skip_overlap_search};
    my $draw_model                = $opts{draw_model};
    my $color_overlapped_features = $opts{color_overlapped_features};
    my $cbc                       = $opts{cbc};
    my $show_cns                  = $opts{cns};
    my $show_ofeat                = $opts{ofeat};
    my $show_gene_space           = $opts{gene_space};
    my $color_anchor              = $opts{color_anchor};

    # TODO: remove the 1. here for testing.
    my $show_contigs = $opts{show_contigs};
    my $feat_labels  = $opts{feat_labels};
    my $font_size    = $opts{font_size} || 13;
    my $accn         = $obj->accn;
    my $track        = 1;
    my %feat_counts;
    unless ( ref $obj ) {
        warn "Possible problem with the object in process_features.  Returning";
        return 0;
    }
    my %prior_feat;
    my $label_location;
    foreach
      my $feat ( sort { $a->strand cmp $b->strand || $a->start <=> $b->start }
        $obj->get_features() )
    {
        if ( !$label_location ) {
            $label_location = "bot";
        }
        elsif ( $label_location eq "bot" ) {
            $label_location = "top";
        }
        else {
            $label_location = "";
        }
        my $anchor;
        if (
            (
                   ref( $feat->qualifiers ) =~ /hash/i
                && $feat->qualifiers->{type}
                && $feat->qualifiers->{type} eq "anchor"
            )
            || $feat->type eq "anchor"
          )
        {
            $anchor = $feat;
        }
        else {
            next if $feat->start > $stop;
            next if $feat->stop < $start;
        }
        my $f;
        my $type = $feat->type;
        my ($name) =
          sort { length($b) <=> length($a) || $a cmp $b }
          @{ $feat->qualifiers->{names} }
          if ref( $feat->qualifiers ) =~ /hash/i;

#	my $anchor = $feat if ref ($feat->qualifiers) =~ /hash/i && $feat->qualifiers->{type} && $feat->qualifiers->{type} eq "anchor";

        if ( $type =~ /pseudogene/i ) {
            next unless $draw_model eq "full";
            $f = CoGe::Graphics::Feature::Gene->new();
            $f->color( [ 255, 33, 0, 50 ] );
            $f->order($track);
            $f->overlay(1);
            $f->mag(0.5);
        }

        if ( $type =~ /gene prediction/i ) {
            next unless $draw_model eq "full";
            $f = CoGe::Graphics::Feature::Gene->new();
            $f->color( [ 50, 50, 255, 100 ] );
            $f->order($track);
            $f->overlay(1);
            $f->mag(1);
        }

        elsif ( $type =~ /Gene$/i ) {
            next
              unless $draw_model eq "full" || $draw_model eq "gene" || $anchor;
            $f = CoGe::Graphics::Feature::Gene->new();

            #$f->color([255,0,0,50]);
            $f->color( [ 219, 219, 219, 50 ] );
            $f->order($track);
            $f->overlay(1);
            $f->mag(0.5);
            $f->label( join( "\t", @{ $feat->qualifiers->{names} } ) )
              if $feat_labels && !$f->label && $feat->qualifiers->{names};
            $f->label_location($label_location)
              if $feat_labels && $feat_labels eq "staggered";
        }

    #Hi Eric, I added this feature type to display exons which were identifed by
    #Co-annotation but aren't part of the official gene models -J
        elsif ( $type =~ /Freeling_predicted_CDS/i ) {
            $f = CoGe::Graphics::Feature::Gene->new();
            $f->color( [ 150, 0, 150, 50 ] );
            $f->order($track);
            $f->overlay(3);
            $f->mag(1);

            #	    $f->label("C");
            $f->label_location($label_location)
              if $feat_labels && $feat_labels eq "staggered";
        }

        elsif ( $type =~ /CDS/i ) {
            next unless $draw_model eq "full" || $draw_model eq "CDS";
            $f = CoGe::Graphics::Feature::Gene->new();
            $f->color( [ 0, 255, 0, 50 ] );
            $f->order($track);
            $f->overlay(3);
            if ($accn) {
                foreach my $name ( @{ $feat->qualifiers->{names} } ) {
                    my $cleaned_name = $name;
                    $cleaned_name =~ s/[\(\)]//g;
                    $cleaned_name =~ s/[\|]/\\\|/g;
                    my $tmp = $accn;
                    $tmp =~ s/\*\*\d+\*\*$//;
                    if ( $tmp =~ /^$cleaned_name\(?\d*\)?$/i ) {
                        $f->color( [ 255, 255, 0 ] ) if $color_anchor;
                        $f->label($name) if $feat_labels;
                    }
                }
            }

            if ($cbc) {

                #my $seq = $feat->genomic_sequence;
                #$f->sequence($seq);
                my $seq;
                my $seq_len = length( $obj->sequence );
                foreach
                  my $block ( sort { $a->[0] <=> $b->[0] } @{ $feat->blocks } )
                {
                    $block->[0] = 1
                      unless $block->[0];    #in case $block is set to 0
                    my $tmp_start = $block->[0];
                    next if $tmp_start > $seq_len;
                    my $tmp_stop = $block->[1];
                    $tmp_stop = $seq_len if $block->[1] > $seq_len;
                    $seq .= substr(
                        $obj->sequence,
                        $tmp_start - 1,
                        $tmp_stop - $tmp_start + 1
                    );
                }
                if ( $feat->location =~ /complement/ ) {
                    $seq = reverse uc($seq);
                    $seq =~ tr /ATCG/TAGC/;
                }
                $f->sequence($seq);
                $f->color_by_codon(1);
            }

        }
        elsif ( $type =~ /repeat/i ) {
            next unless $draw_model eq "full";
            $f = CoGe::Graphics::Feature::Outline->new(
                {
                    start => $feat->blocks->[0][0],
                    stop  => $feat->blocks->[0][1]
                }
            );
            $f->color( [ 0, 0, 255 ] );
            $f->order($track);
            $f->overlay(4);
        }

        elsif ( $type =~ /mrna/i ) {
            next unless $draw_model eq "full" || $draw_model eq "mRNA";
            $f = CoGe::Graphics::Feature::Gene->new();
            $f->color( [ 0, 0, 255, 50 ] );
            $f->order($track);
            $f->overlay(2);
            $f->mag(0.75);
        }

        elsif ( $type =~ /rna/i ) {
            next unless $draw_model eq "full" || $draw_model eq "RNA";
            $f = CoGe::Graphics::Feature::Gene->new();
            $f->color( [ 200, 200, 200, 50 ] );
            $f->order($track);
            $f->overlay(2);
            if ($accn) {
                foreach my $name ( @{ $feat->qualifiers->{names} } ) {
                    if ( $accn =~ /^$name\(?\d*\)?$/i ) {
                        $f->color( [ 255, 255, 0 ] );
                        $f->label($name) if $feat_labels;
                    }
                }
            }

        }
        elsif ( $type =~ /anchor/ ) {
            $f = CoGe::Graphics::Feature::NucTide->new(
                {
                    type   => 'anchor',
                    start  => $feat->blocks->[0][0],
                    strand => 1
                }
            );
            $f->color( [ 0, 0, 255 ] );
            $f->type($type);
            $f->description( $feat->annotation );
            $f->force_draw(1);
            $c->add_feature($f);
            $f = CoGe::Graphics::Feature::NucTide->new(
                {
                    type   => 'anchor',
                    start  => $feat->blocks->[0][0],
                    strand => -1
                }
            );
            $f->color( [ 0, 0, 255 ] );
            $f->type($type);
            $f->description( $feat->annotation );
            $f->force_draw(1);
            $c->add_feature($f);
            next;
        }
        elsif ( $show_cns && $type =~ /cns/i ) {
            my $strand = 1;
            $strand = -1 if $feat->location =~ /complement/;

            $f = CoGe::Graphics::Feature::HSP->new(
                {
                    start  => $feat->blocks->[0][0],
                    stop   => $feat->blocks->[0][1],
                    strand => $strand
                }
            );
            if ( $strand == 1 ) {
                $f->color( [ 99, 0, 99 ] );
            }
            else {
                $f->color( [ 0, 99, 99 ] );
            }

            my $order = 1;

            #	    $order = 0 if $feat->location =~ /complement/;

            $f->order($order);
            $f->overlay(10);
            $f->type($type);
            $f->force_draw(1);
            $f->description( $feat->annotation );
            $c->add_feature($f);
            next;
        }
        elsif ( $show_gene_space && $type =~ /gene_space/i ) {
            $f = CoGe::Graphics::Feature::HSP->new(
                {
                    start => $feat->blocks->[0][0],
                    stop  => $feat->blocks->[0][1]
                }
            );

#$f = CoGe::Graphics::Feature::Line->new({start=>$feat->blocks->[0][0], stop=>$feat->blocks->[0][1]});

            $f->color( [ 255, 255, 102 ] );

            my $order = 1;
            $order = 0 if $feat->location =~ /complement/;

            $f->order($order);

            $f->overlay(-2);
            $f->transparency(0);
            $f->type($type);
            $f->description( $feat->annotation );
            $c->add_feature($f);
            next;
        }
        elsif ( $show_contigs && $type =~ /contig/i ) {
            $f = CoGe::Graphics::Feature::Block->new(
                {
                    start => $feat->blocks->[0][0],
                    stop  => $feat->blocks->[0][1]
                }
            );
            $f->_initialize();
            $f->color( [ 255, 0, 0 ] );

            #$f->bgcolor([255,255,255]);
            $f->order( $feat->location =~ /complement/ ? 0 : 1 );

            $f->overlay(-1);

            #	    $f->type($type);
            #	    $f->description($feat->annotation);
            $f->label( join( "\t", @{ $feat->qualifiers->{names} } ) )
              if $feat->qualifiers->{names};

            #	    $c->add_feature($f);
            #	    next;
        }
        elsif ($show_ofeat) {    #show everything else
            #print STDERR $type, "\n";
            my $strand = 1;
            $strand = -1 if $feat->location =~ /complement/;

            $f = CoGe::Graphics::Feature::HSP->new(
                {
                    start  => $feat->blocks->[0][0],
                    stop   => $feat->blocks->[0][1],
                    strand => $strand
                }
            );
            $f->color( [ 255, 150, 66 ] );
            my $order = 1;

            #	    $order = 0 if $feat->location =~ /complement/;

            $f->order($order);
            $f->overlay(-1);
            $f->type($type);
            $f->force_draw(1);
            $f->description( $feat->annotation );
            $c->add_feature($f);

            #			next;
        }

        else {

            #	    print STDERR "Didn't draw feature_type: ", $type,"\n";
        }

     #need to create an anchor if this feature is an anchor, but not to be drawn
        if ( $anchor && !$f ) {
            $f = CoGe::Graphics::Feature->new(
                {
                    type   => 'anchor',
                    start  => $feat->blocks->[0][0],
                    stop   => $feat->blocks->[0][1],
                    strand => 1
                }
            );
            $f->{anchor} = 1;
            $f->order(0);
            $f->skip_overlap_search(1);
            $f->description("auto generated anchor for GEvo");
            $c->add_feature($f);
            $f->force_draw(1);
            next;
        }
        next unless $f;
        $f->color( [ 255, 0, 255 ] )
          if $color_overlapped_features && $feat->qualifiers->{overlapped_hsp};
        my $strand = 1;
        $strand = -1 if $feat->location =~ /complement/;
        $strand = $f->strand if $f->strand;
        if ( ref($f) =~ /Gene/ ) {
            foreach my $block ( @{ $feat->blocks } ) {
                $block->[0] = 1 unless $block->[0];  #in case $block is set to 0
                $f->add_segment( start => $block->[0], stop => $block->[1] );
                $f->strand($strand);
                print STDERR "\t", join( "-", @$block ), "\n" if $DEBUG;
            }
        }
        print STDERR $name, "\n\n" if $DEBUG;
        $f->font_size($font_size);
        $f->type($type);
        $f->description( $feat->annotation );
        if ( $feat->qualifiers->{id} ) {
            $f->link( "FeatList.pl?fid="
                  . $feat->qualifiers->{id}
                  . ";gstid="
                  . $obj->gstid );
        }
        elsif ( $feat->qualifiers->{names} ) {
            my $names = $feat->qualifiers->{names};
            $f->link( "FeatView.pl?accn=" . $names->[0] )
              if ref($names) =~ /array/i && $names->[0];
        }
        $f->skip_overlap_search($skip_overlap_search);
        if ($anchor) {
            $f->{anchor} =
              1;    #illeagle creation of private variable in feature object!!!!
            $f->force_draw(1);
        }

#	$f->label(join ("\t", @{$feat->qualifiers->{names}})) if !$f->label && $feat->qualifiers->{names};
        $c->add_feature($f);

        if ( $prior_feat{ $feat->type }{$strand} ) {
            if ( $feat->start < $prior_feat{ $feat->type }{$strand}->stop ) {
                $prior_feat{ $feat->type }{$strand} = $feat;
                next;
            }
        }
        $feat_counts{ $feat->type }{count}++;
        $feat_counts{ $feat->type }{overlap}++
          if $feat->qualifiers->{overlapped_hsp};
        $prior_feat{ $feat->type }{$strand} = $feat;
    }
    return ( \%feat_counts );
}

sub process_hsps {
    my %opts       = @_;
    my $c          = $opts{c};
    my $data       = $opts{data};
    my $accn       = $opts{accn};
    my $reverse    = $opts{rev};
    my $seq_len    = $opts{seq_length};
    my $hsp_labels = $opts{hsp_labels};

    #    my $hsp_limit = $opts{hsp_limit};
    #    my $hsp_limit_num = $opts{hsp_limit_num};
    my $gbobj                     = $opts{gbobj};
    my $eval_cutoff               = $opts{eval_cutoff};
    my $bitscore_cutoff           = $opts{bitscore_cutoff};
    my $hsp_top                   = $opts{hsp_top};
    my $hsp_single_color          = $opts{hsp_single_color};
    my $hsp_track                 = $opts{hsp_track};
    my $color_hsp                 = $opts{color_hsp};
    my $colors                    = $opts{colors};
    my $hsp_overlap_limit         = $opts{hsp_overlap_limit};
    my $hsp_size_limit            = $opts{hsp_size_limit};
    my $hsp_overlap_length        = $opts{hsp_overlap_length};
    my $color_overlapped_features = $opts{color_overlapped_features};
    my $show_hsps_with_stop_codon = $opts{show_hsps_with_stop_codon};
    my $skip_overlap_search       = $opts{skip_overlap_search};
    my $font_size                 = $opts{font_size} || 13;
    my $analysis_program          = $opts{analysis_program};

    my $seq_type = $analysis_program =~ /tblastx/i ? " aa" : " nt";

#to reverse hsps when using genomic sequences from CoGe, they need to be drawn on the opposite strand than where blast reports them.  This is because CoGe::Graphics has the option of reverse drawing a region.  However, the sequence fed into blast has already been reverse complemented so that the HSPs are in the correct orientation for the image.  Thus, if the image is reverse, they are drawn on the wrong strand.  This corrects for that problem.   Sorry for the convoluted logic, but it was the simplest way to substantiate this option
    my $i     = 0;
    my $track = scalar @$data + 1;
    $track = 2 if $hsp_track;    #all HSPs to be drawn in same track
    foreach my $item (@$data) {
        my $report = $item->[0];
        my $accn1  = $item->[1];
        my $accn2  = $item->[2];
        my $blast  = $item->[3];
        unless ( $accn eq $accn1 || $accn eq $accn2 ) {
            $track-- unless $hsp_track;
        }
    }
    my @feats;
    my %overlap_count
      ; #count number of features that overlap an HSP if features overlapped by HSPs are to be colored
    foreach my $item (@$data) {
        my $report = $item->[0];
        my $accn1  = $item->[1];
        my $accn2  = $item->[2];
        my $blast  = $item->[3];
        unless ( $accn eq $accn1 || $accn eq $accn2 ) {
            $i++;
            next;
        }
        my ( $accna, $accnb ) =
          $accn1 eq $accn ? ( $accn1, $accn2 ) : ( $accn2, $accn1 );
        my %feats;
        foreach my $hsp ( @{ $blast->hsps } ) {
            next if defined $eval_cutoff && $hsp->eval > $eval_cutoff;
            next
              if defined $bitscore_cutoff
                  && defined $hsp->score
                  && $hsp->score < $bitscore_cutoff;
            my $color = $hsp_single_color ? $colors->[0] : $colors->[$i];
            my $skip = 0;
            next if $hsp_size_limit && $hsp->length < $hsp_size_limit;
            if ( $show_hsps_with_stop_codon
                && ( $hsp->qalign =~ /\*/ || $hsp->salign =~ /\*/ ) )
            {
                for my $i ( 0 .. ( length( $hsp->qalign ) - 1 ) ) {
                    my $chr1 = substr( $hsp->qalign, $i, 1 );
                    my $chr2 = substr( $hsp->salign, $i, 1 );
                    next unless $chr1 eq "*" || $chr2 eq "*";
                    $skip = 1
                      unless ( $chr1 eq "-" && $chr2 eq "*" )
                      || ( $chr2 eq "-" && $chr1 eq "*" );
                }
            }
            next if $skip;
            my ( $start, $stop, $seq );
            if ( $accn1 eq $accn ) {
                $start = $hsp->qb;
                $stop  = $hsp->qe;
                $seq   = $hsp->qalign;
            }
            elsif ( $accn2 eq $accn ) {
                $start = $hsp->sb;
                $stop  = $hsp->se;
                $seq   = $hsp->salign;
            }
            next unless $start && $stop;
            if ( $stop < $start ) {
                my $tmp = $start;
                $start = $stop;
                $stop  = $tmp;
            }
            my $strand = $hsp->strand =~ /-/ ? "-1" : 1;

            #find features in gbobj that overlap hsp
            if (
                $color_overlapped_features
                && (   ( $hsp->qe - $hsp->qb + 1 ) >= $hsp_overlap_length
                    || ( $hsp->se - $hsp->sb + 1 ) >= $hsp_overlap_length )
              )
            {
                foreach my $feat ( sort { $a->start <=> $b->start }
                    $gbobj->get_features( start => $start, stop => $stop ) )
                {
                    next
                      if $feat->stop < $start;   #feature ends before HSP begins
                    next
                      if $feat->start > $stop;   #feature begins after HSP stops
                        #lets skip it if it only has minimal end overlap
                    if (   abs( $start - $feat->stop ) <= 3
                        || abs( $stop - $feat->start ) <= 3 )
                    {
                        next;
                    }
                    $overlap_count{ $feat->type }++
                      unless $feat->qualifiers->{overlapped_hsp};
                    $feat->qualifiers->{overlapped_hsp} = 1;
                }
            }

            my $f;
            if ( $hsp->link_as_gene ) {
                $f =
                    $feats{ $hsp->number }
                  ? $feats{ $hsp->number }
                  : CoGe::Graphics::Feature::Gene->new(
                    { type => 'HSP', no_3D => 1, overlay => 2 } );
                $feats{ $hsp->number } = $f;
                $f->add_segment( start => $start, stop => $stop );

                #add base arrow for looks
                if ( $feats{ $hsp->number . "base" } ) {
                    my $base_arrow = $feats{ $hsp->number . "base" };
                    $base_arrow->start($start) if $start < $base_arrow->start;
                    $base_arrow->stop($stop)   if $stop > $base_arrow->stop;
                }
                else {
                    my $base_arrow = CoGe::Graphics::Feature::Gene->new(
                        {
                            no_3D   => 1,
                            overlay => 1,
                            order   => $track,
                            strand  => $hsp->strand
                        }
                    );
                    $base_arrow->add_segment( start => $start, stop => $stop );
                    $base_arrow->mag(0.5);
                    $base_arrow->description('this is just for looks');
                    $base_arrow->type("base_arrow");
                    $base_arrow->color( [ 200, 200, 200 ] );

                    #$c->add_feature($base_arrow);
                    push @feats, $base_arrow;
                    $feats{ $hsp->number . "base" } = $base_arrow;
                }
            }
            else {
                $f = CoGe::Graphics::Feature::HSP->new(
                    { start => $start, stop => $stop } );
            }

            $f->color($color);
            $f->order($track);
            $f->strand($strand);
            $f->strand(1) if $hsp_top;
            $f->color_matches($color_hsp) if ref($f) =~ /HSP/;

            #	    if ($hsp_limit)
            #	      {
            #		$f->label($hsp->number) if $hsp->number <= $hsp_limit_num;
            #	      }
            #	    else
            #	      {
            $f->label( $hsp->number )
              if $hsp_labels;    #add label to HSP if they are to be drawn

            #	      }
            $f->alt( join( "-", $hsp->number, $accn1, $accn2 ) );
            my ( $ab_start, $ab_stop ) =
              $reverse
              ? ( $gbobj->stop - $stop + 1, $gbobj->stop - $start + 1 )
              : ( $gbobj->start + $start - 1, $gbobj->start + $stop - 1 )
              ;                  #adjust hsp position to real-world coords
            my $desc = "<table><tr>";
            $desc .= "<tr>" . "<td>HSP:<td>" . $hsp->number if $hsp->number;
            if ( $blast->query ) {
                $desc .= qq{  <span class="small">(};
                $desc .= $blast->query;
                $desc .= "-";
                $desc .= $blast->subject;
                $desc .= ")</span>";
            }
            $desc .=
                "<tr><td>Location:<td>"
              . commify($ab_start) . "-"
              . commify($ab_stop) . " ("
              . $hsp->strand . ")"
              if $ab_start && $ab_stop && $hsp->strand;
            $desc .= "<tr><td>Match:<td>" . $hsp->match . $seq_type
              if $hsp->match;
            $desc .= "<tr><td>Length:<td>" . commify( $hsp->length ) . $seq_type
              if $hsp->length;
            $desc .= "<tr><td>Identity:<td>" . $hsp->percent_id . "%"
              if $hsp->percent_id;

            $desc .= "<tr><td>E_val:<td>" . $hsp->pval  if $hsp->pval;
            $desc .= "<tr><td>Score:<td>" . $hsp->score if $hsp->score;
            $desc .=
"<tr><td>Aligned Sequence:<td><span style=\"white-space: nowrap\">"
              . $seq
              . "</span>"
              if $seq;
            $seq =~ s/-//g if $seq;
            $desc .=
                "<tr><td>Sequence:<td><span style=\"white-space: nowrap\">"
              . $seq
              . "</span>"
              if $seq;
            $desc .=
                "<tr><td valign=top>Alignment: <td><pre>"
              . $hsp->alignment
              . "</pre>"
              unless $seq;
            $desc .= qq{<tr><td>cutoff:<td>$eval_cutoff}
              if defined $eval_cutoff;
            $desc .= "</table>";
            $f->description($desc);

            if ($reverse) {
                $strand = $strand =~ /-/ ? "1" : "-1";
                my $tmp = $seq_len - $start + 1;
                $start = $seq_len - $stop + 1;
                $stop  = $tmp;
            }
            $report =~ s/$TEMPDIR/GEvo/;
            my $link =
                "HSPView.pl?report=$report&num="
              . $hsp->number . "&db="
              . $cogeweb->basefilename
              . ".sqlite";
            $link .= join( "&",
                "&qstart=" . ( $gbobj->start + $start - 1 ),
                "qstop=" .   ( $gbobj->start + $stop - 1 ),
                "qchr=" . $gbobj->chromosome,
                "qds=" . $gbobj->dataset,
                "qstrand=" . $strand )
              if $gbobj->dataset;
            $f->link($link) if $link;
            $f->alignment( $hsp->alignment ) if $hsp->alignment;
            push @feats, $f;

            print STDERR $hsp->number, "-", $hsp->strand, $track, ":", $strand,
              "\n"
              if $DEBUG;

        }
        $i++;
        $track-- unless $hsp_track;
    }

    my $label_location = "top";
    my $order;
    @feats = sort {
             $a->strand cmp $b->strand
          || $a->track <=> $b->track
          || $a->start <=> $b->start
    } @feats;

    #    @feats = reverse @feats if $reverse;
    foreach my $f (@feats) {

#	next unless $f->label; #this has been commented out so that if features are labeled, but HSP are not, HSPs will still get drawn.  Not sure why this check was put into place. . .
        $order = $f->track unless $order;
        if ( $order ne $f->track ) {
            $order          = $f->track;
            $label_location = "top";
        }
        $f->label_location($label_location)
          if $hsp_labels && $hsp_labels eq "staggered";
        $f->font_size($font_size);
        $f->skip_overlap_search($skip_overlap_search);
        $c->add_feature($f);

        #	$c->_check_overlap($f);
        if ( !$label_location ) {
            $label_location = "bot";
        }
        elsif ( $label_location eq "top" ) {
            $label_location = "";
        }
        else {
            $label_location = "top";
        }
    }
    if ($hsp_overlap_limit) {
        foreach my $f (@feats) {
            if ( $f->_overlap >= $hsp_overlap_limit ) {
                $c->delete_feature($f);
            }
        }
    }
    if ($color_overlapped_features) {
        my %feat_count;
        map { $feat_count{ $_->type }++ } $gbobj->get_features();
        my $message = "\n" . $gbobj->accn . " contains:\n";
        foreach my $type ( sort keys %feat_count ) {
            $overlap_count{$type} = 0 unless $overlap_count{$type};
            my $percent_overlap = 100 *
              sprintf( "%.4f", $overlap_count{$type} / $feat_count{$type} );
            $message .= "\t"
              . $overlap_count{$type} . "/"
              . $feat_count{$type}
              . " $type ($percent_overlap%) with overlapping HSPs.\n";
        }
        CoGe::Accessory::Web::write_log( $message, $cogeweb->logfile )
          if $message;
    }
}

sub find_hsp_info #sub for retrieving hsp info when there is no query information
{
    my %opts  = @_;
    my $hsp   = $opts{hsp};
    my $start = $opts{start};    #genomic region start position
    my $stop  = $opts{stop};     #genomic region stop position
    my $rev   = $opts{rev};      #do we need to conver for a reverse complement?
    if ( $hsp->query_name =~ /fid:(\d+)/
      ) #we have a feature id in the query name -- HSP was generated by a feature sequence and not a genomic region!
    {
        my $fid  = $1;
        my $feat = $coge->resultset('Feature')->find($fid);
        my ( $new_start, $new_stop );
        if ($rev) {
            $new_start = ( $stop - $start + 1 ) - ( $feat->start - $start + 1 );
            $new_stop  = ( $stop - $start + 1 ) - ( $feat->stop - $start + 1 );
            print STDERR ( $new_stop, "\t", $new_start, "\n" );
            $hsp->strand( $hsp->strand * -1 )
              ; #need to switch strands of HSP is sequence is reverse complemented
        }
        else {
            $new_start = $feat->start - $start + 1;
            $new_stop  = $feat->stop - $start + 1;
        }
        $hsp->query_start($new_start);
        $hsp->query_stop($new_stop);
    }
}

sub generate_obj_from_seq {
    my %opts     = @_;
    my $sequence = $opts{seq};
    my $num      = $opts{num};
    my $rc       = $opts{rc};
    my $start    = $opts{start};                   #used in file name
    my $length   = $opts{length};                  #used in file name
    my $obj      = new CoGe::Accessory::GenBank;
    my $filename =
      $TEMPDIR . "/" . md5_hex( $sequence . $start . $length ) . ".faa";
    $obj->srcfile($filename);
    if ( $sequence =~ /^LOCUS/ ) {
        $obj->parse_genbank( content=>$sequence, rev=>$rc );
    }
    else {
        my $seq;
        if ( $sequence =~ /^>/ ) {

            #fasta sequence
            my ( $header, $seq ) = split( /\n/, $sequence, 2 );

            unless ($seq) {
                $seq = $header unless $seq;
                $header = 'NO_HEADER';
            }
            if ($header) {
                $header =~ s/>//g;
                $header =~ s/\|/_/g;
                $header =~ s/^\s+//;
                $header =~ s/\s+$//;
            }
            if ($seq) {
                $seq =~ s/\n//g;
                $seq =~ s/\s//g;
            }
            return unless $seq;
            my ($accn) = $header =~ /^(\S*)/;
            $obj->accn($header);
            $obj->locus($header);
            $obj->definition($header);
            $seq =~ s/\n|\r//g;
            $sequence = $seq;
        }
        else {

            #just the sequence
            $obj->accn("RAW_SEQUENCE_SUBMISSION_$num");
            $obj->locus("RAW_SEQUENCE_SUBMISSION_$num");
        }
        $sequence =~ s/\n|\r|\s//g;
        $obj->sequence($sequence);
    }
    if ($rc) {
        $obj->sequence(
            CoGeX::Result::Feature->reverse_complement( $obj->sequence ) );
    }
    return $obj;
}

sub get_obj_from_genome_db {
    my %opts   = @_;
    my $accn   = $opts{accn};
    my $featid = $opts{featid};
    my $pos    = $opts{pos};
    my $dsid   = $opts{dsid};
    my $rev    = $opts{rev};
    my $chr    = $opts{chr};
    my $up     = $opts{up} || 0;
    my $down   = $opts{down} || 0;
    my $gstid  = $opts{gstid};       #genomic sequence type database id
    my $mask   = $opts{mask};        #need this to check for seq file
    my $dsgid  = $opts{dsgid};       #dataset group id
                                     #print STDERR Dumper \%opts;
    my $color_anchor = $opts{color_anchor}; #option to color anchor gene yellow: EHL 5/5/15
    $color_anchor = 1 unless defined $color_anchor;
    
    my $gen_prot_sequence = $opts{gen_prot_sequence}
      || 0;    #are we generating a protein sequence file too?
    my $message;
    my $dsg;
    my $ds;

    if ($dsgid) {
        $dsg = $coge->resultset('Genome')->find($dsgid);
        return "", "Can't find entry for $dsgid" unless $dsg;
        return "", "Permission denied"
          if $dsg->restricted && !$USER->has_access_to_genome($dsg);

#($ds) = $dsid ? $coge->resultset('Dataset')->find($dsid) : $dsg->datasets( chr => $chr );
        ($ds) = $dsg->datasets( chr => $chr );

#print STDERR join ("\t", map {$_->id} $dsg->datasets(chromosome=>$chr) ),"::$chr\n";
        return ( "", "Chromosome '$chr' not found for this genome" ) unless $ds;
        $dsid  = $ds->id;
        $gstid = $dsg->type->id;
    }

    #let's get a unique file name for this sequence
    my $seq_file = create_seq_file_name(%opts);
    my $t1       = new Benchmark;
    my ($feat)   = $coge->resultset('Feature')->find($featid);
    my $new_seq  = 0
      ; #flag for whether we got the seq from the database or a previously generated file;
    unless ( ref($feat) =~ /feature/i || ( $pos && $dsid ) ) {
        CoGe::Accessory::Web::write_log(
            "Can't find valid feature database entry for id=$featid",
            $cogeweb->logfile );
        return ( "", "Can;t find valid feature database entry for id $featid" );
    }
    my $t2    = new Benchmark;
    my $start = 0;
    my $stop  = 0;
    my $seq;

    if ($feat) {
        $dsid = $feat->dataset->id unless $dsid;
        if ($rev) {
            $start = $feat->start - $down;
            $stop  = $feat->stop + $up;
        }
        else {
            $start = $feat->start - $up;
            $stop  = $feat->stop + $down;
        }
        $start = 1 if $start < 1;
        $chr = $feat->chr;
    }
    elsif ($pos) {
        if ($rev) {
            $start = $pos - $down;
            $stop  = $pos + $up;
        }
        else {
            $start = $pos - $up;
            $stop  = $pos + $down;
        }
        $start = 1 if $start < 1;
    }

    $ds = $coge->resultset('Dataset')->find($dsid) if $dsid && !$ds;
    return ( "", "unable able to find dataset for id $dsid" ) unless $ds;
    $accn = $dsg->name if !$accn && $dsg;
    $accn = $ds->name  if !$accn && $ds;
    if ( -r $seq_file ) {
        $/ = "\n";
        open( IN, $seq_file );
        while (<IN>) {
            chomp;
            next if /^>/;    #skip header;
            $seq .= $_;
        }
        close IN;

    }
    unless ($seq && length($seq) > 1) {
        ($chr) = $ds->get_chromosomes unless defined $chr;
        my $tmp;
        ( $tmp, $seq_file ) = CoGe::Accessory::Web::check_taint($seq_file);
        unlink($seq_file);
        $seq = $ds->get_genomic_sequence(
            start => $start,
            stop  => $stop,
            chr   => $chr,
            gstid => $gstid
        );
        $new_seq = 1;
        if ( !$seq ) {
            print STDERR "Error retrieving sequence: "
              . join( "\t", $ds->name, $start, $stop, $chr, $gstid ), "\n";
        }
        $seq = CoGeX::Result::Feature->reverse_complement($seq) if $rev;
    }
    if ( $stop - $start + 1 > length($seq) ) {
        my $len = length($seq);
        $stop = $start + $len - 1;
    }
    my $t3       = new Benchmark;
    my $gst      = $coge->resultset('GenomicSequenceType')->find($gstid);
    my $org_name = $ds->organism->name();
    $org_name .= ": " . $dsg->description if $dsg && $dsg->description;
    $org_name .= " ("
      . $ds->data_source->name . " v"
      . $ds->version . ", "
      . $gst->name . ")";
    my $obj = new CoGe::Accessory::GenBank(
        {
            accn                     => $accn,
            locus                    => $accn,
            version                  => $ds->version(),
            data_source              => $ds->data_source->name(),
            dataset                  => $dsid,
            chromosome               => $chr,
            start                    => $start,
            stop                     => $stop,
            organism                 => $org_name,
            seq_length               => length($seq),
            genomic_sequence_type_id => $gstid,
            srcfile                  => $seq_file,
        }
    );

    my %used_names;
    $used_names{$accn} = 1;
    my $t4 = new Benchmark;

    #print STDERR "Region: $chr: $start-$stop\n";# if $DEBUG;
    my %feats = map { $_->id, $_ } $coge->get_features_in_region(
        start      => $start,
        stop       => $stop,
        chr        => $chr,
        dataset_id => $dsid,
        gid        => $dsgid
    );

    #print STDERR join ("\t", $start, $stop, $chr, $dsid),"\n";
    $feats{ $feat->id } = $feat if $feat;

    my $t5 = new Benchmark;
    my $prot_sequence;
    foreach my $f ( values %feats ) {
        my $name;
        my @names = $f->names;
        next if $f->type->name =~ /misc_feat/;
        next if $f->type->name eq "chromosome";
	#new changes to deal with GD failures:  EHL 4/30/15
	next if $f->start == $f->stop;

        #	next if $f->type->name eq "contig";
        foreach my $tmp (@names) {
            $name = $tmp;
            last if ( $tmp =~ /^$accn.+/i );
        }
        unless (@names) {

            #	    next;
            push @names, "No Name Feature: " . $f->type->name, "\n";

#	    print STDERR "Feature has no Name.  Db Id: ",$f->id,"\n" unless $f->type->name =~ /misc/;

        }
        $name = $accn unless $name;
	    print STDERR "\n" if $DEBUG;
        print STDERR $name, " ID: ",$f->id,"\n" if $DEBUG;
        print STDERR "\t", $f->genbank_location_string(), "\n" if $DEBUG;
        print STDERR "\t", $f->genbank_location_string( recalibrate => $start ),
          "\n"
          if $DEBUG;
        my $anno =
          $P->{SERVER} . "FeatAnno.pl?fid=" . $f->id . ";gstid=" . $gstid;
        my $location = $f->genbank_location_string( recalibrate => $start );
        $location = $obj->reverse_genbank_location( loc => $location, ) if $rev;
        print STDERR $name, "\t", $f->type->name, "\t", $location, "\n"
          if $DEBUG;
        my $type = $f->type->name;
        $type = "anchor" if $f->id && $featid && $f->id == $featid;# WORKING ON THIS FEATURE && $color_anchor;
        $obj->add_feature(
            type       => $f->type->name,
            location   => $location,
            qualifiers => {
                names => [@names],
                type  => $type,
                id    => $f->id,
            },
            annotation => $anno,
            force      => 1
            , #force loading of anchor in case it is outside of specified sequence range
        );

        if ( $gen_prot_sequence && $f->type->name eq "CDS" ) {
            $prot_sequence .= $f->fasta( prot => 1, col => 0, add_fid => 1 );
        }
    }
    if ($pos) {
        my $pos_loc = $rev ? $stop - $pos : $pos - $start;
        $obj->add_feature(
            type       => "anchor",
            location   => ($pos_loc),
            annotation => "User specified anchor point",
            qualifiers => {
                names => ['position_anchor'],
                type  => 'anchor',
            },
            force => 1
            , #force loading of anchor in case it is outside of specified sequence range
        );
    }
    if ($gen_prot_sequence) {
        my $prot_file = $seq_file . ".prot";
        my $tmp;
        ( $tmp, $prot_file ) = CoGe::Accessory::Web::check_taint($prot_file);
        open( OUT, ">$prot_file" );
        print OUT $prot_sequence if $prot_sequence;
        close OUT;
        $obj->other_stuff($prot_file);
    }
    if ($new_seq) {

        # mask out cds if required
        $seq = $obj->mask_cds($seq) if ( $mask && $mask eq "cds" );

        # mask out rna if required
        $seq = $obj->mask_rna($seq) if ( $mask && $mask eq "rna" );

        # mask out non coding sequences if required
        $seq = $obj->mask_ncs($seq) if ( $mask && $mask eq "non-cds" );

        # mask out non genic sequences if required
        $seq = $obj->mask_ngene($seq) if ( $mask && $mask eq "non-genic" );
    }
    $obj->sequence($seq);

    my $t6               = new Benchmark;
    my $db_time          = timestr( timediff( $t2, $t1 ) );
    my $seq_time         = timestr( timediff( $t3, $t2 ) );
    my $int_obj_time     = timestr( timediff( $t4, $t3 ) );
    my $feat_region_time = timestr( timediff( $t5, $t4 ) );
    my $pop_obj_time     = timestr( timediff( $t6, $t5 ) );
    print STDERR qq{
Getting feature from DB took:                         $db_time
Getting sequence for region took:                     $seq_time
Initialize obj took:                                  $int_obj_time
Getting feats in region took:                         $feat_region_time
Populating object took:                                $pop_obj_time
Region:         dsid: $dsid $start-$stop($chr)
} if $BENCHMARK && $DEBUG;
    return $obj;
}

sub run_bl2seq {
    my %opts         = @_;
    my $sets         = $opts{sets};
    my $blast_params = $opts{params};
    my $program      = $opts{blast_program};
    my $parser_opts  = $opts{parser_opts};
    my $comp_adj     = $opts{comp_adj};
    my $eval_cutoff  = $opts{eval_cutoff};

    $program = "blastn" unless $program;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count      = 0;
    my $pm         = new Parallel::ForkManager($MAX_PROC);
    for ( my $i = 0 ; $i < scalar @$sets ; $i++ ) {
        for ( my $j = 0 ; $j < scalar @$sets ; $j++ ) {
            next if $comp_adj && !( $i + 1 == $j );
            next unless $j > $i;

            #
            my $seqfile1 = $sets->[$i]->{file};
            my $seqfile2 = $sets->[$j]->{file};
            next
              unless $seqfile1
                  && -r $seqfile1
                  && $seqfile2
                  && -r $seqfile2;    #make sure these files exist

            next
              unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
            $count++;
            my ( $accn1, $accn2 ) = ( $sets->[$i]{accn}, $sets->[$j]{accn} );
            my ($tempfile) =
                $cogeweb->basefile . "_"
              . ( $sets->[$i]->{seq_num} ) . "-"
              . ( $sets->[$j]->{seq_num} )
              . ".bl2seq";
            push @reports, [ $tempfile, $accn1, $accn2 ];
            $pm->start and next;
            my $command = $BL2SEQ;

            # format the bl2seq command
            $command .= " -p $program -o '$tempfile'";
            $command .= " -i '$seqfile1' -j '$seqfile2'";
            $command .= " " . $blast_params;
            my $x = "";
            ( $x, $command ) = CoGe::Accessory::Web::check_taint($command);
            if ($DEBUG) {
                print STDERR "About to execute...\n $command\n";
            }
            unless ($x) {
                next;
            }

            # execute the command
            CoGe::Accessory::Web::write_log(
                "running ($count/$total_runs) " . $command,
                $cogeweb->logfile );
            `$command`;
            system "/bin/chmod +rw '$tempfile'";
            $pm->finish;
        }
    }
    $pm->wait_all_children;
    foreach my $item (@reports) {
        my $tempfile    = $item->[0];
        my $parser_opts = $opts{parser_opts};
        my $blastreport =
          new CoGe::Accessory::bl2seq_report( { file => $tempfile } )
          if -r $tempfile;
        if ($blastreport) {
            push @$item, $blastreport;
        }
        else {
            my $blastreport = new CoGe::Accessory::bl2seq_report();
            push @$item, $blastreport;
            CoGe::Accessory::Web::write_log(
                "WARNING:  no results from blasting "
                  . $item->[1] . " and "
                  . $item->[2]
                  . ".  Possible error\n",
                $cogeweb->logfile
            );
        }
    }
    return ( \@reports );
}

sub run_blastz {
    my %opts        = @_;
    my $sets        = $opts{sets};
    my $params      = $opts{params};
    my $parser_opts = $opts{parser_opts};
    my $comp_adj    = $opts{comp_adj};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count      = 0;
    my $pm         = new Parallel::ForkManager($MAX_PROC);

    for ( my $i = 0 ; $i < scalar @$sets ; $i++ ) {
        for ( my $j = 0 ; $j < scalar @$sets ; $j++ ) {
            next if $comp_adj && !( $i + 1 == $j );
            next unless $j > $i;
            my $seqfile1 = $sets->[$i]->{file};
            my $seqfile2 = $sets->[$j]->{file};
            next
              unless $seqfile1
                  && $seqfile2
                  && -r $seqfile1
                  && -r $seqfile2;    #make sure these files exist
            next
              unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
            $count++;
            my ( $accn1, $accn2 ) = ( $sets->[$i]{accn}, $sets->[$j]{accn} );
            my ($tempfile) =
                $cogeweb->basefile . "_"
              . ( $sets->[$i]->{seq_num} ) . "-"
              . ( $sets->[$j]->{seq_num} )
              . ".blastz";
            push @reports, [ $tempfile, $accn1, $accn2 ];
            $pm->start and next;
            my $command = $BLASTZ;
            $command .= " '$seqfile1'" . "[unmask] '$seqfile2'" . "[unmask]";
            $command .= " " . $params if $params;
            my $x = "";
            ( $x, $command ) = CoGe::Accessory::Web::check_taint($command);

            if ($DEBUG) {
                print STDERR "About to execute...\n $command\n";
            }
            unless ($x) {
                next;
            }
            $command .= " > '$tempfile'";

            # execute the command
            CoGe::Accessory::Web::write_log(
                "running ($count/$total_runs) " . $command,
                $cogeweb->logfile );
            `$command`;
            system "/bin/chmod +rw '$tempfile'";
            $pm->finish;
        }
    }
    $pm->wait_all_children;
    foreach my $item (@reports) {
        my $tempfile = $item->[0];
        my $blastreport =
          new CoGe::Accessory::blastz_report( { file => $tempfile } )
          if -r $tempfile;
        if ($blastreport) {
            push @$item, $blastreport;
        }
        else {
            push @$item,
              "no results from blasting " . $item->[1] . " and " . $item->[2];
        }
    }
    return \@reports;
}

sub run_lagan {
    my %opts        = @_;
    my $sets        = $opts{sets};
    my $params      = $opts{params};
    my $parser_opts = $opts{parser_opts};
    my $comp_adj    = $opts{comp_adj};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count      = 1;

    for ( my $i = 0 ; $i < scalar @$sets ; $i++ ) {
        for ( my $j = 0 ; $j < scalar @$sets ; $j++ ) {
            next if $comp_adj && !( $i + 1 == $j );
            next unless $j > $i;
            my $seqfile1 = $sets->[$i]->{file};
            my $seqfile2 = $sets->[$j]->{file};
            next
              unless -r $seqfile1 && -r $seqfile2;  #make sure these files exist
            next
              unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
            my ( $accn1, $accn2 ) = ( $sets->[$i]{accn}, $sets->[$j]{accn} );
            my ($tempfile) =
                $cogeweb->basefile . "_"
              . ( $sets->[$i]->{seq_num} ) . "-"
              . ( $sets->[$j]->{seq_num} )
              . ".lagan";
            my $command = $LAGAN;
            $command .= " $seqfile1 $seqfile2";
            $command .= " -mfa";
            $command .= " " . $params if $params;
            my $x = "";
            ( $x, $command ) = CoGe::Accessory::Web::check_taint($command);

            if ($DEBUG) {
                print STDERR "About to execute...\n $command\n";
            }
            unless ($x) {
                next;
            }
            $command .= " > '$tempfile'";
            CoGe::Accessory::Web::write_log(
                "running ($count/$total_runs) " . $command,
                $cogeweb->logfile );

            #time for execution
            #	    print STDERR "running $command\n";
            `$command`;
            system "/bin/chmod +rw '$tempfile'";
            my $report = new CoGe::Accessory::lagan_report(
                { file => $tempfile, %$parser_opts } )
              if -r $tempfile;
            my @tmp = ( $tempfile, $accn1, $accn2 );
            if ($report) {
                push @tmp, $report;
            }
            else {
                push @tmp,
                  "no results from comparing $accn1 and $accn2 with LAGAN";
            }
            push @reports, \@tmp;
            $count++;
        }
    }
    return \@reports;
}

sub run_chaos {
    my %opts        = @_;
    my $sets        = $opts{sets};
    my $params      = $opts{params};
    my $parser_opts = $opts{parser_opts};
    my $comp_adj    = $opts{comp_adj};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count      = 1;

    for ( my $i = 0 ; $i < scalar @$sets ; $i++ ) {
        for ( my $j = 0 ; $j < scalar @$sets ; $j++ ) {
            next if $comp_adj && !( $i + 1 == $j );
            next unless $j > $i;
            my $seqfile1 = $sets->[$i]->{file};
            my $seqfile2 = $sets->[$j]->{file};
            next
              unless -r $seqfile1 && -r $seqfile2;  #make sure these files exist
            next
              unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
            my ( $accn1, $accn2 ) = ( $sets->[$i]{accn}, $sets->[$j]{accn} );
            my ($tempfile) =
                $cogeweb->basefile . "_"
              . ( $sets->[$i]->{seq_num} ) . "-"
              . ( $sets->[$j]->{seq_num} )
              . ".chaos";
            my $command = $CHAOS;
            $command .= " $seqfile1 $seqfile2";

            $command .= " " . $params if $params;
            my $x = "";
            ( $x, $command ) = CoGe::Accessory::Web::check_taint($command);
            if ($DEBUG) {
                print STDERR "About to execute...\n $command\n";
            }
            unless ($x) {
                next;
            }
            $command .= " > '$tempfile'";
            CoGe::Accessory::Web::write_log(
                "running ($count/$total_runs) " . $command,
                $cogeweb->logfile );

            #time for execution
            `$command`;
            system "/bin/chmod +rw '$tempfile'";
            my $report = new CoGe::Accessory::chaos_report(
                { file => $tempfile, %$parser_opts } )
              if -r $tempfile;
            my @tmp = ( $tempfile, $accn1, $accn2 );
            if ($report) {
                push @tmp, $report;
            }
            else {
                push @tmp,
                  "no results from comparing $accn1 and $accn2 with Chaos";
            }
            push @reports, \@tmp;
            $count++;
        }
    }
    return \@reports;
}

sub run_dialign {
    my %opts        = @_;
    my $sets        = $opts{sets};
    my $params      = $opts{params};
    my $parser_opts = $opts{parser_opts};
    my $comp_adj    = $opts{comp_adj};
    my $max_length  = $opts{max_length} || 10000;
    my $kill_length = 2 * $max_length;
    my $program_ran = "DIALIGN";
    my $error_message;

    #my $algo_limits = $opts{algo_limits};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count      = 1;
    for ( my $i = 0 ; $i < scalar @$sets ; $i++ ) {
        for ( my $j = 0 ; $j < scalar @$sets ; $j++ ) {
            next if $comp_adj && !( $i + 1 == $j );
            next unless $j > $i;
            my $obj1        = $sets->[$i]->{obj};
            my $obj2        = $sets->[$j]->{obj};
            my $seq_length1 = $obj1->seq_length;
            my $seq_length2 = $obj2->seq_length;
            my $max_seq =
              $seq_length1 > $seq_length2 ? $seq_length1 : $seq_length2;
            my $kb_max_seq = $max_seq / 1000;
            $kb_max_seq =~ s/(\.\d+)$//;

            if ( $max_seq > $max_length ) {
                if ( $max_seq > $kill_length ) {
                    $error_message =
"The sequence you have chosen DIALIGN to align is too long (your longest sequence is $kb_max_seq kb long). To complete the alignment would take a significant amount of time, and use up considerable system resources that others need to run their alignments. Please limit your sequences to no more than "
                      . ( $kill_length / 1000 )
                      . " kb long. Thank you.";
                    return ( '', '', $error_message );
                }
                else {
                    my $changed = $params =~ s/((\s?-cs)|(\s?-nt)|(\s?-ma))//g;
                    $params .= " -n"       unless $params =~ /(-n)/;
                    $params .= " -ds"      unless $params =~ /(-ds)/;
                    $params .= " -thr 2"   unless $params =~ /(-thr)/;
                    $params .= " -lmax 30" unless $params =~ /(-lmax)/;
                    $error_message =
"Your search parameters were altered due to the length of the sequences you requested for alignment. Your alignment was done with the \"Large Sequence Alignment\" option enabled";
                    $error_message .=
                      $changed
                      ? ", and with \"Short Sequence Alignment with Translation\" disabled."
                      : ".";
                    $error_message .=
                        " This is done automatically with sequences over "
                      . ( $max_length / 1000 )
                      . " kb long.\n";

                    #print $error_message;
                }
            }
            my $seqfile1 = $sets->[$i]->{file};
            my $seqfile2 = $sets->[$j]->{file};
            my ( $accn1, $accn2 ) = ( $sets->[$i]{accn}, $sets->[$j]{accn} );
            my $report;    #place to hold report
            my $tempfile;
            if ( $parser_opts->{anchor_params} ) {
                my $anchor_prog;
                if ( $parser_opts->{anchor_params}{prog} =~ /chaos/i ) {
                    $anchor_prog = $CHAOS;
                }
                elsif ( $parser_opts->{anchor_params}{prog} =~ /bl2/i ) {
                    $anchor_prog = $BL2SEQ;
                }
                else {
                    $anchor_prog = $BLASTZ;
                }
                my ( $x, $dialign_opts ) =
                  CoGe::Accessory::Web::check_taint($params);
                $program_ran .=
                  " using anchors from " . $parser_opts->{anchor_params}{prog};
                unless ($x) {
                    next;
                }
                my ( $y, $anchor_opts ) = CoGe::Accessory::Web::check_taint(
                    $parser_opts->{anchor_params}{param_string} );
                unless ($y) {
                    next;
                }
                my $obj = new CoGe::Accessory::dialign_report::anchors(
                    {
                        file1            => $seqfile1,
                        file2            => $seqfile2,
                        base_name        => $cogeweb->basefilename,
                        extension        => $parser_opts->{anchor_params}{prog},
                        run_anchor       => $anchor_prog,
                        run_dialign_opts => $dialign_opts,
                        run_anchor_opts  => $anchor_opts,
                        dialign_report_opts => $parser_opts,
                        anchor_report_opts =>
                          $parser_opts->{anchor_params}{parser_opts},
                        log_file => $cogeweb->logfile,

                        #									DEBUG=>1,
                    }
                );
                $tempfile = $obj->dialign_file;
                $report   = $obj->dialign_report;
            }
            else {
                my $seqfile =
                    $cogeweb->basefile . "_"
                  . ( $sets->[$i]->{seq_num} ) . "-"
                  . ( $sets->[$j]->{seq_num} )
                  . ".fasta";
                next
                  unless -r $seqfile1
                      && -r $seqfile2;    #make sure these files exist
                    #put two fasta files into one for dialign
                my $tmp;
                ( $tmp, $seqfile1 ) =
                  CoGe::Accessory::Web::check_taint($seqfile1);
                ( $tmp, $seqfile2 ) =
                  CoGe::Accessory::Web::check_taint($seqfile2);
                `cat '$seqfile1' > '$seqfile'`;
                `cat '$seqfile2' >> '$seqfile'`;
                next
                  unless $sets->[$i]{reference_seq}
                      || $sets->[$j]{reference_seq};
                ($tempfile) =
                    $cogeweb->basefile . "_"
                  . ( $sets->[$i]->{seq_num} ) . "-"
                  . ( $sets->[$j]->{seq_num} )
                  . ".dialign";
                my $command = $DIALIGN;
                $command .= " " . $params if $params;
                $command .= " -fn '$tempfile'";
                $command .= " '$seqfile'";
                my $x = "";
                ( $x, $command ) = CoGe::Accessory::Web::check_taint($command);

                if ($DEBUG) {
                    print STDERR "About to execute...\n $command\n";
                }
                unless ($x) {
                    next;
                }

                #time for execution
                CoGe::Accessory::Web::write_log(
                    "running ($count/$total_runs) " . $command,
                    $cogeweb->logfile );
                `$command`;
                system "/bin/chmod +rw '$tempfile'";
                $report = new CoGe::Accessory::dialign_report(
                    { file => $tempfile, %$parser_opts } )
                  if -r $tempfile;
            }
            my @tmp = ( $tempfile, $accn1, $accn2 );
            if ($report) {
                push @tmp, $report;
            }
            else {
                push @tmp,
                  "no results from comparing $accn1 and $accn2 with LAGAN";
            }
            push @reports, \@tmp;
            $count++;
        }
    }
    return ( \@reports, $program_ran, $error_message );
}

sub run_genomethreader {
    my %opts        = @_;
    my $sets        = $opts{sets};
    my $params      = $opts{params};
    my $parser_opts = $opts{parser_opts};
    my $comp_adj    = $opts{comp_adj};
    my $matrix      = $opts{parser_opts}{matrix}
      || $P->{BLASTMATRIX} . "aa/BLOSUM62";
    my @reports;
    my $total_runs =
      number_of_runs($sets) * 2
      ; #need to run this in both directions for each genomic and protein sequence
    my $count = 1;
    my $pm    = new Parallel::ForkManager($MAX_PROC);

    for ( my $i = 0 ; $i < scalar @$sets ; $i++ ) {
        for ( my $j = 0 ; $j < scalar @$sets ; $j++ )
        { #need to run this in both directions using each as the genome and protein sequence
            next if $comp_adj && !( $i + 1 == $j );
            next unless $j > $i;
            my $seqfile1  = $sets->[$i]->{file};
            my $protfile1 = $sets->[$i]->{obj}->other_stuff;
            my $seqfile2  = $sets->[$j]->{file};
            my $protfile2 = $sets->[$j]->{obj}->other_stuff;
            next
              unless -r $seqfile1 && -r $seqfile2;  #make sure these files exist
            next
              unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
            my ( $accn1, $accn2 ) = ( $sets->[$i]{accn}, $sets->[$j]{accn} );

#	    foreach my $item ({genomic=>$seqfile1, protein=>$protfile2},{genomic=>$seqfile2, protein=>$protfile1})
            foreach my $item ( { genomic => $seqfile2, protein => $protfile1 } )
            {
                my ($tempfile) =
                    $cogeweb->basefile . "_"
                  . ( $sets->[$i]->{seq_num} ) . "-"
                  . ( $sets->[$j]->{seq_num} )
                  . ".genomethreader$count";
                my $command = $GENOMETHREADER;
                $command .= " -genomic " . $item->{genomic};
                $command .= " -protein " . $item->{protein};
                $command .= " -scorematrix $matrix";
                $command .= " " . $params if $params;
                my $x = "";
                ( $x, $command ) = CoGe::Accessory::Web::check_taint($command);

                if ($DEBUG) {
                    print STDERR "About to execute...\n $command\n";
                }
                unless ($x) {
                    next;
                }
                $command .= " > '$tempfile'";
                CoGe::Accessory::Web::write_log(
                    "running ($count/$total_runs) " . $command,
                    $cogeweb->logfile );

                #time for execution
                `$command`;
                system "/bin/chmod +rw '$tempfile'";
                my $report = new CoGe::Accessory::GenomeThreader_report(
                    { file => $tempfile, %$parser_opts } )
                  if -r $tempfile;

#queries are protein sequences. need to convert their positions to relative genomic coordinates
#is the query a protein sequence?  if so, the query name should contain the fid tag
                foreach my $hsp ( @{ $report->hsps } ) {
                    if ( $hsp->query_name =~ /fid/ ) {
                        find_hsp_info(
                            hsp   => $hsp,
                            start => $sets->[$i]->{obj}->start,
                            stop  => $sets->[$i]->{obj}->stop,
                            rev   => $sets->[$i]->{rev}
                        );
                    }
                }

                my @tmp = ( $tempfile, $accn1, $accn2 );
                if ($report) {
                    push @tmp, $report;
                }
                else {
                    push @tmp,
                      "no results from comparing $accn1 and $accn2 with LAGAN";
                }
                push @reports, \@tmp;
                $count++;
            }
        }

        #	$pm->wait_all_children;
    }
    return \@reports;
}

sub get_substr {
    my %opts  = @_;
    my $seq   = $opts{seq};
    my $start = $opts{start};
    $start =~ s/,|\.//g;
    $start = 1 if !$start || $start < 1;
    my $stop = $opts{stop};
    return $seq if ( $start == 1 && !$stop );
    my $seqlength = length $seq;
    $stop = $seqlength unless $stop;
    $stop =~ s/,|\.//g;
    my $newlength = $stop - $start;

    if ( $start < $seqlength && $newlength ) {
        $seq = substr( $seq, $start - 1, $newlength );
    }
    return $seq;
}

sub write_fasta {
    my %opts     = @_;
    my $gbobj    = $opts{obj};
    my $start    = $opts{start} || $opts{startpos} || 1;
    my $mask     = $opts{mask};
    my $length   = $opts{length};
    my $force    = $opts{force};
    my $fullname = $gbobj->srcfile;
    my ($seq)    = uc( $gbobj->sequence() );
# EHL 5/6/19:  Removing the check for the cached sequence file due to problems generating bad sequences.
#    if ( -r $fullname && !$force ) {
#        CoGe::Accessory::Web::write_log(
#            "sequence file $fullname already exists.",
#            $cogeweb->logfile );
#        return $fullname;
#    }
    my $hdr = $gbobj->get_headerfasta();
    $start = 1 if $start < 1;
    my $stop = $length ? $start + $length - 1 : length($seq);

    # mask out cds if required
    $seq = $gbobj->mask_cds($seq) if ( $mask && $mask eq "cds" );

    # mask out rna if required
    $seq = $gbobj->mask_rna($seq) if ( $mask && $mask eq "rna" );

    # mask out non coding sequences if required
    $seq = $gbobj->mask_ncs($seq) if ( $mask && $mask eq "non-cds" );

    # mask out non genic sequences if required
    $seq = $gbobj->mask_ngene($seq) if ( $mask && $mask eq "non-genic" );
    $seq = substr( $seq, $start - 1, $stop - $start + 1 );
    $gbobj->sequence($seq);    #replace objects sequence with modified sequence
    ($fullname) = CoGe::Accessory::Web::check_filename_taint($fullname);
    open( OUT, ">$fullname" ) or die "Couldn't open '$fullname': $!";
    print OUT "$hdr\n";
    print OUT $seq, "\n";
    close(OUT);
    system "/bin/chmod +rw '$fullname'";
    CoGe::Accessory::Web::write_log( "Created sequence file $fullname",
        $cogeweb->logfile );
    return ($fullname);
}

sub get_bit_score_cutoff {
    my %opts       = @_;
    my $seq_length = $opts{seq_length};
    my $match      = $opts{match};
    my $mismatch   = $opts{mismatch};
    my $feat       = new CoGeX::Result::Feature;
    my $bs         = $feat->blast_bit_score(
        match      => $match,
        mismatch   => $mismatch,
        seq_length => $seq_length
    );
    return $bs;
}

sub generate_annotation {
    my %opts = @_;
    my $obj  = $opts{obj};
    return unless $obj;
    my $rev      = $opts{rev};
    my $seq_num  = $opts{seq_num};
    my $fullname = $cogeweb->basefile . "_" . $seq_num . ".anno";
    my %data;
    my $length = length( $obj->sequence() );

    foreach my $feat ( $obj->get_features() ) {
        my $type = $feat->type;
        my ($name) =
          sort { length($b) <=> length($a) || $a cmp $b }
          @{ $feat->qualifiers->{names} }
          if ref( $feat->qualifiers ) =~ /hash/i && $feat->qualifiers->{names};
        my $dir = $feat->location =~ /complement/ ? "<" : ">";
        foreach my $block ( @{ $feat->blocks } ) {
            $block->[0] = 1 unless $block->[0];
            my $start = $block->[0];
            my $stop  = $block->[1];
            if ($rev) {
                $start = $length - $block->[1];
                $stop  = $length - $block->[0];
                $dir   = $dir eq ">" ? "<" : ">";
            }
            next unless $name;
            push @{ $data{$name}{$type} }, [ $start, $stop, $dir ];
        }
    }
    return unless keys %data;
    open( OUT, ">$fullname" ) || die "Can't open $fullname for writing! $!\n";
    foreach my $name ( keys %data ) {
        my $gene = $data{$name}{gene};
        if ($gene) {
            delete $data{$name}{gene};
            print OUT
              join( " ", $gene->[0][2], $gene->[0][0], $gene->[0][1], $name ),
              "\n";
            foreach my $type ( keys %{ $data{$name} } ) {
                my $outtype;
                $outtype = "utr"  if $type =~ /rna/i;
                $outtype = "exon" if $type =~ /cds/i;
                next unless $outtype;
                foreach my $item ( @{ $data{$name}{$type} } ) {
                    next if $item->[0] == $item->[1];
                    print OUT join( " ", $item->[0], $item->[1], $outtype ),
                      "\n";
                }
            }
        }
    }
    close OUT;
    return $fullname;
}

sub gen_params {
    my $num_seqs = shift || $NUM_SEQS;

    my $params;
    for ( my $i = 1 ; $i <= $num_seqs ; $i++ ) {
        $params .= qq{'args__display_order$i', 'display_order$i',};
        $params .=
qq{'args__draccn$i', 'args__'+encodeURIComponent(document.getElementById('accn$i').value),}
          ;    #qq{'args__draccn$i', 'accn$i',}; # mdb changed 10/9/12 issue #11
        $params .= qq{'args__featid$i', 'featid$i',};
        $params .= qq{'args__drup$i', 'drup$i',};
        $params .= qq{'args__drdown$i', 'drdown$i',};
        $params .= qq{'args__pos$i', 'pos$i',};
        $params .= qq{'args__dsid$i', 'dsid$i',};
        $params .= qq{'args__dsgid$i', 'dsgid$i',};

        #	$params .= qq{'args__posdsid$i', 'posdsid$i',};
        $params .= qq{'args__chr$i', 'chr$i',};
        $params .= qq{'args__gbaccn$i', 'gbaccn$i',};
        $params .= qq{'args__gbstart$i', 'gbstart$i',};
        $params .= qq{'args__gblength$i', 'gblength$i',};
        $params .= qq{'args__dirseq$i', 'dirseq$i',};
        $params .= qq{'args__dirstart$i', 'dirstart$i',};
        $params .= qq{'args__dirlength$i', 'dirlength$i',};
        $params .= qq{'args__ref_seq$i', 'ref_seq$i',};
        $params .= qq{'args__skip_seq$i', 'skip_seq$i',};
        $params .= qq{'args__rev$i', 'rev$i',};
        $params .= qq{'args__mask$i', 'mask$i',};
    }
    for ( my $i = 1 ; $i <= num_colors($num_seqs) ; $i++ ) {
        $params .=
qq{'args__rgb$i', 'args__'+\$('#sample_color$i').css('backgroundColor'),};

        #	$params .= qq{'args__g$i', 'g$i',};
        #	$params .= qq{'args__b$i', 'b$i',};
    }

    $params .= qq{
        'args__pad_gs', 'pad_gs',

	'args__spike', 'spike',
	'args__blast_word', 'blast_wordsize',
	'args__blast_gapo', 'blast_gapopen',
	'args__blast_gape', 'blast_gapextend',
	'args__blast_mmatch', 'blast_mismatch',
	'args__blast_match', 'blast_match',
	'args__blast_eval', 'blast_eval',
	'args__blast_params', 'blastparams',
        'args__blast_filter','blast_filter',

        'args__lagan_min_length','lagan_min_length',
        'args__lagan_max_gap','lagan_max_gap',
        'args__lagan_percent_id','lagan_percent_id',
        'args__lagan_params','lagan_params',

        'args__chaos_word_length','chaos_word_length',
        'args__chaos_score_cutoff','chaos_score_cutoff',
        'args__chaos_rescore_cutoff','chaos_rescore_cutoff',
        'args__chaos_lookback','chaos_lookback',
        'args__chaos_gap_length','chaos_gap_length',
        'args__chaos_gap_start','chaos_gap_start',
        'args__chaos_gap_extension','chaos_gap_extension',
        'args__chaos_params','chaos_params',

        'args__dialign_min_score','dialign_min_score',
        'args__dialign_max_gap','dialign_max_gap',
        'args__dialign_split_score','dialign_split_score',
        'args__dialign_params','dialign_params',
        'args__dialign_motif_motif','dialign_motif_motif',
        'args__dialign_motif_weight','dialign_motif_weight',
        'args__dialign_motif_penalty','dialign_motif_penalty',
        'args__dialign_use_anchor','dialign_use_anchor',
        'args__dialign_anchor_program','dialign_anchor_program',

        'args__blastz_wordsize','blastz_wordsize',
        'args__blastz_chaining','blastz_chaining',
        'args__blastz_threshold','blastz_threshold',
        'args__blastz_mask','blastz_mask',
        'args__blastz_gap_start','blastz_gap_start',
        'args__blastz_gap_extension','blastz_gap_extension',
        'args__blastz_params','blastz_params',

	'args__iw', 'iw',
	'args__fh', 'feat_h',
        'args__gc', 'show_gc',
        'args__nt', 'show_nt',
        'args__cbc', 'show_cbc',
        'args__show_contigs', 'show_contigs',
	'args__color_hsp', 'color_hsp',
	'args__hsp_labels', 'hsp_labels',
	'args__feat_labels', 'feat_labels',
	'args__skip_feat_overlap', 'skip_feat_overlap',
	'args__skip_hsp_overlap', 'skip_hsp_overlap',
	'args__hiqual', 'hiqual',
	'args__email', 'email',
	'args__prog', 'alignment_program',
	'args__showallhsps', 'show_hsps_with_stop_codon',
        'args__padding', 'padding',
        'args__num_seqs','args__$num_seqs',
        'args__draw_model', 'draw_model',
        'args__hsp_overlap_limit', 'hsp_overlap_limit',
        'args__hsp_size_limit', 'hsp_size_limit',
        'args__hsp_top','hsp_top',
        'args__hsp_single_color','hsp_single_color',
        'args__hsp_track','hsp_track',
        'args__comp_adj','comp_adj',
        'args__color_overlapped_features','color_overlapped_features',
        'args__hsp_overlap_length','hsp_overlap_length',
        'args__basefile','args__'+pageObj.basefile,
        'args__show_cns','show_cns',
        'args__show_ofeat','show_ofeat',
        'args__font_size','font_size',
        'args__show_gene_space','show_gene_space',
	'args__ca','ca', 

};

    #deleted params:
    #	'args__hsplim', 'hsp_limit',
    #	'args__hsplimnum', 'hsp_limit_num',
    #
    $params =~ s/\n//g;
    $params =~ s/\s+/ /g;
    return $params;
}

sub gen_go_run {
    my $num_seqs = shift || $NUM_SEQS;
    my $params   = gen_params($num_seqs);
    my $run      = qq!
<SCRIPT language="JavaScript">
function go_run (){
 if (ajax.length)
  {
   setTimeout("go_run()", 100);
   return;
  }
run([$params],[handle_results], 'POST');
setTimeout(" monitor_log()", 5000);
}
</script>!;
    return $run;
}

sub gen_save_settings {
    my $num_seqs      = shift || $NUM_SEQS;
    my $params        = gen_params($num_seqs);
    my $save_settings = qq{save_settings_gevo([$params],[])};
    return $save_settings;
}

sub gen_hsp_colors {
    my %opts     = @_;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my $prefs    = $opts{prefs} || CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GEvo.tmpl' );
    my $hsp_colors;
    my @colors = color_pallet( num_seqs => $num_seqs, prefs => $prefs );
    $template->param( HSP_COLOR_FORM => 1 );
    $template->param( HSP_COLOR_LOOP => \@colors );
    my $count = -2;

    foreach my $line ( split( /\n/, $template->output ) ) {
        next unless $line;
        $count++ if $line =~ /<table>/i;

        if ( $count == 6 ) {
            $line =~ s/(<table>)/<td>$1/i;
            $count = 0;
        }

        $hsp_colors .= $line . "\n";
    }
    return $hsp_colors, scalar @colors;
}

sub color_pallet {
    my %opts     = @_;
    my $start    = $opts{start} || [ 255, 100, 100 ];
    my $offset   = $opts{offset} || 50;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my $prefs    = $opts{prefs} || CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
    $prefs = {} unless $prefs;
    $start = [ $prefs->{r1}, $prefs->{g1}, $prefs->{b1} ]
      if defined $prefs->{r1} && defined $prefs->{g1} && defined $prefs->{b1};

    my @colors;
    my %set = (
        HSP_NUM => 1,
        RED     => $start->[0],
        GREEN   => $start->[1],
        BLUE    => $start->[2]
    );

    #    push @colors, \%set;
    my $temp = [@$start];
    for ( my $i = 1 ; $i <= num_colors($num_seqs) ; $i++ ) {
        my @color;
        @color = @$temp;
        @color =
          ( $prefs->{ "r" . $i }, $prefs->{ "g" . $i }, $prefs->{ "b" . $i } )
          if defined $prefs->{ "r" . $i }
              && defined $prefs->{ "g" . $i }
              && defined $prefs->{ "b" . $i };
        push @colors,
          {
            HSP_NUM => $i,
            RED     => $color[0],
            GREEN   => $color[1],
            BLUE    => $color[2],
          };
        if ( $i % 3 ) {
            $temp = [ map { int( $_ / 1.5 ) } @color ];

        }
        else {
            $start = [ $start->[2], $start->[0], $start->[1] ];
            $temp = [@$start];
        }
        unless ( $i % 9 ) {
            $start = [ $start->[0], int( $start->[0] / 1.5 ), $start->[1] ];
            $temp = [@$start];
        }
        unless ( $i % 18 ) {
            $start = [ $start->[0], $start->[0], int( $start->[0] / 4 ) ];
            $temp = [@$start];
        }
    }
    return wantarray ? @colors : \@colors;
}

sub num_colors {
    my $num_seqs = shift || $NUM_SEQS;
    my $num_colors = 1;
    for ( my $i = $num_seqs - 1 ; $i > 1 ; $i-- ) {
        $num_colors += $i;
    }
    return $num_colors;
}

sub algorithm_list {
    my $program = shift;
    $program = "blastz" unless $program;

    my %progs = (
        "blastn" => "BlastN: Small Regions",
        blastz   => "(B)LastZ: Large Regions",
        "CHAOS"  => "Chaos: Fuzzy Matches",

        #		 "DiAlign_2"=>"DiAlign_2: Glocal Alignment",
        "LAGAN"          => "Lagan: Global Alignment",
        "tblastx"        => "TBlastX: Protein Translation",
        "GenomeThreader" => "GenomeThreader - Spliced Gene Alignments",
        "-=None=-"       => "-=None=-",
    );
    my @programs = sort { lc $a cmp lc $b } keys %progs;

    #    push @programs, "-=NONE=-";
    my $html;
    foreach my $prog (@programs) {
        $html .= "<option value=$prog";
        $html .= $program && $program =~ /$prog/i ? " selected>" : ">";
        $html .= $progs{$prog} . "</option>\n";
    }
    return $html;
}

sub merge_previous    #shabari:for parsing GEvo and tiny urls
{
    my ( $url,     $old_num ) = @_;
    my ( $new_num, $array )   = &parse_url($url);

    my ( $seq_submission, $total_num, $go_run, $ref_seqs ) =
      &add_url_seq( new_num => $new_num, old_num => $old_num, data => $array );

    return ( $seq_submission, $total_num, $go_run, $ref_seqs, $old_num );
}

sub parse_url         #shabari:for parsing GEvo and tiny urls
{
    my $url = shift;
    $url = getlongurl($url) if ( $url =~ /.*tinyurl.*/ || $url =~ /\/r\// );
    my ($numseqs) = $url =~ /num_seqs\=(\d+)/;

    my @array = &makeurlarray( $numseqs, $url );
    return ( $numseqs, \@array );
}

sub getlongurl        #shabari:for parsing GEvo and tiny urls
{
    my $url      = shift;
    my $headHash = head($url);
    my $longurl  = $headHash->{'_previous'}{'_headers'}{'location'};
    return $longurl;
}

sub makeurlarray      #shabari:for parsing GEvo and tiny urls
{
    my ( $numseqs, $longurl ) = @_;
    $numseqs = 2 unless $numseqs;
    my %hash;

    my @array;
    my $i = 1;
    $longurl =~ s/&/;/g;
    while ( $i <= $numseqs ) {
        my $rec = {};
        $rec->{'accn'} = $1 if ( $longurl =~ /accn$i\=(.*?)\;/ );
        $rec->{'accn'} = $1 if ( $longurl =~ /accn$i\=(.*?)\;/ );

        $rec->{'pos'} = $1 if ( $longurl =~ /x$i\=(.*?)\;/ );
        $rec->{'type'} = "position" if ( defined $rec->{'pos'} );

        $rec->{'gbaccn'} = $1 if ( $longurl =~ /gbaccn$i\=(.*?)\;/ );
        $rec->{'type'} = "genbank" if ( defined $rec->{'gbaccn'} );

        $rec->{'fid'} = $1 if ( $longurl =~ /fid$i\=(.*?)\;/ );
        $rec->{'type'} = "feature" if ( defined $rec->{'fid'} );

        $rec->{'dsid'}  = $1 if ( $longurl =~ /dsid$i\=(.*?)\;/ );
        $rec->{'dsgid'} = $1 if ( $longurl =~ /dsgid$i\=(.*?)\;/ );
        $rec->{'gstid'} = $1 if ( $longurl =~ /gstid$i\=(.*?)\;/ );

        $rec->{'chr'} = $1 if ( $longurl =~ /chr$i\=(.*?)\;/ );
        my $up = "dr" . $i . "up";
        $rec->{'drup'} = $1 if ( $longurl =~ /$up\=(.*?)\;/ );
        my $down = "dr" . $i . "down";
        $rec->{'drdown'} = $1 if ( $longurl =~ /$down\=(.*?)\;/ );

        $rec->{'gbstart'}       = $1 if ( $longurl =~ /gbstart$i\=(.*?)\;/ );
        $rec->{'display_order'} = $1 if ( $longurl =~ /do$i\=(.*?)\;/ );
        $rec->{'display_order'} = $i unless $rec->{'display_order'};
        push @array, $rec;

        $i++;
    }
    @array = sort { $a->{display_order} <=> $b->{display_order} } @array;
    return (@array);

}

sub add_url_seq {
    my %opts    = @_;
    my $new_num = $opts{new_num};
    my $old_num = $opts{old_num};
    my $data    = $opts{data};

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GEvo.tmpl' );

    my $total_num = $new_num + $old_num;

    my $prefs = CoGe::Accessory::Web::load_settings(
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );

    my $message;
    my @seq_nums;
    my @seq_sub;
    my $i = $old_num;
    for my $element (@$data) {
        $i++;
        my $draccn = $element->{'accn'} if $element->{'accn'};
        my $pos    = $element->{'pos'}  if $element->{'pos'};
        my $fid    = $element->{'fid'}  if $element->{'fid'};

        if ( $fid && !$draccn ) {
            my $feat = $coge->resultset('Feature')->find($fid);
            if ($feat) {
            	($draccn) = $feat->primary_name;
            }
        }

        my $drup   = $element->{'drup'}   if $element->{'drup'};
        my $drdown = $element->{'drdown'} if $element->{'drdown'};
        $drup   = 10000 unless defined $drup;
        $drdown = 10000 unless defined $drdown;

        my $dsid    = $element->{'dsid'};
        my $dsgid   = $element->{'dsgid'};
        my $gstid   = $element->{'gstid'};
        my $mask    = $element->{'mask'};
        my $gbaccn  = $element->{'gbaccn'} if $element->{'gbaccn'};
        my $gbstart = $element->{'gbstart'} if $element->{'gbstart'};
        $gbstart = 1 unless defined $gbstart;
        my $gblength = $element->{'gblength'} if $element->{gblength};

        my $chr = $element->{'chr'} if defined $element->{'chr'};

        my $revn = "checked";
        my $revy;

        my $refn = " ";
        my $refy = "checked";

        my $org_title;
        my $dsg_menu = qq{<input type="hidden" id="dsgid$i"};
        $dsg_menu .= qq{value ="$dsgid"} if $dsgid;
        $dsg_menu .= qq{>};

        ( $org_title, $dsid, $gstid, $dsgid ) = get_org_info(
            dsid    => $dsid,
            chr     => $chr,
            gstid   => $gstid,
            dsgid   => $dsgid,
            seq_num => $i
        ) if $pos;
        my %opts = (
            SEQ_NUM   => $i,
            REV_YES   => $revy,
            REV_NO    => $revn,
            REF_YES   => $refy,
            REF_NO    => $refn,
            DRUP      => $drup,
            DRDOWN    => $drdown,
            DRACCN    => $draccn,
            DSID      => $dsid,
            DSGID     => $dsgid,
            DSG_STUFF => $dsg_menu,
            GBACCN    => $gbaccn,
            GBSTART   => $gbstart,
            GBLENGTH  => $gblength,
            POS       => $pos,
            ORGINFO   => $org_title,
            CHR       => $chr,

        );

        if ($mask) {
            $opts{MASK_CDS}    = "selected" if $mask eq "cds";
            $opts{MASK_RNA}    = "selected" if $mask eq "rna";
            $opts{MASK_NCDS}   = "selected" if $mask eq "non-cds";
            $opts{MASK_NGENIC} = "selected" if $mask eq "non-genic";
        }
        $opts{COGEPOS} =
qq{<option value="cogepos$i" selected="selected">CoGe Database Position</option>}
          if $pos;
        push @seq_sub, {%opts};
    }

    #generate sequence submission selector
    $template->param( SEQ_SELECT      => 1 );
    $template->param( SEQ_SELECT_LOOP => \@seq_sub );
    my $seq_submission = $template->output;

    $template->param( SEQ_SELECT => 0 );

    $template->param( GEN_REFERENCE_SEQUENCES => 1 );
    $template->param( REF_SEQ => [ { SEQ_NUM => $total_num } ] );
    my $ref_seqs = $template->output;

    $template->param( GEN_REFERENCE_SEQUENCES => 0 );
    my $go_run = gen_go_run($total_num);
    return ( $seq_submission, $total_num, $go_run, $ref_seqs, $old_num );

}

sub add_seq {
    my $num_seq = shift;
    $num_seq++;

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'GEvo.tmpl' );
    my $hsp_colors;
    if ( $num_seq > $MAX_SEQS ) {
        ($hsp_colors) = gen_hsp_colors( num_seqs => $MAX_SEQS );
        my $go_run = gen_go_run($MAX_SEQS);
        return ( '', $MAX_SEQS, $go_run, '', $hsp_colors,
            qq{Exceeded max number of sequences ($MAX_SEQS)} );
    }
    my @seqs = {
        SEQ_NUM   => $num_seq,
        REV_NO    => "checked",
        REF_YES   => "checked",
        DRUP      => 10000,
        DRDOWN    => 10000,
        GBSTART   => 1,
        DSG_STUFF => qq{<input type="hidden" id="dsgid$num_seq" value="">},
    };

    $template->param( SEQ_SELECT      => 1 );
    $template->param( SEQ_SELECT_LOOP => \@seqs );
    my $seq_submission = $template->output;

    $template->param( SEQ_SELECT => 0 );

    $template->param( GEN_REFERENCE_SEQUENCES => 1 );
    $template->param( REF_SEQ => [ { SEQ_NUM => $num_seq } ] );
    my $ref_seqs = $template->output;
    my $go_run   = gen_go_run($num_seq);

    #    $hsp_colors = gen_hsp_colors(num_seqs=>$num_seq);

    return ( $seq_submission, $num_seq, $go_run, $ref_seqs );
}

sub hsp_color_cookie {
    my $cname   = "GEVO_HSP_COLORS";
    my %cookies = fetch CGI::Cookie;
    if ( ref $cookies{$cname} ) {
    }
    else {
    }
}

sub get_algorithm_options {
    my %opts             = @_;
    my $analysis_program = $opts{prog};

    #blast stuff
    my $blast_wordsize  = $opts{blast_word};
    my $blast_gapopen   = $opts{blast_gapo};
    my $blast_gapextend = $opts{blast_gape};
    my $blast_mismatch  = $opts{blast_mmatch};
    my $blast_match     = $opts{blast_match};
    my $blast_eval      = $opts{blast_eval};
    my $blast_params    = $opts{blast_params};
    my $blast_filter    = $opts{blast_filter};

    #blastz stuff
    my $blastz_wordsize      = $opts{blastz_wordsize};
    my $blastz_chaining      = $opts{blastz_chaining};
    my $blastz_threshold     = $opts{blastz_threshold};
    my $blastz_mask          = $opts{blastz_mask};
    my $blastz_gap_start     = $opts{blastz_gap_start};
    my $blastz_gap_extension = $opts{blastz_gap_extension};
    my $blastz_params        = $opts{blastz_params};

    #lagan_stuff
    my $lagan_params     = $opts{lagan_params};
    my $lagan_min_length = $opts{lagan_min_length};
    my $lagan_max_gap    = $opts{lagan_max_gap};
    my $lagan_percent_id = $opts{lagan_percent_id};

    #dialign_stuff
    my $dialign_params         = $opts{dialign_params};
    my $dialign_motif_motif    = uc( $opts{dialign_motif_motif} );
    my $dialign_motif_weight   = $opts{dialign_motif_weight};
    my $dialign_motif_penalty  = $opts{dialign_motif_penalty};
    my $dialign_min_score      = $opts{dialign_min_score};
    my $dialign_max_gap        = $opts{dialign_max_gap};
    my $dialign_split_score    = $opts{dialign_split_score};
    my $dialign_use_anchor     = $opts{dialign_use_anchor};
    my $dialign_anchor_program = $opts{dialign_anchor_program};

    #chaos stuff
    my $chaos_word_length    = $opts{chaos_word_length};
    my $chaos_score_cutoff   = $opts{chaos_score_cutoff};
    my $chaos_rescore_cutoff = $opts{chaos_rescore_cutoff};
    my $chaos_lookback       = $opts{chaos_lookback};
    my $chaos_gap_length     = $opts{chaos_gap_length};
    my $chaos_gap_start      = $opts{chaos_gap_start};
    my $chaos_gap_extension  = $opts{chaos_gap_extension};
    my $chaos_params         = $opts{chaos_params};

    my (
        $blastz_string, $blastz_link,    $blast_string, $blast_link,
        $lagan_string,  $dialign_string, $chaos_string
    );

    #build blastz param string
    $blastz_string = " W=" . $blastz_wordsize if $blastz_wordsize;
    $blastz_string .= " C=" . $blastz_chaining      if $blastz_chaining;
    $blastz_string .= " K=" . $blastz_threshold     if $blastz_threshold;
    $blastz_string .= " M=" . $blastz_mask          if $blastz_mask;
    $blastz_string .= " O=" . $blastz_gap_start     if $blastz_gap_start;
    $blastz_string .= " E=" . $blastz_gap_extension if $blastz_gap_extension;
    $blastz_string .= " " . $blastz_params          if $blastz_params;
    $blastz_link = ";bzW=" . $blastz_wordsize if $blastz_wordsize;
    $blastz_link .= ";bzC=" . $blastz_chaining      if $blastz_chaining;
    $blastz_link .= ";bzK=" . $blastz_threshold     if $blastz_threshold;
    $blastz_link .= ";bzM=" . $blastz_mask          if $blastz_mask;
    $blastz_link .= ";bzO=" . $blastz_gap_start     if $blastz_gap_start;
    $blastz_link .= ";bzE=" . $blastz_gap_extension if $blastz_gap_extension;

    #build blast param string
    $blast_wordsize = 3
      if ( $analysis_program eq "tlastx" && $blast_wordsize > 3 );
    $blast_string = " -W " . $blast_wordsize if defined $blast_wordsize;
    $blast_string .= " -G " . $blast_gapopen
      if defined $blast_gapopen && $analysis_program eq "blastn";
    $blast_string .= " -E " . $blast_gapextend
      if defined $blast_gapextend && $analysis_program eq "blastn";
    $blast_string .= " -q " . $blast_mismatch
      if defined $blast_mismatch && $analysis_program eq "blastn";
    $blast_string .= " -r " . $blast_match
      if defined $blast_match && $analysis_program eq "blastn";
    $blast_string .= " -e " . $blast_eval if defined $blast_eval;
    $blast_string .= " -F " . $blast_filter
      if defined $blast_filter && $analysis_program eq "blastn";
    $blast_string .= " " . $blast_params if defined $blast_params;
    $blast_link = ";bnW=" . $blast_wordsize if defined $blast_wordsize;
    $blast_link .= ";bnG=" . $blast_gapopen
      if defined $blast_gapopen && $analysis_program eq "blastn";
    $blast_link .= ";bnE=" . $blast_gapextend
      if defined $blast_gapextend && $analysis_program eq "blastn";
    $blast_link .= ";bnq=" . $blast_mismatch
      if defined $blast_mismatch && $analysis_program eq "blastn";
    $blast_link .= ";bnr=" . $blast_match
      if defined $blast_match && $analysis_program eq "blastn";
    $blast_link .= ";bne=" . $blast_eval if defined $blast_eval;
    $blast_link .= ";bnF=" . $blast_filter
      if defined $blast_filter && $analysis_program eq "blastn";

    #build lagan param string
    $lagan_string = $lagan_params;

    #build dialign param string
    $dialign_motif_motif =~ s/[^ATCG\[\]]//g if $dialign_motif_motif;
    $dialign_string .=
      "-mot $dialign_motif_motif $dialign_motif_weight $dialign_motif_penalty"
      if $dialign_motif_motif =~ /^[ATGC\[\]]+$/;
    $dialign_string .= " " . $dialign_params;

    #build chaos param string
    $chaos_string .= " -v";    #need verbose mode for the parser
    $chaos_string .= " -b";    #search both strands
    $chaos_string .= " -wl $chaos_word_length" if defined $chaos_word_length;
    $chaos_string .= " -co $chaos_score_cutoff" if defined $chaos_score_cutoff;
    $chaos_string .= " -rsc $chaos_rescore_cutoff"
      if defined $chaos_rescore_cutoff;
    $chaos_string .= " -lb $chaos_lookback"   if defined $chaos_lookback;
    $chaos_string .= " -gl $chaos_gap_length" if defined $chaos_gap_length;
    $chaos_string .= " -gs $chaos_gap_start"  if defined $chaos_gap_start;
    $chaos_string .= " -gc $chaos_gap_extension"
      if defined $chaos_gap_extension;
    $chaos_string .= " $chaos_params" if defined $chaos_params;

    #stuff to be returned
    my $param_string
      ;    #param string to be passed to command line for program execution
    my $link_string;    #params for GEvo link
    my %parser_options
      ; #params used to set the parser object -- key is accessor method, value is value

    if ( $analysis_program eq "blastz" ) {
        $param_string = $blastz_string;
        $link_string  = $blastz_link;
    }
    elsif ( $analysis_program =~ /blast/i )    #match blastn and tblastx
    {
        $param_string = $blast_string;
        $link_string  = $blast_link;
    }
    elsif ( $analysis_program eq "LAGAN" ) {
        $param_string                   = $lagan_string;
        $parser_options{max_gap}        = $lagan_max_gap;
        $parser_options{length_cutoff}  = $lagan_min_length;
        $parser_options{percent_cutoff} = $lagan_percent_id;
    }
    elsif ( $analysis_program eq "DIALIGN" ) {
        $param_string                = $dialign_string;
        $parser_options{min_score}   = $dialign_min_score;
        $parser_options{max_gap}     = $dialign_max_gap;
        $parser_options{split_score} = $dialign_split_score;
        my ( $anchor_program, $anchor_param_string );
        if ( $dialign_anchor_program =~ /chaos/i ) {
            ( $anchor_program, $anchor_param_string ) =
              ( "CHAOS", $chaos_string );
        }
        elsif ( $dialign_anchor_program =~ /bl2/i ) {
            ( $anchor_program, $anchor_param_string ) =
              ( "bl2seq", $blast_string );
        }
        else {
            ( $anchor_program, $anchor_param_string ) =
              ( "blastz", $blastz_string );
        }
        $parser_options{anchor_params} = {
            prog         => $anchor_program,
            param_string => $anchor_param_string,
          }
          if $dialign_use_anchor;
    }
    elsif ( $analysis_program eq "CHAOS" ) {
        $param_string = $chaos_string;
    }
    my ( $x, $clean_param_string ) =
      CoGe::Accessory::Web::check_taint($param_string);
    unless ($x) {
        print STDERR "user "
          . $USER->user_name
          . "'s param_string: $param_string was tainted!\n";
    }
    return ( $analysis_program, $clean_param_string, \%parser_options,
        $link_string );
}

sub number_of_runs {
    my $sets = shift;
    return unless $sets;
    my $refs = 0;
    my $num  = @$sets;
    $num--;
    foreach my $item ( map { $_->{reference_seq} } @$sets ) {
        $refs++ if $item;
    }
    my $total = 0;
    for ( 1 .. $refs ) {
        $total += $num;
        $num--;
    }
    return $total;
}

sub dataset_search {
    my %opts = @_;
    my $accn = $opts{accn};
    my $num  = $opts{num};
    $num = 1 unless $num;
    my $dsid   = $opts{dsid};
    my $featid = $opts{featid};
    my $dsgid  = $opts{dsgid};
    my $gstid  = $opts{gstid};
    $accn =~ s/^\s+// if $accn;
    $accn =~ s/\s+$// if $accn;

    $accn = uri_unescape($accn);    # mdb added 10/9/12 issue #11

    #	print STDERR "dataset_search: accn=$accn\n";

    my $feat = $coge->resultset('Feature')->find($featid) if $featid;
    $dsid = $feat->dataset->id if $feat;
    my $html;
    my %sources;
    my $rs;
    if ($feat) {
        $rs =
          $coge->resultset('Dataset')
          ->search( { dataset_id => $feat->dataset->id } );
    }
    elsif ($accn) {
        $rs = $coge->resultset('Dataset')->search(
            { 'feature_names.name' => $accn, },
            {
                'join' => { 'features' => 'feature_names', },

                #						   'prefetch'=>['datasource', 'organism'],
            }
        );
    }
    return unless $rs;
    while ( my $ds = $rs->next() ) {
        my $skip = 0;
        #next if $ds->deleted;  #This function hasn't been pushed out yet
        foreach my $item ( $ds->genomes ) {
            next if $USER->is_admin;
            $skip = 1 if ($item->deleted || !$USER->has_access_to_genome($item));
        }
        next if $skip;
        my $ver     = $ds->version;
        my $desc    = $ds->description;
        my $sname   = $ds->data_source->name;
        my $ds_name = $ds->name;
        foreach my $seqtype ( $ds->sequence_type ) {
            my $typeid = $seqtype->id if $seqtype;
            unless ($typeid) {
                print STDERR "Error retrieving sequence_type object for ",
                  $ds->name, ": dsid ", $ds->id, "\n";
                next;
            }
            my $title =
                "$ds_name ($sname, "
              . ( $ver ? "v$ver, " : '' ) . "dsid"
              . $ds->id . ")";

            next if $sources{ $ds->id } && $sources{ $ds->id }{typeid} < $typeid;
            if ( $dsgid && !$dsid ) {
                foreach my $item ( $ds->genomes ) {
                    $dsid = $ds->id if $dsgid == $item->id;
                }
            }
            $sources{ $ds->id } = {
                title   => $title,
                version => $ver,
                typeid  => $typeid
            };
        }
    }

    if ( keys %sources ) {

# 	$html .= qq{
# <SELECT name = "dsid$num" id= "dsid$num" onChange="feat_search(['args__accn','accn$num','args__dsid', 'dsid$num', 'args__num','args__$num', 'args__featid'],['feat$num']);">
# };
        $html .= qq{
 <SELECT name = "dsid$num" id= "dsid$num" onChange="genome_search(['args__dsid', 'dsid$num', 'args__dsgid', 'dsgid$num', 'args__gstid', 'gstid$num', 'args__num','args__$num', 'args__featid', 'featid$num'],[feat_search_chain]);">
 };
        foreach my $id (
            sort {
                     $sources{$b}{version} cmp $sources{$a}{version}
                  || $sources{$a}{typeid} <=> $sources{$b}{typeid}
            } keys %sources
          )
        {
            my $val = $sources{$id}{title};
            $html .= qq{  <option value="$id"};
            $html .= qq{ selected } if ( $dsid && $id == $dsid );
            $html .= qq{>$val\n};
        }
        $html .= qq{</SELECT>\n};
        my $count = scalar keys %sources;
        $html .= qq{<font class=small>($count)</font>} if $count > 1;
    }
    else {
        $html .=
qq{<span class="small container">Accession not found for $accn</span><hidden id=dsid$num>}
          if $accn;
    }
    return ( $html, $num, $featid );
}

sub genome_search {
    my %opts   = @_;
    my $dsid   = $opts{dsid};     #dataset database id
    my $num    = $opts{num};      #sequence submission number
    my $featid = $opts{featid};
    my $dsgid  = $opts{dsgid};
    my $gstid  = $opts{gstid};
    my $ds = $coge->resultset('Dataset')->find($dsid);
    unless ($ds) {
        my $html = "<input type=hidden id=dsgid$num>";
        $html .= "No genome found for $dsid" if $dsid;
        return $html, $num;
    }
    my $html =
qq{<SELECT name="dsgid$num" id="dsgid$num" onChange="feat_search(['args__accn','accn$num','args__dsid', 'dsid$num','args__dsgid', 'dsgid$num', 'args__num','args__$num', 'args__featid', 'args__$featid'],['feat$num']);">};
    my $count = 0;
    foreach my $dsg (
        sort {
            versioncmp( $b->version, $a->version )
              || $a->type->id <=> $b->type->id
        } $ds->genomes
      )
    {
        next if $dsg->deleted;
        my $dsgid_tmp = $dsg->id;
        my $title     = $dsg->name;
        $title = $dsg->organism->name unless $title;
        $title .= " (v" . $dsg->version . " " . $dsg->type->name . ")";
        $html  .= "<option value = $dsgid_tmp";
        $html  .= " selected"
          if ( $dsgid && $dsg->id eq $dsgid )
          || ( $gstid && $dsg->type->id eq $gstid );
        $html .= ">$title";
        $count++;
    }
    $html .= "</select>";
    $html .= "<span class=small>($count)</span>" if $count > 1;

    return ( $html, $num, $featid ) if $count > 0;


    $html = qq{<select name="dsgid$num" id="dsgid$num" onChange="feat_search(['args__accn','accn$num','args__dsid', 'dsid$num','args__dsgid', 'dsgid$num', 'args__num','args__$num', 'args__featid', 'args__$featid'],['feat$num']);">};

    $html .= "<option>The genome has been deleted</option></select>";

    return ( $html, $num, $featid);
}

sub save_settings_gevo {
    my %opts = @_;
    my $opts = Dumper \%opts;
    print STDERR $opts;
    my $item = CoGe::Accessory::Web::save_settings(
        opts => $opts,
        user => $USER,
        page => $PAGE_NAME,
        coge => $coge
    );
}

sub reset_settings_gevo {
    my %opts = @_;
    my $item =
      reset_settings( user => $USER, page => $PAGE_NAME, coge => $coge );
}

sub get_opt {
    my %opts   = @_;
    my $params = $opts{params};
    my $form   = $opts{form} || $FORM;
    my $param  = $opts{param};
    my $opt;
    $opt = $form->param($param) if defined $form->param($param);
    $opt = $params->{$param} if ( ref($params) =~ /hash/i & !defined $opt );
    return $opt;

}

sub get_org_info {
    my %opts  = @_;
    my $dsid  = $opts{dsid};
    my $chr   = $opts{chr};
    my $dsgid = $opts{dsgid};
    my $gstid = $opts{gstid};
    my $num   = $opts{seq_num};
    my ( $ds, $dsg, $gst );
    if ($dsgid) {
        $dsg = $coge->resultset('Genome')->find($dsgid);
        $gst = $dsg->type if $dsg;
    }
    elsif ($dsid) {
        my $ds = $coge->resultset('Dataset')->find($dsid);
        if ($ds) {
            foreach my $item ( $ds->genomes( chr => $chr ) ) {
                if ($gstid) {
                    $dsg = $item if $item->type->id eq $gstid;
                }
                else {
                    $dsg = $item;
                    last;
                }
            }
            $gst   = $dsg->type if $dsg;
            $dsgid = $dsg->id   if $dsg;
        }
    }
    return "<span class=\"small alert\">Dataset group was not found</span>"
      . qq{<input type="hidden" id="dsgid$num">}
      unless $dsg;
    my $dsg_menu =
qq{<span class="small">Genome: </span><SELECT name="dsgid$num" id="dsgid$num">};
    foreach my $item (
        sort {
            versioncmp( $b->version, $a->version )
              || $a->type->id <=> $b->type->id
        } $dsg->organism->genomes
      )
    {
        if ( $dsg->restricted && !$USER->has_access_to_genome($dsg) ) {
            $dsg_menu .= "<option>Restricted</option>";
            next;
        }

        if ($dsg->deleted) {
            $dsg_menu .= "<option>The genome has been been deleted</option>";
            next;
        }

        my $dsgid_tmp = $item->id;
        my $title;
        $title = $item->name . " " if $item->name;

        #	$title = $dsg->organism->name unless $title;
        $title    .= "v" . $item->version . " " . $item->type->name;
        $title    .= " (dsgid" . $item->id . ")";
        $dsg_menu .= "<option value = $dsgid_tmp";
        $dsg_menu .= " selected" if ( $dsgid && $item->id eq $dsgid );
        $dsg_menu .= ">$title";
    }
    $dsg_menu .= "</select>";
    $gstid = $gst->id if $gst;

    my $ver  = $dsg->version;
    my $org  = $dsg->organism->name;
    my $oid  = $dsg->organism->id;
    my $type = $gst->name if $gst;
    $chr = join( ", ", $dsg->get_chromosomes ) unless defined $chr;

#    my $title = qq{<span class="link"><a class="link" href="OrganismView.pl?oid=$oid;dsgid=$dsgid" target=_new>$org $type (v$ver):</a></span> chr. $chr};
    my $title =
qq{<span class="link"><a class="link" href="OrganismView.pl?oid=$oid" target=_new>$org</a></span>};
    return ( $title, $dsid, $gstid, $dsgid, $dsg_menu );
}

sub feat_search {
    my %opts   = @_;
    my $accn   = $opts{accn};
    my $dsid   = $opts{dsid};
    my $dsgid  = $opts{dsgid};
    my $gstid  = $opts{gstid};
    my $num    = $opts{num};
    my $featid = $opts{featid};
    $accn = uri_unescape($accn);    # mdb added 10/9/12 issue #11

    #print STDERR "feat_search: accn=$accn\n";
    $accn =~ s/^\s+// if $accn;
    $accn =~ s/\s+$// if $accn;
    my $feat = $coge->resultset('Feature')->find($featid) if $featid;
    if ($feat) {
        $dsid = $feat->dataset->id unless $dsid;
    }
    unless ($gstid) {
        if ( $featid =~ /_/ ) {
            ( $featid, $gstid ) = split( /_/, $featid );
        }
        if ($dsgid) {
            my $dsg = $coge->resultset('Genome')->find($dsgid);
            $gstid = $dsg->type->id if $dsg;
        }
    }
    $gstid = 1 unless $gstid;    #need a default; go unmasked
    return qq{<input type="hidden" id="featid$num">\n} unless $dsid;
    my @feats;
    my $rs = $coge->resultset('Feature')->search(
        {
            'feature_names.name' => $accn,
            'dataset.dataset_id' => "$dsid",
        },
        {
            'join'     => [ 'feature_type', 'dataset', 'feature_names' ],
            'prefetch' => [ 'feature_type', 'dataset' ],
        }
    );
    my %seen;
    my @genes;    #get genes on top of list!
    while ( my $f = $rs->next() ) {
        next unless $f->dataset->id == $dsid;
        unless ( $seen{ $f->id } ) {
            if ( $f->type->name eq "gene" ) {
                push @genes, $f;
            }
            else {
                push @feats, $f;
            }
        }
        $seen{ $f->id } = 1;
    }
    unless (@feats) {
        push @feats, $feat if $feat;
    }
    @feats = sort { $a->type->name cmp $b->type->name } @feats;
    unshift @feats, @genes if @genes;

    my $feat_count;
    my $html;
    if (@feats) {
        $html .= qq{<SELECT name="featid$num" id="featid$num">};
        foreach my $feat (@feats) {
            next unless $feat->dataset->id == $dsid; # mdb added 2/12/15 COGE-588
            my $loc = "(" . $feat->type->name . ") Chr:"
              . $feat->chromosome
              . " " . commify( $feat->start ) . "-" . commify( $feat->stop );
            $loc =~ s/(complement)|(join)//g;
            my $fid = $feat->id;
            $fid .= "_" . $gstid if $gstid;
            $html .= qq{  <option value="$fid"};
            $html .= qq{ selected } if $featid && $featid == $feat->id;
            $html .= qq{>$loc \n};
            $feat_count++;
        }
        $html .= qq{</SELECT>\n};
        $html .= qq{<font class="small">($feat_count)</font>} if $feat_count > 1;
    }
    else {
        $html .= qq{<input type="hidden" id="featid$num">\n};
    }
    return $html;
}

sub email_results {
    my %opts          = @_;
    my $email_address = $opts{email};
    my $basefilename  = $opts{basefile};
    my $server        = $P->{SERVER};
    my $full_gevo_url = $opts{full_gevo_url};
    return
      unless $email_address =~
/^[_a-zA-Z0-9-]+(\.[_a-zA-Z0-9-]+)*@[a-zA-Z0-9-]+(\.[a-zA-Z0-9-]+)*\.(([0-9]{1,3})|([a-zA-Z]{2,3})|(aero|coop|info|museum|name))$/;
    my $mailer = Mail::Mailer->new("sendmail");
    $mailer->open(
        {
            From    => 'GEvo <gevo_results@genomevolution.org>',
            To      => $email_address,
            Subject => 'GEvo Analysis Results Ready',
        }
    ) or die "Can't open: $!\n";
    my $name = $USER->first_name;
    $name .= " " . $USER->last_name if $USER->last_name;    #user_name;
    my $body = qq{Dear $name,

Thank you for using the Genome Evolution Analysis Tool! The results from your latest analysis are ready, and can be viewed here:
} . $server . qq{/GEvo_direct.pl?name=$basefilename

To contact us or to cite CoGe please visit:
http://genomevolution.org/wiki/index.php/Contact_Page

You can use this URL for regenerating your results:
$full_gevo_url

Thank you for using the CoGe Software Package.

- The CoGe Team
};

    print $mailer $body;
    $mailer->close();
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

sub create_seq_file_name {
    my %opts = @_;
    my $name;
    $name .= "f" . $opts{featid}  if $opts{featid};
    $name .= "p" . $opts{pos}     if $opts{pos};
    $name .= "ds" . $opts{dsid}   if $opts{dsid};
    $name .= "r" . $opts{rev}     if $opts{rev};
    $name .= "c" . $opts{chr}     if $opts{chr};
    $name .= "u" . $opts{up}      if $opts{up};
    $name .= "d" . $opts{down}    if $opts{down};
    $name .= "g" . $opts{gstid}   if $opts{gstid};
    $name .= "m" . $opts{mask}    if $opts{mask};
    $name .= "dsg" . $opts{dsgid} if $opts{dsgid};
    return $TEMPDIR . "/" . $name . ".faa";
}

sub get_tiny_url {
    my %opts = @_;
    my $url  = $opts{url};
    my $tiny = CoGe::Accessory::Web::get_tiny_link(
        db      => $coge,
        user_id => $USER->id,
        page    => $PAGE_NAME,
        url     => $url,
        log_msg => "GEvo link",
    );
    my $html .= qq{<a href="$tiny" onclick="window.open('$tiny');" target=_new>$tiny</a>};
    return $html;
}

sub add_to_user_history {
    my %opts = @_;
    my ($link) = $opts{url} =~ />(http.+)<br/;
    $USER->add_to_works(
        {
            'name'        => $opts{work_name},
            'archive'     => $opts{archive},
            'page'        => $PAGE_NAME,
            'parameter'   => $link,
            'description' => $opts{description},
            'note'        => $opts{note},
        }
    );
}

sub get_image_info {
    my %opts         = @_;
    my $idx          = $opts{id};
    my $basefilename = $opts{basename};
    return ("no opts specified") unless ( $idx && $basefilename );
    $cogeweb = CoGe::Accessory::Web::initialize_basefile(
        basename => $basefilename,
        tempdir  => $TEMPDIR
    );
    my $dbh =
      DBI->connect( "dbi:SQLite:dbname=" . $cogeweb->sqlitefile, "", "" );
    my $query = qq{select * from image_info where id = $idx;};
    my (
        $id,         $display_id, $name,  $title,
        $px_witdth,  $bpmin,      $bpmax, $dsid,
        $chromosome, $rc,         $px_height
    ) = $dbh->selectrow_array($query);
    return ("$bpmin||$bpmax");
}
