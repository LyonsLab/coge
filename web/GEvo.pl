#!/usr/bin/perl -w
use strict;
use warnings;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Cookie;
use CGI::Ajax;
use HTML::Template;
use Data::Dumper;
use File::Basename;
use File::Temp;
use CoGe::Accessory::GenBank;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::Restricted_orgs;
use CoGe::Accessory::bl2seq_report;
use CoGe::Accessory::blastz_report;
use CoGe::Accessory::lagan_report;
use CoGe::Accessory::chaos_report;
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
use CoGeX;
use CoGeX::Feature;
use DBIxProfiler;
use DBI;
use LWP::Simple;
use Parallel::ForkManager;
use Statistics::Basic;
use Benchmark qw(:all);

# for security purposes

$ENV{PATH} = "/opt/apache/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
#for chaos
$ENV{'LAGAN_DIR'} = '/opt/apache/CoGe/bin/lagan/';
#for dialign
$ENV{'DIALIGN2_DIR'} = '/opt/apache/CoGe/bin/dialign2_dir/';
use vars qw( $PAGE_NAME $DATE $DEBUG $BL2SEQ $BLASTZ $LAGAN $CHAOS $DIALIGN $TEMPDIR $TEMPURL $USER $FORM $cogeweb $BENCHMARK $coge $NUM_SEQS $MAX_SEQS $REPEATMASKER $MAX_PROC);
$PAGE_NAME = "GEvo.pl";
$BL2SEQ = "/usr/bin/bl2seq ";
$BLASTZ = "/usr/bin/blastz ";
$LAGAN = "/opt/apache/CoGe/bin/lagan/lagan.pl";
$CHAOS = "/opt/apache/CoGe/bin/lagan/chaos_coge";
$DIALIGN = "/opt/apache/CoGe/bin/dialign2_dir/dialign2-2_coge";
$REPEATMASKER = "/opt/apache/CoGe/bin/RepeatMasker/RepeatMasker";
$TEMPDIR = "/opt/apache/CoGe/tmp/GEvo";
$TEMPURL = "/CoGe/tmp/GEvo";
$MAX_PROC=8;
# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$BENCHMARK = 1;
$NUM_SEQS = 3;
$MAX_SEQS = 21;
$| = 1; # turn off buffering 
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$CGI::POST_MAX= 60 * 1024 * 1024; # 24MB
$CGI::DISABLE_UPLOADS = 0; 
($USER) = CoGe::Accessory::LogUser->get_user();

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

my %ajax = CoGe::Accessory::Web::ajax_func();
$ajax{dataset_search} = \&dataset_search; #override this method from Accessory::Web for restricted organisms
$ajax{feat_search} = \&feat_search; 
my $pj = new CGI::Ajax(
		       run=>\&run,
		       loading=>\&loading,
		       add_seq=>\&add_seq,
		       get_file=>\&get_file,
		       gen_go_run=>\&gen_go_run,
		       gen_hsp_colors =>\&gen_hsp_colors,
		       save_settings_gevo=>\&save_settings_gevo,
		       reset_settings_gevo=>\&reset_settings_gevo,
		       %ajax,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
#$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);

#$USER="elyons";print gen_html();

sub loading
  {
    my $message = shift || "Generating results. . .";
    return qq{<div class="dna"><div id="loading">$message</div></div>}; 
  }

sub gen_html
  {
    my $html;# =  "Content-Type: text/html\n\n";
    unless ($USER)
      {
	$html = login();
      }
    else
      {
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	$template->param(LOGO_PNG=>"GEvo-logo.png");
	$template->param(TITLE=>'Genome Evolution Analysis');
	$template->param(HELP=>'GEvo');
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(NO_BOX=>1);
	$template->param(BODY=>gen_body());
	my $prebox = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
	$prebox->param(RESULTS_DIV=>1);
	$template->param(PREBOX=>$prebox->output);
	
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $form = $FORM;
    my $prefs = load_settings(user=>$USER, page=>$PAGE_NAME);
    my $num_seqs = get_opt(params=>$prefs, form=>$form, param=>'num_seqs');
    $num_seqs = $NUM_SEQS unless defined $num_seqs;
    $MAX_SEQS = 10 if $form->param('override');
    my $message;
    if (! ($num_seqs =~ /^\d+$/) )
      {
	$message .= "Problem with requested number of sequences: '$num_seqs'.  Defaulting to $NUM_SEQS input sequences.";
	$num_seqs = $NUM_SEQS;
      }
    elsif ($num_seqs < 2)
      {
	$message .= "Minimum number of sequences to compare is two.";
	$num_seqs = 2;
      }
    elsif ($num_seqs > $MAX_SEQS)
      {
	$message .= "Maximum number of sequence is $MAX_SEQS.";
	$num_seqs = $MAX_SEQS;
      }
    my @seq_nums;
    my @seq_sub;
    my $autosearch_string;
    for (my $i = 1; $i <= $num_seqs; $i++)
      {
	my $draccn;
	$draccn= $form->param("accn".$i) if $form->param("accn".$i);
	my $pos = $form->param("x".$i) if $form->param("x".$i);
	my $chr = $form->param("chr".$i) if $form->param("chr".$i);
	my $drfid = $form->param("fid".$i) if $form->param("fid".$i);
	if ($drfid && !$draccn)
	  {
	    my $feat = $coge->resultset('Feature')->find($drfid);
	    ($draccn) = $feat->primary_name if $feat;
	  }
	my $drup = $form->param('dr'.$i.'up') if defined $form->param('dr'.$i.'up');
	my $drdown = $form->param('dr'.$i.'down') if defined $form->param('dr'.$i.'down');
	$drup = $form->param('drup'.$i) if defined $form->param('drup'.$i);
	$drdown = $form->param('drdown'.$i) if defined $form->param('drdown'.$i);
	$drup = 10000 unless defined $drup;
	$drdown = 10000 unless defined $drdown;
#	$drup += $pad_gs if $pad_gs;
#	$drdown += $pad_gs if $pad_gs;
	my $dsid = $form->param('dsid'.$i) if $form->param('dsid'.$i);
	my $gbaccn = $form->param("gbaccn".$i) if $form->param("gbaccn".$i);
	my $gbstart = $form->param("gbstart".$i) if $form->param("gbstart".$i);
	$gbstart = 1 unless defined $gbstart;
	my $gblength = $form->param("gblength".$i) if $form->param("gblength".$i);
	my $revy = "checked" if $form->param('rev'.$i);
	my $revn = "checked" unless $revy;
	my $maskexonon = "checked" if $form->param('maskexon'.$i);
	my $maskexonoff = "checked" unless $maskexonon;
	my $refn = "checked" if $form->param('nref'.$i);
	my $refy = "checked" unless $refn;
	$autosearch_string .= 'if ($'.qq!('#accn$i').val()) {dataset_search(['args__accn','accn$i','args__num', 'args__$i'!;
	$autosearch_string .= ", 'args__dsid', 'args__$dsid'" if $dsid;
	$autosearch_string .= ", 'args__featid', 'args__$drfid'" if $drfid;
	$autosearch_string .=qq!],[feat_search_chain]);}!;
	my $dsinfo = $dsid ? qq{<input type="hidden" id="posdsid$i" value="$dsid">} : qq{<input type="hidden" id="posdsid$i">};
	$dsinfo .= $chr ? qq{<input type="hidden" id="chr$i" value="$chr">} : qq{<input type="hidden" id="chr$i">};
	$dsinfo .=get_dataset_info(dsid=>$dsid, chr=>$chr) if $dsid;
	push @seq_nums, {
			  SEQ_NUM=>$i,
			 };
	my %opts = (
		    SEQ_NUM=>$i,
		    REV_YES=>$revy,
		    REV_NO=>$revn,
		    REF_YES=>$refy,
		    REF_NO=>$refn,
		    DRUP=>$drup,
		    DRDOWN=>$drdown,
		    DRACCN=>$draccn,
		    DSID=>$dsid,
		    GBACCN=>$gbaccn,
		    GBSTART=>$gbstart,
		    GBLENGTH=>$gblength,
		    POS=>$pos,
		    DSINFO=>$dsinfo,
		    EXON_MASK_ON=>$maskexonon,
		    EXON_MASK_OFF=>$maskexonoff,

		   );
#	print STDERR "pos $i: $pos\n";
	$opts{COGEPOS} = qq{<option value="cogepos$i" selected="selected">CoGe Database Position</option>} if $pos;
	push @seq_sub, {%opts}
	  
      }
    #page preferences
    my $pad_gs = $form->param("pad_gs") if $form->param("pad_gs");
    $pad_gs = 0 unless $pad_gs;
    my $prog = get_opt(params=>$prefs, form=>$form, param=>'prog');
    $prog = "blastz" unless $prog;
    my $image_width = get_opt(params=>$prefs, form=>$form, param=>'iw');
    $image_width = 1000 unless $image_width;
    my $image_height = get_opt(params=>$prefs, form=>$form, param=>'ih');
    $image_height = 100 unless $image_height;
    my $feature_height = get_opt(params=>$prefs, form=>$form, param=>'fh');
    $feature_height = 10 unless $feature_height;
    my $padding  = get_opt(params=>$prefs, form=>$form, param=>'padding');
    $padding = 2 unless defined $padding;
    my $gc_color = get_opt(params=>$prefs, form=>$form, param=>'gc');
    $gc_color = 0 unless $gc_color;
    my $nt_color = get_opt(params=>$prefs, form=>$form, param=>'nt');
    $nt_color = 1 unless $nt_color;
    my $auto_adjust_feats = get_opt(params=>$prefs, form=>$form, param=>'overlap');
    $auto_adjust_feats = 0 unless defined $auto_adjust_feats;
    my $hiqual = get_opt(params=>$prefs, form=>$form, param=>'hiqual');
    $hiqual = 0 unless $hiqual;
    my $color_hsp = get_opt(params=>$prefs, form=>$form, param=>'colorhsp');
    $color_hsp = 0 unless $color_hsp;
    my $hsp_label = get_opt(params=>$prefs, form=>$form, param=>'hsplabel');
    $hsp_label = undef unless defined $hsp_label;
    my $hsp_limit = get_opt(params=>$prefs, form=>$form, param=>'hsplim');
    $hsp_limit = 0 unless $hsp_limit;
    my $hsp_limit_num = get_opt(params=>$prefs, form=>$form, param=>'hsplimnum');
    $hsp_limit_num = 20 unless defined $hsp_limit_num;
    my $draw_model = get_opt(params=>$prefs, form=>$form, param=>'draw_model');
    $draw_model = "full" unless $draw_model;
    my $hsp_overlap_limit = get_opt(params=>$prefs, form=>$form, param=>'hsp_overlap_limit');
    $hsp_overlap_limit = 0 unless $hsp_overlap_limit;

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
    $template->param(PAD_GS=>$pad_gs);
    $template->param(IMAGE_WIDTH=>$image_width);
    $template->param(IMAGE_HEIGHT=>$image_height);
    $template->param(FEAT_HEIGHT=>$feature_height);
    $template->param(PADDING=>$padding);
    $template->param(HSP_OVERLAP_LIMIT=>$hsp_overlap_limit);
    if ($gc_color) {$template->param(GC_COLOR_YES=>"checked");}
    else {$template->param(GC_COLOR_NO=>"checked");}
    if ($nt_color) {$template->param(NT_COLOR_YES=>"checked");}
    else {$template->param(NT_COLOR_NO=>"checked");}
    if ($auto_adjust_feats) {$template->param(OVERLAP_YES=>"checked");}
    else {$template->param(OVERLAP_NO=>"checked");}
    if ($hiqual) {$template->param(HIQUAL_YES=>"checked");}
    else {$template->param(HIQUAL_NO=>"checked");}
    if ($color_hsp) {$template->param(COLOR_HSP_YES=>"checked");}
    else {$template->param(COLOR_HSP_NO=>"checked");}
    if ($hsp_label eq "staggered") {$template->param(HSP_LABELS_STAG=>"selected");}
    elsif ($hsp_label eq "linear") {$template->param(HSP_LABELS_LIN=>"selected");}
    else {$template->param(HSP_LABELS_NO=>"selected");}
    if ($hsp_limit) {$template->param(HSP_LIMIT_YES=>"checked");}
    else {$template->param(HSP_LIMIT_NO=>"checked");}
    $template->param(HSP_LIMIT_NUM=>$hsp_limit_num);
    if ($draw_model eq "full") {$template->param(DRAW_MODEL_FULL=>"selected");}
    elsif ($draw_model eq "gene") {$template->param(DRAW_MODEL_GENE=>"selected");}
    elsif ($draw_model eq "mRNA") {$template->param(DRAW_MODEL_mRNA=>"selected");}
    elsif ($draw_model eq "CDS") {$template->param(DRAW_MODEL_CDS=>"selected");}
    elsif ($draw_model eq "RNA") {$template->param(DRAW_MODEL_RNA=>"selected");}
    else {$template->param(DRAW_MODEL_NO=>"selected");}

    my $box = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
    #generate sequence submission selector
    $template->param(SEQ_SELECT=>1);
    $template->param(SEQ_SELECT_LOOP=>\@seq_sub);
    my $seq_submission = $template->output;
    $template->param(SEQ_SELECT=>0);

    #generate the hsp color option
    my $hsp_colors = gen_hsp_colors(num_seqs=>$num_seqs, prefs=>$prefs);
    my $html;
    my $spike_len = spike_filter_select();

    $template->param(SPIKE_LEN=>$spike_len);
    $template->param(SEQ_RETRIEVAL=>1);
    $template->param(NUM_SEQS=>$num_seqs);
    $message .= "<BR/>" if $message;
    $template->param(MESSAGE=>$message);
    $template->param(SEQ_SUB=>$seq_submission);
    $template->param(HSP_COLOR=>$hsp_colors);
    $template->param(GO_RUN=>gen_go_run($num_seqs));
    $template->param(AUTOSEARCH=>$autosearch_string);
    my $gobe_version = `svnversion /opt/apache/CoGe/gobe/flash`;
    $gobe_version =~ s/\n//g;;
    $template->param(GOBE_VERSION=>$gobe_version);
    $box->param(BOX_NAME=>"Options:");
    $template->param(OPTIONS=>1);
    $template->param(ALIGNMENT_PROGRAMS=>algorithm_list($prog));
    $template->param(SAVE_SETTINGS=>gen_save_settings($num_seqs)) unless !$USER || $USER =~ /public/i;
    $box->param(BODY=>$template->output);
    $html .= $box->output;
    return $html;
  }

sub run
  {
    my %opts = @_;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my $spike_len = $opts{spike};
    my $iw = $opts{iw};
    my $ih = $opts{ih};
    my $feat_h = $opts{fh};
    my $show_gc = $opts{gc};
    my $show_nt = $opts{nt};
    my $color_hsp = $opts{colorhsp};
    my $hsp_label = $opts{hsplabel};
    my $draw_model = $opts{draw_model};
    my $hsp_overlap_limit = $opts{hsp_overlap_limit};
    my $overlap_adjustment = $opts{overlap};
    my $hiqual = $opts{hiqual};
    my $hsp_limit = $opts{hsplim};
    my $hsp_limit_num = $opts{hsplimnum};
    my $show_hsps_with_stop_codon = $opts{showallhsps};
    my $padding = $opts{padding};
    my ($analysis_program, $param_string, $parser_opts) = get_algorithm_options(%opts);
    my $basefilename = $opts{basefile};
    my $show_spike = $opts{show_spike} || 0;
    my $pad_gs = $opts{pad_gs} || 0;
    my $color_overlapped_features = $opts{color_overlapped_features};
    my $hsp_overlap_length = $opts{hsp_overlap_length};
    my $message;
    $cogeweb = initialize_basefile(basename=>$basefilename, prog=>"GEvo");
    my @hsp_colors;
    for (my $i = 1; $i <= num_colors($num_seqs); $i++)
      {
	my $r = $opts{"r$i"};
	my $g = $opts{"g$i"};
	my $b = $opts{"b$i"};
	my @tmp;
	foreach my $color ($r, $g, $b)
	  {
	    $color = 0 unless $color && $color =~ /^\d+$/;
	    $color = 0 if $color < 0;
	    $color = 255 if $color > 255;
	    push @tmp, $color;
	  }
	push @hsp_colors,\@tmp;
      }



    my $stagger_label = $hsp_label && $hsp_label =~ /staggered/i ? 1 : 0;
    my $feature_labels = !$hsp_label ? 0 : 1;
    my $form = $FORM;
    my $gevo_link = $form->url."?prog=$analysis_program";
    $gevo_link .= ";spike_len=$spike_len";
    
    my @sets;
    my $html;
    my $t1 = new Benchmark;
    my $spike_seq;
    my $seqcount = 1;
    for (my $i = 1; $i <= $num_seqs; $i++)
      {
	my $skip_seq =$opts{"skip_seq$i"};
	next if $skip_seq;
	my $accn = $opts{"draccn$i"};
	my $featid = $opts{"featid$i"};
	my $feat = $coge->resultset('Feature')->find($featid) if $featid;
	my $dsid;
	$dsid = $feat->dataset_id if $feat;
	$dsid = $opts{"dsid$i"} unless $dsid;
	$dsid = $opts{"posdsid$i"} unless $dsid;
	my $chr;
	$chr = $feat->chromosome if $feat;
	$chr = $opts{"chr$i"} unless $chr;

	my $drup = $opts{"drup$i"};
	my $drdown = $opts{"drdown$i"};
	$drup += $pad_gs;
	$drdown += $pad_gs;
	my $pos = $opts{"pos$i"};
	my $gbaccn = $opts{"gbaccn$i"};
	my $gbstart = $opts{"gbstart$i"};
	$gbstart = 1 unless defined $gbstart;
	my $gblength = $opts{"gblength$i"};

	my $dirseq = $opts{"dirseq$i"};
	my $dirstart = $opts{"dirstart$i"};
	my $dirlength = $opts{"dirlength$i"};

	my $rev = $opts{"rev$i"};
	my $mask_cds_flag = $opts{"maskcds$i"};
	my $mask_ncs_flag = $opts{"maskncs$i"};
	my ($up, $down, $seq);
	my ($file, $obj);
	my $reference_seq =$opts{"ref_seq$i"};
	my $repeat_mask =$opts{"repmask$i"};
	next unless $accn || $featid || $gbaccn || $dirseq|| $pos;
	$gevo_link .= ";accn$seqcount=".CGI::escape($accn) if $accn;
	$gevo_link .= ";x$seqcount=".CGI::escape($pos) if $pos;
	$gevo_link .= ";fid$seqcount=".CGI::escape($featid) if $featid;
	$gevo_link .= ";dsid$seqcount=".CGI::escape($dsid) if $dsid;
	$gevo_link .= ";chr$seqcount=".CGI::escape($chr) if $chr;
	$gevo_link .= ";dr$seqcount"."up=$drup" if defined $drup;
	$gevo_link .= ";dr$seqcount"."down=$drdown" if defined $drdown;
	$gevo_link .= ";gbaccn$seqcount=".CGI::escape($gbaccn) if $gbaccn;
	$gevo_link .= ";gbstart$seqcount=$gbstart" if $gbstart;
	$gevo_link .= ";gblength$seqcount=$gblength" if $gblength;
	$gevo_link .= ";rev$seqcount=1" if $rev;
	$gevo_link .= ";nref$seqcount=1" unless $reference_seq;
	$seqcount++;
#	print STDERR "pos: $pos, dsid: $dsid, chr: $chr\n";
	if ($featid || $pos)
	  {
	    $obj = get_obj_from_genome_db( accn=>$accn, featid=>$featid, pos=>$pos, dsid=>$dsid, rev=>$rev, up=>$drup, down=>$drdown, chr=>$chr );
	    if ($obj)
	      {
		($file, $spike_seq, $seq) = 
		  generate_seq_file(obj=>$obj,
				    mask_cds=>$mask_cds_flag,
				    mask_ncs=>$mask_ncs_flag,
				    spike_len=>$spike_len,
				    seq_num=>$i,
				    repeat_mask=>$repeat_mask,
				   );
		$up = $drup;
		$down = $drdown;
	      }
	    else
	      {
		$message .=  "Unable to generate sequence for sequence $i.  (Skipped.)";
	      }
	  }
 	elsif ($dirseq )
 	  {
 	    ($obj) = generate_obj_from_seq($dirseq, $i, $rev);
	    $dirlength = length($obj->sequence)-$dirstart+1 unless $dirlength;
 	    if ($obj)
	      {
		#add an anchor
		$obj->add_feature(
			      type=>"direct sequence submission",
			      location=>(1-$dirstart+1)."..".(length($obj->sequence)-$dirstart+1),
			      strand=>1,
			      qualifiers=>{
					   type=>"anchor",
					   names=>[$obj->accn],
					  }
			     );

		($file, $spike_seq, $seq) = 
		  generate_seq_file (
				     obj=>$obj,
				     mask_cds=>$mask_cds_flag,
				     mask_ncs=>$mask_ncs_flag,
				     startpos=>$dirstart,
				     length=>$dirlength,
				     spike_len=>$spike_len, 
				     seq_num=>$i,
				     repeat_mask=>$repeat_mask,
				    );
		$obj->start($dirstart);
		$obj->stop($dirstart+length($seq)-1);
		$obj->chromosome(1);
#		$obj->dataset("NA");
		$up = $dirstart;
		$down = $dirlength;
	      }
	    else
	      {
		$message .= "Problem with direct sequence submission\n";
	      }
	  }
	elsif ($gbaccn )
	  {
	    $obj = new CoGe::Accessory::GenBank;
	    $obj->add_gene_models(1); #may want to make a user selectable option
 	    my ($res, $error) = $obj->get_genbank_from_ncbi($gbaccn, $rev);
	    $message .= $error."\n" unless $res;
	    
	    if ($obj->accn)
	      {
		#add an anchor
		$obj->add_feature(
				  type=>"genbank entry",
				  location=>(1-$gbstart+1)."..".(length($obj->sequence)-$dirstart+1),
				  strand=>1,
				  qualifiers=>{
					       type=>"anchor",
					       names=>[$gbaccn],
					      }
				 );

		($file, $spike_seq, $seq) = 
		  generate_seq_file (
				     obj=>$obj,
				     mask_cds=>$mask_cds_flag,
				     mask_ncs=>$mask_ncs_flag,
				     startpos=>$gbstart,
				     length=>$gblength,
				     spike_len=>$spike_len,
				     seq_num=>$i,
				     repeat_mask=>$repeat_mask,
				    );
		$obj->start($gbstart);
		$obj->stop($gbstart + length($seq)-1);
		$obj->chromosome("?") unless $obj->chromosome;
		$up = $gbstart;
		$down = $gblength;
	      }
	    else
	      {
		$message .= "No GenBank entry found for $gbaccn\n";
	      }
 	  }
	next unless $obj;
	unless ($show_spike)
	  {
	    $seq =~ s/N*$spike_seq$//;
	    $seq =~ s/N*$spike_seq$//;
	  }
	$obj->sequence($seq);

	if ($obj && $obj->sequence)
	  {
	    #need to check for duplicate accession names -- sometimes happens and major pain in the ass for other parts of the code
	    my $accn = $obj->accn;
	    my $count = 0;
	    foreach my $accn2 (map {$_->{obj}->accn()} @sets)
	      {
		$accn2 =~ s/\*\*\d+\*\*$//;
		$count++ if $accn eq $accn2;
	      }
	    $accn .= "**".($count+1)."**" if $count;
	    $obj->accn($accn);
	    push @sets, {
			 obj=>$obj,
			 file=>$file,
#			 file_begin=>$file_begin,
#			 file_end=>$file_end,
			 accn=>$accn,
			 rev=>$rev,
			 up=>$up,
			 down=>$down,
			 spike_seq=>$spike_seq,
			 reference_seq=>$reference_seq,
			 seq_num=>$i,
			};
	  }
	else
	  {
#	    push @sets, {seq_num=>$i};
	  }
      }
    $seqcount--;

    $gevo_link .= ";num_seqs=".$seqcount;
    $gevo_link .= ";hsp_overlap_limit=".$hsp_overlap_limit if defined $hsp_overlap_limit;
    unless (@sets >1)
      {
	$message .= "Problem retrieving information.  Please check submissions.\n";
	return '', '', '', '',0,$message;
      }
    my $t2 = new Benchmark;
    # set up output page
    
    # run bl2seq
    
    my $analysis_reports;

    if ($analysis_program eq "blastz")
      {
	$analysis_reports = run_blastz(sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts);
      }
    elsif ($analysis_program eq "LAGAN")
      {
	$analysis_reports = run_lagan (sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts);
      }
    elsif ($analysis_program eq "CHAOS")
      {
	$analysis_reports = run_chaos (sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts);
      }
    elsif ($analysis_program eq "DiAlign_2")
      {
	($analysis_reports, $analysis_program,$message) = run_dialign (sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts);
      }
    else
      {
	$analysis_reports = run_bl2seq(sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts, blast_program=>$analysis_program, spike_seq=>$spike_seq);
      }
    $analysis_reports = [] unless ref($analysis_reports) =~ /ARRAY/i;
    write_log($message, $cogeweb->logfile) if $message;
   #sets => array or data for blast
   #blast_reports => array of arrays (report, accn1, accn2, parsed data)
    my $t3 = new Benchmark;
    my $count = 1;
    my $frame_height;
    initialize_sqlite();
    my @gfxs;
    foreach my $item (@sets)
      {
	my $obj = $item->{obj};
	next unless $obj->sequence;
	my $rev = $item->{rev};
	my $spike_seq = $item->{spike_seq};
	if ($obj)
	  {
	    write_log("generating image ($count/".scalar @sets.")for ".$obj->accn, $cogeweb->logfile);
	    $count++;
#	    my ($image, $map, $mapname, $gfx, $eval_cutoff, $feat_count, $overlap_count) = 
	    my ($gfx, $stats) = 
	      generate_image(
			     set=>$item,
			     start=>1,
			     stop=>length($obj->sequence),
			     data=>$analysis_reports,
			     iw=>$iw,
			     ih=>$ih,
			     fh=>$feat_h,
			     show_gc=>$show_gc,
			     show_nt=>$show_nt,
			     stagger_label=>$stagger_label,
			     overlap_adjustment=>$overlap_adjustment,
			     feature_labels=>$feature_labels,
			     hsp_limit=>$hsp_limit,
			     hsp_limit_num=>$hsp_limit_num,
			     color_hsp=>$color_hsp,
			     hsp_colors=>\@hsp_colors,
			     spike_sequence=>$spike_seq,
			     show_hsps_with_stop_codon => $show_hsps_with_stop_codon,
			     hiqual=>$hiqual,
			     padding=>$padding,
#			     seq_num=>$item->{seq_num},
			     reverse_image=>$rev,
			     draw_model=>$draw_model,
			     hsp_overlap_limit=>$hsp_overlap_limit,
			     color_overlapped_features=>$color_overlapped_features,
			     hsp_overlap_length=>$hsp_overlap_length,
			    );
#	    $item->{feat_count} = $feat_count;
#	    $item->{overlap_count} = $overlap_count;
	    my $filename = $cogeweb->basefile."_".$item->{seq_num}.".png";
	    $filename = check_filename_taint($filename);
	    my $image = basename($filename);
	    $item->{image}=$image;
	    $item->{stats}=$stats;
	    $gfx->generate_png(file=>$filename);
	    generate_image_db(set=>$item, gfx=>$gfx);
	    $frame_height += $gfx->ih + $gfx->ih*.1;
	  }
      }
    my $t3_5 = new Benchmark;
    image_db_create_hsp_pairs();

    my $t4 = new Benchmark;
    $html .= qq{<DIV id=flash_viewer></DIV>};
    $html .= qq{<table>};
    $html .= qq{<tr valign=top><td class = small>Alignment reports};
    my $stats_file = $cogeweb->basefile."_stats.txt";
    $stats_file = check_filename_taint($stats_file);

    if ($analysis_reports && @$analysis_reports)
      {
	foreach my $item (@$analysis_reports)
	  {
	    my $report = $item->[0];
	    my $accn1 = $item->[1];
	    my $accn2 = $item->[2];
	    my $basereportname = basename( $report );
	    $basereportname = $TEMPURL . "/$basereportname\n";
	    $html .= "<div><font class=small><A HREF=\"$basereportname\" target=_new>View alignment output for $accn1 versus $accn2</A></font></DIV>\n";
	  }
      }
    else
      {
	$html .= "<div class=small>No alignment reports were generated</DIV>\n";
      }
    $html .= qq{<td class = small>Fasta files};
    foreach my $item (@sets)
      {
	next unless $item->{file};
	my $basename = $TEMPURL."/".basename ($item->{file});
	print STDERR "basename is undefined: $basename\n" if $basename =~ /defined/i;
	my $accn = $item->{accn};
	$html .= "<div><font class=small><A HREF=\"$basename\" target=_new>Fasta file for $accn</A></font></DIV>\n";
      }
    $html .= qq{<td class = small><a href = "http://baboon.math.berkeley.edu/mavid/gaf.html">GAF</a> annotation files};
    write_log("Features overlapped by HSPs:", $stats_file);
    foreach my $item (@sets)
      {
	my $anno_file = generate_annotation(%$item);
	next unless $anno_file;
	my $basename = $TEMPURL."/".basename ($anno_file);
	my $accn = $item->{accn};
	$html .= "<div><font class=small><A HREF=\"$basename\" target=_new>Annotation file for $accn</A></font></DIV>\n";
	my $stats = $item->{stats};
	foreach my $accn1 (sort keys %{$stats->{data}})
	  {
	    foreach my $accn2 (sort keys %{$stats->{data}{$accn1}})
	      {
		foreach my $type (sort keys %{$stats->{data}{$accn1}{$accn2}})
		  {
		    next if $type eq "hsp";
		    write_log(join ("\t", $accn1, $accn2, $type, $stats->{data}{$accn1}{$accn2}{$type}{overlapped_hsps}), $stats_file);
		  }
	      }
	  }
      }
    $html .= qq{<td class = small>SQLite db};
    my $dbname = $TEMPURL."/".basename($cogeweb->sqlitefile);
    
    $html .= "<div class=small><A HREF=\"$dbname\" target=_new>SQLite DB file</A></DIV>\n";
    $html .= qq{<td class = small>Log File};
    my $logfile = $TEMPURL."/".basename($cogeweb->logfile);
    $html .= "<div class=small><A HREF=\"$logfile\" target=_new>Log</A></DIV>\n";
    my $tiny = get("http://tinyurl.com/create.php?url=$gevo_link");
    ($tiny) = $tiny =~ /<b>(http:\/\/tinyurl.com\/\w+)<\/b>/;
    $html .= qq{<td class = small>GEvo Link<div class=small><a href=$tiny target=_new>$tiny<br>(See log file for full link)</a></div>};
    $html .=qq{<a href="GEvo_direct.pl?name=$basefilename" target=_new>Results only</a>};
    $html .= qq{<td class = small>Overlap Feature Stats:<table>};
    $html .= qq{<tr class=small>}.join ("<th align=left>", qw(Dataset ));
    write_log("\nAverage percent identity for all HSPs (normalized to HSP length):", $stats_file);
    foreach my $item (@sets)
      {
	$html .= "<tr class=small>";
	$html .= "<td>".$item->{obj}->accn."<td nowrap>".commify(length $item->{obj}->sequence)."bp"."<td nowrap>".$item->{stats}{overlap_count}."/".$item->{stats}{feat_count};
	$html .= "<td nowrap>".sprintf("%.2f", $item->{stats}{overlap_count}/$item->{stats}{feat_count}*100)."%" if $item->{stats}{feat_count};
	my $stats = $item->{stats};
	foreach my $accn1 (sort keys %{$stats->{data}})
	  {
	    foreach my $accn2 (sort keys %{$stats->{data}{$accn1}})
	      {
		my $len= 0;
		map {$len+=$_->[1]} @{$stats->{data}{$accn1}{$accn2}{hsp}};
		my $pid = 0;
		map {$pid+=$_->[1]*$_->[0]} @{$stats->{data}{$accn1}{$accn2}{hsp}};
		write_log("$accn1\t$accn2\t".sprintf("%.2f",$pid/$len), $stats_file);
	      }
	  }
      }
    my $stats_url = $TEMPURL."/".basename($stats_file);
    $html .= qq{<tr class=small><td><a href=$stats_url target=_new>Stats file</a>};
    $html .= qq{</table>};
    $html .= qq{</table>};

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
    my $results_name = "Results: $analysis_program";
    $results_name .= qq{ <span class="small">(spike sequence filter length: $spike_len)</span>} if $spike_len;
    $template->param(BOX_NAME=>$results_name);
    $template->param(BODY=>$html);
    my $outhtml = $template->output;
    my $t5 = new Benchmark;
    my $db_time = timestr(timediff($t2,$t1));
    my $blast_time = timestr(timediff($t3,$t2));
    my $image_time = timestr(timediff($t3_5,$t3));
    my $db_hsp_time = timestr(timediff($t4,$t3_5));
    my $html_time = timestr(timediff($t5,$t4));
    my $bench =  qq{
GEvo Benchmark: $DATE
Time to get sequence                              : $db_time
Time to run $analysis_program                     : $blast_time
Time to generate images, maps, and sqlite database: $image_time
Time to find and update sqlite database for HSPs  : $db_hsp_time
Time to process html                              : $html_time
};
    print STDERR $bench if $BENCHMARK;
    write_log($bench, $cogeweb->logfile);
    write_log("Finished!", $cogeweb->logfile);
    write_log("\nGEvo link: $gevo_link\n", $cogeweb->logfile);
    write_log("Tiny url: $tiny", $cogeweb->logfile);
    $count--;
    return $outhtml, $iw+400, $frame_height, $cogeweb->basefilename,$count,$message;

}


sub generate_image
  {
    my %opts = @_;
    my $set = $opts{set};
    my $gbobj = $set->{obj};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $data = $opts{data};
    my $masked_exons = $opts{mask};
    my $masked_ncs = $opts{mask_ncs};
    my $spike_seq = $opts{spike_sequence};
    my $iw = $opts{iw} || 1600;
    my $ih = $opts{ih} || 200;
    my $fh = $opts{fh} || 25;
    my $show_gc = $opts{show_gc};
    my $show_nt = $opts{show_nt};
    my $reverse_image = $opts{reverse_image};
    my $stagger_label = $opts{stagger_label};
    my $overlap_adjustment = $opts{overlap_adjustment};
    my $feature_labels = $opts{feature_labels};
    my $hsp_limit = $opts{hsp_limit};
    my $color_hsp = $opts{color_hsp};
    my $hsp_limit_num = $opts{hsp_limit_num};
    my $eval_cutoff = $opts{eval_cutoff};
    my $hsp_colors = $opts{hsp_colors};
    my $show_hsps_with_stop_codon = $opts{show_hsps_with_stop_codon};
    my $hiqual = $opts{hiqual};
    my $padding = $opts{padding} || 5;
#    my $seq_num = $opts{seq_num};
    my $draw_model = $opts{draw_model};
    my $hsp_overlap_limit = $opts{hsp_overlap_limit};
    my $color_overlapped_features = $opts{color_overlapped_features};
    my $hsp_overlap_length = $opts{hsp_overlap_length};
    my $graphic = new CoGe::Graphics;
    my $gfx = new CoGe::Graphics::Chromosome;
#    print STDERR "$start"."::"."$stop"." -- ".length($gbobj->sequence)."\n";
    $graphic->initialize_c (
			    c=>$gfx,
			    iw=>$iw,
			    start=> $start,
			    stop => $stop,
			    draw_chr=>1,
			    draw_ruler=>1,
			    draw_chr_end=>0,
			    chr_start_height=>$ih,
			    chr_mag_height=>5,
			    feature_start_height=>$fh,
			    mag=>0,
			    mag_off=>1,
			    chr_length => length($gbobj->sequence),
			    fill_labels=>1,
			    forcefit=>1,
			    minor_tick_labels=>1,
			    feature_labels=>$feature_labels,
			    draw_hi_qual=>$hiqual,
			    padding=>$padding,
			   );
 #   print STDERR $gfx->_region_start ."--".$gfx->_region_stop."\n";
    $gfx->overlap_adjustment(1);
    $gfx->top_padding(15);
    $gfx->skip_duplicate_features(1);
    $gfx->DEBUG(0);
    $gfx->major_tick_labels(0);
    $graphic->process_nucleotides(c=>$gfx, seq=>$gbobj->sequence, layers=>{gc=>$show_gc, nt=>$show_nt});
    my $stats = process_hsps(
			  c=>$gfx, 
			  data=>$data, 
			  accn=>$gbobj->accn, 
			  rev=>$reverse_image, 
			  seq_length=> length($gbobj->sequence), 
			  stagger_label=>$stagger_label, 
			  hsp_limit=>$hsp_limit,
			  hsp_limit_num=>$hsp_limit_num,
			  gbobj=>$gbobj,
			  spike_seq=>$spike_seq,
			  color_hsp=>$color_hsp,
			  colors=>$hsp_colors,
			  show_hsps_with_stop_codon=>$show_hsps_with_stop_codon,
			  hsp_overlap_limit=>$hsp_overlap_limit,
			  hsp_overlap_length=>$hsp_overlap_length,
			 );
    my ($feat_count, $overlap_count) = process_features(c=>$gfx, obj=>$gbobj, start=>$start, stop=>$stop, overlap_adjustment=>$overlap_adjustment, draw_model=>$draw_model, color_overlapped_features=>$color_overlapped_features);
    $stats->{feat_count} = $feat_count;
    $stats->{overlap_count} = $overlap_count;
    return ($gfx, $stats);
  }

sub initialize_sqlite
  {
    my %opts = @_;
    my $dbfile = $cogeweb->sqlitefile;
    return if -r $dbfile;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
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
color varchar
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
iname varchar,
title varchar,
px_width integer,
bpmin integer,
bpmax integer,
dsid integer,
chromosome varchar
)
};
#TODO: make sure to populate the bpmin and pbmax and image width!!!!!!!!!!!!!!!!!
    $dbh->do($create);
    $dbh->do('CREATE INDEX id ON image_info (id)');
    system "chmod +rw $dbfile";
  }

sub generate_image_db
  {
    my %args = @_;
    my $gfx = $args{gfx};
    my $set = $args{set};
    next unless $gfx;
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");
    my $image = $set->{image};
    my $accn = $set->{accn};
    my $title;
    $title = $accn;
    $title .= " ".$set->{obj}->organism() if $set->{obj}->organism();
    $title .= "(chr: ".$set->{obj}->chromosome." ". $set->{obj}->start."-".$set->{obj}->stop.")" if defined $set->{up};
    $title .= qq! Reverse Complement! if $set->{rev};
    my $width = $gfx->image_width;
    my $dsid = $set->{obj}->dataset;
    my ($chr) =  $set->{obj}->chromosome =~ /(\d+)/;
    $chr = "NULL" unless $chr;
    $dsid = "NULL" unless $dsid;
    my $image_start = $set->{obj}->start;
    my $image_stop = $set->{obj}->stop;
    my $image_id = $set->{seq_num};
    my $statement = qq{
INSERT INTO image_info (id, iname, title, px_width,dsid, chromosome, bpmin, bpmax) values ($image_id, "$image", "$title", $width, "$dsid", $chr, $image_start, $image_stop)
};
    print STDERR $statement unless $dbh->do($statement);
    foreach my $feat ($gfx->get_feats)
      {
	if ($feat->fill)
	  {
	    next unless $feat->type eq "anchor";
	  }
#	print STDERR Dumper $feat;
#	print STDERR Dumper $feat if $feat->{anchor};
	next unless $feat->image_coordinates;
	my $type = $feat->type;
	my $pair_id = "-99";
	my $coords = $feat->image_coordinates;
	$coords =~ s/\s//g;
	next if $feat->type eq "unknown";
	my $name = $feat->type =~ /HSP/i ? $feat->alt : $feat->label;
	$name .= "_".$feat->type;# unless $name;

	my $color = "NULL";
	if ($feat->type =~ /HSP/)
	  {
	    $color ="#";
	    foreach my $c (@{$feat->color})
	      {
		$color .= sprintf("%X",$c);
	      }
	  }
	#generate link
	my $link = $feat->link;
	$link = " " unless $link;
	$link =~ s/'//g if $link;
	#generate image track
	my $image_track = $feat->track;
	$image_track = "-".$image_track if $feat->strand =~ /-/;
	
	my ($xmin, $ymin, $xmax, $ymax) = split /,/, $coords;
	$xmin++;
	$xmax++;
	my $anno = $feat->description;
	#	    $anno =~ s/'|"//g;
	$anno =~ s/'//g if $anno;
	$anno =~ s/<br\/?>/&#10;/ig if $anno;
	$anno =~ s/\n/&#10;/g if $anno;
	$anno =~ s/[\[\]\(\)]/ /g if $anno;
	$anno = " " unless $anno;
	my $start = $feat->start;
	my $stop = $feat->stop;
	my $length_nt = $stop-$start+1;
	my $length_pix = $xmax-$xmin;
	next if !$feat->{anchor} && ($length_nt == 0 );
	$type = "anchor" if $feat->{anchor};
	my $bpmin = $set->{obj}->start+$start-1;
	my $bpmax = $set->{obj}->start+$stop-1;
	$statement = qq{
INSERT INTO image_data (name, type, xmin, xmax, ymin, ymax, bpmin,bpmax,image_id, image_track,pair_id, link, annotation, color) values ("$name", "$type", $xmin, $xmax, $ymin, $ymax, $bpmin, $bpmax,$image_id, "$image_track",$pair_id, '$link', '$anno', '$color')
};
	print STDERR $statement unless $dbh->do($statement);
      }
  }

sub image_db_create_hsp_pairs
  {
    my $dbh = DBI->connect("dbi:SQLite:dbname=".$cogeweb->sqlitefile,"","");
    my $query = qq{select id,name from image_data where type = "HSP"};
    my %ids;
    foreach my $item (@{$dbh->selectall_arrayref($query)})
      {
	push @{$ids{$item->[1]}},$item->[0];
      }
    foreach my $name (keys %ids)
      {
	if (@{$ids{$name}} != 2)
	  {
	    print STDERR "Problem with $name.  Retrieved", Dumper $ids{$name} if $DEBUG;
	    next;
	  }
	my ($id1, $id2) = @{$ids{$name}};
	my $statement = "update image_data set pair_id = $id1 where id = $id2";
	$dbh->do($statement);
	$statement = "update image_data set pair_id = $id2 where id = $id1";
	$dbh->do($statement);
      }
  }


sub process_features
  {
    #process features
    my %opts = @_;
    my $c = $opts{c};
    my $obj=$opts{obj};
    my $start=$opts{start};
    my $stop = $opts{stop};
    my $overlap = $opts{overlap_adjustment};
    my $draw_model = $opts{draw_model};
    my $color_overlapped_features = $opts{color_overlapped_features};
#    return unless $draw_model;
    my $accn = $obj->accn;
    my $track = 1;
    my $feat_count = 0;
    my $overlap_count = 0;
    my %feat_count;
    my @opts = ($start, $stop) if $start && $stop;
    unless (ref $obj)
      {
	warn "Possible problem with the object in process_features.  Returning";
	return 0;
      }
    
    foreach my $feat($obj->get_features())
      {
        my $f;
	my $type = $feat->type;
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} @{$feat->qualifiers->{names}} if ref ($feat->qualifiers) =~ /hash/i;
	my $anchor = $feat if ref ($feat->qualifiers) =~ /hash/i && $feat->qualifiers->{type} && $feat->qualifiers->{type} eq "anchor";
        if ($type =~ /pseudogene/i)
          {
	    next unless $draw_model eq "full";
	    $f = CoGe::Graphics::Feature::Gene->new();
	    $f->color([255,33,0,50]);
	    $f->order($track);
	    $f->overlay(1);
	    $f->mag(0.5);
          }
        elsif ($type =~ /Gene/i)
          {
	    next unless $draw_model eq "full" || $draw_model eq "gene" || $anchor;
	    $f = CoGe::Graphics::Feature::Gene->new();
	    #$f->color([255,0,0,50]);
        $f->color([219,219,219,50]);
	    $f->order($track);
	    $f->overlay(1);
	    $f->mag(0.5);
          }
        elsif ($type =~ /CDS/i)
          {
	    next unless $draw_model eq "full" || $draw_model eq "CDS";
	    $f = CoGe::Graphics::Feature::Gene->new();
	    $f->color([0,255,0, 50]);
	    $f->order($track);
	    $f->overlay(3);
	    if ($accn)
	      {
		foreach my $name (@{$feat->qualifiers->{names}})
		  {
            my $cleaned_name = $name;
            $cleaned_name =~ s/[\(\)]//g;
		    if ($accn =~ /^$cleaned_name\(?\d*\)?$/i)
		      {
			$f->color([255,255,0]) ;
			$f->label($name);
		      }
		  }
	      }
	    
		
          }
        elsif ($type =~ /mrna/i)
          {
	    next unless $draw_model eq "full" || $draw_model eq "mRNA";
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,0,255, 50]);
        	$f->order($track);
		$f->overlay(2);
		$f->mag(0.75);
          }
        elsif ($type =~ /rna/i)
          {
	    next unless $draw_model eq "full" || $draw_model eq "RNA";
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([200,200,200, 50]);
        	$f->order($track);
		$f->overlay(2);
		if ($accn)
		  {
		    foreach my $name (@{$feat->qualifiers->{names}})
		      {
			if ($accn =~ /^$name\(?\d*\)?$/i)
			  {
			    $f->color([255,255,0]) ;
			    $f->label($name);
			  }
		      }
		  }

          }
	elsif ($type =~ /anchor/)
	  {
	    $f = CoGe::Graphics::Feature::NucTide->new({type=>'anchor',start=>$feat->blocks->[0][0], strand=>1});
	    $f->color([0,0,255]);
	    $f->type($type);
	    $f->description($feat->annotation);
	    $c->add_feature($f);
	    $f = CoGe::Graphics::Feature::NucTide->new({type=>'anchor',start=>$feat->blocks->[0][0], strand=>-1});
	    $f->color([0,0,255]);
	    $f->type($type);
	    $f->description($feat->annotation);
	    $c->add_feature($f);
	    next;
	  }
	elsif ($type =~ /cns/i)
	  {
	    $f = CoGe::Graphics::Feature::HSP->new({start=>$feat->blocks->[0][0], stop=>$feat->blocks->[0][1]});
	    $f->color([255,51,0]);
	    $f->order($track);
	    $f->overlay(4);
	    $f->type($type);
	    $f->force_draw(1);
	    $c->add_feature($f);
	    next;
	  }
	#need to create an anchor if this feature is an anchor, but not to be drawn
	if ($anchor && !$f)
	  {
	    $f = CoGe::Graphics::Feature->new({type=>'anchor',start=>$feat->blocks->[0][0], stop=>$feat->blocks->[0][1], strand=>1});
	    $f->{anchor}=1;
	    $f->order(0);
	    $f->skip_overlap_search(1);
	    $f->description("auto generated anchor for GEvo");
	    $c->add_feature($f);
	    next;
	  }
        next unless $f;
	$f->color([255,0,255]) if $color_overlapped_features && $feat->qualifiers->{overlapped_hsp};
	my $strand = 1;
 	$strand = -1 if $feat->location =~ /complement/;
	foreach my $block (@{$feat->blocks})
	  {
	    $block->[0] =1 unless $block->[0]; #in case $block is set to 0
	    $f->add_segment(start=>$block->[0], stop=>$block->[1]);
	    $f->strand($strand);
	    print STDERR "\t", join ("-", @$block),"\n" if $DEBUG;
	  }

	print STDERR $name,"\n\n" if $DEBUG;
        $f->type($type);
	$f->description($feat->annotation);
	$f->link("FeatList.pl?fid=".$feat->qualifiers->{id}) if $feat->qualifiers->{id};
	my $foverlap = $overlap ? 0 : 1; #I think I need this to get this to work as expected
	$f->skip_overlap_search($foverlap);
	$f->{anchor}=1 if $anchor;
        $c->add_feature($f);
	if ($feat->type =~ /gene/ &! $feat_count{$feat->start}{$feat->stop})
	  {
	    $feat_count++;
	    $overlap_count++ if $feat->qualifiers->{overlapped_hsp};
	    $feat_count{$feat->start}{$feat->stop}=1;
	  }

    }
    return ($feat_count, $overlap_count);
  }

sub process_hsps
  {
    my %opts = @_;
    my $c = $opts{c};
    my $data = $opts{data};
    my $accn = $opts{accn};
    my $reverse = $opts{rev};
    my $seq_len = $opts{seq_length};
    my $stagger_label = $opts{stagger_label};
    my $hsp_limit = $opts{hsp_limit};
    my $hsp_limit_num = $opts{hsp_limit_num};
    my $spike_seq = $opts{spike_seq};
    my $gbobj = $opts{gbobj};
    my $eval_cutoff = $opts{eval_cutoff};
    my $score_cutoff = $opts{score_cutoff};
    my $color_hsp = $opts{color_hsp};
    my $colors = $opts{colors};
    my $hsp_overlap_limit = $opts{hsp_overlap_limit};
    my $hsp_overlap_length = $opts{hsp_overlap_length};
    my $show_hsps_with_stop_codon = $opts{show_hsps_with_stop_codon};
    #to reverse hsps when using genomic sequences from CoGe, they need to be drawn on the opposite strand than where blast reports them.  This is because CoGe::Graphics has the option of reverse drawing a region.  However, the sequence fed into blast has already been reverse complemented so that the HSPs are in the correct orientation for the image.  Thus, if the image is reverse, they are drawn on the wrong strand.  This corrects for that problem.   Sorry for the convoluted logic, but it was the simplest way to substantiate this option
    my %stats;
    my $i = 0;
    my $track = 2;
    my @feats;
    my %overlapped_feats;
    foreach my $item (@$data)
      {
	my $report = $item->[0];
	my $accn1 = $item->[1];
	my $accn2 = $item->[2];
	my $blast = $item->[3];
#	next unless ref($blast) =~ /bl2seq/i || ref($blast) =~ /blastz/i;
	unless ($accn eq $accn1 || $accn eq $accn2)
	  {
	    $i++;
	    next;
	  }
	my ($accna, $accnb) = $accn1 eq $accn ? ($accn1,$accn2) : ($accn2,$accn1);
	if ($spike_seq)
	  {
	    my $found_spike = 0;
	    foreach my $hsp (@{$blast->hsps})
	      {
		if ($hsp->qalign =~ /^$spike_seq$/i  || $hsp->salign =~ /^$spike_seq$/i)
		  {
		    $hsp->contains_spike(1);
		    if (defined $hsp->eval && $hsp->eval ne "N/A")
		      {
			$eval_cutoff = $hsp->eval;
		      }
		    elsif (defined $hsp->score)
		      {
			$score_cutoff=$hsp->score;
#			print STDERR "score cutoff: $score_cutoff\n";
		      }
              # ERIC DO NOT REMOVE THESE, I NEED THEM FOR THE PAIR TRACKING STUFF
		    write_log("Found spike sequence for $accn1 and $accn2: eval cutoff set to $eval_cutoff", $cogeweb->logfile) if defined $eval_cutoff;
		    write_log("Found spike sequence for $accn1 and $accn2: score cutoff set to $score_cutoff", $cogeweb->logfile) if defined $score_cutoff;
		    $found_spike =1;
#		    last;
		  }
	      }
	    write_log("WARNING:  Did not find spike sequence: $spike_seq in any HSP", $cogeweb->logfile) unless $found_spike;
	  }
	print STDERR "\t",$blast->query," ", $blast->subject,"\n" if $DEBUG;
	foreach my $hsp (@{$blast->hsps})
	  {
	    next if defined $eval_cutoff && $hsp->eval > $eval_cutoff;
	    next if defined $score_cutoff && $hsp->score < $score_cutoff;
	    my $color = $colors->[$i];
	    my $skip = 0;

	    if ($show_hsps_with_stop_codon && ($hsp->qalign =~ /\*/ || $hsp->salign =~ /\*/))
	      {
		for my $i (0..(length ($hsp->qalign)-1))
		  {
		    my $chr1 = substr($hsp->qalign, $i, 1);
		    my $chr2 = substr($hsp->salign, $i, 1);
		    next unless $chr1 eq "*" || $chr2 eq "*"; 
		    $skip = 1 unless ($chr1 eq "-" && $chr2 eq "*") || ($chr2 eq "-" && $chr1 eq "*");
		  }
	      }
	    next if $skip;
	    my ($start, $stop, $seq);
	    if ($accn1 eq $accn)
	      {
		$start = $hsp->qb;
		$stop = $hsp->qe;
		$seq = $hsp->qalign;
	      }
	    elsif ($accn2 eq $accn)
	      {
		$start = $hsp->sb;
		$stop = $hsp->se;
		$seq = $hsp->salign;
	      }
	    if ($stop < $start)
	      {
		my $tmp = $start;
		$start = $stop;
		$stop = $tmp;
	      }
	    print STDERR "\t",$hsp->number,": $start-$stop\n" if $DEBUG;
	    my $strand = $hsp->strand =~ /-/ ? "-1" : 1;
	    push @{$stats{data}{$accna}{$accnb}{hsp}},[$hsp->percent_id, $hsp->length];

	    #find overlapping features in gbobj
	    if (($hsp->qe-$hsp->qb+1)>=$hsp_overlap_length || ($hsp->se-$hsp->sb+1)>=$hsp_overlap_length)
	      {
		foreach my $feat ($gbobj->get_features(start=>$start, stop=>$stop))
		  {
		    #lets skip it if it only has minimal end overlap
		    next if abs($start - $feat->stop)<=3 || abs($stop - $feat->start) <=3;
		    $feat->qualifiers->{overlapped_hsp}=1;
		    next if $overlapped_feats{$feat->type}{$feat->start}{$feat->stop};
		    $overlapped_feats{$feat->type}{$feat->start}{$feat->stop}=1;
		    $stats{data}{$accna}{$accnb}{$feat->type}{overlapped_hsps}++;
		  }
	      }

	    my $f = CoGe::Graphics::Feature::HSP->new({start=>$start, stop=>$stop});
	    $color = [100,100,100] if $spike_seq && $hsp->contains_spike;
	    $f->color($color);
	    $f->order($track);
	    $f->strand($strand);
	    $f->color_matches($color_hsp);
	    if ($hsp_limit)
	      {
		$f->label($hsp->number) if $hsp->number <= $hsp_limit_num;
	      }
	    else
	      {
		$f->label($hsp->number);
	      }
	    $f->alt(join ("-",$hsp->number,$accn1,$accn2));
	    my $desc = join ("<br/>", "HSP: ".$hsp->number. qq{  <span class="small">(}.$blast->query."-". $blast->subject.")</span>", "Location: ".$start."-".$stop." (".$hsp->strand.")","Match: ".$hsp->match,"Length: ".$hsp->length,"Identity: ".$hsp->percent_id);
	    
	    $desc .= "<br/>E_val: ".$hsp->pval if $hsp->pval;
	    $desc .= "<br/>Score: ".$hsp->score if $hsp->score;
	    $desc .= "<br/>Sequence: ".$seq;
	    $desc .= qq{<br/><span class="small">(cutoff: $eval_cutoff)</span>} if defined $eval_cutoff;
	    $desc .= qq{<br/>(contains spike sequence) } if $hsp->contains_spike;
	    $f->description($desc);
	    if ($reverse)
	      {
		$strand = $strand =~ /-/ ? "1" : "-1";
		my $tmp = $seq_len - $start+1;
		$start = $seq_len - $stop+1;
		$stop = $tmp;
	      }
	    my $link = "HSPView.pl?report=$report&num=".$hsp->number."&db=".$cogeweb->basefilename.".sqlite";
#	    print STDERR $link,"\n";# if $cogeweb->basefilename =~ /_$/;
	    $link .= join ("&","&qstart=".($gbobj->start+$start-1), "qstop=".($gbobj->start+$stop-1),"qchr=".$gbobj->chromosome, "qds=". $gbobj->dataset,"qstrand=".$strand) if $gbobj->dataset;
	    $f->link($link) if $link;
	    $f->alignment($hsp->alignment);
	    push @feats, $f;
	    print STDERR $hsp->number,"-", $hsp->strand, $track,":", $strand,"\n" if $DEBUG;

	  }
	$i++;
	$track++;
      }
    my $label_location = "top";
    my $order;
    @feats = sort{$a->strand cmp $b->strand || $a->track <=> $b->track || $a->start <=> $b->start} @feats;
#    @feats = reverse @feats if $reverse;
    foreach my $f (@feats)
      {
	next unless $f->label;
	$order = $f->track unless $order;
	if ($order ne $f->track)
	  {
	    $order = $f->track;
	    $label_location = "top";
	  }
	$f->label_location($label_location) if $stagger_label;
	$c->add_feature($f);
	$c->_check_overlap($f);
	if (!$label_location)
	  {
	    $label_location = "bot";
	  }
	elsif ($label_location eq "top")
	  {
	    $label_location = "";
	  }
	else
	  {
	    $label_location = "top";
	  }
      }
    if ($hsp_overlap_limit)
      {
	foreach my $f (@feats)
	  {
	    if ($f->_overlap >=$hsp_overlap_limit)
	      {
		$c->delete_feature($f);
	      }
	  }
      }
    $stats{eval_cutoff} = $eval_cutoff;
    return \%stats;;
  }

sub generate_output
  {
    my %opts = @_;
    my $file = $opts{file};
    my $c = $opts{c};
    if ($file)
      {
        $c->generate_png(file=>$file);
      }
    else
      {
        print "Content-type: image/png\n\n";
        $c->generate_png();
      }
  } 

sub generate_obj_from_seq
  {
    my $sequence = shift;
    my $num = shift;
    my $rc = shift;
    
#    my ($obj, $file, $file_begin, $file_end, $spike_seq);
    my $obj = new CoGe::Accessory::GenBank;
    if ($sequence =~ /^LOCUS/)
      {
	#genbank sequence
	$obj->parse_genbank($sequence, $rc);
      }
   elsif ($sequence =~ /^>/)
      {
	#fasta sequence
	my ($header, $seq) = split /\n/, $sequence, 2;

	$header =~ s/>//g;
	$header =~ s/\|/_/g;
	$header =~ s/^\s+//;
	$header =~ s/\s+$//;
	my ($accn) = $header=~/^(\S*)/;
	$obj->accn($accn);
	$obj->locus($header);
	$obj->definition($header);
	$seq =~ s/\n|\r//g;
	$obj->sequence($seq);
      }
    else
      {
	#just the sequence
	$obj->accn("RAW_SEQUENCE_SUBMISSION_$num");
	$obj->locus("RAW_SEQUENCE_SUBMISSION_$num");
	$sequence =~ s/\n|\r//g;
	$obj->sequence($sequence);
      }
    if ($rc)
      {
	$obj->sequence(CoGeX::Feature->reverse_complement($obj->sequence));
      }
    return $obj
  }

sub generate_seq_file
  {
    my %opts = @_;
    my $obj = $opts{obj} || $opts{gbobj};
    my $start = $opts{start} || $opts{startpos} || 1;
    my $length = $opts{length} || 0;
    my $spike_len = $opts{spike_len} || 0;
    my $mask = $opts{mask_cds};
    my $mask_ncs = $opts{mask_ncs};
    my $seq_num = $opts{seq_num};
    my $repeat_mask = $opts{repeat_mask};
    my $t1 = new Benchmark;
    my ($file, $spike_seq, $seq) = 
      write_fasta(
		  obj=>$obj,
		  accn=>$obj->accn,
		  mask=>$mask,
		  mask_ncs=>$mask_ncs,
		  start=>$start,
		  length=>$length,
		  spike=>$spike_len,
		  seq_num=>$seq_num,
		 );
    if ($repeat_mask)
      {
	write_log("repeat masking $file", $cogeweb->logfile);
	`$REPEATMASKER $file`;
	`mv $file.masked $file` if -r "$file.masked";
      }
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2,$t1));
    print STDERR "Time to generate Sequence file:   $time\n" if $BENCHMARK && $DEBUG;
    
    return ($file, $spike_seq, $seq);
  }

sub get_obj_from_genome_db
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $featid = $opts{featid};
    my $pos = $opts{pos};
    my $dsid = $opts{dsid};
    my $rev = $opts{rev};
    my $chr = $opts{chr};
    my $up = $opts{up} || 0;
    my $down = $opts{down} || 0;
    my $t1 = new Benchmark;
    my ($feat) = $coge->resultset('Feature')->esearch({"me.feature_id"=>$featid})->next;
    unless (ref ($feat)=~ /feature/i || ($pos && $dsid))
      {
	write_log("Can't find valid feature database entry for id=$featid", $cogeweb->logfile);
	return;
      }
    my $t2 = new Benchmark;
    my $start = 0;
    my $stop = 0;
    my $seq;
    
    if ($feat)
      {
	$dsid = $feat->dataset->id;
	if ($rev)
	  {
	    $start = $feat->start-$down;
	    $stop = $feat->stop+$up;
	  }
	else
	  {
	    $start = $feat->start-$up;
	    $stop = $feat->stop+$down;
	  }
	$start = 1 if $start < 1;
	$chr = $feat->chr;
      }
    elsif ($pos)
      {
	if ($rev)
	  {
	    $start = $pos-$down;
	    $stop = $pos+$up;
	  }
	else
	  {
	    $start = $pos-$up;
	    $stop = $pos+$down;
	  }
	$start = 1 if $start < 1;
	$accn = get_dataset_info(dsid=>$dsid, chr=>$chr);
	$accn =~ s/<br>/, /g;
	$accn =~ s/<.*?>//g;
      }
    my $ds = $coge->resultset('Dataset')->find($dsid);
    return unless $ds;
    ($chr) = $ds->get_chromosomes unless $chr;
    $seq = $ds->get_genomic_sequence(
				     start => $start,
				     stop => $stop,
				     chr => $chr,
				    );
    if ($stop-$start+1 > length($seq))
      {
	my $len = length($seq);
	$stop = $start+$len-1;
      }
    $seq = CoGeX::Feature->reverse_complement($seq) if $rev;

    my $t3 = new Benchmark;
    my $obj= new CoGe::Accessory::GenBank({
					   accn=>$accn,
					   locus=>$accn,
					   version=>$ds->version(),
					   data_source=>$ds->datasource->name(),
					   dataset=>$dsid,
					   chromosome=>$chr,
					   start=>$start,
					   stop=>$stop,
					   organism=>$ds->organism->name()."(v".$ds->version.")",
					   seq_length=>length($seq),
					   sequence=>$seq,
					 });

    my %used_names;
    $used_names{$accn} = 1;
    my $t4 = new Benchmark;

    print STDERR "Region: $chr: $start-$stop\n" if $DEBUG;
#    print STDERR "Region: $chr: ",$start-$start+1,"-",$stop-$start,"\n";
    my %feats = map{$_->id,$_} $coge->get_features_in_region(start=>$start, stop=>$stop, chr=>$chr, dataset_id=>$dsid);
    $feats{$feat->id}=$feat if $feat;
    my $t5 = new Benchmark;
    foreach my $f (values %feats)
      {
	my $name;
	my @names = $f->names;
	next if $f->type->name =~ /misc_feat/;
	foreach my $tmp (@names)
	  {
	    $name = $tmp;
	    last if ($tmp =~ /^$accn.+/i);
	  }
	unless (@names)
	  {
	    next;
#	    print STDERR "No Name: ",$f->id,"\n" unless $f->type->name =~ /misc/;

	  }
	$name = $accn unless $name;
	print STDERR $name,"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(),"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(recalibrate=>$start),"\n\n" if $DEBUG;
	my $anno = $f->annotation_pretty_print();
	$anno =~ s/\n/<br\/>/ig;
	my $location = $f->genbank_location_string(recalibrate=>$start);
	$location = $obj->reverse_genbank_location(loc=>$location, ) if $rev;
	print STDERR $name, "\t",$f->type->name ,"\t",$location,"\n" if $DEBUG;
	my $type = $f->type->name;
	$type = "anchor" if $f->id && $featid && $f->id == $featid;
	$obj->add_feature (
			   type=>$f->type->name,
			   location=> $location,
			   qualifiers=>{
                                        names=> [@names],
					type=>$type,
					id=>$f->id,
                                       },
			   annotation=>$anno,
			  );
      }
    if ($pos)
      {
	my $pos_loc = $rev ? $stop-$pos: $pos-$start;
	$obj->add_feature(type=>"anchor", 
			  location=>($pos_loc), 
			  annotation=>"User specified anchor point",
			  qualifiers=>{names=>[]},
			 );
      }
    my $t6 = new Benchmark;
    my $db_time = timestr(timediff($t2,$t1));
    my $seq_time = timestr(timediff($t3,$t2));
    my $int_obj_time = timestr(timediff($t4,$t3));
    my $feat_region_time = timestr(timediff($t5,$t4));
    my $pop_obj_time = timestr(timediff($t6,$t5));
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
  my %opts = @_;
  my $sets = $opts{sets};
  my $blast_params = $opts{params};
  my $program = $opts{blast_program};
  my $parser_opts = $opts{parser_opts};
  my $eval_cutoff = $opts{eval_cutoff};
  my $spike_seq = $opts{spike_seq};
  
  $program = "blastn" unless $program;
  my @reports;
  my $total_runs = number_of_runs($sets);
  my $count = 0;
  my $pm = new Parallel::ForkManager($MAX_PROC);
  for (my $i=0; $i<scalar @$sets; $i++)
    {
      for (my $j=0; $j<scalar @$sets; $j++)
	{
	  next unless $j > $i;
	  #
	  my $seqfile1 = $sets->[$i]->{file};
	  my $seqfile2 = $sets->[$j]->{file};
	  #check_sequence_files_spike($spike_seq, $seqfile1, $seqfile2) if ($spike_seq);
	  next unless $seqfile1 && -r $seqfile1 && $seqfile2 && -r $seqfile2; #make sure these files exist
	  
	  next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
	  $count++;
	  my ($accn1, $accn2) = ($sets->[$i]{accn}, $sets->[$j]{accn});
	  my ($tempfile) = $cogeweb->basefile."_".($sets->[$i]->{seq_num})."-".($sets->[$j]->{seq_num}).".bl2seq";;
	  push @reports,[$tempfile, $accn1, $accn2];
	  $pm->start and next;
	  my $command = $BL2SEQ;
	  
	  
	# need to create a temp filename here

	  
	  # format the bl2seq command
	  $command .= "-p $program -o $tempfile ";
	  $command .= "-i $seqfile1 -j $seqfile2";
	  $command .= " " . $blast_params;
	  my $x = "";
	  ($x,$command) = check_taint( $command );
	  if ( $DEBUG ) {
	    print STDERR "About to execute...\n $command\n";
	  }
	  unless ($x)
	    {
	      next;
	    }
	  # execute the command
	  write_log("running ($count/$total_runs) ".$command, $cogeweb->logfile);
	  `$command`;
	  system "chmod +rw $tempfile";
	  $pm->finish;
	}
    }
  $pm->wait_all_children;
  foreach my $item (@reports)
    {
      my $tempfile = $item->[0];
      my $parser_opts = $opts{parser_opts};
      my $blastreport = new CoGe::Accessory::bl2seq_report({file=>$tempfile}) if -r $tempfile;
      if ($blastreport)
	{
	  push @$item, $blastreport
	}
      else
	{
	  push @$item, "no results from blasting ".$item->[1]." and ".$item->[2];
	}
    }
  return( \@reports );
}

sub run_blastz
  {
    my %opts = @_;
    my $sets = $opts{sets};
    my $params= $opts{params};
    my $parser_opts = $opts{parser_opts};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count = 0;
    my $pm = new Parallel::ForkManager($MAX_PROC);
    for (my $i=0; $i<scalar @$sets; $i++)
      {
	for (my $j=0; $j<scalar @$sets; $j++)
	  {
	    next unless $j > $i;
	    my $seqfile1 = $sets->[$i]->{file};
	    my $seqfile2 = $sets->[$j]->{file};
	    next unless $seqfile1 && $seqfile2 && -r $seqfile1 && -r $seqfile2; #make sure these files exist
	    next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
	    $count++;
	    my ($accn1, $accn2) = ($sets->[$i]{accn}, $sets->[$j]{accn});
	    my ($tempfile) = $cogeweb->basefile."_".($sets->[$i]->{seq_num})."-".($sets->[$j]->{seq_num}).".blastz";
	    push @reports, [$tempfile, $accn1, $accn2];
	    $pm->start and next;
	    my $command = $BLASTZ;
	    $command .= " $seqfile1 $seqfile2";
	    $command .= " ".$params if $params;
	    my $x = "";
	    ($x,$command) = check_taint( $command );
	    if ( $DEBUG ) {
	      print STDERR "About to execute...\n $command\n";
	    }
	    unless ($x)
	      {
		next;
	      }
	    $command .= " > ".$tempfile;
	    	  # execute the command
	    write_log("running ($count/$total_runs) ".$command, $cogeweb->logfile);
#	    print STDERR $command,"\n";
	    `$command`;
	    system "chmod +rw $tempfile";
	    $pm->finish;
	  }
      }
    $pm->wait_all_children;
    foreach my $item (@reports)
      {
	my $tempfile = $item->[0];
	my $blastreport = new CoGe::Accessory::blastz_report({file=>$tempfile}) if -r $tempfile;
	if ($blastreport)
	  {
	    push @$item, $blastreport
	  }
	else
	  {
	    push @$item, "no results from blasting ".$item->[1]." and ".$item->[2];
	  }
      }
    return \@reports;
  }


sub run_lagan
  {
    my %opts = @_;
    my $sets = $opts{sets};
    my $params= $opts{params};
    my $parser_opts = $opts{parser_opts};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count = 1;
    for (my $i=0; $i<scalar @$sets; $i++)
      {
	for (my $j=0; $j<scalar @$sets; $j++)
	  {
	    next unless $j > $i;
	    my $seqfile1 = $sets->[$i]->{file};
	    my $seqfile2 = $sets->[$j]->{file};
	    next unless -r $seqfile1 && -r $seqfile2; #make sure these files exist
	    next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
	    my ($accn1, $accn2) = ($sets->[$i]{accn}, $sets->[$j]{accn});
	    my ($tempfile) = $cogeweb->basefile."_".($sets->[$i]->{seq_num})."-".($sets->[$j]->{seq_num}).".lagan";
	    my $command = $LAGAN;
	    $command .= " $seqfile1 $seqfile2";
	    $command .= " -mfa";
	    $command .= " ".$params if $params;
	    my $x = "";
	    ($x,$command) = check_taint( $command );
	    if ( $DEBUG ) {
	      print STDERR "About to execute...\n $command\n";
	    }
	    unless ($x)
	      {
		next;
	      }
	    $command .= " > ".$tempfile;
	    write_log("running ($count/$total_runs) ".$command, $cogeweb->logfile);
	    #time for execution
	    `$command`;
	    system "chmod +rw $tempfile";
	    my $report = new CoGe::Accessory::lagan_report({file=>$tempfile, %$parser_opts}) if -r $tempfile;
	    my @tmp = ($tempfile, $accn1, $accn2);
	    if ($report)
	    {
	      push @tmp, $report
	    }
	    else
	      {
		push @tmp, "no results from comparing $accn1 and $accn2 with LAGAN";
	      }
	    push @reports, \@tmp;
	    $count++;
	  }
      }
    return \@reports;
  }

sub run_chaos
  {
    my %opts = @_;
    my $sets = $opts{sets};
    my $params= $opts{params};
    my $parser_opts = $opts{parser_opts};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count = 1;
    for (my $i=0; $i<scalar @$sets; $i++)
      {
	for (my $j=0; $j<scalar @$sets; $j++)
	  {
	    next unless $j > $i;
	    my $seqfile1 = $sets->[$i]->{file};
	    my $seqfile2 = $sets->[$j]->{file};
	    next unless -r $seqfile1 && -r $seqfile2; #make sure these files exist
	    next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
	    my ($accn1, $accn2) = ($sets->[$i]{accn}, $sets->[$j]{accn});
	    my ($tempfile) = $cogeweb->basefile."_".($sets->[$i]->{seq_num})."-".($sets->[$j]->{seq_num}).".chaos";
	    my $command = $CHAOS;
	    $command .= " $seqfile1 $seqfile2";
	    
	    $command .= " ".$params if $params;
	    my $x = "";
	    ($x,$command) = check_taint( $command );
	    if ( $DEBUG ) {
	      print STDERR "About to execute...\n $command\n";
	    }
	    unless ($x)
	      {
		next;
	      }
	    $command .= " > ".$tempfile;
	    write_log("running ($count/$total_runs) ".$command, $cogeweb->logfile);
	    #time for execution
	    `$command`;
	    system "chmod +rw $tempfile";
	    my $report = new CoGe::Accessory::chaos_report({file=>$tempfile, %$parser_opts}) if -r $tempfile;
	    my @tmp = ($tempfile, $accn1, $accn2);
	    if ($report)
	    {
	      push @tmp, $report
	    }
	    else
	      {
		push @tmp, "no results from comparing $accn1 and $accn2 with Chaos";
	      }
	    push @reports, \@tmp;
	    $count++;
	  }
      }
    return \@reports;
  }

sub run_dialign
  {
    my %opts = @_;
    my $sets = $opts{sets};
    my $params= $opts{params};
    my $parser_opts = $opts{parser_opts};
    my $max_length = $opts{max_length} || 10000;
    my $kill_length = 2 * $max_length;
    my $program_ran = "DIALIGN";
    my $error_message;
    #my $algo_limits = $opts{algo_limits};
    my @files;
    my @reports;
    my $total_runs = number_of_runs($sets);
    my $count = 1;
    for (my $i=0; $i<scalar @$sets; $i++)
      {
	for (my $j=0; $j<scalar @$sets; $j++)
	  {
	    next unless $j > $i;
	    my $obj1 = $sets->[$i]->{obj};
	    my $obj2 = $sets->[$j]->{obj};
	    my $seq_length1 = $obj1->seq_length;
	    my $seq_length2 = $obj2->seq_length;
	    my $max_seq = $seq_length1 > $seq_length2 ? $seq_length1 : $seq_length2;
	    my $kb_max_seq = $max_seq / 1000;
	    $kb_max_seq =~s/(\.\d+)$//;
	    if ($max_seq > $max_length)
	    {
	      if ($max_seq > $kill_length)
	      {
	        $error_message = "The sequence you have chosen DIALIGN to align is too long (your longest sequence is $kb_max_seq kb long). To complete the alignment would take a significant amount of time, and use up considerable system resources that others need to run their alignments. Please limit your sequences to no more than ".($kill_length/1000)." kb long. Thank you.";
	       return ('','',$error_message);
	      }
	      else
	      {
	       my $changed = $params=~s/((\s?-cs)|(\s?-nt)|(\s?-ma))//g;
	       $params .= " -n" unless $params =~ /(-n)/;
	       $params .= " -ds" unless $params =~ /(-ds)/;
	       $params .= " -thr 2" unless $params =~ /(-thr)/;
	       $params .= " -lmax 30" unless $params =~ /(-lmax)/;
	       $error_message = "Your search parameters were altered due to the length of the sequences you requested for alignment. Your alignment was done with the \"Large Sequence Alignment\" option enabled";
	       $error_message .= $changed ? ", and with \"Short Sequence Alignment with Translation\" disabled." : ".";
	       $error_message.= " This is done automatically with sequences over ".($max_length/1000)." kb long.\n";
	       #print $error_message;
	      }
	     }
	    my $seqfile1 = $sets->[$i]->{file};
	    my $seqfile2 = $sets->[$j]->{file};
	    my ($accn1, $accn2) = ($sets->[$i]{accn}, $sets->[$j]{accn});
	    my $report; #place to hold report
	    my $tempfile;
	    if ($parser_opts->{anchor_params})
	      {
		my $anchor_prog;
		if ($parser_opts->{anchor_params}{prog} =~ /chaos/i)
		{
		   $anchor_prog = $CHAOS;
		}
		elsif ($parser_opts->{anchor_params}{prog} =~ /bl2/i)
		{
		   $anchor_prog = $BL2SEQ;
		}
		else
		{
		   $anchor_prog = $BLASTZ;
		}
		my ($x,$dialign_opts) = check_taint( $params);
		$program_ran .= " using anchors from ".$parser_opts->{anchor_params}{prog};
		unless ($x)
		  {
		    next;
		  }
		my ($y,$anchor_opts) = check_taint($parser_opts->{anchor_params}{param_string});
		unless ($y)
		  {
		    next;
		  }
		my $obj = new CoGe::Accessory::dialign_report::anchors({
									file1=>$seqfile1,
									file2=>$seqfile2,
									base_name=>$cogeweb->basefilename,
									extension=>$parser_opts->{anchor_params}{prog},
									run_anchor=>$anchor_prog,
									run_dialign_opts=>$dialign_opts,
									run_anchor_opts=>$anchor_opts,
									dialign_report_opts=>$parser_opts,
									anchor_report_opts=>$parser_opts->{anchor_params}{parser_opts},
									log_file=>$cogeweb->logfile,
#									DEBUG=>1,
								       });
		$tempfile = $obj->dialign_file;
		$report = $obj->dialign_report;
	      }
	    else
	      {
		my $seqfile = $cogeweb->basefile."_".($sets->[$i]->{seq_num})."-".($sets->[$j]->{seq_num}).".fasta";
		next unless -r $seqfile1 && -r $seqfile2; #make sure these files exist
		#put two fasta files into one for dialign
		`cat $seqfile1 > $seqfile`;
		`cat $seqfile2 >> $seqfile`;
		next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
		($tempfile) = $cogeweb->basefile."_".($sets->[$i]->{seq_num})."-".($sets->[$j]->{seq_num}).".dialign";
		my $command = $DIALIGN;
		$command .= " ".$params if $params;
		$command .= " -fn ".$tempfile;
		$command .= " $seqfile";
		my $x = "";
		($x,$command) = check_taint( $command );
		if ( $DEBUG ) {
		  print STDERR "About to execute...\n $command\n";
		}
		unless ($x)
		  {
		    next;
		  }
		#time for execution
		write_log("running ($count/$total_runs) ".$command, $cogeweb->logfile);
		`$command`;
		system "chmod +rw $tempfile";
		$report = new CoGe::Accessory::dialign_report({file=>$tempfile, %$parser_opts}) if -r $tempfile;
	      }
	    my @tmp = ($tempfile, $accn1, $accn2);
	    if ($report)
	      {
		push @tmp, $report
	      }
	    else
	      {
		push @tmp, "no results from comparing $accn1 and $accn2 with LAGAN";
	      }
	    push @reports, \@tmp;
	    $count++;
	  }
      }
    return (\@reports,$program_ran,$error_message);
  }
  
sub get_substr
  {
    my %opts = @_;
    my $seq = $opts{seq};
    my $start = $opts{start};
    $start =~ s/,|\.//g;
    $start = 1 if !$start || $start < 1;
    my $stop = $opts{stop};
    return $seq if ($start == 1 && !$stop);
    my $seqlength = length $seq;
    $stop = $seqlength unless $stop;
    $stop =~ s/,|\.//g;
    my $newlength = $stop-$start;
    if ($start < $seqlength && $newlength)
      {
	$seq = substr($seq, $start-1, $newlength);
      }
    return $seq;
  }

sub write_fasta 
  {
    my %opts = @_;
    my $gbobj = $opts{obj};
#    print STDERR Dumper $gbobj;
    my $start = $opts{start};
    my $mask = $opts{mask};
    my $mask_ncs = $opts{mask_ncs};
    my $length = $opts{length};
    my $spike = $opts{spike};
    my $seq_num = $opts{seq_num};
    my $fullname = $cogeweb->basefile."_".$seq_num.".faa";
    my $hdr = $gbobj->get_headerfasta( );
    $start = 1 if $start < 1;
    my ($seq) = uc($gbobj->sequence());

    my $stop = $length? $start+$length-1 :length($seq);
    # mask out exons if required
    $seq = $gbobj->mask_exons( $seq ) if ( $mask );
    # mask out non coding sequences if required
    $seq = $gbobj->mask_ncs( $seq ) if ( $mask_ncs );
    ($seq) = $gbobj->subsequence( $start, $stop, $seq );
    my $spike_seq = "";
    if ($spike)
      {
	($spike_seq) = generate_spike_seq($spike);
	$seq = $seq."N" x length($spike_seq).$spike_seq."N" x length($spike_seq).$spike_seq;#."N" x length($spike_seq);
      }
    ($fullname) = check_filename_taint( $fullname );
    $length = length($seq);
    open(OUT, ">$fullname") or die "Couldn't open $fullname!: $!\n";
    print OUT "$hdr\n";
    my $max = 100;
    my $i = 0;
    print OUT $seq,"\n";
    close(OUT);
    system "chmod +rw $fullname";
    write_log("Created sequence file $seq_num for $hdr.  Length:". $length, $cogeweb->logfile);
    write_log('spike:' . $spike_seq, $cogeweb->logfile) if $spike_seq;
    return($fullname, $spike_seq, $seq);
  }

sub generate_spike_seq { 
	my $spikesize = shift;
	return unless $spikesize;
	my %nt = (0=>'A', 3=>'T', 2=>'C', 1=>'G');
	my $spike_seq = "";
	my $idx=1;
	my $i = 1;
	while (length($spike_seq) <=$spikesize)
	  {
	    $spike_seq .= $nt{($i+$idx)%4}x($idx);
	    $idx++ unless $i%4;
	    $i+=$idx;
	  }
	$spike_seq = substr($spike_seq, 0, $spikesize);

	return($spike_seq);
}

sub spike_filter_select
  {
    my $form = shift || $FORM;
    my $match = $form->param('spike_len') ? $form->param('spike_len') : 0; ;
    return $match;
  }

sub generate_annotation
  {
    my %opts = @_;
    my $obj = $opts{obj};
    return unless $obj;
    my $rev = $opts{rev};
    my $seq_num = $opts{seq_num};
    my $fullname = $cogeweb->basefile."_".$seq_num.".anno";
    my %data;
    my $length = length($obj->sequence());
    foreach my $feat($obj->get_features())
      {
	my $type = $feat->type;
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} @{$feat->qualifiers->{names}} if ref($feat->qualifiers) =~ /hash/i && $feat->qualifiers->{names};
	my $dir = $feat->location =~ /complement/ ? "<" : ">";
	foreach my $block (@{$feat->blocks})
	  {
	    $block->[0] =1 unless $block->[0];
	    my $start = $block->[0];
	    my $stop = $block->[1];
	    if ($rev)
	      {
		$start = $length - $block->[1];
		$stop = $length - $block->[0];
		$dir = $dir eq ">" ? "<" : ">";
	      }
	    next unless $name;
	    push @{$data{$name}{$type}},[$start, $stop, $dir];
	  }
      }
    return unless keys %data;
    open (OUT, ">$fullname") || die "Can't open $fullname for writing! $!\n";
    foreach my $name (keys %data)
      {
	my $gene = $data{$name}{gene};
	if ($gene)
	  {
	    delete $data{$name}{gene};
	    print OUT join (" ", $gene->[0][2], $gene->[0][0], $gene->[0][1], $name),"\n";
	    foreach my $type (keys %{$data{$name}})
	      {
		my $outtype;
		$outtype = "utr" if $type =~ /rna/i;
		$outtype = "exon" if $type =~ /cds/i;
		next unless $outtype;
		foreach my $item (@{$data{$name}{$type}})
		  {
		    next if $item->[0] == $item->[1];
		    print OUT join (" ", $item->[0], $item->[1],$outtype),"\n";
		  }
	      }
	  }
      }
    close OUT;
    return $fullname;
  }

sub gen_params
  {
    my $num_seqs = shift || $NUM_SEQS;
    my $params;
    for (my $i = 1; $i <=$num_seqs; $i++)
      {
	$params .= qq{'args__draccn$i', 'accn$i',};
	$params .= qq{'args__featid$i', 'featid$i',};
	$params .= qq{'args__drup$i', 'drup$i',};
	$params .= qq{'args__drdown$i', 'drdown$i',};
	$params .= qq{'args__pos$i', 'pos$i',};
	$params .= qq{'args__dsid$i', 'dsid$i',};
	$params .= qq{'args__posdsid$i', 'posdsid$i',};
	$params .= qq{'args__chr$i', 'chr$i',};
	$params .= qq{'args__gbaccn$i', 'gbaccn$i',};
	$params .= qq{'args__gbstart$i', 'gbstart$i',};
	$params .= qq{'args__gblength$i', 'gblength$i',};
	$params .= qq{'args__dirseq$i', 'dirseq$i',};
	$params .= qq{'args__dirstart$i', 'dirstart$i',};
	$params .= qq{'args__dirlength$i', 'dirlength$i',};
	$params .= qq{'args__ref_seq$i', 'ref_seq$i',};
	$params .= qq{'args__skip_seq$i', 'skip_seq$i',};
	$params .= qq{'args__repmask$i', 'repmask$i',};
	$params .= qq{'args__rev$i', 'rev$i',};
	$params .= qq{'args__maskcds$i', 'mask_cds$i',};
	$params .= qq{'args__maskncs$i', 'mask_ncs$i',};

      }
    for (my $i = 1; $i <=num_colors($num_seqs); $i++)
      {
	$params .= qq{'args__r$i', 'r$i',};
	$params .= qq{'args__g$i', 'g$i',};
	$params .= qq{'args__b$i', 'b$i',};
      }


    $params .= qq{
        'args__pad_gs', 'pad_gs',

	'args__spike', 'spike', 
        'args__show_spike', 'show_spike',

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
	'args__ih', 'ih', 
	'args__fh', 'feat_h',
        'args__gc', 'show_gc',
        'args__nt', 'show_nt',
	'args__colorhsp', 'color_hsp',
	'args__hsplabel', 'hsp_labels',
	'args__overlap', 'overlap_adjustment',
	'args__hiqual', 'hiqual',
	'args__hsplim', 'hsp_limit',
	'args__hsplimnum', 'hsp_limit_num',
	'args__prog', 'alignment_program',
	'args__showallhsps', 'show_hsps_with_stop_codon',
        'args__padding', 'padding',
        'args__num_seqs','args__$num_seqs',
        'args__draw_model', 'draw_model',
        'args__hsp_overlap_limit', 'hsp_overlap_limit',
        'args__color_overlapped_features','color_overlapped_features',
        'args__hsp_overlap_length','hsp_overlap_length',
        'args__basefile','args__'+pageObj.basefile

};
    $params =~ s/\n//g;
    $params =~ s/\s+/ /g;
    return $params;
  }
sub gen_go_run
  {
    my $num_seqs = shift || $NUM_SEQS;
    my $params = gen_params($num_seqs);
    my $run = qq!
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

sub gen_save_settings
  {
    my $num_seqs = shift || $NUM_SEQS;
    my $params = gen_params($num_seqs);
    my $save_settings = qq{save_settings_gevo([$params],[])};
    return $save_settings;
  }

sub gen_hsp_colors
  {
    my %opts = @_;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my $prefs = $opts{prefs} || load_settings(user=>$USER, page=>$PAGE_NAME);
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
    my $hsp_colors;
    my @colors = color_pallet(num_seqs=>$num_seqs, prefs=>$prefs);
    $template->param(HSP_COLOR_FORM=>1);
    $template->param(HSP_COLOR_LOOP=>\@colors);
    my $count = -2;
    foreach my $line (split /\n/, $template->output)
      {
	next unless $line;
	$count ++ if $line =~/<table>/i;
	
	if ($count == 6)
	  {
	    $line =~ s/(<table>)/<td>$1/i;
	    $count = 0;
	  }
	
	$hsp_colors .= $line."\n";
      }
    return $hsp_colors;
  }

sub color_pallet
  {
    my %opts = @_;
    my $start = $opts{start} || [255,100,100];
    my $offset = $opts{offset} || 50;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my $prefs = $opts{prefs} || load_settings(user=>$USER, page=>$PAGE_NAME);
    $prefs = {} unless $prefs;
    $start = [$prefs->{r1}, $prefs->{g1}, $prefs->{b1}] if defined $prefs->{r1} && defined $prefs->{g1} && defined $prefs->{b1};

    my @colors;
    my %set = (HSP_NUM=>1,
	       RED=>$start->[0],
	       GREEN=>$start->[1],
	       BLUE=>$start->[2]);
#    push @colors, \%set;
    my $temp = [@$start];
    for (my $i = 1; $i <= num_colors($num_seqs); $i++)
      {
	my @color;
	@color = @$temp;
	@color = ($prefs->{"r".$i}, $prefs->{"g".$i}, $prefs->{"b".$i}) if defined $prefs->{"r".$i} && defined $prefs->{"g".$i} && defined $prefs->{"b".$i};
	push @colors, {
		       HSP_NUM=>$i,
		       RED=>$color[0],
		       GREEN=>$color[1],
		       BLUE=>$color[2],
		      };
	if ($i%3)
	  {
	    $temp = [map {int($_/1.5)} @color];
	    
	  }
	else
	  {
	    $start = [$start->[2], $start->[0], $start->[1]];
	    $temp =[@$start];
	  }
	unless ($i%9)
	  {
	    $start = [$start->[0], int($start->[0]/1.5), $start->[1]];
	    $temp =[@$start];
	  }
	unless ($i%18)
	  {
	    $start = [$start->[0], $start->[0], int($start->[0]/4)];
	    $temp =[@$start];
	  }
      }
    return wantarray ? @colors : \@colors;
  }

sub num_colors
  {
    my $num_seqs = shift || $NUM_SEQS;
    my $num_colors = 1;
    for (my $i = $num_seqs-1; $i > 1; $i--)
      {
	$num_colors += $i;
      }
    return $num_colors
  }

sub algorithm_list
  {
    my $program = shift;
    $program = "blastz" unless $program;
    my @programs = sort {lc $a cmp lc $b} qw(blastn blastz CHAOS DiAlign_2 LAGAN tblastx);
    my $html;
    foreach my $prog (@programs)
      {
	$html .= "<option";
	$html .= $program && $program =~ /$prog/i ? " selected>" : ">";
	$html .= $prog."</option>\n";
      }
    return $html;
  }

sub add_seq
  {
    my $num_seq = shift;
    $num_seq ++;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
    my $hsp_colors;
    if ($num_seq > $MAX_SEQS)
      {
	$hsp_colors = gen_hsp_colors(num_seqs=>$MAX_SEQS);
	my $go_run = gen_go_run($MAX_SEQS);
	return ('',$MAX_SEQS, $go_run, '', $hsp_colors,qq{Exceeded max number of sequences ($MAX_SEQS)});
      }
    my $dsinfo =qq{<input type="hidden" id="posdsid$num_seq">};
    $dsinfo .= qq{<input type="hidden" id="chr$num_seq">};
    my @seqs = {
		SEQ_NUM=>$num_seq,
		REV_NO=>"checked",
		EXON_MASK_OFF=>"checked",
		REF_YES=>"checked",

		DRUP=>10000,
		DRDOWN=>10000,
		DSINFO=>$dsinfo,
		GBSTART=>1,
	       };

    $template->param(SEQ_SELECT=>1);
    $template->param(SEQ_SELECT_LOOP=>\@seqs);
    my $seq_submission = $template->output;
    $template->param(SEQ_SELECT=>0);

    $template->param(GEN_REFERENCE_SEQUENCES=>1);
    $template->param(REF_SEQ=>[{SEQ_NUM=>$num_seq}]);	
    my $ref_seqs = $template->output;
    $template->param(GEN_REFERENCE_SEQUENCES=>0);
    my $go_run = gen_go_run($num_seq);
    $hsp_colors = gen_hsp_colors(num_seqs=>$num_seq);
    return ($seq_submission, $num_seq, $go_run, $ref_seqs, $hsp_colors);
  }

sub hsp_color_cookie
  {
    my $cname = "GEVO_HSP_COLORS";
    my %cookies = fetch CGI::Cookie;
    if (ref $cookies{$cname})
      {
      }
    else
      {
      }
  }

sub get_algorithm_options
  {
    my %opts = @_;
    my $analysis_program = $opts{prog};
    #blast stuff
    my $blast_wordsize = $opts{blast_word};
    my $blast_gapopen = $opts{blast_gapo};
    my $blast_gapextend = $opts{blast_gape};
    my $blast_mismatch = $opts{blast_mmatch};
    my $blast_match = $opts{blast_match};
    my $blast_eval = $opts{blast_eval};
    my $blast_params = $opts{blast_params};
    my $blast_filter = $opts{blast_filter};

    #blastz stuff
    my $blastz_wordsize = $opts{blastz_wordsize};
    my $blastz_chaining = $opts{blastz_chaining};
    my $blastz_threshold = $opts{blastz_threshold};
    my $blastz_mask = $opts{blastz_mask};
    my $blastz_gap_start = $opts{blastz_gap_start};
    my $blastz_gap_extension = $opts{blastz_gap_extension};
    my $blastz_params = $opts{blastz_params};

    #lagan_stuff
    my $lagan_params = $opts{lagan_params};
    my $lagan_min_length = $opts{lagan_min_length};
    my $lagan_max_gap = $opts{lagan_max_gap};
    my $lagan_percent_id = $opts{lagan_percent_id};
    
    #dialign_stuff
    my $dialign_params = $opts{dialign_params};
    my $dialign_motif_motif = uc($opts{dialign_motif_motif});
    my $dialign_motif_weight = $opts{dialign_motif_weight};
    my $dialign_motif_penalty = $opts{dialign_motif_penalty};
    my $dialign_min_score = $opts{dialign_min_score};
    my $dialign_max_gap = $opts{dialign_max_gap};
    my $dialign_split_score = $opts{dialign_split_score};
    my $dialign_use_anchor = $opts{dialign_use_anchor};
    my $dialign_anchor_program = $opts{dialign_anchor_program};
    
    #chaos stuff
    my $chaos_word_length = $opts{chaos_word_length};
    my $chaos_score_cutoff = $opts{chaos_score_cutoff};
    my $chaos_rescore_cutoff = $opts{chaos_rescore_cutoff};
    my $chaos_lookback = $opts{chaos_lookback};
    my $chaos_gap_length = $opts{chaos_gap_length};
    my $chaos_gap_start = $opts{chaos_gap_start};
    my $chaos_gap_extension = $opts{chaos_gap_extension};
    my $chaos_params = $opts{chaos_params};

    my ($blastz_string, $blast_string, $lagan_string, $dialign_string, $chaos_string);
    #build blastz param string
    $blastz_string = " W=" .$blastz_wordsize if $blastz_wordsize;
    $blastz_string .= " C=" .$blastz_chaining if $blastz_chaining;
    $blastz_string .= " K=" .$blastz_threshold if $blastz_threshold;
    $blastz_string .= " M=" .$blastz_mask if $blastz_mask;
    $blastz_string .= " O=" .$blastz_gap_start if $blastz_gap_start;
    $blastz_string .= " E=" .$blastz_gap_extension if $blastz_gap_extension;
    $blastz_string .= " " .$blastz_params if $blastz_params;
    #build blast param string
    $blast_wordsize = 3 if ($analysis_program eq "tlastx" && $blast_wordsize > 3);
    $blast_string = " -W " . $blast_wordsize if defined $blast_wordsize;
    $blast_string .= " -G " . $blast_gapopen if defined $blast_gapopen;
    $blast_string .= " -E " . $blast_gapextend if defined $blast_gapextend;
    $blast_string .= " -q " . $blast_mismatch if defined $blast_mismatch;
    $blast_string .= " -r " . $blast_match if defined $blast_match;
    $blast_string .= " -e " . $blast_eval if defined $blast_eval;
    $blast_string .= " -F " . $blast_filter if defined $blast_filter;
    $blast_string .= " "    . $blast_params if defined $blast_params;
    #build lagan param string
    $lagan_string = $lagan_params;
    #build dialign param string
    $dialign_motif_motif =~ s/[^ATCG\[\]]//g if $dialign_motif_motif;
    $dialign_string .= "-mot $dialign_motif_motif $dialign_motif_weight $dialign_motif_penalty" if $dialign_motif_motif =~/^[ATGC\[\]]+$/;
    $dialign_string .= " ".$dialign_params;
    #build chaos param string
    $chaos_string .= " -v"; #need verbose mode for the parser
    $chaos_string .= " -b"; #search both strands
    $chaos_string .= " -wl $chaos_word_length" if defined $chaos_word_length;
    $chaos_string .= " -co $chaos_score_cutoff" if defined $chaos_score_cutoff;
    $chaos_string .= " -rsc $chaos_rescore_cutoff" if defined $chaos_rescore_cutoff;
    $chaos_string .= " -lb $chaos_lookback" if defined $chaos_lookback;
    $chaos_string .= " -gl $chaos_gap_length" if defined $chaos_gap_length;
    $chaos_string .= " -gs $chaos_gap_start" if defined $chaos_gap_start;
    $chaos_string .= " -gc $chaos_gap_extension" if defined $chaos_gap_extension;
    $chaos_string .= " $chaos_params" if defined $chaos_params;


    #stuff to be returned
    my $param_string;  #param string to be passed to command line for program execution
    my %parser_options; #params used to set the parser object -- key is accessor method, value is value

    if ($analysis_program eq "blastz")
      {
	$param_string = $blastz_string;
      }
    elsif ($analysis_program =~ /blast/i) #match blastn and tblastx
      {
	$param_string = $blast_string;
      }
    elsif ($analysis_program eq "LAGAN")
      {
	$param_string = $lagan_string;
	$parser_options{max_gap} = $lagan_max_gap;
	$parser_options{length_cutoff} = $lagan_min_length;
	$parser_options{percent_cutoff} = $lagan_percent_id;
      }
    elsif ($analysis_program eq "DIALIGN")
      {
	$param_string = $dialign_string;
        $parser_options{min_score} = $dialign_min_score;
        $parser_options{max_gap} = $dialign_max_gap;
        $parser_options{split_score} = $dialign_split_score;
#	print STDERR "use anchors: $dialign_use_anchor\n";
        my ($anchor_program, $anchor_param_string);
        if ($dialign_anchor_program =~ /chaos/i)
		{
		   ($anchor_program, $anchor_param_string) = ("CHAOS", $chaos_string);
		}
	elsif ($dialign_anchor_program =~ /bl2/i)
		{
		   ($anchor_program, $anchor_param_string) = ("bl2seq", $blast_string);
		}
	else
		{
		   ($anchor_program, $anchor_param_string) = ("blastz", $blastz_string);
		}
	$parser_options{anchor_params} = {
					  prog=>$anchor_program,
					  param_string=>$anchor_param_string,
					 } if $dialign_use_anchor;
      }
    elsif ($analysis_program eq "CHAOS")
      {
	$param_string = $chaos_string;
      }
    my ($x, $clean_param_string) = check_taint($param_string);
    unless ($x)
      {
	print STDERR "user $USER's param_string: $param_string was tainted!\n";
      }
    return ($analysis_program, $clean_param_string, \%parser_options);
  }

sub check_sequence_files_spike
  {
    my $spike = shift;
    return unless $spike;
    my $start = length($spike) * -1;
    my @files = @_;
    my $check_nt;
    foreach my $file (@files)
      {
	next unless $file && -r $file;
	$/ = "\n>";
	open (IN, $file);
	my ($name, $seq);
	while (<IN>)
	  {
	    ($name, $seq) = split /\n/, $_,2;
	  }
	close IN;
	$/ = "\n";
	$seq =~ s/>//g;
	$name =~ s/>//g;
	$seq =~ s/\n//g;
	my $nt = substr($seq, $start-1,1);
	unless ($check_nt)
	  {
	    $check_nt = $nt;
	    next;
	  }
	if ($nt =~ /$check_nt/i)
	  {
	    write_log(join ("-", @files)." had similar ends.  Additional sequence added between before spike sequence: "."N" x length($spike), $cogeweb->logfile);
	    substr($seq, $start, 0) = "N" x length($spike);
	    open (OUT, ">$file");
	    print OUT ">$name\n";
	    print OUT $seq,"\n";
	    close OUT;
	  }
      }
  }

sub number_of_runs
  {
    my $sets = shift;
    return unless $sets;
    my $refs = 0;
    my $num = @$sets;
    $num--;
    foreach my $item (map {$_->{reference_seq}} @$sets)
      {
	$refs++ if $item;
      }
    my $total =0;
    for (1..$refs)
      {
	$total+=$num;
	$num--;
      }
    return $total;
  }

sub dataset_search
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $num = $opts{num};
    $num = 1 unless $num;
    my $dsid = $opts{dsid};
    my $featid = $opts{featid};
    my $feat = $coge->resultset('Feature')->find($featid) if $featid;
    $dsid = $feat->dataset->id if $feat;
    return ( qq{<input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">}, $num )unless $accn;
    my $html;
    my %sources;
    ($USER) = CoGe::Accessory::LogUser->get_user();
    my $restricted_orgs = restricted_orgs(user=>$USER);
    my $rs = $coge->resultset('Dataset')->search(
						  {
						   'feature_names.name'=> $accn,
						  },
						  {
						   'join'=>{
							    'features' => 'feature_names',
							   },
							    
						   'prefetch'=>['datasource', 'organism'],
						  }
						 );
    while (my $ds = $rs->next())
      {
	my $name = $ds->name;
	my $ver = $ds->version;
	my $desc = $ds->description;
	my $sname = $ds->datasource->name;
	my $ds_name = $ds->name;
	my $org = $ds->organism->name;
	my $title = "$org: $ds_name ($sname, v$ver)";
	next if $restricted_orgs->{$org};
	$sources{$ds->id} = {
			     title=>$title,
			     version=>$ver,
			    };
      }
     if (keys %sources)
       {
 	$html .= qq{
 <SELECT name = "dsid$num" id= "dsid$num" onChange="feat_search(['args__accn','accn$num','args__dsid', 'dsid$num', 'args__num','args__$num', 'args__featid'],['feat$num']);">
 };
 	foreach my $id (sort {$sources{$b}{version} <=> $sources{$a}{version}} keys %sources)
 	  {
 	    my $val = $sources{$id}{title};
 	    $html  .= qq{  <option value="$id"};
	    $html .= qq{ selected } if $dsid && $id == $dsid;
	    $html .= qq{>$val\n};
 	  }
 	$html .= qq{</SELECT>\n};
 	my $count = scalar keys %sources;
 	$html .= qq{<font class=small>($count)</font>};
       }
     else
       {
 	$html .= qq{<span class=container>Accession not found.</span> <input type="hidden" id="dsid$num">\n<input type="hidden" id="featid$num">\n};	
       }    
    return ($html,$num, $featid);
  }


sub save_settings_gevo
  {
    my %opts = @_;
    my $opts = Dumper \%opts;
    my $item = save_settings(opts=>$opts, user=>$USER, page=>$PAGE_NAME);
  }

sub reset_settings_gevo
  {
    my %opts = @_;
    my $item = reset_settings(user=>$USER, page=>$PAGE_NAME);
  }



sub get_opt
  {
    my %opts = @_;
    my $params = $opts{params};
    my $form = $opts{form} || $FORM;
    my $param = $opts{param};
    my $opt;
    $opt = $form->param($param) if defined $form->param($param);
    $opt = $params->{$param} if (ref ($params) =~ /hash/i &! defined $opt);
    return $opt;
      
  }
sub get_dataset_info
  {
    my %opts = @_;
    my $dsid = $opts{dsid};
    my $chr = $opts{chr};
    
    my ($ds) = $coge->resultset('Dataset')->find($dsid);
    return "<span class=\"small alert\">Dataset id $dsid was not found</span>"unless $ds;
    my $name = $ds->name;
    my $ver = $ds->version;
    my $desc = $ds->description;
    my $sname = $ds->datasource->name;
    my $ds_name = $ds->name;
    my $org = $ds->organism->name;
    $chr = join (", ",$ds->get_chromosomes) unless $chr;
    my $title = qq{<span class="species"><a href=GenomeView.pl?org_name=$org target=_new>$org v$ver:</a></span> <a href=GenomeView.pl?dsid=$dsid target=_new>chr, $chr</a>};
    return $title;
  }  

sub feat_search
  {
    my %opts = @_;
    my $accn = $opts{accn};
    my $dsid=$opts{dsid};
    my $num = $opts{num};
    my $featid = $opts{featid};
    return qq{<input type="hidden" id="featid$num">\n} unless $dsid;
    my @feats;
    my $rs = $coge->resultset('Feature')->search(
						  {
						   'feature_names.name'=> $accn,
						   'dataset.dataset_id' => "$dsid",
						  },
						  {
						   'join'=>['feature_type','dataset', 'feature_names'],
						   'prefetch'=>['feature_type', 'dataset'],
						  }
						 );
    my %seen;
    my @genes; #get genes on top of list!
    while( my $f =$rs->next())
      {
	next unless $f->dataset->id == $dsid;
	unless ($seen{$f->id})
	  {
	    if ($f->type->name eq "gene")
	      {
		push @genes, $f;
	      }
	    else
	      {
		push @feats, $f;
	      }
	  }
	$seen{$f->id}=1;
      }
    @feats = sort {$a->type->name cmp $b->type->name} @feats;
    unshift @feats, @genes if @genes;
    my $html;
    if (@feats)
      {
	
	$html .= qq{
<SELECT name = "featid$num" id = "featid$num" >
  };
	foreach my $feat (@feats)
	  {
	    my $loc = "(".$feat->type->name.") Chr:".$feat->locations->next->chromosome." ".$feat->start."-".$feat->stop;
	    $loc =~ s/(complement)|(join)//g;
	    my $fid = $feat->id;
	    $html .= qq {  <option value="$fid"};
	    $html .= qq { selected } if $featid && $featid == $fid;
	    $html .= qq{>$loc \n};
	  }
	$html .= qq{</SELECT>\n};
	my $count = scalar @feats;
	$html .= qq{<font class=small>($count)</font>};
      }
    else
      {
	$html .=  qq{<input type="hidden" id="featid$num">\n}
      }
    return $html;
  }
sub commify {
        my $input = shift;
        $input = reverse $input;
        $input =~ s<(\d\d\d)(?=\d)(?!\d*\.)><$1,>g;
        return scalar reverse $input;
}
