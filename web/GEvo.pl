#!/usr/bin/perl -w
use strict;
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
#use Text::Wrap qw($columns &wrap);
use Benchmark qw(:all);

# for security purposes

$ENV{PATH} = "/opt/apache/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };
#for chaos
$ENV{'LAGAN_DIR'} = '/opt/apache/CoGe/bin/lagan/';
#for dialign
$ENV{'DIALIGN2_DIR'} = '/opt/apache/CoGe/bin/dialign2_dir/';
use vars qw( $DATE $DEBUG $BL2SEQ $BLASTZ $LAGAN $CHAOS $DIALIGN $TEMPDIR $TEMPURL $USER $FORM $cogeweb $BENCHMARK $coge $NUM_SEQS $MAX_SEQS $BASEFILE $BASEFILENAME);
$BL2SEQ = "/opt/bin/bio/bl2seq ";
$BLASTZ = "/usr/bin/blastz ";
$LAGAN = "/opt/apache/CoGe/bin/lagan/lagan.pl";
$CHAOS = "/opt/apache/CoGe/bin/lagan/chaos_coge";
$DIALIGN = "/opt/apache/CoGe/bin/dialign2_dir/dialign2-2_coge";
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$BENCHMARK = 1;
$NUM_SEQS = 3;
$MAX_SEQS = 6;
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
#print STDERR Dumper $USER;

my %ajax = CoGe::Accessory::Web::ajax_func();
#$ajax{dataset_search} = \&dataset_search_for_feat_name; #override this method from Accessory::Web
my $pj = new CGI::Ajax(
		       run=>\&run,
		       loading=>\&loading,
		       add_seq=>\&add_seq,
		       get_file=>\&get_file,
		       gen_go_run=>\&gen_go_run,
		       gen_hsp_colors =>\&gen_hsp_colors,
		       %ajax,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);

#$USER=1;print gen_html();

sub loading
  {
    return qq{<font class="loading">Generating results. . .</font>}; 
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
	my $num_seqs = $FORM->param("num_seqs") || $NUM_SEQS;
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	$template->param(LOGO_PNG=>"GEvo-logo.png");
	$template->param(TITLE=>'Genome Evolution Analysis');
	$template->param(HELP=>'GEvo');
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(NO_BOX=>1);
	$template->param(BODY=>gen_body($num_seqs));
	my $prebox = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
	$prebox->param(RESULTS_DIV=>1);
	$template->param(PREBOX=>$prebox->output);
	
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $num_seqs = shift;;
    my $form = $FORM;
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
	my $draccn = $form->param("accn".$i) if $form->param("accn".$i);
	my $drup = $form->param('dr'.$i.'up') if $form->param('dr'.$i.'up');
	my $drdown = $form->param('dr'.$i.'down') if $form->param('dr'.$i.'down');
	$drup = $form->param('drup'.$i) if $form->param('drup'.$i);
	$drdown = $form->param('drdown'.$i) if $form->param('drdown'.$i);
	$drup = 10000 unless defined $drup;
	$drdown = 10000 unless defined $drdown;
	my $dsid = $form->param('dsid'.$i) if $form->param('dsid'.$i);
	my $gbaccn = $form->param("gbaccn".$i) if $form->param("gbaccn".$i);
	my $gbstart = $form->param("gbstart".$i) if $form->param("gbstart".$i);
	$gbstart = 1 unless defined $gbstart;
	my $gblength = $form->param("gblen".$i) if $form->param("gblen".$i);
	my $revy = "checked" if $form->param('rev'.$i);
	my $revn = "checked" unless $revy;
	$autosearch_string .= 'if ($'.qq!('#accn$i').val()) {dataset_search(['accn$i','args__$i', 'args__!;
	$autosearch_string .= $dsid if $dsid;
	$autosearch_string .=qq!'],[feat_search_chain]);}!;
	push @seq_nums, {
			  SEQ_NUM=>$i,
			 };
	push @seq_sub, {
			SEQ_NUM=>$i,
			REV_YES=>$revy,
			REV_NO=>$revn,
			DRUP=>$drup,
			DRDOWN=>$drdown,
			DRACCN=>$draccn,
			DSID=>$dsid,
			GBACCN=>$gbaccn,
			GBSTART=>$gbstart,
			GBLENGTH=>$gblength,
		       }
      }

    
    my $prog = lc($form->param('prog')) if $form->param('prog');
    my $exon_mask = $form->param('exon_mask');

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
    my $box = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
    #generate sequence submission selector
    $template->param(SEQ_SELECT=>1);
    $template->param(SEQ_SELECT_LOOP=>\@seq_sub);
    my $seq_submission = $template->output;
    $template->param(SEQ_SELECT=>0);

    #generate the hsp color option
    my $hsp_colors = gen_hsp_colors($num_seqs);
    my $html;
    my $spike_len = spike_filter_select();
    $template->param(GEN_REFERENCE_SEQUENCES=>1);
    $template->param(REF_SEQ=>\@seq_nums);
    my $ref_seqs = $template->output;
    $template->param(GEN_REFERENCE_SEQUENCES=>0);

    $template->param(SPIKE_LEN=>$spike_len);
    $template->param(SEQ_RETRIEVAL=>1);
    $template->param(NUM_SEQS=>$num_seqs);
    $message .= "<BR/>" if $message;
    $template->param(MESSAGE=>$message);
    $template->param(SEQ_SUB=>$seq_submission);
    $template->param(HSP_COLOR=>$hsp_colors);
    $template->param(GO_RUN=>gen_go_run($num_seqs));
    $template->param(AUTOSEARCH=>$autosearch_string);
    $box->param(BOX_NAME=>"Options:");
#    $box->param(BOX_NAME=>"Sequence Retrieval:");
#    $box->param(BODY=>$template->output);
    
#    $html .= $box->output;
#    $template->param(SEQ_RETRIEVAL=>0);
    $template->param(OPTIONS=>1);
    if ($exon_mask)
      {
	$template->param(EXON_MASK_ON=>"checked");
      }
    else
      {
	$template->param(EXON_MASK_OFF=>"checked");
      }
    $template->param(REFERENCE_SEQUENCES=>$ref_seqs);
    $template->param(ALIGNMENT_PROGRAMS=>algorithm_list($prog));
    $box->param(BODY=>$template->output);
    $html .= $box->output;
    return $html;
  }

sub run
  {
    my %opts = @_;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my $spike_len = $opts{spike};
    my $mask_cds_flag = $opts{maskcds};
    my $mask_ncs_flag = $opts{maskncs};
    my $iw = $opts{iw};
    my $ih = $opts{ih};
    my $feat_h = $opts{fh};
    my $show_gc = $opts{gc};
    my $show_nt = $opts{nt};
    my $color_hsp = $opts{colorhsp};
    my $hsp_label = $opts{hsplabel};
    my $overlap_adjustment = $opts{overlap};
    my $hiqual = $opts{hiqual};
    my $hsp_limit = $opts{hsplim};
    my $hsp_limit_num = $opts{hsplimnum};
    my $show_hsps_with_stop_codon = $opts{showallhsps};
    my $padding = $opts{padding};
    my $viewer = $opts{viewer};
    my ($analysis_program, $param_string, $parser_opts) = get_algorithm_options(%opts);
    my $message;
    initialize_basefile();

    my @hsp_colors;
    for (my $i = 1; $i <= num_colors($num_seqs); $i++)
      {
	my $r = $opts{"r$i"};
	my $g = $opts{"g$i"};
	my $b = $opts{"b$i"};
	my @tmp;
	foreach my $color ($r, $g, $b)
	  {
	    $color = 0 unless $color =~ /^\d+$/;
	    $color = 0 if $color < 0;
	    $color = 255 if $color > 255;
	    push @tmp, $color;
	  }
	push @hsp_colors,\@tmp;
      }


#    print Dumper \@_;

    my $stagger_label = $hsp_label =~ /staggered/i ? 1 : 0;
    my $feature_labels = $hsp_label eq "0" ? 0 : 1;
    my $form = $FORM;

    my @sets;
    my $html;
    my $t1 = new Benchmark;
    my $spike_seq;
    for (my $i = 1; $i <= $num_seqs; $i++)
      {
	my $accn = $opts{"draccn$i"};
	my $featid = $opts{"featid$i"};
	my $drup = $opts{"drup$i"};
	my $drdown = $opts{"drdown$i"};
	my $drrev = $opts{"drrev$i"};

	my $gbaccn = $opts{"gbaccn$i"};
	my $gbstart = $opts{"gbstart$i"};
	my $gblength = $opts{"gblength$i"};
	my $gbrev = $opts{"gbrev$i"};

	my $dirseq = $opts{"dirseq$i"};
	my $dirrev = $opts{"dirrev$i"};
	my $dirstart = $opts{"dirstart$i"};
	my $dirstop = $opts{"dirstop$i"};
	my $rev = 0;
	my ($up, $down);
	my ($file, $file_begin, $file_end, $obj);
	my $reference_seq =$opts{"ref_seq$i"};
	
	if ($featid)
	  {
	    $obj = get_obj_from_genome_db( $accn, $featid, $drrev, $drup, $drdown );
	    if ($obj)
	      {
		($file, $file_begin, $file_end,$spike_seq) = 
		  generate_seq_file(obj=>$obj,
				    mask_cds=>$mask_cds_flag,
				    mask_ncs=>$mask_ncs_flag,
				    spike_len=>$spike_len,
				    seq_num=>$i,
				   );
		$up = $drup;
		$down = $drdown;
		$rev = 1 if $drrev;
	      }
	    else
	      {
		$message .=  "No entry found for $featid.\n";
	      }

	  }
 	elsif ($dirseq )
 	  {
	    $dirseq = CoGeX::Feature->reverse_complement($dirseq) if $dirrev;
 	    my $seq = get_substr(seq=>$dirseq, start=>$dirstart, stop=>$dirstop);
	    $dirstop = length($seq) unless $dirstop;
 	    ($obj) = generate_obj_from_seq($seq, $i);
 	    
 	    if ($obj)
	      {
		($file, $file_begin, $file_end, $spike_seq) = 
		  generate_seq_file (
				     obj=>$obj,
				     mask_cds=>$mask_cds_flag,
				     mask_ncs=>$mask_ncs_flag,
				     spike_len=>$spike_len, 
				     seq_num=>$i,
				    );
		$up = $dirstart;
		$down = $dirstop;
		$rev = 1 if $dirrev;
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
 	    my $res = $obj->get_genbank_from_nbci($gbaccn, $gbrev);	    
	    if ($obj)
	      {
		$obj->sequence(CoGeX::Feature->reverse_complement($obj->sequence)) if $gbrev;
		my $seq;
		($file, $file_begin, $file_end,$spike_seq, $seq) = 
		  generate_seq_file (
				     obj=>$obj,
				     mask_cds=>$mask_cds_flag,
				     mask_ncs=>$mask_ncs_flag,
				     startpos=>$gbstart,
				     downstream=>$gblength,
				     spike_len=>$spike_len,
				     seq_num=>$i,
				    );
		$obj->sequence($seq);
		$rev = 1 if ($gbrev);
		$up = $gbstart;
		$down = $gblength;
	      }
	    else
	      {
		$message .= "No GenBank entry found for $gbaccn\n";
	      }
 	  }
	if ($obj)
	  {
	    push @sets, {
			 obj=>$obj,
			 file=>$file,
			 file_begin=>$file_begin,
			 file_end=>$file_end,
			 accn=>$obj->accn,
			 rev=>$rev,
			 up=>$up,
			 down=>$down,
			 spike_seq=>$spike_seq,
			 reference_seq=>$reference_seq,
			 seq_num=>$i,
			};
	  }
      }

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
    elsif ($analysis_program eq "DIALIGN")
      {
	($analysis_reports, $analysis_program,$message) = run_dialign (sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts);
      }
    else
      {
	$analysis_reports = run_bl2seq(sets=>\@sets, params=>$param_string, parser_opts=>$parser_opts, blast_program=>$analysis_program, spike_seq=>$spike_seq);
      }
    $analysis_reports = [] unless ref($analysis_reports) =~ /ARRAY/i;
   #sets => array or data for blast
   #blast_reports => array of arrays (report, accn1, accn2, parsed data)
#    print STDERR Dumper $analysis_reports;
    my $t3 = new Benchmark;
    my $count = 1;
    my $frame_height;
    my $html_viewer;
    foreach my $item (@sets)
      {
	my $accn = $item->{accn};
	my $obj = $item->{obj};
	my $file = $item->{file};
	my $file_begin = $item->{file_begin};
	my $file_end = $item->{file_end};
	my $rev = $item->{rev};
	my $up = $item->{up};
	my $down = $item->{down};
	$down = $obj->seq_length unless $down;
	my $spike_seq = $item->{spike_seq};
	if ($obj)
	  {
	    my ($image, $map, $mapname, $gfx, $eval_cutoff) = generate_image(
									     gbobj=>$obj, 
									     start=>$file_begin,
									     stop => $file_end,
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
									     seq_num=>$count,
									     reverse_image=>$rev,
									    );
	    $frame_height += $gfx->ih + $gfx->ih*.1;
	    $html_viewer .= qq!<div>$accn!;
	    $html_viewer .= qq!(<font class=species>!.$obj->organism.qq!</font>)! if $obj->organism;
	    $html_viewer .= "(".$up."::".$down.")" if defined $up;
	    $html_viewer .= qq!<font class=small> Reverse Complement</font>! if $rev;
#	    $html_viewer .= qq!<font class=small> (eval cutoff: $eval_cutoff)</font>! if defined $eval_cutoff;
	    $html_viewer .= qq!</DIV>\n!;
	    $html_viewer .= qq!<IMG SRC="$TEMPURL/$image" !;
	    my $tmp_count = $count -1;
	    $html_viewer .= qq!BORDER=0 ismap usemap="#$mapname">\n!;
	    $html_viewer .= "$map\n";
	    $item->{image} = $image;
	    $item->{gfx} = $gfx;
	  }
	$count++;
      }
    $count--;

    my $str = qq{/bpederse/gobe?imgdir=/CoGe/tmp/&img=$BASEFILENAME&n=}.($count-1);
    my $t4 = new Benchmark;
    $html_viewer .= qq!<br/>!;
    $html_viewer .= qq!<FORM NAME=\"info\">\n!;
    $html_viewer .= $form->br();
    $html_viewer .= "Information: ";
    $html_viewer .= $form->br();
    $html_viewer .= qq!<DIV id="info"></DIV>!;
    $html_viewer .= "</FORM>\n";
    my $flash_viewer = qq{<DIV id=flash_viewer></DIV>};
#    $html .= $viewer eq "flash" ? $flash_viewer : $html_viewer;
    $html .= qq{<div id=html_viewer>$html_viewer</div>};
    $html .= $flash_viewer;
    $html .= qq{<table>};
    $html .= qq{<tr valign=top><td class = small>Alignment reports};
    if ($analysis_reports && @$analysis_reports)
      {
	foreach my $item (@$analysis_reports)
	  {
	    my $report = $item->[0];
	    my $accn1 = $item->[1];
	    my $accn2 = $item->[2];
	    my $basereportname = basename( $report );
	    $basereportname = $TEMPURL . "/$basereportname\n";
	    $html .= "<div><font class=xsmall><A HREF=\"$basereportname\">View alignment output for $accn1 versus $accn2</A></font></DIV>\n";
	  }
      }
    else
      {
	$html .= "<div class=xsmall>No alignment reports were generated</DIV>\n";
      }
    $html .= qq{<td class = small>Fasta files};
	foreach my $item (@sets)
	  {
	    my $basename = $TEMPURL."/".basename ($item->{file});
	    my $accn = $item->{accn};
	    $html .= "<div><font class=xsmall><A HREF=\"$basename\">Fasta file for $accn</A></font></DIV>\n";
	  }
    $html .= qq{<td class = small><a href = "http://baboon.math.berkeley.edu/mavid/gaf.html">GAF</a> annotation files};
    foreach my $item (@sets)
      {
	my $basename = $TEMPURL."/".basename (generate_annotation(%$item));
	my $accn = $item->{accn};
	$html .= "<div><font class=xsmall><A HREF=\"$basename\">Annotation file for $accn</A></font></DIV>\n";
      }
    $html .= qq{<td class = small>SQLite db};
    my $dbname = generate_imagemap_db(\@sets);
    $html .= "<div class=xsmall><A HREF=\"$dbname\">SQLite DB file</A></DIV>\n";
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
    my $image_time = timestr(timediff($t4,$t3));
    my $html_time = timestr(timediff($t5,$t4));
    print STDERR qq{
GEvo Benchmark: $DATE
Time to get DB info             : $db_time
Time to run $analysis_program   : $blast_time
Time to generate images and maps: $image_time
Time to process html            : $html_time
} if $BENCHMARK;

    return $outhtml, $iw+400, $frame_height, $BASEFILENAME,$count,$message;

}


sub generate_image
  {
    my %opts = @_;
    my $gbobj = $opts{gbobj};
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
    my $seq_num = $opts{seq_num};
    my $graphic = new CoGe::Graphics;
    my $gfx = new CoGe::Graphics::Chromosome;
    $gfx->overlap_adjustment(1);
    $gfx->skip_duplicate_features(1);
    $gfx->DEBUG(0);
#    print STDERR "$start"."::"."$stop"." -- ".length($gbobj->sequence)."\n";
    $graphic->initialize_c (
			    c=>$gfx,
			    iw=>$iw,
			    start=> 1,
			    stop => length($gbobj->sequence),
			    draw_chr=>1,
			    draw_ruler=>1,
			    draw_chr_end=>0,
			    chr_start_height=>$ih,
			    chr_mag_height=>5,
			    feature_start_height=>$fh,
			    mag=>0,
			    mag_off=>1,
			    chr_length => length($gbobj->sequence),
			    feature_labels=>1,
			    fill_labels=>1,
			    forcefit=>1,
			    minor_tick_labels=>1,
#			    overlap_adjustment=>$overlap_adjustment,
			    feature_labels=>$feature_labels,
			    draw_hi_qual=>$hiqual,
			    padding=>$padding,
			   );
    $gfx->major_tick_labels(0);
    my $f1= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => 1});
    $f1->merge_percent(0);
    $gfx->add_feature($f1);
    my $f2= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => -1});
    $f2->merge_percent(0);
    $gfx->add_feature($f2);
    $graphic->process_nucleotides(c=>$gfx, seq=>$gbobj->sequence, layers=>{gc=>$show_gc, nt=>$show_nt});
#    print STDERR join ("\t", $gbobj->accn, $start, $stop),"\n";
    process_features(c=>$gfx, obj=>$gbobj, start=>$start, stop=>$stop, overlap_adjustment=>$overlap_adjustment);
    $eval_cutoff = process_hsps(
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
				eval_cutoff=>$eval_cutoff,
				color_hsp=>$color_hsp,
				colors=>$hsp_colors,
				show_hsps_with_stop_codon=>$show_hsps_with_stop_codon,
			       );
#    $gfx->DEBUG(1);
#    print STDERR Dumper $gfx;
    my $filename = $BASEFILE."_".$seq_num.".png";
    $filename = check_filename_taint($filename);
    $gfx->generate_png(file=>$filename);
    my $mapname = "map_".$seq_num;
    my ($map)=$gfx->generate_imagemap(name=>$mapname);
    return (basename($filename), $map, $mapname, $gfx, $eval_cutoff);
  }

sub generate_imagemap_db
  {
    my ($sets) = @_;
    my $tempfile = $BASEFILE.".sqlite";
    $tempfile = $TEMPDIR."/".$tempfile unless $tempfile =~ /$TEMPDIR/;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$tempfile","","");
    my $create = qq{
CREATE TABLE image_data
(
id integer(10) PRIMARY KEY,
name varchar(50),
xmin integer(10),
xmax integer(10),
ymin integer(10),
ymax integer(10),
image varchar(50),
image_track varchar(10),
pair_id integer(10),
link varchar(50),
annotation blob,
color varchar(10)
)
};
    $dbh->do($create);
     my $index = qq{
 CREATE INDEX name ON image_data (name)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX image ON image_data (image)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX xmin ON image_data (xmin)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX xmax ON image_data (xmax)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX ymin ON image_data (ymin)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX ymax ON image_data (ymax)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX image_track ON image_data (image_track)
 };
     $dbh->do($index);
     $index = qq{
 CREATE INDEX pair_id ON image_data (pair_id)
 };
     $dbh->do($index);
    $create = qq{
CREATE TABLE image_info
(
id integer(10) PRIMARY KEY,
iname varchar(50),
title varchar(1024)
)
};
    $dbh->do($create);
     $index = qq{
 CREATE INDEX iname ON image_info (iname)
 };
     $dbh->do($index);
    my $i = 1;
    my $j = 1;
#    my $sth_insert = $dbh->prepare("INSERT INTO image_data VALUES(?,?,?,?,?,?,?,?,?,?)");
    foreach my $item (@$sets)
      {
	my $image = $item->{image};
	my $gfx = $item->{gfx};
	my $accn = $item->{accn};
	my $title;
	$title = $accn;
	$title .= " ".$item->{obj}->organism() if $item->{obj}->organism();
	$title .= "(".$item->{up}."::".$item->{down}.")" if defined $item->{up};
	$title .= qq! Reverse Complement! if $item->{rev};
	my $statement = qq{
INSERT INTO image_info (id, iname, title) values ($j, "$image", "$title")
};
	print STDERR $statement unless $dbh->do($statement);
	$j++;
	foreach my $feat ($gfx->get_feats)
	  {
	    next if $feat->fill; #skip backgroup images;
	    next unless $feat->image_coordinates;
	    my $pair_id = "NULL";
	    my $coords = $feat->image_coordinates;
	    $coords =~ s/\s//g;
	    #next unless $feat->type =~ /HSP/i;
	    next if $feat->type eq "unknown";
	    my $name = $feat->type =~ /HSP/i ? $feat->alt : $feat->label;
	    $name = $feat->type unless $name;
	    my $query = qq{select id from image_data where name = "$name"};
	    my $sth = $dbh->prepare($query);
	    $sth->execute;
	    my $res = $sth->fetchrow_array();
	    my $color = "NULL";
	    if ($res && $feat->type =~ /HSP/)
	      {
		$color ="#";
		$pair_id = $res;
		my $statement = "update image_data set pair_id = $i where id = $res";
		$dbh->do($statement);
		foreach my $c (@{$feat->color})
		  {
		    $color .= sprintf("%X",$c);
		  }
	      }
	    #generate link
	    my $link = $feat->link;
	    $link =~ s/'//g;
	    #generate image track
	    my $image_track = $feat->track;
	    $image_track = "-".$image_track if $feat->strand =~ /-/;

	    my ($xmin, $ymin, $xmax, $ymax) = split /,/, $coords;
	    my $anno = $feat->description;
#	    $anno =~ s/'|"//g;
	    $anno =~ s/<br\/?>/&#10;/ig;
	    $anno =~ s/\n/&#10;/g;
#	    print STDERR $anno if $anno =~ /Location/;
#	    print STDERR $anno,"\n" if $anno =~ /01020/;
	    $statement = qq{
INSERT INTO image_data (id, name, xmin, xmax, ymin, ymax, image, image_track,pair_id, link, annotation, color) values ($i, "$name", $xmin, $xmax, $ymin, $ymax, "$image", "$image_track",$pair_id, '$link', '$anno', '$color')
};
	    print STDERR $statement unless $dbh->do($statement);
#	    $sth_insert->execute($i,$name, $xmin, $xmax, $ymin, $ymax, $image, $pair_id, $link, $anno);
	    $i++;
	  }
      }
#     my $sth = $dbh->prepare("SELECT * from image_data");
#     $sth->execute();
#     while(my $item = $sth->fetchrow_arrayref)
#       {
# 	print STDERR Dumper $item;
#       }
    system "chmod +rw $tempfile";
    return $tempfile;
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
    my $accn = $obj->accn;
    my $track = 1;
    my @opts = ($start, $stop) if $start && $stop;
    unless (ref $obj)
      {
	warn "Possible problem with the object in process_features.  Returning";
	return 0;
      }
    foreach my $feat($obj->get_features(@opts))
      {

        my $f;
	my $type = $feat->type;
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} @{$feat->qualifiers->{names}} if ref ($feat->qualifiers) =~ /hash/i ;
        if ($type =~ /pseudogene/i)
          {
#	    next;
	    $f = CoGe::Graphics::Feature::Gene->new();
	    $f->color([255,33,0,50]);
	    $f->order($track);
	    $f->overlay(1);
	    $f->mag(0.5);
          }
        elsif ($type =~ /Gene/i)
          {
#	    next;
	    $f = CoGe::Graphics::Feature::Gene->new();
	    $f->color([255,0,0,50]);
	    $f->order($track);
	    $f->overlay(1);
	    $f->mag(0.5);
          }
        elsif ($type =~ /CDS/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,255,0, 50]);
        	$f->order($track);
		$f->overlay(3);
		if ($accn)
		  {
		    foreach my $name (@{$feat->qualifiers->{names}})
		      {
			$f->color([255,255,0]) if $name =~ /^$accn/i;
			$f->label($name) if $name =~ /^$accn/i;
		      }
		  }

		
          }
        elsif ($type =~ /mrna/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,0,255, 50]);
        	$f->order($track);
		$f->overlay(2);
		$f->mag(0.75);
          }
        elsif ($type =~ /rna/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([200,200,200, 50]);
        	$f->order($track);		$f->overlay(2);
		if ($accn)
		  {
		    foreach my $name (@{$feat->qualifiers->{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
		      }
		  }

          }

        next unless $f;
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
	$f->link("FeatView.pl?accn=$name") if $name;
	$f->skip_overlap_search($overlap);
        $c->add_feature($f);
    }
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
    my $color_hsp = $opts{color_hsp};
    my $colors = $opts{colors};
    my $show_hsps_with_stop_codon = $opts{show_hsps_with_stop_codon};
    #to reverse hsps when using genomic sequences from CoGe, they need to be drawn on the opposite strand than where blast reports them.  This is because CoGe::Graphics has the option of reverse drawing a region.  However, the sequence fed into blast has already been reverse complemented so that the HSPs are in the correct orientation for the image.  Thus, if the image is reverse, they are drawn on the wrong strand.  This corrects for that problem.   Sorry for the convoluted logic, but it was the simplest way to substantiate this option
    my $i = 0;
    my $track = 2;
    my @feats;
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
	if ($spike_seq)
	  {
	    foreach my $hsp (@{$blast->hsps})
	      {
		if ($hsp->qalign =~ /^$spike_seq$/i)
		  {
		    $eval_cutoff = $hsp->eval;
		    last;
		  }
	      }
	  }
	print STDERR "\t",$blast->query," ", $blast->subject,"\n" if $DEBUG;
	foreach my $hsp (@{$blast->hsps})
	  {
	    next if defined $eval_cutoff && $hsp->eval > $eval_cutoff;
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
	    my $f = CoGe::Graphics::Feature::HSP->new({start=>$start, stop=>$stop});
	    $color = [100,100,100] if $spike_seq && $hsp->qalign =~ /$spike_seq$/i;
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
	    my $desc = join ("<br/>", "HSP: ".$hsp->number. qq{  <span class="small">(}.$blast->query."-". $blast->subject.")</span>", "Location: ".$start."-".$stop." (".$hsp->strand.")", "Sequence: ".$seq,"Match: ".$hsp->match,"Length: ".$hsp->length,"Identity: ".$hsp->percent_id);
	    $desc .= "<br/>E_val: ".$hsp->pval if $hsp->pval;
	    $desc .= "<br/>Score: ".$hsp->score if $hsp->score;
	    $desc .= qq{<span class="small"> (cutoff: $eval_cutoff)</span>} if defined $eval_cutoff;
	    $f->description($desc);
	    if ($reverse)
	      {
		$strand = $strand =~ /-/ ? "1" : "-1";
		my $tmp = $seq_len - $start+1;
		$start = $seq_len - $stop+1;
		$stop = $tmp;
	      }

	    my $link = "HSPView.pl?report=$report&num=".$hsp->number."&db=".$BASEFILENAME.".sqlite";
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
    return $eval_cutoff;
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
    
    my ($obj, $file, $file_begin, $file_end, $spike_seq);
    $obj = new CoGe::Accessory::GenBank;
    if ($sequence =~ /^LOCUS/)
      {
	#genbank sequence
	$obj->parse_genbank($sequence);
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
	$obj->accn("RAW_SEQUENCE_SUBMISSION $num");
	$obj->locus("RAW_SEQUENCE_SUBMISSION $num");
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
    my %options = @_;
    my $obj = $options{obj} || $options{gbobj};
    my $start = $options{start} || $options{startpos} || 1;
    my $up = $options{up} || $options{upstream} || 0;
    my $down = $options{down} || $options{downstream} || 0;
    my $spike_len = $options{spike_len} || 0;
    my $mask = $options{mask_cds};
    my $mask_ncs = $options{mask_ncs};
    my $seq_num = $options{seq_num};

    my $t1 = new Benchmark;
    my ($file, $file_begin, $file_end, $spike_seq, $seq) = 
      write_fasta(
		  obj=>$obj,
		  accn=>$obj->accn,
		  mask=>$mask,
		  mask_ncs=>$mask_ncs,
		  start=>$start,
		  upstream=>$up,
		  downstream=>$down,
		  spike=>$spike_len,
		  seq_num=>$seq_num,
		 );
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2,$t1));
    print STDERR "Time to generate Sequence file:   $time\n" if $BENCHMARK;
    return ($file, $file_begin, $file_end, $spike_seq, $seq);
  }

sub get_obj_from_genome_db
  {
    my $accn = shift;
    my $featid = shift;
    my $rev = shift;
    my $up = shift || 0;
    my $down = shift || 0;
    my $t1 = new Benchmark;
    my ($feat) = $coge->resultset('Feature')->esearch({"me.feature_id"=>$featid})->next;
    my $t2 = new Benchmark;
    my $start = 0;
    my $stop = 0;
    my $chr;
    my $seq;
    my $ds_id = $feat->dataset->id;
    if ($feat)
      {
	
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

	$seq = $coge->resultset('Dataset')->find($ds_id)->get_genome_sequence(
									      start => $start,
									      stop => $stop,
									      chr => $chr,
									     );
	$seq = CoGeX::Feature->reverse_complement($seq) if $rev;
      }
    my $t3 = new Benchmark;

    my $obj= new CoGe::Accessory::GenBank({
					   accn=>$accn,
					   locus=>$accn,
					   version=>$feat->dataset->version(),
					   data_source=>$feat->dataset->datasource->name(),
					   dataset=>$ds_id,
					   chromosome=>$chr,
					   start=>$start,
					   stop=>$stop,
					   organism=>$feat->org->name()."(v".$feat->dataset->version.")",
					   seq_length=>length($seq),
					   sequence=>$seq,
					 });

    my $fnum = 1;
    my %used_names;
    $used_names{$accn} = 1;
    my $t4 = new Benchmark;

    print STDERR "Region: $chr: $start-$stop\n" if $DEBUG;
#    print STDERR "Region: $chr: ",$start-$start+1,"-",$stop-$start,"\n";
    my @feats = $coge->get_features_in_region(start=>$start, stop=>$stop, chr=>$chr, dataset_id=>$ds_id);
    my $t5 = new Benchmark;
    foreach my $f (@feats)
      {

	my $name;
	my @names = $f->names;
	foreach my $tmp (@names)
	  {
	    $name = $tmp;
	    last if ($tmp =~ /^$accn.+/i);
	  }
	unless (@names)
	  {
	    print STDERR "No Name: ",$f->id,"\n";;
	  }
	$name = $accn unless $name;
	print STDERR $name,"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(),"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(recalibrate=>$start),"\n\n" if $DEBUG;
	my $anno = $f->annotation_pretty_print_html;
	$anno =~ s/\n/<br\/>/ig;
	my $location = $f->genbank_location_string(recalibrate=>$start);
	$location = $obj->reverse_genbank_location(loc=>$location, ) if $rev;
	print STDERR $name, "\t",$f->type->name ,"\t",$location,"\n" if $DEBUG;
	$obj->add_feature (
			   number=>$fnum,
			   type=>$f->type->name,
			   location=> $location,
			   qualifiers=>{
                                        names=> [@names],
                                       },
			   annotation=>$anno,
			  );
	$fnum++;
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
Region:         ds_id: $ds_id $start-$stop($chr)
} if $BENCHMARK;
    return $obj;
  }

sub check_filename_taint {
	my $v = shift;
	if ($v =~ /^([A-Za-z0-9\-\.=\/_]*)$/) {
		my $v1 = $1;
		return($v1);
	} else {
		return(0);
	}
}

sub check_taint {
	my $v = shift;
	if ($v =~ /^([-\w._=\s+\/]+)$/) {
			$v = $1;
			# $v now untainted
			return(1,$v);
	} else {
	# data should be thrown out
	  carp "$v failed taint check\n";
			return(0);
	}
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
  for (my $i=0; $i<scalar @$sets; $i++)
    {
      for (my $j=0; $j<scalar @$sets; $j++)
	{
	  next unless $j > $i;
	  #
	  my $seqfile1 = $sets->[$i]->{file};
	  my $seqfile2 = $sets->[$j]->{file};
	  check_sequence_files_spike($spike_seq, $seqfile1, $seqfile2) if ($spike_seq);
	  next unless -r $seqfile1 && -r $seqfile2; #make sure these files exist
	  
	  next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
	  my ($accn1, $accn2) = ($sets->[$i]{accn}, $sets->[$j]{accn});
	  my $command = $BL2SEQ;
	  
	  
	# need to create a temp filename here
	  my ($tempfile) = $BASEFILE."_".($i+1)."-".($j+1).".bl2seq";;
	  
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
	  print STDERR $command,"\n";
	  `$command`;
	  system "chmod +rw $tempfile";
	  my $blastreport = new CoGe::Accessory::bl2seq_report({file=>$tempfile}) if -r $tempfile;
	  my @tmp = ($tempfile, $accn1, $accn2);
	  if ($blastreport)
	    {
	      $blastreport->eval_cutoff($eval_cutoff);
	      push @tmp, $blastreport
	    }
	  else
	    {
	      push @tmp, "no results from blasting $accn1 and $accn2";
	    }
	  push @reports, \@tmp;
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
	    my ($tempfile) = $BASEFILE."_".($i+1)."-".($j+1).".blastz";
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
#	    print STDERR $command,"\n";
	    `$command`;
	    system "chmod +rw $tempfile";
	    my $blastreport = new CoGe::Accessory::blastz_report({file=>$tempfile}) if -r $tempfile;
	    my @tmp = ($tempfile, $accn1, $accn2);
	    if ($blastreport)
	    {
	      push @tmp, $blastreport
	    }
	  else
	    {
	      push @tmp, "no results from blasting $accn1 and $accn2";
	    }
	  push @reports, \@tmp;
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
	    my ($tempfile) = $BASEFILE."_".($i+1)."-".($j+1).".lagan";
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
	    my ($tempfile) = $BASEFILE."_".($i+1)."-".($j+1).".chaos";
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
	    print STDERR $command,"\n";
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
	      push @tmp, "no results from comparing $accn1 and $accn2 with LAGAN";
	    }
	  push @reports, \@tmp;
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
									base_name=>$BASEFILENAME,
									extension=>$parser_opts->{anchor_params}{prog},
									run_anchor=>$anchor_prog,
									run_dialign_opts=>$dialign_opts,
									run_anchor_opts=>$anchor_opts,
									dialign_report_opts=>$parser_opts,
									anchor_report_opts=>$parser_opts->{anchor_params}{parser_opts},
#									DEBUG=>1,
								       });
		$tempfile = $obj->dialign_file;
		$report = $obj->dialign_report;
	      }
	    else
	      {
		my $seqfile = $TEMPDIR."/".$BASEFILENAME."_".($i+1)."-".($j+1).".fasta";
		next unless -r $seqfile1 && -r $seqfile2; #make sure these files exist
		#put two fasta files into one for dialign
		`cat $seqfile1 > $seqfile`;
		`cat $seqfile2 >> $seqfile`;
		next unless $sets->[$i]{reference_seq} || $sets->[$j]{reference_seq};
		($tempfile) = $BASEFILE."_".($i+1)."-".($j+1).".dialign";
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
#		print STDERR $command,"\n";
		#time for execution
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
    my $start = $opts{start};
    my $mask = $opts{mask};
    my $mask_ncs = $opts{mask_ncs};
    my $downstream = $opts{downstream};
    my $upstream = $opts{upstream};
    my $spike = $opts{spike};
    my $seq_num = $opts{seq_num};
    # vars
    my($seq,$seq_begin,$seq_end,$db_begin,$db_end,$chr);
    my($gene_end,$gene_begin);
    my $fullname = $BASEFILE."_".$seq_num.".fa";
    my $hdr = $gbobj->get_headerfasta( );
    $seq_begin = $start - $upstream;
    if ( $seq_begin < 1 ) {
      $seq_begin = 1;
    }
    ($seq) = uc($gbobj->sequence());

    if ($downstream)
      {
	$seq_end = $start + $downstream;
      } else {
	$seq_end = length($seq);
      }
    $seq_end = length($seq) if ( $seq_end > length( $seq ));
    # mask out exons if required
    $seq = $gbobj->mask_exons( $seq ) if ( $mask );
    # mask out non coding sequences if required
    $seq = $gbobj->mask_ncs( $seq ) if ( $mask_ncs );
    ($seq) = $gbobj->subsequence( $seq_begin, $seq_end, $seq );
    my $full_seq = $seq;
    my $spike_seq = "";
    if ($spike)
      {
	($spike_seq) = generate_spike_seq($spike);
	$seq .= $spike_seq;
      }
    ($fullname) = check_filename_taint( $fullname );
    my $length = length($seq);
    open(OUT, ">$fullname") or die "Couldn't open $fullname!\n";
    print OUT "$hdr\n";
    my $max = 100;
    my $i = 0;
    while ($i < $length && length($seq) > $max)
      {
	print OUT substr ($seq, 0, $max),"\n";;
	substr ($seq, 0, $max) = "";
	$i+=$max;
      }
    print OUT $seq,"\n" if $seq;
    close(OUT);
    system "chmod +rw $fullname";
    return($fullname,$seq_begin,$seq_end,$spike_seq, $full_seq);
  }

sub generate_spike_seq { 
	my $spikesize = shift;
	return unless $spikesize;
	my %nt = (1=>'A', 2=>'T', 3=>'C', 4=>'G');
	my $spike_seq;
	my $idx = 1;
	for (my $i =1; $i<=$spikesize; $i+=3)
	  {
	    $spike_seq .= $nt{$idx}x3;
	    $idx++;
	    $idx = 1 if $idx == 5;
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
    my $start = $opts{file_begin};
    my $stop = $opts{file_end};
    my $rev = $opts{rev};
    my $seq_num = $opts{seq_num};
    my @opts = ($start, $stop);
    my $fullname = $BASEFILE."_".$seq_num.".anno";
    my %data;
    my $length = length($obj->sequence());
    foreach my $feat($obj->get_features(@opts))
      {
	my $type = $feat->type;
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} @{$feat->qualifiers->{names}} if ref($feat->qualifiers) =~ /hash/i;
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
	    push @{$data{$name}{$type}},[$start, $stop, $dir];
	  }
      }
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

sub gen_go_run
  {
    my $num_seqs = shift || $NUM_SEQS;
    my $params;
    for (my $i = 1; $i <=$num_seqs; $i++)
      {
	$params .= qq{'args__draccn$i', 'accn$i',};
	$params .= qq{'args__featid$i', 'featid$i',};
	$params .= qq{'args__drup$i', 'drup$i',};
	$params .= qq{'args__drdown$i', 'drdown$i',};
	$params .= qq{'args__drrev$i', 'drrev$i',};
	$params .= qq{'args__gbaccn$i', 'gbaccn$i',};
	$params .= qq{'args__gbstart$i', 'gbstart$i',};
	$params .= qq{'args__gblength$i', 'gblength$i',};
	$params .= qq{'args__gbrev$i', 'gbrev$i',};
	$params .= qq{'args__dirseq$i', 'dirseq$i',};
	$params .= qq{'args__dirrev$i', 'dirrev$i',};
	$params .= qq{'args__dirstart$i', 'dirstart$i',};
	$params .= qq{'args__dirstop$i', 'dirstop$i',};
	$params .= qq{'args__ref_seq$i', 'ref_seq$i',};
      }
    for (my $i = 1; $i <=num_colors($num_seqs); $i++)
      {
	$params .= qq{'args__r$i', 'r$i',};
	$params .= qq{'args__g$i', 'g$i',};
	$params .= qq{'args__b$i', 'b$i',};
      }

# old params
#	'args__matchfilter', 'matchfilter', 

    $params .= qq{
	'args__spike', 'spike', 
	'args__maskcds', 'mask_cds', 
	'args__maskncs', 'mask_ncs', 

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

        'args__viewer', 'viewer',
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
};
    $params =~ s/\n//g;
    $params =~ s/\s+/ /g;
#    my $run = "<SCRIPT language=\"JavaScript\">alert('hi');</script>";
    my $run = qq!<SCRIPT language=\"JavaScript\">function go_run (){ setTimeout("run([$params],[handle_results], 'POST')",500)}</script>!;
#
#    my $button = qq{<a href="#" onMouseOver="activeButton('go')" onMouseOut="inactiveButton('go')" onMouseDown="pressedButton('go')" onMouseUp="inactiveButton('go')" onClick="loading([],['results']); setTimeout('go_run()', 1000);"><img name="go" src="/CoGe/picts/buttons/GEvo/go/inactive.png" width="47" height="27" border="0" alt="javascript button"></a>};
    return $run;
  }

sub gen_hsp_colors
  {
    my $num_seqs = shift || $NUM_SEQS;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
    my $hsp_colors;
    my @colors = color_pallet(num_seqs=>$num_seqs);
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
    my $offset = $opts{offset} || 75;
    my $num_seqs = $opts{num_seqs} || $NUM_SEQS;
    my @colors;
    push @colors, {
		   HSP_NUM=>1,
		   RED=>$start->[0],
		   GREEN=>$start->[1],
		   BLUE=>$start->[2],
		  };
    for (my $i = 2; $i <= num_colors($num_seqs); $i++)
      {
	my @color;
	unless ($i%5)
	  {
	    $start = [$start->[1], $start->[2], $start->[0]];
	  }

	foreach my $color (@$start)
	  {
	    $color -= $offset;
	    $color += 255 if $color < 0;
	    push @color, $color;
	  }
	
	push @colors, {
		       HSP_NUM=>$i,
		       RED=>$color[0],
		       GREEN=>$color[1],
		       BLUE=>$color[2],
		      };
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
    my @programs = sort {lc $a cmp lc $b} qw(blastn blastz CHAOS DIALIGN LAGAN tblastx);
    my $html;
    foreach my $prog (@programs)
      {
	$html .= "<option";
	$html .= $program && $program =~ /$prog/i ? " selected>" : ">";
	$html .= $prog."</option>\n";
      }
    return $html;
  }

sub initialize_basefile
  {
    my $file = new File::Temp ( TEMPLATE=>'GEvo_XXXXXXXX',
				   DIR=>$TEMPDIR,
				    #SUFFIX=>'.png',
				    UNLINK=>1);
    ($BASEFILE)= $file->filename;
    
    ($BASEFILENAME) = $file->filename =~ /([^\/]*$)/;
  }

sub add_seq
  {
    my $num_seq = shift;
    $num_seq ++;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/GEvo.tmpl');
    my $hsp_colors;
    if ($num_seq > $MAX_SEQS)
      {
	$hsp_colors = gen_hsp_colors($MAX_SEQS);
	my $go_run = gen_go_run($MAX_SEQS);
	return ('',$MAX_SEQS, $go_run, '', $hsp_colors,qq{Exceeded max number of sequences ($MAX_SEQS)});
      }
    my @seqs = {
		SEQ_NUM=>$num_seq,
		REV_NO=>"checked",
		DRUP=>10000,
		DRDOWN=>10000,
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
    $hsp_colors = gen_hsp_colors($num_seq);
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
#    print STDERR Dumper [map {"$_ => $opts{$_}" } sort keys %opts];
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
	    $nt =~ tr/ATCGatcg/TAGCtagc/;
	    substr($seq, $start, 0) = $nt x length($spike);
	    open (OUT, ">$file");
	    print OUT ">$name\n";
	    print OUT $seq,"\n";
	    close OUT;
	  }
      }
  }
