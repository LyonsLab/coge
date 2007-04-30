#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use HTML::Template;
use GS::MyDB;                  # to hook into local gb db
use Data::Dumper;
use File::Basename;
use File::Temp;
use GS::MyDB::GBDB::GBObject;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGe::Accessory::bl2seq_report;
use CoGe::Accessory::blastz_report;
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
#use Text::Wrap qw($columns &wrap);
use Benchmark qw(:all);

# for security purposes
$ENV{PATH} = "/opt/apache2/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw( $DATE $DEBUG $BL2SEQ $BLASTZ $TEMPDIR $TEMPURL $USER $FORM $cogeweb $BENCHMARK $coge);
$BL2SEQ = "/opt/bin/bio/bl2seq ";
$BLASTZ = "/usr/bin/blastz ";
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$BENCHMARK = 1;

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
		       run=>\&Show_Summary,
		       loading=>\&loading,
		       %ajax,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);

#$USER=1;print gen_html();

#print Show_Summary('contig_16224','2665176','10000','10000','502','At1g07300','124379','10000','10000','6','At2g29640','134498','10000','10000','7','','1','0','','','1','0','','','1','0','','','0','1','','','0','1','','','0','1','','1','15','0','0','7','5','2','-2','','1000','150','20');
#Show_Summary('supercontig_226','3686634','-126900','-232900','513','0','At1g07370','122228','10','10','6','1','At2g29570','142412','10','10','7','0','','1','0','','','1','0','','','1','0','','','0','1','','','0','1','','','0','1','','1','15','0','0','7','5','2','-2','10','','1000','100','20','1','linear','1','0','20','blastn');

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
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	$template->param(LOGO_PNG=>"SynView-logo.png");
	$template->param(TITLE=>'Synteny Viewer');
	$template->param(HELP=>'SynView');
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(NO_BOX=>1);
	$template->param(BODY=>gen_body());
#	$template->param(BODY_ONLOAD=>"setup();");
	my $prebox = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynView.tmpl');
	$prebox->param(RESULTS_DIV=>1);
	$template->param(PREBOX=>$prebox->output);
	
	$html .= $template->output;
      }
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $accn1 = $form->param('accn1') if $form->param('accn1');
    my $accn2 = $form->param('accn2') if $form->param('accn2');
    my $accn3 = $form->param('accn3') if $form->param('accn3');
    my $dr1up = $form->param('dr1up') if $form->param('dr1up');
    my $dr1down = $form->param('dr1down') if $form->param('dr1down');
    my $dr2up = $form->param('dr2up') if $form->param('dr2up');
    my $dr2down = $form->param('dr2down') if $form->param('dr2down');
    my $dr3up = $form->param('dr3up') if $form->param('dr3up');
    my $dr3down = $form->param('dr3down') if $form->param('dr3down');
    $dr1up = 10000 unless defined $dr1up;
    $dr1down = 10000 unless defined $dr1down;
    $dr2up = 10000 unless defined $dr2up;
    $dr2down = 10000 unless defined $dr2down;
    $dr3up = 10000 unless defined $dr3up;
    $dr3down = 10000 unless defined $dr3down;
    my $rev1y = "checked" if $form->param('rev1');
    my $rev1n = "checked" unless $rev1y;
    my $rev2y = "checked" if $form->param('rev2');
    my $rev2n = "checked" unless $rev2y;
    my $rev3y = "checked" if $form->param('rev3');
    my $rev3n = "checked" unless $rev3y;

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynView.tmpl');
    my $box = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
    $template->param(ACCN1=>$accn1);
    $template->param(ACCN2=>$accn2);
    $template->param(ACCN3=>$accn3);
    $template->param(DR1UP=>$dr1up);
    $template->param(DR1DOWN=>$dr1down);
    $template->param(DR2UP=>$dr2up);
    $template->param(DR2DOWN=>$dr2down);
    $template->param(DR3UP=>$dr3up);
    $template->param(DR3DOWN=>$dr3down);
    $template->param(REV1_YES=>$rev1y);
    $template->param(REV1_NO=>$rev1n);
    $template->param(REV2_YES=>$rev2y);
    $template->param(REV2_NO=>$rev2n);
    $template->param(REV3_YES=>$rev3y);
    $template->param(REV3_NO=>$rev3n);

    my $html;
    my $spike_len = match_filter_select();
    $template->param(SPIKE_LEN=>$spike_len);
    $template->param(SEQ_RETRIEVAL=>1);
    $box->param(BOX_NAME=>"Sequence Retrieval:");
    $box->param(BODY=>$template->output);
    $html .= $box->output;
    $template->param(SEQ_RETRIEVAL=>0);
    $template->param(OPTIONS=>1);
    $box->param(BOX_NAME=>"Options:");
    $box->param(BODY=>$template->output);
    $html .= $box->output;
    return $html;
  }

sub Show_Summary 
  {
    my %opts = @_;

    my $accn1 = $opts{accn1};
    my $featid1 = $opts{featid1};
    my $dr1up = $opts{up1};
    my $dr1down = $opts{down1};
    my $rev1 = $opts{rev1};
    
    my $accn2 = $opts{accn2};
    my $featid2 = $opts{featid2};
    my $dr2up = $opts{up2};
    my $dr2down = $opts{down2};
    my $rev2 = $opts{rev2};
    
    my $accn3 = $opts{accn3};
    my $featid3 = $opts{featid3};
    my $dr3up = $opts{up3};
    my $dr3down = $opts{down3};
    my $rev3 = $opts{rev13};
    
    my $gbaccn1 = $opts{gbaccn1};
    my $gb1start = $opts{gb1start};
    my $gb1up = $opts{gb1up};
    my $gb1down = $opts{gb1down};
    
    
    my $gbaccn2 = $opts{gbaccn2};
    my $gb2start = $opts{gb2start};
    my $gb2up = $opts{gb2up};
    my $gb2down = $opts{gb2down};
    
    
    my $gbaccn3 = $opts{gbaccn3};
    my $gb3start = $opts{gb3start};
    my $gb3up = $opts{gb3up};
    my $gb3down = $opts{gb3down};
    
    
    my $seq1 = $opts{seq1};
    my $dscomp1 = $opts{dscomp1};
    my $dsstart1 = $opts{dsstart1};
    my $dsstop1 = $opts{dsstop1};
    
    my $seq2 = $opts{seq2};
    my $dscomp2 = $opts{dscomp2};
    my $dsstart2 = $opts{dsstart2};
    my $dsstop2 = $opts{dsstop2};
    
    my $seq3 = $opts{seq3};
    my $dscomp3 = $opts{dscomp3};
    my $dsstart3 = $opts{dsstart3};
    my $dsstop3 = $opts{dsstop3};
    
    my $match_filter = $opts{matchfilter};
    my $spike_len = $opts{spike};
    my $mask_flag = $opts{mask};
    my $mask_ncs_flag = $opts{maskcns};
    my $wordsize = $opts{word};
    my $gapopen = $opts{gapo};
    my $gapextend = $opts{gape};
    my $mismatch = $opts{mmatch};
    my $eval = $opts{eval};
    my $blastparams = $opts{bp};
    my $iw = $opts{iw};
    my $ih = $opts{ih};
    my $feat_h = $opts{fh};
    my $show_gc = $opts{gc};
    my $color_hsp = $opts{colorhsp};
    my $hsp_label = $opts{hsplabel};
    my $overlap_adjustment = $opts{overlap};
    my $hiqual = $opts{hiqual};
    my $hsp_limit = $opts{hsplim};
    my $hsp_limit_num = $opts{hsplimnum};
    my $blast_program = $opts{prog};
    my $show_hsps_with_stop_codon = $opts{showallhsps};
    my $hsp1r = $opts{h1r};
    my $hsp1g = $opts{h1g};
    my $hsp1b = $opts{h1b};
    my $hsp2r = $opts{h2r};
    my $hsp2g = $opts{h2g};
    my $hsp2b = $opts{h2b};
    my $hsp3r = $opts{h3r};
    my $hsp3g = $opts{h3g};
    my $hsp3b = $opts{h3b};


    $spike_len = 0 unless $match_filter;
    my @hsp_colors;
    foreach my $set ([$hsp1r,$hsp1g,$hsp1b],[$hsp2r,$hsp2g,$hsp2b],[$hsp3r,$hsp3g,$hsp3b])
      {
	my @tmp;
	foreach my $color (@$set)
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

    my ($seq_file1, $seq_file2, $seq_file3);
    my ($f1start, $f1up, $f1down);
    my ($f2start, $f2up, $f2down);
    my ($f3start, $f3up, $f3down);
    my ($obj1, $obj2, $obj3);


    my ($file1, $file1_begin, $file1_end);
    my ($file2, $file2_begin, $file2_end);
    my ($file3, $file3_begin, $file3_end);


    my @sets;
    my $html;
    my $spike_seq;
    my $db = new GS::MyDB;
    my $t1 = new Benchmark;
    if ($featid1)
      {
	$obj1 = get_obj_from_genome_db( $accn1, $featid1, $rev1, $dr1up, $dr1down );
	return "<font class=error>No entry found for $featid1</font>" unless ($obj1);
	($file1, $file1_begin, $file1_end,$spike_seq) = 
	  generate_seq_file(obj=>$obj1,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"Q",
			    spike_len=>$spike_len,
			   );
      }
    elsif ($seq1 )
      {
	$seq1 = get_substr(seq=>$seq1, start=>$dsstart1, stop=>$dsstop1);
	($obj1) = generate_obj_from_seq($seq1, 1, $dscomp1);
	return "<font class=error>Problem with direct sequence submission</font>" unless ($obj1);
	($file1, $file1_begin, $file1_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj1,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     spike_type=>"Q",
			     spike_len=>$spike_len, 
			    );
      }
    elsif ($gbaccn1 )
      {
	$obj1 = $db->{GBSyntenyViewer}->get_genbank_from_nbci($gbaccn1);
	return "<font class=error>No entry found for $gbaccn1</font>" unless ($obj1);
	($file1, $file1_begin, $file1_end,$spike_seq) = 
	      generate_seq_file (
				 obj=>$obj1,
				 mask=>$mask_flag,
				 mask_ncs=>$mask_ncs_flag,
				 startpos=>$gb1start,
				 upstream=>$gb1up,
				 downstream=>$gb1down,
				 spike_type=>"Q",
				 spike_len=>$spike_len, 
			    );
      }
    if ($obj1)
      {
	push @sets, {
		     obj=>$obj1,
		     file=>$file1,
		     file_begin=>$file1_begin,
		     file_end=>$file1_end,
		     accn=>$obj1->{ACCN},
		    };
      }
    #create object 2
    if ($featid2)
      {
	$obj2 = get_obj_from_genome_db( $accn2, $featid2, $rev2, $dr2up, $dr2down);
	return "<font class=error>No entry found for $featid2</font>" unless ($obj2);
	($file2, $file2_begin, $file2_end,$spike_seq) =
	  generate_seq_file(obj=>$obj2,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"S",
			    spike_len=>$spike_len,
			   );
      }
   elsif ($seq2 )
      {
	$seq2 = get_substr(seq=>$seq2, start=>$dsstart2, stop=>$dsstop2);
	($obj2) = generate_obj_from_seq($seq2, 2, $dscomp2);
	return "<font class=error>Problem with direct sequence submission</font>" unless ($obj2);
	($file2, $file2_begin, $file2_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj2,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     spike_type=>"Q",
			     spike_len=>$spike_len, 
			    );
      }
    elsif ($gbaccn2)
      {
	$obj2 = $db->{GBSyntenyViewer}->get_genbank_from_nbci($gbaccn2);	
	return "<font class=error>No entry found for $gbaccn2</font>" unless ($obj2);
	($file2, $file2_begin, $file2_end,$spike_seq) = 
	  generate_seq_file (
			     obj=>$obj2,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     startpos=>$gb2start,
			     upstream=>$gb2up,
			     downstream=>$gb2down,
			     spike_type=>"Q",
			     spike_len=>$spike_len, 
			    );
       }
    if ($obj2)
      {
	push @sets, {
		     obj=>$obj2,
		     file=>$file2,
		     file_begin=>$file2_begin,
		     file_end=>$file2_end,
		     accn=>$obj2->{ACCN},
		    };
      }
    #create object 3
    if ($featid3)
      {
	$obj3 = get_obj_from_genome_db( $accn3, $featid3, $rev3, $dr3up, $dr3down );
	return "<font class=error>No entry found for $featid3</font>" unless ($obj3);
	($file3, $file3_begin, $file3_end,$spike_seq) =
	  generate_seq_file(obj=>$obj3,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"S",
			    spike_len=>$spike_len,
			   );
      }
   elsif ($seq3 )
      {
	$seq3 = get_substr(seq=>$seq3, start=>$dsstart3, stop=>$dsstop3);
	($obj3) = generate_obj_from_seq($seq3, 3, $dscomp3);
	return "<font class=error>Problem with direct sequence submission</font>" unless ($obj3);
	($file3, $file3_begin, $file3_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj3,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     spike_type=>"Q",
			     spike_len=>$spike_len, 
			    );
      }
    elsif ($gbaccn3)
      {
	$obj3 = $db->{GBSyntenyViewer}->get_genbank_from_nbci($gbaccn3);	
	return "<font class=error>No entry found for $gbaccn3</font>" unless ($obj3);
	($file3, $file3_begin, $file3_end,$spike_seq) = 
	  generate_seq_file (
			     obj=>$obj3,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     startpos=>$gb3start,
			     upstream=>$gb3up,
			     downstream=>$gb3down,
			     spike_type=>"Q",
			     spike_len=>$spike_len, 
			    );
       }

    if ($obj3)
      {
	push @sets, {
		     obj=>$obj3,
		     file=>$file3,
		     file_begin=>$file3_begin,
		     file_end=>$file3_end,
		     accn=>$obj3->{ACCN},
		    };
      }
#    print STDERR "\n";

    my @reverse_image;
    push @reverse_image, $rev1 if $obj1;
    push @reverse_image, $rev2 if $obj2;
    push @reverse_image, $rev3 if $obj3;
    my @locs;
    push @locs, [$dr1up,$dr1down] if $obj1;
    push @locs, [$dr2up,$dr2down] if $obj2;
    push @locs, [$dr3up,$dr3down] if $obj3;

    my $obj_count = 0;
    foreach ($obj1, $obj2, $obj3)
      {
	$obj_count++ if ref($_) =~ /GBObject/ || ref($_) =~ /hash/i;
      }
    unless ($obj_count >1)
      {
	return "<h3><font color = red>Problem retrieving information.  Please try again.</font></h3>";
      }
    my $t2 = new Benchmark;
    # set up output page
    
    # run bl2seq
    $wordsize = 3 if ($blast_program eq "tblastx" && $wordsize > 3);
    my $bl2seq_params = " -W " . $wordsize if defined $wordsize;
    $bl2seq_params .= " -G " . $gapopen if defined $gapopen;
    $bl2seq_params .= " -X " . $gapextend if defined $gapextend;
    $bl2seq_params .= " -q " . $mismatch if defined $mismatch;
    $bl2seq_params .= " -e " . $eval if defined $eval;
    $bl2seq_params .= " "    . $blastparams if defined $blastparams;
    my $blast_reports;
    if ($blast_program eq "blastz")
      {
	$blast_reports = run_blastz( files=>[map {$_->{file}} @sets], accns=>[map {$_->{accn}}@sets], params=>undef);
      }
    else
      {
	$blast_reports = run_bl2seq( files=>[map {$_->{file}} @sets], accns=>[map {$_->{accn}}@sets], blast_params=>$bl2seq_params, blast_program=>$blast_program );
      }
#    print Dumper $blast_reports;
   #sets => array or data for blast
   #blast_reports => array of arrays (report, accn1, accn2, parsed data)

    my $i = 0;
    my $t3 = new Benchmark;
    foreach my $item (@sets)
      {
	my $accn = $item->{accn};
	my $obj = $item->{obj};
	my $file_begin = $item->{file_begin};
	my $file_end = $item->{file_end};
	if ($obj)
	  {
	    my ($image, $map, $mapname) = generate_image(
							 gbobj=>$obj, 
							 start=>$file_begin,
							 stop => $file_end,
							 data=>$blast_reports,
							 iw=>$iw,
							 ih=>$ih,
							 fh=>$feat_h,
							 show_gc=>$show_gc,
							 reverse_image => $reverse_image[$i],
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
							);
#	    print STDERR $map;
	    $html .= qq!<div>$accn!;
	    $html .= qq!(<font class=species>!.$obj->{ORGANISM}.qq!</font>)! if $obj->{ORGANISM};
	    $html .= qq!($locs[$i][0]::$locs[$i][1])! if defined $locs[$i][0];
	    $html .= qq!<font class=small> Image Inverted</font>! if $reverse_image[$i];
	    $html .= qq!</DIV>\n!;
	    $html .= qq!<IMG SRC="$TEMPURL/$image" !;
	    $html .= qq!BORDER=0 ismap usemap="#$mapname">\n!;
	    $html .= "$map\n";
	  }
	$i++;
      }
    my $t4 = new Benchmark;
    $html .= qq!<br>!;
    $html .= qq!<FORM NAME=\"info\">\n!;
    $html .= $form->br();
    $html .= "Information: ";
    $html .= $form->br();
    $html .= qq!<DIV id="info"></DIV>!;
    $html .= "</FORM>\n";
    $html .= qq{<table>};
    if ($blast_reports && @$blast_reports)
      {
	$html .= qq{<tr><td class = small>Blast reports};
	foreach my $item (@$blast_reports)
	  {
	    my $report = $item->[0];
	    my $accn1 = $item->[1];
	    my $accn2 = $item->[2];
	    my $basereportname = basename( $report );
	    $basereportname = $TEMPURL . "/$basereportname\n";
	    $html .= "<div><font class=xsmall><A HREF=\"$basereportname\">View blast output for $accn1 versus $accn2</A></font></DIV>\n";
	  }
	$html .= qq{<td class = small>Fasta files};
	foreach my $item (@sets)
	  {
	    my $basename = $TEMPURL."/".basename ($item->{file});
	    my $accn = $item->{accn};
	    $html .= "<div><font class=xsmall><A HREF=\"$basename\">Fasta file for $accn</A></font></DIV>\n";
	  }
	$html .= qq{<td class = small>Annotation files};
	foreach my $item (@sets)
	  {
	    my $basename = $TEMPURL."/".basename (generate_annotation(%$item));
	    my $accn = $item->{accn};
	    $html .= "<div><font class=xsmall><A HREF=\"$basename\">Annotation file for $accn</A></font></DIV>\n";
	  }
      }
    $html .= qq{</table>};

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
#    print STDERR $html;
    $template->param(BOX_NAME=>"Results: $blast_program");
    $template->param(BODY=>$html);
    my $outhtml = $template->output;
#    open (OUT, ">$TEMPDIR/results.html");
#    print STDERR $outhtml;
#    print OUT $outhtml;
#    close OUT;
    my $t5 = new Benchmark;
    my $db_time = timestr(timediff($t2,$t1));
    my $blast_time = timestr(timediff($t3,$t2));
    my $image_time = timestr(timediff($t4,$t3));
    my $html_time = timestr(timediff($t5,$t4));
    print STDERR qq{
SynView Benchmark: $DATE
Time to get DB info             : $db_time
Time to run blast               : $blast_time
Time to generate images and maps: $image_time
Time to process html            : $html_time
} if $BENCHMARK;

    return $outhtml;;

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
    my $graphic = new CoGe::Graphics;
    my $gfx = new CoGe::Graphics::Chromosome;
    $gfx->overlap_adjustment(1);
    $gfx->skip_duplicate_features(1);
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
			    chr_length => length($gbobj->{SEQUENCE}),
			    feature_labels=>1,
			    fill_labels=>1,
			    forcefit=>1,
			    invert_chromosome=>$reverse_image,
			    minor_tick_labels=>-1,
#			    overlap_adjustment=>$overlap_adjustment,
			    feature_labels=>$feature_labels,
			    draw_hi_qual=>$hiqual,
			   );
    $gfx->major_tick_labels(0);
    my $f1= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => 1});
    $f1->merge_percent(0);
    $gfx->add_feature($f1);
    my $f2= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => -1});
    $f2->merge_percent(0);
    $gfx->add_feature($f2);
    $graphic->process_nucleotides(c=>$gfx, seq=>$gbobj->{SEQUENCE}, layers=>{gc=>$show_gc});
    process_features(c=>$gfx, obj=>$gbobj, start=>$start, stop=>$stop, overlap_adjustment=>$overlap_adjustment);
    process_hsps(c=>$gfx, data=>$data, accn=>$gbobj->{ACCN}, rev=>$reverse_image, seq_length=> length($gbobj->{SEQUENCE}), stagger_label=>$stagger_label, hsp_limit=>$hsp_limit, hsp_limit_num=>$hsp_limit_num, gbobj=>$gbobj, spike_seq=>$spike_seq, eval_cutoff=>$eval_cutoff, color_hsp=>$color_hsp, colors=>$hsp_colors, show_hsps_with_stop_codon=>$show_hsps_with_stop_codon);
    my $file = new File::Temp ( TEMPLATE=>'SynView__XXXXX',
				   DIR=>$TEMPDIR,
				    SUFFIX=>'.png',
				    UNLINK=>0);
    my ($filename) = $file->filename =~ /([^\/]*$)/;;
    $gfx->generate_png(file=>$file->filename);
    close($file);
    system "chmod +rw ".$file->filename;
    my $mapname = $filename."map";
    my ($map)=$gfx->generate_imagemap(name=>$mapname);
    return ($filename, $map, $mapname);
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
    my $accn = $obj->{ACCN};
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
	my $type = $feat->{F_KEY};
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} @{$feat->{QUALIFIERS}{names}};
        if ($type =~ /Gene/i)
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
		    foreach my $name (@{$feat->{QUALIFIERS}{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
			$f->label($name) if $name =~ /$accn/i;
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
        	$f->order($track);
		$f->overlay(2);
		if ($accn)
		  {
		    foreach my $name (@{$feat->{QUALIFIERS}{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
		      }
		  }

          }

        next unless $f;
	my $strand = 1;
 	$strand = -1 if $feat->location =~ /complement/;
#	print STDERR $type,"\n";
	foreach my $block (@{$feat->{'blocks'}})
	  {
#	    print STDERR $feat->{QUALIFIERS}{names}[0],", ",  $feat->{F_KEY},": ", $block->[0],"-", $block->[1],"\n";;
	    $block->[0] =1 unless $block->[0]; #in case $block is set to 0
	    $f->add_segment(start=>$block->[0], stop=>$block->[1]);
	    $f->strand($strand);
	    print STDERR "\t", join ("-", @$block),"\n" if $DEBUG;
	  }

	print STDERR $name,"\n\n" if $DEBUG;
        $f->type($type);
	$f->description($feat->{QUALIFIERS}{annotation});
	$f->link("FeatView.pl?accn=$name\" target=\"_new");
	$f->skip_overlap_search($overlap);
        $c->add_feature($f);
#	print STDERR Dumper ($f) if $f->{description} =~ /at2g29570/i;#unless ($f->start && $f->stop);
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
#	my $set = $data->{$accn}{$accn2};
#	next unless ($set->[1] eq $accn || $set->[2] eq $accn);
	my $report = $item->[0];
	my $accn1 = $item->[1];
	my $accn2 = $item->[2];
	my $blast = $item->[3];
#	print STDERR "!",join ("!\t!", $accn1, $accn2, $accn),"!\n";
	next unless ref($blast) =~ /bl2seq/i || ref($blast) =~ /blastz/i;
	unless ($accn eq $accn1 || $accn eq $accn2)
	  {
	    $i++;
	    next;
	  }
#	print Dumper $blast;
	if ($spike_seq)
	  {
	    foreach my $hsp (@{$blast->hsps})
	      {
		if ($hsp->qalign =~ /$spike_seq$/i)
		  {
		    $eval_cutoff = $hsp->eval;
		    last;
		  }
	      }
	  }
#	print STDERR Dumper $blast;
	print STDERR $blast->query," ", $blast->subject,"\n" if $DEBUG;
	foreach my $hsp (@{$blast->hsps})
	  #	while (my $hsp = $blast->nextHSP)
	  {
#	    print STDERR Dumper $hsp;
#	    print STDERR $hsp->qalign,"\n";
	    next if defined $eval_cutoff && $hsp->eval > $eval_cutoff;
	    my $color = $colors->[$i];
	    my $skip = 0;

	    if (!$show_hsps_with_stop_codon && ($hsp->qalign =~ /\*/ || $hsp->salign =~ /\*/))
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
#	    print STDERR $hsp->{'orientation'},"\n";
#	    print STDERR Dumper 
	    if ($reverse)
	      {
		$strand = $strand =~ /-/ ? "1" : "-1";
		my $tmp = $seq_len - $start+1;
		$start = $seq_len - $stop+1;
		$stop = $tmp;
	      }
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
	    my $desc = join ("<br>", "HSP: ".$hsp->number. "  <font class=small>(".$blast->query."-". $blast->subject.")</font>", $start."-".$stop." (".$hsp->strand.")", $seq,"Match: ".$hsp->match,"Length: ".$hsp->length,"Identity: ".$hsp->percent_id,"E_val: ".$hsp->pval);
	    $f->description($desc);
#	    my $link = "bl2seq_summary.pl?".join("&", "blast_report=".$report, "accnq=$accn1", "accns=$accn2", "qbegin=".($gbobj->{start}+$start-1), "qend=".($gbobj->{start}+$stop-1),"qchr=".$gbobj->{chr}, "qds=". $gbobj->{ds}, "sbegin=","send=","submit=GO") if $gbobj->{ds};
	    my $link = "HSPView.pl?blast_report=$report&hsp_num=".$hsp->number;
	    $link .= join ("&","&qstart=".($gbobj->{start}+$start-1), "qstop=".($gbobj->{start}+$stop-1),"qchr=".$gbobj->{chr}, "qds=". $gbobj->{ds},"qstrand=".$strand) if $gbobj->{ds};
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
    @feats = reverse @feats if $reverse;
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
    if ($sequence =~ /^LOCUS/)
      {
	#genbank sequence
	my $db = new GS::MyDB;
	$obj = $db->{GBSyntenyViewer}->parse_genbank($sequence);
      }
   elsif ($sequence =~ /^>/)
      {
	#fasta sequence
	my ($header, $seq) = split /\n/, $sequence, 2;
	my ($accn) = $header=~/>(\S*)/;
	$accn =~ s/\|/_/g;
	$obj = new GS::MyDB::GBDB::GBObject(
					    ACCN=>$accn,
					    LOCUS=>$accn,
					    DEFINITION=>$header,
					    );
	$seq =~ s/\n|\r//g;
	$obj->sequence($seq);
      }
    else
      {
	#just the sequence
	$obj = new GS::MyDB::GBDB::GBObject(
					     ACCN=>"RAW_SEQUENCE_SUBMISSION $num",
					     LOCUS=>"RAW_SEQUENCE_SUBMISSION $num",
					    );
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
    my $spike_type = $options{spike_type};
    my $mask = $options{mask} || $options{mask_flag};
    my $mask_ncs = $options{mask_ncs};
    my $t1 = new Benchmark;
    my ($file, $file_begin, $file_end, $spike_seq) = 
      write_fasta(
		  GBOBJ=>$obj,
		  accn=>$obj->{ACCN},
		  mask=>$mask,
		  mask_ncs=>$mask_ncs,
		  startpos=>$start,
		  upstream=>$up,
		  downstream=>$down,
		  spike=>[$spike_type,$spike_len],
		  path=>$TEMPDIR,
		  
		 );
    my %tmp = (
	       start=>$start,
	       upstream=>$up,
	       downstream=>$down,
	       seq_len=>length($obj->{SEQUENCE}),
	       file_begin=>$file_begin,
	       file_end=>$file_end,
	       );
    my $t2 = new Benchmark;
    my $time = timestr(timediff($t2,$t1));
    print STDERR "Time to generate Sequence file:   $time\n" if $BENCHMARK;
    return ($file, $file_begin, $file_end, $spike_seq);
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
#	print STDERR join ("\t", $accn, $chr, $start, $stop, $feat->start, $feat->stop, $up, $down),"\n";

	$seq = $coge->resultset('Dataset')->find($ds_id)->get_genome_sequence(
									      start => $start,
									      stop => $stop,
									      chr => $chr,
									     );
	$seq = CoGeX::Feature->reverse_complement($seq) if $rev;
      }
    my $t3 = new Benchmark;

    my $obj= new GS::MyDB::GBDB::GBObject(
					  ACCN=>$accn,
					  LOCUS=>$accn,
					  VERSION=>$feat->dataset->version(),
					  SOURCE=>$feat->dataset->datasource->name(),
					  ORGANISM=>$feat->org->name(),
					  					 );

    $obj->sequence($seq);
    $obj->{start} = $start;
    $obj->{stop} = $stop;
    $obj->{chr} = $chr;
    $obj->{ds} = $ds_id;
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
	$name = $accn unless $name;
	print STDERR $name,"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(),"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(recalibrate=>$start),"\n\n" if $DEBUG;
	my $anno = $f->annotation_pretty_print;
	$anno =~ s/\n/<br>/g;
	$anno =~ s/\t/&nbsp;&nbsp;/g;
	my $location = $f->genbank_location_string(recalibrate=>$start);
	print STDERR $name, "\t",$f->type->name ,"\t",$location,"\n" if $DEBUG;
	
	$obj->add_feature (
			   F_NUMBER=>$fnum,
			   F_KEY=>$f->type->name,
			   LOCATION=> $location,
			   QUALIFIERS=>[
                                        [annotation=>$anno],
                                        [names=> [@names]],
                                       ],
			   ACCN=>$name,
			   LOCUS=>$name,
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
#	print STDERR "check_taint received $v\n";
	if ($v =~ /^([A-Za-z0-9\-\.=\/_]*)$/) {
		$v = $1;
#		print STDERR "check_taint turned $v into $1\n";
		return(1,$v);
	} else {
		return(0,0);
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
  my $files = $opts{files};
  my $accns = $opts{accns};
  my $blast_params = $opts{blast_params};
  my $program = $opts{blast_program};
  my $eval_cutoff = $opts{eval_cutoff};
  $program = "blastn" unless $program;
  my @files;
  foreach my $item (@$files)
    {
      next unless $item;
      if (-r $item)
	{
	  push @files, $item
	}
    }
  my @reports;
  for (my $i=0; $i<scalar @files; $i++)
    {
      for (my $j=0; $j<scalar @files; $j++)
	{
	  next unless $j > $i;
	  my $seqfile1 = $files[$i];
	  my $seqfile2 = $files[$j];
	  my ($accn1, $accn2) = ($accns->[$i], $accns->[$j]);
	  my $command = $BL2SEQ;
	  
	  
	# need to create a temp filename here
	  my $tmp_file = new File::Temp ( TEMPLATE=>'Bl2Seq__XXXXX',
					  DIR=>$TEMPDIR,
					  SUFFIX=>'.txt',
					  UNLINK=>0);
	  my ($tempfile) = $tmp_file->filename;# =~ /([^\/]*$)/;
	  
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
	  `$command`;
	  system "chmod +rw $tempfile";
	  #$reports{$accns->[$i]}{$accns->[$j]} = $tempfile;
	  my $blastreport = new CoGe::Accessory::bl2seq_report($tempfile) if -r $tempfile;
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
    my $files = $opts{files};
    my $params = $opts{params};
    my $accns = $opts{accns};
     my @files;
  foreach my $item (@$files)
    {
      next unless $item;
      if (-r $item)
	{
	  push @files, $item
	}
    }
    my @reports;
    for (my $i=0; $i<scalar @files; $i++)
      {
	for (my $j=0; $j<scalar @files; $j++)
	  {
	    next unless $j > $i;
	    my $seqfile1 = $files[$i];
	    my $seqfile2 = $files[$j];
	    my ($accn1, $accn2) = ($accns->[$i], $accns->[$j]);
	    my $tmp_file = new File::Temp ( TEMPLATE=>'blastz__XXXXX',
					  DIR=>$TEMPDIR,
					  SUFFIX=>'.txt',
					  UNLINK=>0);
	    my ($tempfile) = $tmp_file->filename;# =~ /([^\/]*$)/;
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
	    `$command`;
	    system "chmod +rw $tempfile";
	    my $blastreport = new CoGe::Accessory::blastz_report($tempfile) if -r $tempfile;
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

sub write_fasta {
	my $options =  {};
	# main input is ACCN, startpos, upstream, downstream, mask and path
	# return value for this function is the sequence and the start/end
	# position of the sequence.
	# (10/10/03):  Now check if the spike_flag is true.  If so,
	# construct and attach an appropriately sized spike sequence at the
	# end of each sequence.  the options parameter spike is set to the
	# size of the spike, or 0 if not set on the original form.

	# vars
	my($seq,$startpos,$seq_begin,$seq_end,$db_begin,$db_end,$chr);
	my($gene_end,$gene_begin);
	if ( @_ ) {
		while ( @_ ) {
			my($key,$value) = splice(@_,0,2);;
			$options->{$key} = $value;
		}
	} else {
		return( 0 );
	}

	my $gbobj = $options->{GBOBJ};

	if ( not defined $options->{path} ) {
		$options->{path} = "";
	} else {
		# does the path end in a "/"?  add if no
		if ( $options->{path} !~ /(.*)\/$/ ) {
			$options->{path} .= "/";
		}
		# other checks here?
	}
	my $tmp_file = new File::Temp ( TEMPLATE=>'SEQ__XXXXX',
					DIR=>$options->{path},
					SUFFIX=>'.fa',
					UNLINK=>0);
	my $fullname = $tmp_file->filename;

	my $hdr = $gbobj->get_headerfasta( );
	
	$seq_begin = $options->{startpos} - $options->{upstream};
	if ( $seq_begin < 1 ) {
		$seq_begin = 1;
	}
	($seq) = uc($gbobj->sequence());

	if ($options->{downstream})
	{
	    $seq_end = $options->{startpos} + $options->{downstream};
	} else {
	  $seq_end = length($seq);
	}
	if ( $seq_end > length( $seq )) {
		$seq_end = length($seq);
	}
	# mask out exons if required
	if ( defined $options->{mask} and $options->{mask} == 1 ) {
		$seq = $gbobj->mask_exons( $seq );
	}
	# mask out non coding sequences if required
	if ( defined $options->{mask_ncs} and $options->{mask_ncs} == 1 ) {
		$seq = $gbobj->mask_ncs( $seq );
	}
	($seq) = $gbobj->subsequence( $seq_begin, $seq_end, $seq );
	my $spike_seq = "";
	my($pair,$spikesize);
	if ( defined $options->{spike} ) {
		($pair,$spikesize) = @{ $options->{spike} };
		if ( $spikesize > 0 ) {
			($seq,$spike_seq) = spike( $seq, $spikesize, $pair );
		}
	}
	my $error = 0;
	($error,$fullname) = check_filename_taint( $fullname );
	if ( $error ) 
	  {
	    my $length = length($seq);
	    open(OUT, ">$fullname") or die "Couldn't open $fullname!\n";
	    print OUT "$hdr\n";
#	    my $max = ceil ($length/10000);
	    my $max = 100;# if $max < 100;
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
	    return($fullname,$seq_begin,$seq_end,$spike_seq);
	  } 
	else 
	  {
	    return(0,0,0,"");
	  }
}

sub spike { 
	my $seq = shift;
	my $spikesize = shift;
	return ($seq) unless $spikesize;
	$spikesize -= 12; # to account for the "CCCAAAGGGTTT" tag
	my $pair = shift;
	my $spike_seq = "CCCAAAGGGTTT";
	my @chars = split( //, $seq );

	if ( $pair eq "Q" ) {
		if ( $chars[-1] eq "A" ) {
			$seq .= "T";
		} else {
			$seq .= "ANNNNNNNNNNN";
		}
	} else {
		if ( $chars[-1] eq "G" ) {
			$seq .= "C";
		} else {
			$seq .= "G";
		}
	}

	for ( my $n = 0; $n < $spikesize; $n++ ) {
		$spike_seq .= "C";
	}
	$seq .= $spike_seq;
	return($seq,$spike_seq);
}

sub match_filter_select
  {
    my $match = shift;
    my $form = shift || $FORM;
    $match = $form->param('spike_len') if $form->param('spike_len');
    $match = 15 unless $match;
    my $html;
    for (my $i = 13; $i<=18; $i++)
      {
	my $item = qq{<label><input class="backbox" type="radio" name="spike" id="spike" value="$i"};
	$item .= " checked=\"checked\"" if $match && $match == $i;
	$item .= qq{ />$i/$i</label>};
	$html .= $item;
      }
    $html .= " nucleotides";
    return $html;
  }

sub generate_annotation
  {
    my %opts = @_;
    my $obj = $opts{obj};
    my $start = $opts{file_begin};
    my $stop = $opts{file_end};
    my @opts = ($start, $stop);# if $start && $stop;
    my $tmp_file = new File::Temp ( TEMPLATE=>'Anno__XXXXX',
				    DIR=>$TEMPDIR,
				    SUFFIX=>'.txt',
				    UNLINK=>0);
    my $fullname = $tmp_file->filename;
    my %data;
    foreach my $feat($obj->get_features(@opts))
      {
	my $type = $feat->{F_KEY};
	my ($name) = sort { length ($b) <=> length ($a) || $a cmp $b} @{$feat->{QUALIFIERS}{names}};
	my $dir = $feat->location =~ /complement/ ? "<" : ">";
	foreach my $block (@{$feat->{'blocks'}})
	  {
	    $block->[0] =1 unless $block->[0];
	    push @{$data{$name}{$type}},[$block->[0], $block->[1], $dir];
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
		$type = "utr" if $type =~ /rna/i;
		foreach my $item (@{$data{$name}{$type}})
		  {
		    print OUT join (" ", $item->[0], $item->[1],$type),"\n";
		  }
	      }
	  }
      }
    close OUT;
    return $fullname;
  }
