#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use GS::LogUser;
use HTML::Template;
use GS::MyDB;                  # to hook into local gb db
use Data::Dumper;
use File::Basename;
use File::Temp;
use GS::MyDB::GBDB::GBObject;
use CoGe::Genome;
use CoGe::Accessory::Web;
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature;
use CoGe::Graphics::Feature::Gene;
use CoGe::Graphics::Feature::NucTide;
use CoGe::Graphics::Feature::GAGA;
use CoGe::Graphics::Feature::Exon_motifs;
use CoGe::Graphics::Feature::AminoAcid;
use CoGe::Graphics::Feature::Domain;
# for security purposes
$ENV{PATH} = "/opt/apache2/CoGe/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw( $DATE $DEBUG $BL2SEQ $TEMPDIR $TEMPURL $USER $FORM $cogeweb);
$BL2SEQ = "/opt/bin/bio/bl2seq ";
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
# set this to 1 to print verbose messages to logs
$DEBUG = 0;

$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$CGI::POST_MAX= 60 * 1024 * 1024; # 24MB
$CGI::DISABLE_UPLOADS = 0; 
($USER) = GS::LogUser->get_user();
my $pj = new CGI::Ajax(
		       rset=>\&Rset,
		       run=>\&Show_Summary,
		       CoGe::Accessory::Web::ajax_func(),
		       loading=>\&loading,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print gen_html();

#print Show_Summary('contig_16224','2665176','10000','10000','502','At1g07300','124379','10000','10000','6','At2g29640','134498','10000','10000','7','','1','0','','','1','0','','','1','0','','','0','1','','','0','1','','','0','1','','1','15','0','0','7','5','2','-2','','1000','150','20');

sub Rset
  {
    return " ";
  }
sub loading
  {
    return qq{<font class="loading">Generating results. . .</font>};
  }

sub gen_html
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(LOGO_PNG=>"SynView-logo.png");
    $template->param(TITLE=>'Synteny Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(NO_BOX=>1);
    $template->param(BODY=>gen_body());
    my $prebox = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynView.tmpl');
    $prebox->param(RESULTS_DIV=>1);
    $template->param(PREBOX=>$prebox->output);
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = shift || $FORM;
    my $accn1 = $form->param('accn1') if $form->param('accn1');
    my $accn2 = $form->param('accn2') if $form->param('accn2');
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynView.tmpl');
    my $box = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
    $template->param(ACCN1=>$accn1);
    $template->param(ACCN2=>$accn2);
    my $html;
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
    my (
	$accn1,
	$featid1,
	$dr1up,
	$dr1down,
	$infoid1,
	$rev1,
	
	$accn2,
	$featid2,
	$dr2up,
	$dr2down,
	$infoid2,
	$rev2,
	
	$accn3,
	$featid3,
	$dr3up,
	$dr3down,
	$infoid3,
	$rev3,
	
	$gbaccn1,
	$gb1start,
	$gb1up,
	$gb1down,
	

	$gbaccn2,
	$gb2start,
	$gb2up,
	$gb2down,
	

	$gbaccn3,
	$gb3start,
	$gb3up,
	$gb3down,
	

	$seq1,
	$dscomp1,
	$dsstart1,
	$dsstop1,

	$seq2,
	$dscomp2,
	$dsstart2,
	$dsstop2,

	$seq3,
	$dscomp3,
	$dsstart3,
	$dsstop3,

	$match_filter,
	$spike_flag,
	$mask_flag,
	$mask_ncs_flag,
	$wordsize,
	$gapopen,
	$gapextend,
	$mismatch,
	$blastparams,
	$iw,
	$ih,
	$feat_h,
	$show_gc,
	$hsp_label,
	$overlap_adjustment,
	$hsp_limit,
	$hsp_limit_num,
	$blast_program,
	#	$seq_file1, $seq_file2,
#	$f1start, $f1up, $f1down,
#	$f2start, $f2up, $f2down,

       ) = @_;
#    print STDERR join "\n",$accn1,$featid1,$dr1up,$dr1down,$infoid1,$rev1,"\n";
	

#    print STDERR Dumper \@_;

    my $stagger_label = $hsp_label =~ /staggered/i ? 1 : 0;
    my $feature_labels = $hsp_label eq "0" ? 0 : 1;
    my @reverse_image = ($rev1, $rev2, $rev3);
    my @locs = ([$dr1up,$dr1down],
		[$dr2up,$dr2down],
		[$dr3up,$dr3down],);
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

    if ($featid1)
      {
	$obj1 = get_obj_from_genome_db( $accn1, $featid1, $infoid1, $rev1, $dr1up, $dr1down );
	return "<font class=error>No entry found for $featid1</font>" unless ($obj1);
	($file1, $file1_begin, $file1_end,$spike_seq) = 
	  generate_seq_file(obj=>$obj1,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"Q",
			    spike_flag=>$spike_flag,
			   );
      }
    elsif ($seq_file1)
      {
	my $sequence = "";
#	print STDERR Dumper $form;
#	print STDERR Dumper \%ENV;
#	print STDERR $seq_file1,"\n";
#	$form->read_multipart(undef, $ENV{'CONTENT_TYPE'});
#	print STDERR Dumper $form->uploadInfo('seq_file1');
	my $fh = $form->upload( 'seq_file1' );
	while ( <$fh> ) {
	  $sequence .= $_;
	}
	return "<font class=error>Problem uploading file $seq_file1</font>" unless ($sequence);
	($obj1) = generate_obj_from_seq($sequence, 1);
	
	($file1, $file1_begin, $file1_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj1,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
#			     revcomp=>$dscomp1,
			     startpos=>$f1start,
			     upstream=>$f1up,
			     downstream=>$f1down,
			     spike_type=>"Q",
			     spike_flag=>$spike_flag, 
			    );
      }
    elsif ($seq1 )
      {
	$seq1 = get_substr(seq=>$seq1, start=>$dsstart1, stop=>$dsstop1);
	($obj1) = generate_obj_from_seq($seq1, 1);
	return "<font class=error>Problem with direct sequence submission</font>" unless ($obj1);
	($file1, $file1_begin, $file1_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj1,
			     revcomp=>$dscomp1,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     spike_type=>"Q",
			     spike_flag=>$spike_flag, 
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
				 spike_flag=>$spike_flag, 
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
	$obj2 = get_obj_from_genome_db( $accn2, $featid2, $infoid2, $rev2, $dr2up, $dr2down);
	return "<font class=error>No entry found for $featid2</font>" unless ($obj2);
	($file2, $file2_begin, $file2_end,$spike_seq) =
	  generate_seq_file(obj=>$obj2,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"S",
			    spike_flag=>$spike_flag,
			   );
      }
   elsif ($seq2 )
      {
	$seq2 = get_substr(seq=>$seq2, start=>$dsstart2, stop=>$dsstop2);
	($obj2) = generate_obj_from_seq($seq2, 2);
	return "<font class=error>Problem with direct sequence submission</font>" unless ($obj2);
	($file2, $file2_begin, $file2_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj2,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     revcomp=>$dscomp2,
			     spike_type=>"Q",
			     spike_flag=>$spike_flag, 
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
			     spike_flag=>$spike_flag, 
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
	$obj3 = get_obj_from_genome_db( $accn3, $featid3, $infoid3, $rev3, $dr3up, $dr3down );
	return "<font class=error>No entry found for $featid3</font>" unless ($obj3);
	($file3, $file3_begin, $file3_end,$spike_seq) =
	  generate_seq_file(obj=>$obj3,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"S",
			    spike_flag=>$spike_flag,
			   );
      }
   elsif ($seq3 )
      {
	$seq3 = get_substr(seq=>$seq3, start=>$dsstart3, stop=>$dsstop3);
	($obj3) = generate_obj_from_seq($seq3, 3);
	return "<font class=error>Problem with direct sequence submission</font>" unless ($obj3);
	($file3, $file3_begin, $file3_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj3,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
			     revcomp=>$dscomp3,
			     spike_type=>"Q",
			     spike_flag=>$spike_flag, 
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
			     spike_flag=>$spike_flag, 
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
    unless ((ref ($obj1) =~ /GBObject/ || ref ($obj1) =~ /hash/i)  && ( ( (ref ($obj2) =~ /GBObject/) || ref ($obj2) =~ /hash/i ) || ( (ref ($obj3) =~ /GBObject/) || ref ($obj3) =~ /hash/i ) ) )
      {
	return "<h3><font color = red>Problem retrieving information.  Please try again.</font></h3>";
	#return 0;
      }
#    return "<pre>".Dumper($obj1, $obj2)."</pre>";
    # set up output page
    
    # run bl2seq
    my $wordsize = $form->param('wordsize');
    $wordsize = 3 if ($blast_program eq "tblastx" && $wordsize > 3);
    my $bl2seq_params = " -W " . $wordsize;
    $bl2seq_params .= " -G " . $form->param('gapopen');
    $bl2seq_params .= " -X " . $form->param('gapextend');
    $bl2seq_params .= " -q " . $form->param('mismatch');
    $bl2seq_params .= " ".join(" ", $form->param('blastparams'));
    my $blast_reports = run_bl2seq( files=>[map {$_->{file}} @sets], accns=>[map {$_->{accn}}@sets], blast_params=>$bl2seq_params, spike_seq=>$spike_seq, match_filter=>$match_filter, blast_program=>$blast_program );

   #sets => array or data for blast
   #blast_reports => array of arrays (report, accn1, accn2, parsed data)

    my $i = 0;
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
							 spike_len=>length($spike_seq),
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
							);
	    $html .= qq!<div>$accn (<font class=species>!.$obj->{ORGANISM}.qq!)</font>!;
	    $html .= qq!($locs[$i][0]::$locs[$i][1])! if defined $locs[$i][0];
	    $html .= qq!<font class=small> Image Inverted</font>! if $reverse_image[$i];
	    $html .= qq!</DIV>\n!;
	    $html .= qq!<IMG SRC="$TEMPURL/$image" !;
	    $html .= qq!BORDER=0 ismap usemap="#$mapname">\n!;
	    $html .= "$map\n";
	  }
	$i++;
      }
    
    $html .= qq!<br>!;
    $html .= qq!<FORM NAME=\"info\">\n!;
    $html .= $form->br();
    $html .= "Information: ";
    $html .= $form->br();
    $html .= qq!<DIV id="info"></DIV>!;
    $html .= "</FORM>\n";
#    $html .= "<P><B>ERROR - NO HITS WERE RETURNED FROM THE BL2SEQ REPORT!</B><P>\n";
    if ($blast_reports && @$blast_reports)
      {
	foreach my $item (@$blast_reports)
	  {
	    my $report = $item->[0];
	    my $accn1 = $item->[1];
	    my $accn2 = $item->[2];
	    my $basereportname = basename( $report );
	    $basereportname = $TEMPURL . "/$basereportname\n";
	    $html .= "<div><font class=xsmall><A HREF=\"$basereportname\">View bl2seq output for $accn1 versus $accn2</A></font></DIV>\n";
	  }
      }


    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
#    print STDERR $html;
    $template->param(BOX_NAME=>"Results:");
    $template->param(BODY=>$html);
    return $template->output;

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
    my $spike_len = $opts{spike_len};
    my $reports = $opts{reports};
    my $iw = $opts{iw} || 1600;
    my $ih = $opts{ih} || 200;
    my $fh = $opts{fh} || 25;
    my $show_gc = $opts{show_gc};
    my $reverse_image = $opts{reverse_image};
    my $stagger_label = $opts{stagger_label};
    my $overlap_adjustment = $opts{overlap_adjustment};
    my $feature_labels = $opts{feature_labels};
    my $hsp_limit = $opts{hsp_limit};
    my $hsp_limit_num = $opts{hsp_limit_num};
    my $graphic = new CoGe::Graphics;
    my $gfx = new CoGe::Graphics::Chromosome;
    $gfx->overlap_adjustment(0);
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
			    overlap_adjustment=>$overlap_adjustment,
			    feature_labels=>$feature_labels,
			   );
    my $f1= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => 1});
    $f1->merge_percent(0);
    $gfx->add_feature($f1);
    my $f2= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => -1});
    $f2->merge_percent(0);
    $gfx->add_feature($f2);
    $graphic->process_nucleotides(c=>$gfx, seq=>$gbobj->{SEQUENCE}, layers=>{gc=>$show_gc});
    process_features(c=>$gfx, obj=>$gbobj, start=>$start, stop=>$stop);
    process_hsps(c=>$gfx, data=>$data, reports=>$reports, accn=>$gbobj->{ACCN}, rev=>$reverse_image, seq_length=> length($gbobj->{SEQUENCE}), stagger_label=>$stagger_label, hsp_limit=>$hsp_limit, hsp_limit_num=>$hsp_limit_num);
    my $file = new File::Temp ( TEMPLATE=>'SynView__XXXXX',
				   DIR=>$TEMPDIR,
				    SUFFIX=>'.png',
				    UNLINK=>0);
    my ($filename) = $file->filename =~ /([^\/]*$)/;;
    $gfx->generate_png(file=>$file->filename);
    close($file);
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
        $c->add_feature($f);
#	print STDERR Dumper ($f) if $f->{description} =~ /at2g29570/i;#unless ($f->start && $f->stop);
    }
  }

sub process_hsps
  {
    my %opts = @_;
    my $c = $opts{c};
    my $data = $opts{data};
    my $reports = $opts{reports};
    my $accn = $opts{accn};
    my $reverse = $opts{rev};
    my $seq_len = $opts{seq_length};
    my $stagger_label = $opts{stagger_label};
    my $hsp_limit = $opts{hsp_limit};
    my $hsp_limit_num = $opts{hsp_limit_num};
    #to reverse hsps when using genomic sequences from CoGe, they need to be drawn on the opposite strand than where blast reports them.  This is because CoGe::Graphics has the option of reverse drawing a region.  However, the sequence fed into blast has already been reverse complemented so that the HSPs are in the correct orientation for the image.  Thus, if the image is reverse, they are drawn on the wrong strand.  This corrects for that problem.   Sorry for the convoluted logic, but it was the simplest way to substantiate this option
    my @colors = (
		  [ 100, 100, 255],
		  [ 0, 255, 0],
		  [ 255, 0, 0],
		 );
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
	my $set = $item->[3];
	unless ($accn eq $accn1 || $accn eq $accn2)
	  {
	    $i++;
	    next;
	  }
	foreach my $item (@$set)
	  {
	    my $color = $colors[$i];

	    my ($start, $stop, $seq);
	    if ($accn1 eq $accn)
	      {
		$start = $item->{qb};
		$stop = $item->{qe};
		$seq = $item->{qmatchseq};
	      }
	    elsif ($accn2 eq $accn)
	      {
		$start = $item->{sb};
		$stop = $item->{se};
		$seq = $item->{smatchseq};
	      }
	    if ($stop < $start)
	      {
		my $tmp = $start;
		$start = $stop;
		$stop = $tmp;
	      }
	    print STDERR "\t",$item->{number},": $start-$stop\n"  if $DEBUG;
	    my $strand = $item->{'orientation'} =~ /-/ ? "-1" : 1;
	    if ($reverse)
	      {
		$strand = $strand =~ /-/ ? "1" : "-1";
		my $tmp = $seq_len - $start+1;
		$start = $seq_len - $stop+1;
		$stop = $tmp;
	      }
	    my $f = CoGe::Graphics::Feature->new({start=>$start, stop=>$stop});
	    $color = [100,100,100] if $item->{'spike_flag'};
	    $f->iw(5);
	    $f->ih(5);
	    $f->skip_duplicate_search(1);
	    $f->gd->fill(0,0,$f->get_color(@$color));
	    $f->color($color);
	    $f->mag(1);
	    $f->order($track);
	    $f->strand($strand);
	    if ($hsp_limit)
	      {
		$f->label($item->{'number'}) if $item->{'number'} <= $hsp_limit_num;
	      }
	    else
	      {
		$f->label($item->{'number'});
	      }
	    $f->force_label(1);
	    my $desc = join ("<br>", "HSP: ".$item->{number}, $start."-".$stop." (".$item->{orientation}.")", $seq,"Match: ".$item->{hspmatch},"Length: ".$item->{length},"Identity: ".$item->{identity},"E_val: ".$item->{eval});
	    $f->description($desc);
	    my $link = "bl2seq_summary.pl?".join("&", "blast_report=".$report, "accnq=", "accns=", "qbegin=", "qend=", "sbegin=","send=","submit=GO");
	    $f->link($link."&"."hsp=".$item->{number});
	    push @feats, $f;
	    print STDERR $item->{number},"-", $item->{orientation}, $track,":", $strand,"\n" if $DEBUG;
	    $f->font_size(10);
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
    return $obj
  }

sub generate_seq_file
  {
    my %options = @_;
    my $obj = $options{obj} || $options{gbobj};
    my $start = $options{start} || $options{startpos} || 1;
    my $up = $options{up} || $options{upstream} || 0;
    my $down = $options{down} || $options{downstream} || 0;
    my $spike_flag = $options{spike_flag} || 0;
    my $spike_type = $options{spike_type};
    my $sequence = $options{sequence} || $options{seq};
    my $revcomp = $options{revcomp} || $options{comp};
    my $mask = $options{mask} || $options{mask_flag};
    my $mask_ncs = $options{mask_ncs},
     my $db = new GS::MyDB;
    my ($file, $file_begin, $file_end, $spike_seq) = $sequence ?
	      $db->{GBSyntenyViewer}->write_fasta_from_sequence(
						  sequence=>$sequence,
						  revcomp=>$revcomp,
						  startpos=>$start,
						  upstream=>$up,
						  downstream=>$down,
						  spike=>[$spike_type,$spike_flag], 
						  path=>$TEMPDIR ) :
	      $db->{GBSyntenyViewer}->write_fasta(
						  GBOBJ=>$obj,
						  accn=>$obj->{ACCN},
						  mask=>$mask,
						  mask_ncs=>$mask_ncs,
						  startpos=>$start,
						  upstream=>$up,
						  downstream=>$down,
						  spike=>[$spike_type,$spike_flag], 
						  path=>$TEMPDIR );
    my %tmp = (
	       start=>$start,
	       upstream=>$up,
	       downstream=>$down,
	       seq_len=>length($obj->{SEQUENCE}),
	       file_begin=>$file_begin,
	       file_end=>$file_end,
	       );

	
    return ($file, $file_begin, $file_end, $spike_seq);
  }

sub get_obj_from_genome_db
  {
    my $accn = shift;
    my $featid = shift;
    my $ds_id = shift;
    my $rev = shift;
    my $up = shift || 0;
    my $down = shift || 0;
    my $db = new CoGe::Genome;
    my $feat = $db->get_feature_obj->retrieve($featid);
    my $start = 0;
    my $stop = 0;
    my $feat_start = 0;
    my $chr;
    my $seq;
    if ($feat)
      {
	$start = $feat->begin_location-$up;
	$start = 1 if $start < 1;
	$stop = $feat->end_location+$down;
#	print STDERR "$start-$stop getting genomic sequence\n";
	$feat_start = $feat->begin_location;
	$chr = $feat->chr;
	$seq = $db->get_genomic_seq_obj->get_sequence(
						      start => $start,
						      stop => $stop,
						      chr => $chr,
						      org_id => $feat->org->id,
						      dataset_id => $ds_id,
						     );
#	print STDERR "seq_len: ", length($seq), ", start: $start, stop: $stop, chr: $chr, org_id: ".$feat->org->id.", dataset_id: $ds_id\n";
	$seq = $db->get_feature_obj->reverse_complement($seq) if $rev;
	
      }
    my $obj= new GS::MyDB::GBDB::GBObject(
					  ACCN=>$accn,
					  LOCUS=>$accn,
					  VERSION=>$feat->dataset->version(),
					  SOURCE=>$feat->dataset->data_source->name(),
					  ORGANISM=>$feat->org->name(),
					  					 );
    $obj->sequence($seq);
    my $fnum = 1;
    my %used_names;
    $used_names{$accn} = 1;

    print STDERR "Region: $chr: $start-$stop\n" if $DEBUG;
#    print STDERR "Region: $chr: ",$start-$start+1,"-",$stop-$start,"\n";

    foreach my $f ($db->get_feature_obj->get_features_in_region(start=>$start, stop=>$stop, chr=>$chr, dataset_id=>$ds_id))
      {
	my $name;
	foreach my $tmp (sort {$b->name cmp $a->name} $f->names)
	  {
	    $name = $tmp->name;
	    last if ($tmp->name =~ /^$accn.+/i);
	  }
	$name = $accn unless $name;
	print STDERR $name,"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(),"\n" if $DEBUG;
	print STDERR "\t", $f->genbank_location_string(recalibrate=>$start),"\n\n" if $DEBUG;
	my $anno = $f->annotation_pretty_print;
	$anno =~ s/\n/<br>/g;
	$anno =~ s/\t/&nbsp;&nbsp;/g;
	my $location = $f->genbank_location_string(recalibrate=>$start);
	print STDERR $f->type->name ,"\t",$location,"\n" if $DEBUG;
	
	$obj->add_feature (
			   F_NUMBER=>$fnum,
			   F_KEY=>$f->type->name,
			   LOCATION=> $location,
			   QUALIFIERS=>[
                                        [annotation=>$anno],
                                        [names=> [map {$_->name} sort {$a->name cmp $b->name} $f->names]],
                                       ],
			   ACCN=>$name,
			   LOCUS=>$name,
			  );
	$fnum++;
      }
    return $obj;
  }


sub check_taint {
	my $v = shift;
	if ($v =~ /^([-\w._=\s+\/]+)$/) {
			$v = $1;
			# $v now untainted
			return(1,$v);
	} else {
	# data should be thrown out
			return(0);
	}
}


sub run_bl2seq {
  my %opts = @_;
  my $files = $opts{files};
  my $accns = $opts{accns};
  my $blast_params = $opts{blast_params};
  my $match_filter = $opts{match_filter};
  my $spike_seq = $opts{spike_seq};
  my $program = $opts{blast_program};
  $program = "blastn" unless $program;
  my $db = new GS::MyDB;
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
	  my $tmp_file = new File::Temp ( TEMPLATE=>'CNS__XXXXX',
					  DIR=>$TEMPDIR,
					  SUFFIX=>'.blast',
					  UNLINK=>0);
	  my ($tempfile) = $tmp_file->filename;# =~ /([^\/]*$)/;
	  
	  # format the bl2seq command
	  $command .= "-p $program -o $tempfile ";
	  $command .= "-i $seqfile1 -j $seqfile2";
	  $command .= " " . $blast_params;
	  my $x = "";
	  ($x,$command) = check_taint( $command );
	  print STDERR $command,"\n" if $STDERR;
	  if ( $DEBUG ) {
	    print STDERR "About to execute...\n $command\n";
	  }
	  # execute the command
	  `$command`;
	  #$reports{$accns->[$i]}{$accns->[$j]} = $tempfile;
	  my $rc;
	  my $data;

	  if (  $match_filter ) 
	    { 
	      ($rc,$data) = $db->{GBSyntenyViewer}->parse_bl2seq( $tempfile, $accn1, $accn2, $spike_seq);
	    } 
	  else 
	    {
	      ($rc,$data) = $db->{GBSyntenyViewer}->parse_bl2seq( $tempfile, $accn1, $accn2);
	    }
	  my @tmp = ($tempfile, $accn1, $accn2);
	  if ($rc)
	    {
	      push @tmp, $data
	    }
	  else
	    {
	      push @tmp, "no results form blasting $accn1 and $accn2";
	    }
	  push @reports, \@tmp;
	}
    }
  return( \@reports );
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
