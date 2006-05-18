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

use vars qw( $DATE $DEBUG $BL2SEQ $TEMPDIR $TEMPURL $USER);
$BL2SEQ = "/opt/bin/bio/bl2seq ";
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
# set this to 1 to print verbose messages to logs
$DEBUG = 0;

$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

my $form = new CGI;
$CGI::POST_MAX= 60 * 1024 * 1024; # 24MB
$CGI::DISABLE_UPLOADS = 0; 
$USER = GS::LogUser->get_user();
my $pj = new CGI::Ajax(
		       rset=>\&Rset,
		       run=>\&Show_Summary,
		       source_search=>\&get_data_source_info_for_accn,
		       loading=>\&loading,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
$pj->js_encode_function('escape');
print $pj->build_html($form, \&gen_html);
#print gen_html();


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
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynView.tmpl');
    $template->param(TITLE=>'Synteny Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    my $accn1 = $form->param('accn1') if $form->param('accn1');
    my $accn2 = $form->param('accn2') if $form->param('accn2');
    $template->param(ACCN1=>$accn1);
    $template->param(ACCN2=>$accn2);
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub Show_Summary 
  {
    my (
	$accn1, $accn2,
	$dr1up, $dr1down,
	$dr2up, $dr2down,
	$infoid1, $infoid2,
#	$seq_file1, $seq_file2,
#	$f1start, $f1up, $f1down,
#	$f2start, $f2up, $f2down,
	$gbaccn1, $gbaccn2,
	$gb1start, $gb1up, $gb1down,
	$gb2start, $gb2up, $gb2down,
	$seq1, $seq2,
	$dscomp1, $dscomp2,
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
       ) = @_;
    my ($seq_file1, $seq_file2);
    my ($f1start, $f1up, $f1down);
    my ($f2start, $f2up, $f2down);

#    print STDERR"<pre>".Dumper(\@_),"</pre>";
    my $form = new CGI;
    my $html;
    my ($obj1, $obj2);
    my ($file1, $file1_begin, $file1_end);
    my ($file2, $file2_begin, $file2_end);
    my $spike_seq;
    my $db = new GS::MyDB;
    if ($accn1)
      {
	$obj1 = get_obj_from_genome_db( $accn1, $infoid1,$dr1up, $dr1down );
	return "<font class=error>No entry found for $accn1</font>" unless ($obj1);
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
	print STDERR Dumper \%ENV;
	print STDERR $seq_file1,"\n";
#	$form->read_multipart(undef, $ENV{'CONTENT_TYPE'});
	print STDERR Dumper $form->uploadInfo('seq_file1');
	my $fh = $form->upload( 'seq_file1' );
	while ( <$fh> ) {
	  $sequence .= $_;
	}
	return "<font class=error>Problem uploading file $seq_file1</font>" unless ($sequence);
	($obj1) = generate_obj_from_seq($sequence);
	
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
	($obj1) = generate_obj_from_seq($seq1);
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
    #create object 2
    if ($accn2)
      {
	$obj2 = get_obj_from_genome_db( $accn2, $infoid2, $dr2up, $dr2down );
	return "<font class=error>No entry found for $accn1</font>" unless ($obj1);
	($file2, $file2_begin, $file2_end,$spike_seq) =
	  generate_seq_file(obj=>$obj2,
			    mask=>$mask_flag,
			    mask_ncs=>$mask_ncs_flag,
			    spike_type=>"S",
			    spike_flag=>$spike_flag,
			   );
      }
     elsif ($seq_file2)
      {
	my $sequence = "";
#	print STDERR Dumper $form;
	print STDERR Dumper \%ENV;
	print STDERR $seq_file2,"\n";
	$form->read_multipart(undef, $ENV{'CONTENT_TYPE'});
	print STDERR Dumper $form->uploadInfo('seq_file2');
	my $fh = $form->upload( 'seq_file2' );
	while ( <$fh> ) {
	  $sequence .= $_;
	}
	return "<font class=error>Problem uploading file $seq_file2</font>" unless ($sequence);
	($obj2) = generate_obj_from_seq($sequence);

	($file2, $file2_begin, $file2_end, $spike_seq) = 
	  generate_seq_file (
			     obj=>$obj2,
			     mask=>$mask_flag,
			     mask_ncs=>$mask_ncs_flag,
#			     revcomp=>$dscomp2,
			     startpos=>$f2start,
			     upstream=>$f2up,
			     downstream=>$f2down,
			     spike_type=>"Q",
			     spike_flag=>$spike_flag, 
			    );
      }
    elsif ($seq2 )
      {
	($obj2) = generate_obj_from_seq($seq2);
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
    unless ((ref ($obj1) =~ /GBObject/ || ref ($obj1) =~ /hash/i)  && (ref ($obj2) =~ /GBObject/) || ref ($obj2) =~ /hash/i)
      {
	return "<h3><font color = red>Problem retrieving information.  Please try again.</font></h3>";
	#return 0;
      }
    my $rc = 1;
#    return "<pre>".Dumper($obj1, $obj2)."</pre>";
    # set up output page
    
    # run bl2seq
    my $bl2seq_params = " -W " . $form->param('wordsize');
    $bl2seq_params .= " -G " . $form->param('gapopen');
    $bl2seq_params .= " -X " . $form->param('gapextend');
    $bl2seq_params .= " -q " . $form->param('mismatch');
    $bl2seq_params .= join(" ", $form->param('blastparams'));
    my $blast_report = run_bl2seq( $file1, $file2, $bl2seq_params );
    
    # it doesn't matter which $dbo[12], cuz they both point to this
    # blastreport
    my $data = "";
    if (  $match_filter ) 
      { 
	($rc,$data) = $db->{GBSyntenyViewer}->parse_bl2seq( $blast_report, $obj1->{ACCN}, $obj2->{ACCN}, $spike_seq);
      } 
    else 
      {
	($rc,$data) = $db->{GBSyntenyViewer}->parse_bl2seq( $blast_report, $obj1->{ACCN}, $obj2->{ACCN});
      }
    if ( $rc ) 
      {
 	my($filename,$maphtml,$mapname);# =
# 	  $db->{GBSyntenyViewer}->draw_image(
# 					     ACCN_Q_OBJ=>$obj1,
# 					     Q_BEGIN=>$file1_begin, 
# 					     Q_END=>$file1_end,
# 					     ACCN_S_OBJ=>$obj2,
# 					     S_BEGIN=>$file2_begin, 
# 					     S_END=>$file2_end, 
# 					     HSPS=>$data, 
# 					     DIRNAME=>$TEMPDIR,
# 					     BLAST_REPORT=>$blast_report,
# 					     MASKED_EXONS=>$mask_flag,
# 					     SPIKE_LEN=>length($spike_seq),
# 					    );

	$accn1 = 0 unless $accn1;
	$accn2 = 0 unless $accn2;
	my ($qfile, $qmap, $qmapname) = generate_image(
						       gbobj=>$obj1, 
						       start=>$file1_begin,
						       stop => $file1_end,
						       hsps=>[map {{
							 match=>$_->{hspmatch},
							 number=>$_->{number},
							 start=>$_->{qb},
							 stop=>$_->{qe},
							 seq=>$_->{qmatchseq},
							 orientation=>$_->{orientation},
							 length=>$_->{length},
							 identity=>$_->{identity},
							 eval=>$_->{eval},
							 spike_flag=>$_->{spike_flag},
						       }} @$data],
						       spike_len=>length($spike_seq),
						       report=>$blast_report,
						       iw=>$iw,
						       ih=>$ih,
						       fh=>$feat_h,
						       );
	my ($sfile, $smap, $smapname) = generate_image(
						       gbobj=>$obj2, 
						       start=>$file2_begin,
						       stop => $file2_end,
						       hsps=>[map {{
							 match=>$_->{hspmatch},
							 number=>$_->{number},
							 start=>$_->{sb},
							 stop=>$_->{se},
							 seq=>$_->{smatchseq},
							 orientation=>$_->{orientation},
							 length=>$_->{length},
							 identity=>$_->{identity},
							 eval=>$_->{eval},
							 spike_flag=>$_->{spike_flag},
						       }} @$data],
						       spike_len=>length($spike_seq),
						       report=>$blast_report,
						       iw=>$iw,
						       ih=>$ih,
						       fh=>$feat_h,
						       );
	$html .= qq!<IMG SRC="$TEMPURL/$qfile" !;
	$html .= qq!BORDER=0 ismap usemap="#$qmapname">\n!;
	$html .= "$qmap\n";
	$html .= qq!<br>!;
	$html .= qq!<IMG SRC="$TEMPURL/$sfile" !;
	$html .= qq!BORDER=0 ismap usemap="#$smapname">\n!;
	$html .= "$smap\n";
	
	$html .= qq!<br>!;
	$html .= qq!<FORM NAME=\"info\">\n!;
	$html .= $form->br();
	$html .= "Information: ";
	$html .= $form->br();
	$html .= qq!<DIV id="info"></DIV>!;
	$html .= "</FORM>\n";
      } 
    else 
      {
	$html .= "<P><HR><B>ERROR - NO HITS WERE RETURNED FROM THE BL2SEQ REPORT!</B><P>\n";
      }
    my $basereportname = basename( $blast_report );
    $basereportname = $TEMPURL . "/$basereportname\n";
    $html .= "<font class=xsmall><A HREF=\"$basereportname\">View bl2seq output: $basereportname</A></font>\n";

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SynView_results.tmpl');
#    print STDERR $html;
    $template->param(HEADING=>"Results");
    $template->param(RESULTS=>$html);
    return $template->output;

}


sub generate_image
  {
    my %opts = @_;
    my $gbobj = $opts{gbobj};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $hsps = $opts{hsps};
    my $masked_exons = $opts{mask};
    my $masked_ncs = $opts{mask_ncs};
    my $spike_len = $opts{spike_len};
    my $report = $opts{report};
    my $iw = $opts{iw} || 1600;
    my $ih = $opts{ih} || 200;
    my $fh = $opts{fh} || 25;
#    print STDERR Dumper (\%opts);
#    print STDERR Dumper $hsps;
    my $c = new CoGe::Graphics::Chromosome;
    $c->chr_length(length($gbobj->{SEQUENCE}));
#    $c->mag_scale_type("constant_power");
    $c->iw($iw);
    $c->max_mag(10);
    $c->feature_labels(1);
    $c->fill_labels(1);
    $c->draw_chromosome(1);
    $c->draw_ruler(1);
    $c->draw_chr_end(0);
    $c->chr_start_height($ih);
    $c->feature_height($fh);
    $c->chr_mag_height(5);
    $c->set_region(start=>$start, stop=>$stop);
    $c->mag(1);
#    $c->start_picture('left');
    my $f1= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => 1});
    $f1->merge_percent(0);
    $c->add_feature($f1);
    my $f2= CoGe::Graphics::Feature->new({start=>1, order => 2, strand => -1});
    $f2->merge_percent(0);
    $c->add_feature($f2);
    my $link = "bl2seq_summary.pl?".join("&", "blast_report=$report", "accnq=", "accns=", "qbegin=", "qend=", "sbegin=","send=","submit=GO");
    process_nucleotides(c=>$c, seq=>$gbobj->{SEQUENCE});
    process_features(c=>$c, obj=>$gbobj, start=>$start, stop=>$stop);
    process_hsps(c=>$c, hsps=>$hsps, link=>$link);
    my $file = new File::Temp ( TEMPLATE=>'SynView__XXXXX',
				   DIR=>$TEMPDIR,
				    SUFFIX=>'.png',
				    UNLINK=>0);
    my ($filename) = $file->filename =~ /([^\/]*$)/;;
    $c->generate_png(file=>$file->filename);
    close($file);
    my $mapname = $filename."map";
    my ($map)=$c->generate_imagemap(name=>$mapname);
    return ($filename, $map, $mapname);
  }

sub process_nucleotides
  {
    my %opts = @_;
    my $c = $opts{c};
    #process nucleotides
    my $seq = uc($opts{seq});
    
    my $seq_len = length $seq;
    my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw);
    $chrs = 1 if $chrs < 1;
    my $pos = 0;
    my $start = 1;# if $start < 1;    
    while ($pos < $seq_len)
      {
        my $subseq = substr ($seq, $pos, $chrs);
        my $rcseq = substr ($seq, $pos, $chrs);
        $rcseq =~ tr/ATCG/TAGC/;
        next unless $subseq && $rcseq;
        my $f1 = CoGe::Graphics::Feature::NucTide->new({nt=>$subseq, strand=>1, start =>$pos+$start});
	my $f2 = CoGe::Graphics::Feature::NucTide->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
	#my $f2 = CoGe::Graphics::Feature::Exon_motifs->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});
        #my $f2 = CoGe::Graphics::Feature::GAGA->new({nt=>$rcseq, strand=>-1, start =>$pos+$start});        
        $c->add_feature($f1, $f2);
        $pos+=$chrs;
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
    my $accn = $obj->{ACCN};
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
	    next;
	    $f = CoGe::Graphics::Feature::Gene->new();
	    $f->color([255,0,0,50]);
	    if ($accn)
	      {
		foreach my $name (@{$feat->{QUALIFIERS}{names}})
		  {
		    $f->color([255,255,0]) if $name =~ /$accn/i;
		  }
	      }
	    $f->order(2);
	    $f->overlay(1);
	    $f->mag(0.5)
          }
        elsif ($type =~ /CDS/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,255,0, 50]);
        	$f->order(2);
		$f->overlay(3);
		if ($accn)
		  {
		    foreach my $name (@{$feat->{QUALIFIERS}{names}})
		      {
			$f->color([255,255,0]) if $name =~ /$accn/i;
		      }
		  }

#		$f->label($name);
		
#        	draw_prots(genomic_feat=>$feat, c=>$c, chrom_feat=>$f);
          }
        elsif ($type =~ /mrna/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([0,0,255, 50]);
        	$f->order(2);
		$f->overlay(2);
		$f->mag(0.75);
          }
        elsif ($type =~ /rna/i)
          {
        	$f = CoGe::Graphics::Feature::Gene->new();
        	$f->color([200,200,200, 50]);
        	$f->order(2);
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
	print STDERR $type,"\n" if $DEBUG;
	foreach my $block (@{$feat->{'blocks'}})
	  {
	    $block->[0] =1 unless $block->[0]; #in case $block is set to 0
	    $f->add_segment(start=>$block->[0], stop=>$block->[1]);
	    $f->strand($strand);
	    print STDERR "\t", join ("-", @$block),"\n" if $DEBUG;
	  }

	print STDERR $name,"\n\n" if $DEBUG;
        $f->type($type);
	$f->description($feat->{QUALIFIERS}{annotation});
	$f->link("FeatView.pl?accn=$name");
        $c->add_feature($f);
    }
  }

sub process_hsps
  {
    my %opts = @_;
    my $c = $opts{c};
    my $hsps = $opts{hsps};
    my $link = $opts{link};
   
    print STDERR "HSPS:\n" if $DEBUG;
    foreach my $hsp (@$hsps)
      {
	my $start = $hsp->{start};
	my $stop = $hsp->{stop};
	if ($start > $stop)
	  {
	    $start = $hsp->{stop};
	    $stop = $hsp->{start};
	  }
	print STDERR "\t",$hsp->{number},": $start-$stop\n" if $DEBUG;
	my $f = CoGe::Graphics::Feature->new({start=>$start, stop=>$stop});
	#my $color = $hsp->{'orientation'} =~ /-/ ? [100,100,255]: [ 255, 100, 255];
	my $color = [ 255, 100, 255];
	my $strand = $hsp->{'orientation'} =~ /-/ ? "-1" : 1;
	$color = [100,100,100] if $hsp->{'spike_flag'};
	$f->iw(5);
	$f->ih(5);
	$f->gd->fill(0,0,$f->get_color(@$color));
	$f->color($color);
	$f->mag(1.5);
	$f->order(1);
#	$f->add_segment();
	$f->strand($strand);
	$f->label($hsp->{'number'});
	$f->force_label(1);
	my $desc = join ("<br>", "HSP: ".$hsp->{number}, $hsp->{start}."-".$hsp->{stop}." (".$hsp->{orientation}.")", $hsp->{seq},"Match: ".$hsp->{match},"Length: ".$hsp->{length},"Identity: ".$hsp->{identity},"E_val: ".$hsp->{eval});
	$f->description($desc);
	$f->link($link."&"."hsp=".$hsp->{number});
#	print STDERR Dumper $f;
	$c->add_feature($f);
	
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

#foreach my $i (1..10)
# { 
#   $c->mag($i);
#   $c->generate_png(file=>"tmp/test$i.png");
# }

sub draw_prots
  {
    my %opts = @_;
    my $feat = $opts{genomic_feat};
    my $c = $opts{c};
    my $f = $opts{chrom_feat};
    #Do we have any protein sequence we can use?
    foreach my $seq ($feat->sequences)
      {
	next unless $seq->seq_type->name =~ /prot/i;
	my ($pseq) = $seq->sequence_data;
	my $chrs = int (($c->_region_stop-$c->_region_start)/$c->iw)/3;
	$chrs = 1 if $chrs < 1;
	my $pos = 0;
	while ($pos <= length $pseq)
	  {
	    my $aseq = substr($pseq, $pos, $chrs);
	    foreach my $loc ($seq->get_genomic_locations(start=>$pos+1, stop=>$pos+$chrs))
	      {
		my $ao = CoGe::Graphics::Feature::AminoAcid->new({aa=>$aseq, start=>$loc->start, stop=>$loc->stop, strand => $f->strand, order=>4});
		$ao->skip_overlap_search(1);
		$c->add_feature($ao);
		delete $loc->{__Changed}; #silence the warning from Class::DBI
	      }
	    
	    $pos+=$chrs;
	  }
      }
  }



sub get_data_source_info_for_accn
  {
    my $accn = shift;
    return unless $accn;
    my $num = shift;
    my $db = new CoGe::Genome;
    my @feats = $db->get_feats_by_name($accn);
    my %sources;
    foreach my $feat (@feats)
      {
	$sources{$feat->data_info->id} = $feat->data_info;
      }
    my $html;
    if (keys %sources)
      {
	$html .= qq{
<SELECT name = "infoid$num" id= "infoid$num">
};
	my $count = 0;
	foreach my $id (sort {$b <=> $a} keys %sources)
	  {
	    my $val = $sources{$id};
	    my $name = $val->name;
	    my $ver = $val->version;
	    my $desc = $val->description;
	    my $sname = $val->data_source->name;
	    my $org = $val->org->name;
	    $html  .= qq{  <option value="$id">$org: $sname, version $ver\n};
	    $html =~ s/option/option selected/ unless $count;
	    $count++;
	  }
	$html .= qq{</SELECT>\n};
      }
    else
      {
	$html .= qq{Accession not found <input type="hidden" id="infoid$num">\n};
      }
    
    return $html;
  }

sub generate_obj_from_seq
  {
    my $sequence = shift;
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
					     ACCN=>"RAW_SEQUENCE_SUBMISSION",
					     LOCUS=>"RAW_SEQUENCE_SUBMISSION",
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

    return ($file, $file_begin, $file_end, $spike_seq);
  }

sub get_obj_from_genome_db
  {
    my $accn = shift;
    my $info_id = shift;
    my $up = shift || 0;
    my $down = shift || 0;
    my $db = new CoGe::Genome;
    my @feats;
    my $start = 0;
    my $stop = 0;
    my $feat_start = 0;
    my $chr;
    my $seq;
    foreach my $feat ($db->get_feats_by_name($accn))
      {
	if ($feat->data_information->id() eq $info_id)
	  {
	    push @feats, $feat;
	    if ($feat->type->name eq "gene")
	      {
		$start = $feat->begin_location-$up;
		$start = 1 if $start < 1;
		$stop = $feat->end_location+$down;
		$feat_start = $feat->begin_location;
		$chr = $feat->chr;
		$seq = $db->get_genomic_seq_obj->get_sequence(
							      start => $start,
							      stop => $stop,
							      chr => $chr,
							      org_id => $feat->info->org->id,
							      info_id => $info_id,
							     );
	      }
	  }
      }
    my $obj= new GS::MyDB::GBDB::GBObject(
					  ACCN=>$accn,
					  LOCUS=>$accn,
					  VERSION=>$feats[0]->data_information->version(),
					  SOURCE=>$feats[0]->data_information->data_source->name(),
					  ORGANISM=>$feats[0]->org->name(),
					  					 );
    $obj->sequence($seq);
    my $fnum = 1;
    my %used_names;
    $used_names{$accn} = 1;

    print STDERR "Region: $chr: $start-$stop\n" if $DEBUG;

    foreach my $feat ($db->get_feature_obj->get_features_in_region(start=>$start, stop=>$stop, chr=>$chr, info_id=>$info_id))
      {
	my $name;
	foreach my $tmp (sort {$b->name cmp $a->name} $feat->names)
	  {
	    $name = $tmp->name;
	    last if ($tmp->name =~ /^$accn.+/i);
	  }
	$name = $accn unless $name;
	print STDERR $name,"\n" if $DEBUG;
	print STDERR "\t", $feat->genbank_location_string(),"\n" if $DEBUG;
	print STDERR "\t", $feat->genbank_location_string(recalibrate=>$start),"\n\n" if $DEBUG;
	my $anno = $feat->annotation_pretty_print;
	$anno =~ s/\n/<br>/g;
	$obj->add_feature (
			   F_NUMBER=>$fnum,
			   F_KEY=>$feat->type->name,
			   LOCATION=>$feat->genbank_location_string(recalibrate=>$start),
			   QUALIFIERS=>[
                                        [annotation=>$anno],
                                        [names=> [map {$_->name} sort {$a->name cmp $b->name} $feat->names]],
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
	my $seqfile1 = shift;
	my $seqfile2 = shift;
	my(@param) = @_;
	my $command = $BL2SEQ;
	my $program = "blastn";

	# need to create a temp filename here
	my $tmp_file = new File::Temp ( TEMPLATE=>'CNS__XXXXX',
					DIR=>$TEMPDIR,
					SUFFIX=>'.blast',
					UNLINK=>0);
	my ($tempfile) = $tmp_file->filename;# =~ /([^\/]*$)/;

	# format the bl2seq command
	$command .= "-p $program -o $tempfile ";
	$command .= "-i $seqfile1 -j $seqfile2";

	if ( @param > 0 ) {
		$command .= " " . join( " ", @param );
	}

	my $x = "";
	($x,$command) = check_taint( $command );
	if ( $DEBUG ) {
		print STDERR "About to execute...\n $command\n";
	}
	# execute the command
	`$command`;

	return( $tempfile );
}

