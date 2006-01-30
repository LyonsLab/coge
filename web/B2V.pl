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
use GS::MyDB::GBDB::GBObject;
use CoGe::Genome;
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
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/B2V.tmpl');
    $template->param(TITLE=>'Bl2Seq Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
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
	$wordsize,
	$gapopen,
	$gapextend,
	$mismatch,
	$blastparams,
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
 	my($filename,$maphtml,$mapname) =
 	  $db->{GBSyntenyViewer}->draw_image(
 					     ACCN_Q_OBJ=>$obj1,
 					     Q_BEGIN=>$file1_begin, 
 					     Q_END=>$file1_end,
 					     ACCN_S_OBJ=>$obj2,
 					     S_BEGIN=>$file2_begin, 
 					     S_END=>$file2_end, 
 					     HSPS=>$data, 
 					     DIRNAME=>$TEMPDIR,
 					     BLAST_REPORT=>$blast_report,
 					     MASKED_EXONS=>$mask_flag,
 					     SPIKE_LEN=>length($spike_seq),
 					    );
	$html .= qq!<IMG SRC="$TEMPURL/$filename" !;
	$html .= qq!BORDER=1 ismap usemap="#$mapname">\n!;
	$html .= "$maphtml\n";
	
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

    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/B2V_results.tmpl');
#    print STDERR $html;
    $template->param(HEADING=>"Results");
    $template->param(RESULTS=>$html);
    return $template->output;

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
					  DEFINITION=>$feats[0]->annotation_pretty_print(),
					 );
    $obj->sequence($seq);
    my $fnum = 1;
    my %used_names;
    $used_names{$accn} = 1;
    foreach my $feat ($db->get_feature_obj->get_features_in_region(start=>$start, stop=>$stop, chr=>$chr, info_id=>$info_id))
      {
	my $name;
	foreach my $tmp (sort {$b->name cmp $a->name} $feat->names)
	  {
	    $name = $tmp->name;
	    last if ($tmp->name =~ /^$accn.+/i);
	  }
	$name = $accn unless $name;
	$obj->add_feature (
			   F_NUMBER=>$fnum,
			   F_KEY=>$feat->type->name,
			   LOCATION=>$feat->genbank_location_string(recalibrate=>$start),
			   QUALIFIERS=>[
                                        [annotation=>$feat->annotation_pretty_print],
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

