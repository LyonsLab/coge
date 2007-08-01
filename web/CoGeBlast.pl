#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use LWP::Simple;
use LWP::Simple::Post qw(post post_xml);
use URI::Escape;
use CoGeX;
use CoGeX::Feature;
use POSIX;
use Digest::MD5 qw(md5_hex);
use DBIxProfiler;
use File::Temp;
use CoGe::Accessory::blast_report;
use CoGe::Graphics::GenomeView;

$ENV{PATH} = "/opt/apache/CoGe/";
$ENV{BLASTDB}="/opt/apache/CoGe/data/blast/db/";
$ENV{BLASTMAT}="/opt/apache/CoGe/data/blast/matrix/";
use vars qw( $TEMPDIR $TEMPURL $DATADIR $FASTADIR $BLASTDBDIR $FORMATDB $BLAST $FORM $USER $DATE $BASEFILENAME $BASEFILE $LOGFILE %restricted_orgs $coge);

$TEMPDIR = "/opt/apache/CoGe/tmp/";
$DATADIR = "/opt/apache/CoGe/data/";
$FASTADIR = $DATADIR.'/fasta/';
$BLASTDBDIR = $DATADIR.'/blast/db/';
$TEMPURL = "/CoGe/tmp/";
$FORMATDB = "/usr/bin/formatdb";
$BLAST = "/usr/bin/blast";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
if (!$USER || $USER =~ /public/i)
  {
    $restricted_orgs{papaya} = 1;
  }
$FORM = new CGI;

my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

#my %ajax = CoGe::Accessory::Web::ajax_func();

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seq=>\&get_seq,
		       cogesearch=>\&cogesearch,
		       source_search=>\&get_data_source_info_for_accn,
		       get_types=>\&get_types,
		       get_sequence=>\&get_sequence,
		       get_url=>\&get_url,
		       check_seq=>\&check_seq,
		       set_seq=>\&set_seq,
		       get_pretty_print=>\&get_pretty_print,
		       blast_param=>\&blast_param,
		       database_param=>\&database_param,
		       get_orgs=>\&get_orgs,
		       get_from_id=>\&get_from_id,
		       blastoff_search=>\&blastoff_search,
		       #%ajax,
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#$USER=1;print gen_html();


sub gen_html
  {
    my $html;
    unless ($USER)
      {
	$html = login();
      }
    else
     {
    my ($body) = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'CoGe: BLAST');
    $template->param(HELP=>'BLAST');
   # print STDERR "user is: ",$USER,"\n";
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"CoGeBlast-logo.png");
    $template->param(BOX_NAME=>'CoGe: Blast');
    $template->param(BODY=>$body);
    my $html;
    $html .= $template->output;
    }
  }

sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    my $form = $FORM;
    my $featid = $form->param('featid') || 0;
    my $chr = $form->param('chr') || 0;
    my $upstream = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $dsid = $form->param('dsid') || 0;
    #my $feat_name = $form->param('featname');
    my $rc = $form->param('rc') || 0;
    $template->param(JAVASCRIPT=>1);
    $template->param(BLAST_FRONT_PAGE=>1);
    $template->param(UPSTREAM=>$upstream);
    $template->param(DOWNSTREAM=>$downstream);
    $template->param(DSID=>$dsid);
    $template->param(ORG_LIST=>get_orgs());
    $template->param(RC=>$rc);
    $template->param(FEATID=>$featid);
    $template->param(CHR=>$chr);
    if ($featid)
    {
    	$template->param(DISPLAY_FEAT=>1);
    	$template->param(SEQVIEW=>1);
    	my $seq = get_sequence($featid, $dsid, 0, 1, $upstream, $downstream,$rc);
        $template->param(SEQUENCE=>$seq);
    }
    elsif ($chr)
    {
    	$template->param(SEQVIEW=>2);
    	my $seq = get_sequence($chr, $dsid, 0, 2, $upstream, $downstream,$rc);
    	$template->param(SEQUENCE=>$seq);
    }
    else{
        $template->param(SEQVIEW=>0);
    }
    $template->param(USER_NAME=>$USER);
    #$template->param(DEFAULT_PARAM=>$param);
    $template->param(REST=>1);
    $template->output;
  }
  
sub get_sequence
  {
    #my $accn = shift;
    my $fid = shift;
    my $dsid = shift;
    my $blast_type = shift || 0;
    my $seqview = shift || 0;
    my $upstream = shift || 0;
    my $downstream = shift || 0;
    my $rc = shift || 0;
    my $seq;
    my $featid;
    #print STDERR "type: ", $fid, ", dsid: ", $dsid, ", blast_type: ", $blast_type, "\n";
    my $fasta;
    if ($seqview == 2)
    {
    	$seq = $coge->resultset('Dataset')->find($dsid)->get_genome_sequence
	  (
	   start=>$upstream,
	   stop=>$downstream,
	   chr=>$fid,
	  );
       $seq = CoGeX::Feature->reverse_complement($seq) if $rc;
       $fasta = generate_fasta_without_featid(chr=>$fid, dsid=>$dsid, start=>$upstream, stop=>$downstream);
       if ($blast_type eq  "blast_type_p") {
       my $key;
       my $sixframe;
       my $sequence = CoGeX::Feature->frame6_trans(seq=>$seq);
       foreach $key (sort {abs($a) <=> abs($b) || $b <=> $a} keys %$sequence)
           {
      	     $seq = join ("\n", wrap('','',$sequence->{$key}));
      	     $sixframe .= qq/$fasta Frame $key\n$seq\n/;
           }
         $seq = $sixframe;
        }
        else {
	$seq = join ("\n", wrap('','',$seq));
	$seq = ($fasta. $seq);
	}
	return $seq;
    }
    else{
    $featid = $fid;
    my $feat = $coge->resultset('Feature')->find($featid);
    #print STDERR "dsid: ", $dsid, "\n";
    #print STDERR "\nfeatid: ", $featid, ", up: ", $upstream, ", down: ", $downstream, "\n";
    $fasta = generate_fasta_with_featid(featid=>$featid, dsid=>$dsid, rc=>$rc, blast_type=>$blast_type);
    
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $fid";
    }
    else 
    {	  
      if ($blast_type eq  "blast_type_p") {
       ($seq) = $feat->protein_sequence;
       $seq = "No sequence available" unless $seq;
      }
      else {
      $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);
      $seq = reverse_complement($seq) if $rc;
    }
      $seq = join ("\n", wrap('','',$seq));
      $seq = ($fasta. $seq);
    return $seq;
    #print STDERR "accn select: ",$featid, ", type name: ",$fid, ", dsid: ", $dsid, "\n";
  }
 }
}

sub get_url
  {
    my $url = shift;
    #print STDERR "expect: ", $expect, "\n";
    if ($url eq "blastn") {
      $url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=on&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    elsif ($url eq "blastp") {
      $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
      elsif ($url eq "blastx") {
        $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=blastx&BLAST_PROGRAMS=blastx&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    elsif ($url eq "tblastn") {
      $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastn&BLAST_PROGRAMS=tblastn&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    elsif ($url eq "tblastx") {
      $url = "http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?PAGE=Translations&PROGRAM=tblastx&BLAST_PROGRAMS=tblastx&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on";}
    else {
      $url = 1;}
    return $url;
  }

sub generate_fasta_with_featid
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my $dsid = $opts{dsid};
    my $rc = $opts{rc};
    my $blast_type = $opts{blast_type};
    my $ds = $coge->resultset("Dataset")->find($dsid);
    my ($feat) = $coge->resultset("Feature")->find($featid);
    my ($strand) = $feat->strand;
    if ($rc)
    	{$strand *= -1 unless ($blast_type eq "blast_type_p");}
    my $fasta = ">".$ds->organism->name."(v.".$feat->version.")".", Type: ".$feat->type->name.", Location: ".$feat->genbank_location_string.", Chromosome: ".$feat->chr.", Strand: ".$strand."\n";
    return $fasta;
  }
  
sub generate_fasta_without_featid
  {
    my %opts = @_;
    my $chr = $opts{chr};
    my $dsid = $opts{dsid};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my ($ds) = $coge->resultset("Dataset")->find($dsid);
    my $fasta = ">".$ds->organism->name.", Location: ".$start."-".$stop.", Chromosome: ".$chr."\n";
    return $fasta;
  }


sub get_pretty_print
{
    my $featid = shift;
    my $dsid = shift;
    my ($feat) = $coge->resultset("Dataset")->find($dsid);
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $featid";
    }
    my $html = $feat->annotation_pretty_print_html();
    return $html;
}

sub blast_param
{
    my $seq_type = shift || "blast_type_n";
    my $pro = shift;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    if ($seq_type eq "blast_type_n") {
        $template->param(BLAST_NU=>1);}
    else {
	$template->param(BLAST_PRO=>1);
	unless ($pro)
	 {$template->param(BLAST_PRO_COMP=>1);}
	 }
    my $html = $template->output;
    return $html;
}

sub reverse_complement
  {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATCG/TAGC/;
    return $seq;
  }
  
sub database_param
  {
    my $program = shift || "blastn";
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
    if ($program eq "blastn")
      {$template->param(NU_DB=>1);}
    elsif (($program eq "blastp") || ($program eq "blastx"))
      {$template->param(PRO_DB=>1);}
    else
      {$template->param(T_DB=>1);}
    my $html = $template->output;
    return $html;
  }

sub get_orgs
  {
    my $name = shift;
    my @db = $name ? $coge->resultset('Organism')->search({name=>{like=>"%".$name."%"}})
      : $coge->resultset('Organism')->all();
    my %restricted;
    if (!$USER || $USER =~ /public/i)
      {
	$restricted{papaya} = 1;
      }
    my @opts;
    foreach my $item (sort {uc($a->name) cmp uc($b->name)} @db)
      {
	next if $restricted{$item->name};
	push @opts, "<OPTION value=\"".$item->id."\" id=\"o".$item->id."\">".$item->name."</OPTION>";
      }
    my $html;
    $html .= qq{<FONT CLASS ="small">Organism count: }.scalar @opts.qq{</FONT>\n<BR>\n};
    unless (@opts) 
      {
	$html .=  qq{<input type = hidden name="org_id" id="org_id">};
	return $html;
      }

    $html .= qq{<SELECT id="org_id" SIZE="8" MULTIPLE onClick="\$('#remove').hide(0);\$('#add').show(0);" ondblclick="get_from_id(['org_id'],[add_to_list]);">\n};
    $html .= join ("\n", @opts);
    $html .= "\n</SELECT>\n";
    $html =~ s/OPTION/OPTION SELECTED/;
    return $html;
  }
  
sub get_from_id
  {
    my $id = shift;
    #print STDERR "search for $id\n";
    my ($obj) = $coge->resultset('Organism')->find($id);
    my $org = $obj->name;
    return ($id,$org);
  }

sub blastoff_search
  {
    my %opts = @_;
    my $program = $opts{program};
    my $expect = $opts{expect};
    my $job_title = $opts{job_title};
    my $wordsize = $opts{wordsize};
    my $comp = $opts{comp};
    my $matrix = $opts{matrix};
    my $gapcost = $opts{gapcost};
    my $match_score = $opts{matchscore};
    my $seq = $opts{seq};
    my $blastable = $opts{blastable};
    initialize_basefile();
    my @org_ids = split(/,/,$blastable);
    my $fasta_file = create_fasta_file($seq);
    my ($nuc_penalty,$nuc_reward,$exist,$extent);
    if ($gapcost =~/^(\d)\s+(\d)/) {($nuc_penalty,$nuc_reward) = ($1,$2);}
    
    if ($match_score=~/^(\d)\,(-\d)/) {($exist,$extent) = ($1,$2);}
    my $pre_command = "$BLAST -p $program -i $fasta_file";
    if ($program =~ /blastn/i)
      {
	$pre_command .= " -q $nuc_penalty -r $nuc_reward";
      }
    else
      {
	$pre_command .= " -M $matrix";
      }
    $pre_command .=" -W $wordsize";
    $pre_command .= " -G $exist -E $extent";
    $pre_command .= " -e $expect";
    $pre_command .= " -C $comp" if $program =~ /tblastn/i;
    my $x;
    ($x, $pre_command) = check_taint($pre_command);
    my @results;
    my $count =1;
    foreach my $orgid (@org_ids)
      {
	my ($db, $org) = get_blast_db($orgid);
	next unless $db;
	my $command = $pre_command." -d $db";
	my $outfile = "$BASEFILE-$count.blast";
	write_log("running $command" ,$LOGFILE);
	`$command > $outfile`;
	my $report = new CoGe::Accessory::blast_report({file=>$outfile}) if -r $outfile;
	my $file = $report->file();
	$file =~ s/$TEMPDIR//;
	$file = $TEMPURL."/".$file;
	push @results, {
			command=>$command,
			file=>$outfile,
			report=>$report,
			link=>$file,
			organism=>$org,
		       };
	$count++;
      }
    my $html = gen_results_page(results=>\@results);
    return $html;
    #return qq{<a href =$TEMPURL/$BASEFILENAME.log target=_new>$LOGFILE</a>};
  }
 
 
sub gen_results_page
   {
     my %opts = @_;
     my $results = $opts{results};
     my @table;
     my $length;
     my $flag;
     my @check;
     foreach my $set (@$results)
     {
      if (@{$set->{report}->hsps()})
      {
       foreach my $hsp (@{$set->{report}->hsps()})
       {
        my ($dsid) = $hsp->subject_name =~ /id: (\d+)/;
        my ($chr) = $hsp->subject_name =~ /chromosome: (\d+)/;
        next unless $dsid && $chr;
        my @feat = $coge->get_features_in_region(start=>$hsp->subject_start,  
				     stop=>$hsp->subject_stop,
				     chr=>$chr, 
				     dataset_id=>$dsid,
				     );
	foreach my $feature (@feat)
	{
	 #print STDERR  $feature->type->name,"\n";
	 next unless $feature->type->name =~ /gene/i;
	 $length = (($feature->stop) - ($feature->start));		     
	 my ($name) = sort $feature->names;
	 foreach my $data (@check)
	 {
	   $flag = 1 if (($data->{name} eq $name) && ($data->{score} == $hsp->score));
	 }
	 unless ($flag) {
	   push @table, {FID=>$feature->id,FEATURE_NAME=>$name,FEATURE_LENGTH=>$length,FEATURE_EVAL=>$hsp->pval,FEATURE_PID=>$hsp->percent_id,FEATURE_SCORE=>$hsp->score,};
	   push @check,{name=>$name,score=>$hsp->score};
	  }
	}
       }
      }
    }     
#     my $chromosome_data = [{DB_NAME=>"AT",CHR_IMAGE=>"Pretty Pictures",},{DB_NAME=>"OS",CHR_IMAGE=>"Goody Pictures",}];
     my ($chromosome_data) = generate_chromosome_images(results=>$results);
     unless (@table) 
     {
     push @table,{FID=>'N/A',FEATURE_NAME=>'N/A',FEATURE_LENGTH=>'N/A',FEATURE_EVAL=>'N/A',FEATURE_PID=>'N/A',FEATURE_SCORE=>'N/A',};
     }
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
     $template->param(OVERLAP_FEATURE_IF=>1);
     $template->param(FEATURE_TABLE=>\@table);
     my $overlap_feature_element = $template->output;
     $template->param(OVERLAP_FEATURE_IF=>0);
     $template->param(OVERLAP_FEATURE_LIST=>$overlap_feature_element);
     $template->param(CHROMOSOMES_IF=>1);
     $template->param(CHROMOSOME_LOOP=>$chromosome_data);
     my $chromosome_element = $template->output;
     $template->param(CHROMOSOMES_IF=>0);
     $template->param(CHROMOSOMES=>$chromosome_element);
     $template->param(BLAST_RESULTS=>1);
     my $html = $template->output;
     return $html;
   }
 

sub generate_chromosome_images
  {
    my %opts = @_;
    my $results = $opts{results};
    my %data;
    my $color;
    my (@data, @no_data);
    foreach my $set (@$results)
      {
	my $org = $set->{organism};
	$color = get_color_scheme($set);
#	print STDERR Dumper $set->{report}->query;
#	print STDERR Dumper $set->{report}->subject;
	$data{$org}{file}=$set->{link};

	if (@{$set->{report}->hsps()})
	  {
	    foreach my $hsp (@{$set->{report}->hsps()})
	      {
		#first, initialize graphic
#		my ($org) = $hsp->subject_name =~ /^(.*?)\(v/;
		$org =~ s/\s+$//;
		my ($chr) = $hsp->subject_name =~ /chromosome: (\d+)/;
		$data{$org}{image} =new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>300, chromosome_height=>25}) unless $data{$org}{image};
		#add chromosome to graphic
		unless ($data{$org}{chr}{$chr})
		  {
		    my ($dsid) = $hsp->subject_name =~ /id: (\d+)/;
		    next unless $dsid;
		    my $ds = $coge->resultset('Dataset')->find($dsid);
		    my $last_pos = $ds->last_chromosome_position($chr);
		    $data{$org}{image}->add_chromosome(name=>"Chr: $chr",
						       end=>$last_pos,
						      );
		  }
		$data{$org}{chr}{$chr}=1;
		my $up = $hsp->strand eq "++" ? 1 : 0;
		$data{$org}{image}->add_feature(name=>$hsp->number,
						start=>$hsp->sstart,
						stop=>$hsp->sstop,
						chr=>"Chr: $chr",
						up=>$up,
						color=>[$color->{$hsp->number}{r},$color->{$hsp->number}{g},$color->{$hsp->number}{b}],
					       );
	      }
	  }
	else
	  {
	  }
      }
    my $count = 1;
    foreach my $org (sort keys %data)
      {
	if ($data{$org}{image})
	  {
	    my $image_file = $BASEFILE."_$count.png";
	    $data{$org}{image}->generate_png(filename=>$image_file);
	    $image_file =~ s/$TEMPDIR/$TEMPURL/;
	    push @data,  {DB_NAME=>"<a href=".$data{$org}{file}. " target=_new>$org</a><br>", CHR_IMAGE=>"<img src=$image_file>"};
	    $count++;
	  }
	else
	  {
	    push @no_data,  {DB_NAME=>"No Hits: <a href=".$data{$org}{file}. " target=_new>$org</a>"};
	  }
      }


    return [@data,@no_data];
  }
sub create_fasta_file
  {
    my $seq = shift;
    write_log("creating user's fasta file",$LOGFILE);
    open(NEW,"> $BASEFILE.fasta");
    print NEW $seq;
    close NEW;
    return "$BASEFILE.fasta";
  }
    

sub get_color_scheme
  {
    my $set = shift;
    my %hsps;
    my $length = 0;
    my %colors;
    my @sorted;
    if (@{$set->{report}->hsps()})
    {
	foreach my $hsp (@{$set->{report}->hsps()})
	{
	  $hsps{$hsp->number} = $hsp->percent_id;
	  @sorted = sort {$hsps{$b} <=> $hsps{$a}} keys %hsps;
	  $length++;
	}
  	print STDERR Dumper \%hsps;
  	print STDERR Dumper \@sorted;
    	my $color = generate_colors($length);
    	my $j = 1;
    	foreach my $hsp (@sorted)
    	{
    	  $colors{$j} = $color->[$hsp];
    	  $j++;
    	}
    }
    return \%colors;
   }

sub generate_colors
  {
    my $length = shift;
    my ($r,$g,$b) = (255,0,0);
    my @colors = ({r=>255,b=>255,g=>255},
    		 {r=>255,b=>0,g=>0},
    		 {r=>0,b=>255,g=>0},
    		 {r=>0,b=>0,g=>255},
    		 {r=>25,b=>150,g=>50});
    return \@colors;
  }
	
sub get_blast_db
  {
    my $orgid = shift;
    my $org = $coge->resultset('Organism')->find($orgid);
    my @ds = $coge->get_current_datasets_for_org(org=>$orgid);
    @ds = sort {$a->id <=> $b->id }@ds;
    return unless @ds;
    my $org_name = $ds[0]->organism->name;
    my $title = $org_name ." (";
    my %vers = map {$_->version,1} @ds;
    if (keys %vers > 1)
      {
	my @chrs;
	foreach my $ds (@ds)
	  {
	    push @chrs, join ("", map {"chr".$_." v".$ds->version." ds:".$ds->id} $ds->get_chromosomes);
	  }
	$title .= join (", ", @chrs);
      }
    else
      {
	$title .= "v".join ("",keys %vers)." ds:".$ds[0]->id();
      }
    $title .= ")";

    my $md5 = md5_hex($title);
    my $file = $FASTADIR."/$md5.fasta";
    my $res;
    if (-r $file)
      {
	write_log("fasta file for $title exists", $LOGFILE);
	$res = 1;
      }
    else
      {
	$res = generate_fasta(dslist=>\@ds, file=>$file) unless -r $file;
      }
    my $blastdb = "$BLASTDBDIR/$md5";
    if (-r $blastdb.".nsq")
      {
	write_log("blastdb file for $title exists", $LOGFILE);
	$res = 1;
      }
    else
      {
	$res = generate_blast_db(fasta=>$file, blastdb=>$blastdb, title=>$title);
      }
    return $blastdb ,$org_name if $res;
    return 0;
    
  }


sub generate_fasta
  {
    my %opts = @_;
    my $dslist = $opts{dslist};
    my $file = $opts{file};
    $file = $FASTADIR."/$file" unless $file =~ /$FASTADIR/;
    write_log("creating fasta file.", $LOGFILE);
    open (OUT, ">$file");
    foreach my $ds (@$dslist)
      {
	foreach my $chr (sort $ds->get_chromosomes)
	  {
	    my $title =  $ds->organism->name." (v". $ds->version.") "."chromosome: $chr".", CoGe database id: ".$ds->id;
	    $title =~ s/^>+/>/;
	    #write_log("adding sequence $title to $file");
	    write_log("adding sequence $title", $LOGFILE);
	    print OUT ">".$title."\n";
	    print OUT $ds->get_genomic_sequence(chr=>$chr),"\n";
	  }
      }
    close OUT;
    write_log("Completed fasta creation", $LOGFILE);
    return 1 if -r $file;
    write_log("Error with fasta file creation", $LOGFILE);
    return 0;
  }

sub generate_blast_db
  {
    my %opts = @_;
    my $fasta = $opts{fasta};
    my $blastdb = $opts{blastdb};
    my $title= $opts{title};
    my $command = $FORMATDB." -p F";
    $command .= " -i '$fasta'";
    $command .= " -t '$title'";
    $command .= " -n '$blastdb'";
    write_log("creating blastdb for $title",$LOGFILE);
    `$command`;
    return 1 if -r "$blastdb.nsq";
    write_log("error creating blastdb for $title",$LOGFILE);
    return 0;
  }

sub initialize_basefile
  {
    my $basename = shift;
    if ($basename)
      {
	my ($x, $cleanname) = check_taint($basename);
	$BASEFILENAME = $cleanname;
	$BASEFILE = $TEMPDIR."/".$cleanname;
	$LOGFILE = $BASEFILE.".log";
      }
    else
      {
	my $file = new File::Temp ( TEMPLATE=>'CoGeBlast_XXXXXXXX',
				    DIR=>$TEMPDIR,
				    #SUFFIX=>'.png',
				    UNLINK=>1);
	($BASEFILE)= $file->filename;
	$LOGFILE = $BASEFILE.".log";
	($BASEFILENAME) = $file->filename =~ /([^\/]*$)/;
      }
    return $BASEFILENAME;
  }
