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
use CoGe::Graphics;
use CoGe::Graphics::Chromosome;
use CoGe::Graphics::Feature::HSP;
use CoGe::Genome;


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
		       blast_param=>\&blast_param,
		       database_param=>\&database_param,
		       get_orgs=>\&get_orgs,
		       get_from_id=>\&get_from_id,
		       blastoff_search=>\&blastoff_search,
		       generate_feat_info=>\&generate_feat_info,
		       get_hsp_info=>\&get_hsp_info,
		       generate_overview_image=>\&generate_overview_image,
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
        $template->param(SEQUENCE=>'Enter a fasta sequence here');
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
    my $width = $opts{width};
    my $type = $opts{type};
    #print STDERR Dumper \%opts;
    #print STDERR $width,"\n";
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
    my $html = gen_results_page(results=>\@results,width=>$width,type=>$type);
    return $html,$BASEFILENAME;
    #return qq{<a href =$TEMPURL/$BASEFILENAME.log target=_new>$LOGFILE</a>};
  }
 
 
sub gen_results_page
   {
     my %opts = @_;
     my $results = $opts{results};
     my $width = $opts{width};
     my $type = $opts{type};
     #print STDERR $width,"\n";
     my @table;
     my $length;
     my $flag;
     my @check;
     my @no_feat;
     foreach my $set (@$results)
     {
      if (@{$set->{report}->hsps()})
      {
       foreach my $hsp (@{$set->{report}->hsps()})
       {
#        print STDERR Dumper \$hsp;
        my ($dsid) = $hsp->subject_name =~ /id: (\d+)/;
        my ($chr) = $hsp->subject_name =~ /chromosome: (\d+)/;
#	print STDERR $hsp->subject_name,"\n";
        my ($org) = $hsp->subject_name =~ /^\s*(.*?)\s*\(/;
        next unless $dsid && $chr;
        my @feat = $coge->get_features_in_region(start=>$hsp->subject_start,  
				     stop=>$hsp->subject_stop,
				     chr=>$chr, 
				     dataset_id=>$dsid,
				     );
	   if (@feat) {
		foreach my $feature (@feat)
		{
		  next unless $feature->type->name =~ /gene/i;

	 	  $length = (($feature->stop) - ($feature->start));		     
		  my ($name) = sort $feature->names;
	 	  foreach my $data (@check)
	 	  {
	  		 $flag = 1 if (($data->{name} eq $name) && ($data->{score} == $hsp->score));
		  }
	 	  unless ($flag) {
	   	    my $fid = $feature->id."_".$hsp->number;
	   	    my $pid = $hsp->percent_id =~ /\./ ? $hsp->percent_id : $hsp->percent_id.".0";
	   #print STDERR $fid,"\n";
	  	    push @table, {FID=>$fid,FEATURE_NAME=>$name,FEATURE_HSP=>$hsp->number,FEATURE_EVAL=>$hsp->pval,FEATURE_PID=>$hsp->percent_id,FEATURE_SCORE=>$hsp->score,FEATURE_LENGTH=>$length,FEATURE_ORG=>$org,};
	   	    push @check,{name=>$name,score=>$hsp->score};
	  	  }
	 	  $flag=0;
		}
	  }
	  else{
	    push @no_feat, {NO_FEAT=>$hsp->number};
	  }
      }
     }
    }     
#     my $chromosome_data = [{DB_NAME=>"AT",CHR_IMAGE=>"Pretty Pictures",},{DB_NAME=>"OS",CHR_IMAGE=>"Goody Pictures",}];
     my ($chromosome_data, $chromosome_data_large) = generate_chromosome_images(results=>$results,large_width=>$width,hsp_type=>$type);
     #print STDERR Dumper \$chromosome_data;
#     my ($chromosome_data_large) = generate_chromosome_images(results=>$results,width=>$width);
#     for(my $i=0;$i < scalar (@$chromosome_data);$i++)
#     {
#       delete $chromosome_data_large->[$i] if $chromosome_data->[$i]{DB_NAME} =~ /No\s+Hits/i;
#     }
    # print STDERR Dumper \$chromosome_data_large;
     #print STDERR Dumper \@no_feat;
     unless (@table) 
     {
     push @table,{FID=>'N/A',FEATURE_NAME=>'N/A',FEATURE_HSP=>'N/A',FEATURE_ORG=>'N/A',FEATURE_LENGTH=>'N/A',FEATURE_EVAL=>'N/A',FEATURE_PID=>'N/A',FEATURE_SCORE=>'N/A',};
     }
     #table sort!
     @table = sort {$a->{FEATURE_ORG} cmp $b->{FEATURE_ORG} || $a->{FEATURE_NAME} cmp $b->{FEATURE_NAME} || $a->{FEATURE_EVAL} <=> $b->{FEATURE_EVAL} } @table;
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
     $template->param(OVERLAP_FEATURE_IF=>1);
     $template->param(FEATURE_TABLE=>\@table);
     if (@no_feat)
     {
       @no_feat = sort {$a->{NO_FEAT} <=> $b->{NO_FEAT}} @no_feat if @no_feat;
       $template->param(NO_FEAT_IF=>1);
       $template->param(NO_FEATS=>\@no_feat);
     }
     my $overlap_feature_element = $template->output;
     $template->param(NO_FEAT_IF=>0);
     $template->param(OVERLAP_FEATURE_IF=>0);
     $template->param(OVERLAP_FEATURE_LIST=>$overlap_feature_element);
     $template->param(CHROMOSOMES_IF=>1);
     $template->param(CHROMOSOME_LOOP=>$chromosome_data);
     my $chromosome_element = $template->output;
     $template->param(CHROMOSOMES_IF=>0);
     if (scalar(@$chromosome_data_large) > 0)
     {
       $template->param(CHROMOSOMES_LARGE_IF=>1);
       $template->param(CHROMOSOME_LOOP_LARGE=>$chromosome_data_large);
       my $chr_large_element = $template->output;
       $template->param(CHROMOSOMES_LARGE_IF=>0);
       $template->param(CHR_LARGE=>$chr_large_element);
     }
     $template->param(CHROMOSOMES=>$chromosome_element);
     $template->param(BLAST_RESULTS=>1);
     my $html = $template->output;
     return $html;
   }
 

sub generate_chromosome_images
  {
    my %opts = @_;
    my $results = $opts{results};
    my $hsp_type = $opts{hsp_type};
    my $width = $opts{width} || 400;
    my $large_width = $opts{large_width} || 3*$width;
    my $imagefile_name = $opts{filename} || "null";
    my $height = ($width / 16);
    my $large_height = ($large_width / 16) <= 64 ? ($large_width / 16) : 64;
    #print STDERR $width,"\n";
    my $scale = $opts{scale} || 'linear';
    my %data;
    my $filename;
    my ($hsp_info,$max,$min,$length);
    my (@data, @large_data,@no_data);
    foreach my $set (@$results)
    {
	my $org = $set->{organism};
	$data{$org}{file}=$set->{link};
	$filename = $imagefile_name."_*.png";
	if ($imagefile_name ne "null")
	{
	  if (`ls $filename 2>>/dev/null`)
	  {
	    $data{$org}{image} = $imagefile_name;
	    $data{$org}{skip} = 1;
	    next;
	  }
	}
	$filename = $set->{link};
	($hsp_info,$max,$min,$length) = get_color_scheme($set, $hsp_type);
	if (@{$set->{report}->hsps()})
	  {
	    my @hsps;
	    if ($hsp_type eq "eval")
	      {
		@hsps = sort {$b->eval <=> $a->eval} @{$set->{report}->hsps()}
	      }
	    elsif ($hsp_type eq "length")
	      {
		@hsps = sort {$a->length <=> $b->length} @{$set->{report}->hsps()}
	      }
	    elsif ($hsp_type eq "score")
	      {
		@hsps = sort {$a->score <=> $b->score} @{$set->{report}->hsps()}
	      }
	    else
	      {
		@hsps = sort {$a->percent_id <=> $b->percent_id} @{$set->{report}->hsps()}
	      }


	    
	    foreach my $hsp (@hsps)
	      {
		#first, initialize graphic
#		my ($org) = $hsp->subject_name =~ /^(.*?)\(v/;
		$org =~ s/\s+$//;
		my ($chr) = $hsp->subject_name =~ /chromosome: (\d+)/;
		$data{$org}{image} =new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>$width, chromosome_height=>$height}) unless $data{$org}{image};
		$data{$org}{large_image} =new CoGe::Graphics::GenomeView({color_band_flag=>1, image_width=>$large_width, chromosome_height=>$large_height}) unless $data{$org}{large_image};
		
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
		my $num = $hsp->number;
		my $up = $hsp->strand eq "++" ? 1 : 0;
		my $r = generate_colors(max=>$max,
					min=>$min,
					length=>$length,
					scale=>$scale,
					val=>$hsp_info->{$hsp->number},
					);
		my $b = generate_colors(max=>$max,
					min=>$min,
					length=>$length,
					scale=>$scale,
					val=>$hsp_info->{$hsp->number},
					reverse_flag=>1,
					);
		#Reverse color scheme for eval, as less is more
		($r, $b) = ($b, $r) if $hsp_type eq "eval";
		$data{$org}{image}->add_feature(name=>$hsp->number,
						start=>$hsp->sstart,
						stop=>$hsp->sstop,
						chr=>"Chr: $chr",
						imagemap=>qq/class="imagemaplink" title="HSP No. /.$hsp->number.qq/" onclick="\$('#big_picture').slideToggle(pageObj.speed);show_hsp_div();loading('image_info','Information');loading('query_image','Image');loading('subject_image','Image');get_hsp_info(['args__blastfile','args__$filename','args__num','args__$num'],['image_info','query_image','subject_image']);"/,
						up=>$up,
						color=>[$r,0,$b],
					       );
	      }
	  }
      }
    my $count = 1;
    foreach my $org (sort keys %data)
    {
	if ($data{$org}{image})
	  {
	    my $image_file;
	    my $image_map;
	    my $large_image_file;
	    my $image_map_large;
	    unless ($data{$org}{skip})	    
	      {
		my $x;
		$large_image_file = $BASEFILE."_".$hsp_type."_$count"."_large.png";
		($x, $large_image_file) = check_taint($large_image_file);
		$image_file = $BASEFILE."_".$hsp_type."_$count.png";
		($x, $image_file) = check_taint($image_file);
		$data{$org}{image}->generate_png(filename=>$image_file);
		$image_map = $data{$org}{image}->generate_imagemap(mapname=>$BASEFILENAME."_".$count);

		my $map_file = "$TEMPDIR/$BASEFILENAME"."_$count.$hsp_type.map";
		($x, $map_file) = check_taint($map_file);
		open (MAP, ">$map_file");
		print MAP $image_map;
		close MAP;
		$data{$org}{image}->image_width($large_width);
		$data{$org}{image}->chromosome_height($large_height);
		$data{$org}{image}->generate_png(filename=>$large_image_file);
		$image_map_large = $data{$org}{image}->generate_imagemap(mapname=>$BASEFILENAME."_".$count."_large");
		$map_file = "$TEMPDIR/$BASEFILENAME"."_$count.$hsp_type.large.map";
		($x, $map_file) = check_taint($map_file);
		open (MAP, ">$map_file");
		print MAP $image_map_large;
		close MAP;

	    }
	    else
	      {
		my $x;
		$image_file = $data{$org}{image}."_$count.png";
		$image_map = get_map("$TEMPDIR/$BASEFILENAME"."_$count.$hsp_type.map");
		$large_image_file = $data{$org}{image}."_$count"."_large.png";
		$image_map_large = get_map("$TEMPDIR/$BASEFILENAME"."_$count.$hsp_type.large.map");
		print STDERR $image_map_large,"\n";
		($x, $image_file) = check_taint($image_file);
		($x, $large_image_file) = check_taint($large_image_file);
	      }
	    
	    $image_file =~ s/$TEMPDIR/$TEMPURL/;
	    $large_image_file =~ s/$TEMPDIR/$TEMPURL/;
	    
	    push @large_data,  {DB_NAME_LARGE=>"<a href=".$data{$org}{file}. " target=_new>$org</a><br>", CHR_IMAGE_LARGE=>"<img src=$large_image_file ismap usemap='$BASEFILENAME"."_"."$count"."_large' border=0>$image_map_large",IMAGE_ID_LARGE=>$count,};
	    push @data,  {DB_NAME=>"<a href=".$data{$org}{file}. " target=_new>$org</a><br>", CHR_IMAGE=>"<img src=$image_file ismap usemap='$BASEFILENAME"."_"."$count' border=0>$image_map",HIT=>1,IMAGE_ID=>$count,};
	    $count++;
	  }
	else
	  {
	    push @no_data,  {DB_NAME=>"No Hits: <a href=".$data{$org}{file}. " target=_new>$org</a>"};
	  }
      }


    return [@data,@no_data], \@large_data;
  }

sub get_map
  {
    my $file = shift;
    my $map;
    open (IN, $file) || die "$!";;
    while (<IN>)
      {
	$map .= $_;
      }
    close IN;
    return $map;
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
    my $type = shift || 0;
    my %hsps;
    my %info;
    my @sorted_vals;
    my $max;
    my $min;
    my $code_ref;
    if ($type eq "eval")
      {
	$code_ref = sub {$_->eval;};
      }
    elsif ($type eq "score")
      {
	$code_ref = sub {$_->score;};
      }
    elsif ($type eq "length")
      {
	$code_ref = sub {$_->query_length;};
      }
    else
      {
	$code_ref = sub {$_->percent_id;};
      }
    if (@{$set->{report}->hsps()})
      {
	%hsps = map {$_->number, &$code_ref} @{$set->{report}->hsps()};
#	print STDERR Dumper \%hsps;
	@sorted_vals = sort{$a<=> $b} values %hsps;
#	while($sorted_vals[0] == 0) {shift @sorted_vals;}
	$min = $sorted_vals[0];
	$max = $sorted_vals[-1];
      }
    my $length = keys %hsps;

    return \%hsps,$max,$min,$length;
   }


sub generate_colors
  {
    my %opts = @_;
    my $max = $opts{max};
    my $min = $opts{min};
    my $length = $opts{length};
    my $scale = $opts{scale} || 'linear';
    my $val = $opts{val};
    my $color_max = $opts{color_max} || 255;
    my $color_min = $opts{color_min} || 50;
    my $flag = $opts{reverse_flag} || 0;
#    print STDERR Dumper \%opts;
    my $color;
    ($color_max,$color_min) = ($color_min,$color_max) if $flag;
    return $color_max if ($val >= $max);
    return $color_min if ($val <= $min);
    if ($scale =~ /log/)
      {
        if ($val == 0) {$color = $color_max;}
        else{
          $color = log($val/$min)/log($max/$min)*($color_max-$color_min)+$color_min;
          }
      }
    else #linear
      {
	$color = ($val-$min)/($max-$min)*($color_max-$color_min)+$color_min;
      }
   # print STDERR "Color: $color\n";
    return $color;
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

sub generate_feat_info 
  {
    my $featid = shift;
    $featid =~ s/^table_row//;
    $featid =~ s/_\d+$//;
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $featid";
    }
    my $html = qq{<a href="#" onClick="\$('#overlap_box').slideToggle(pageObj.speed);" style="float: right;"><img src='/CoGe/picts/delete.png' width='16' height='16' border='0'></a>};
    $html .= $feat->annotation_pretty_print_html();
    return $html;
  }

sub get_hsp_info
  {
    my %opts = @_;
    my $hsp_num = $opts{num};
    my $filename = $opts{blastfile};
    
    $filename = "/opt/apache".$filename;
    ($BASEFILENAME) = $filename=~ /$TEMPDIR\/*(CoGeBlast_\w+)-/;
    
    my $report = new CoGe::Accessory::blast_report({file=>$filename});# if -r $filename;
    my $hsp = $report->hsps->[$hsp_num-1];
    #print STDERR Dumper \$hsp;
    
    #my $query_loc = $hsp->query_start."-".$hsp->query_stop;
    #$query_loc = $hsp->query_start."-<br>".$hsp->query_stop if (length $query_loc > 8);
   # my $query_start = $hsp->query_start;
   # my $query_stop = $hsp->query_stop;
    my $query_length = $hsp->query_stop > $hsp->query_start ? (($hsp->query_stop) - ($hsp->query_start) + 1) : (($hsp->query_start) - ($hsp->query_stop) + 1);
    my $query_loc = $hsp->query_stop > $hsp->query_start ? $hsp->query_start."-".$hsp->query_stop : $hsp->query_stop."-".$hsp->query_start;
    my $query_mismatch = $query_length - $hsp->match;
    
    my $query_name = "<pre>".$hsp->query_name."</pre>";
    $query_name = wrap('','',$query_name);
    $query_name =~ s/\n/<br>/g;
    
    #my $subject_start = $hsp->subject_start;
   # my $subject_stop = $hsp->subject_stop;
   # my $subject_loc = $hsp->subject_start."-".$hsp->subject_stop;
   # $subject_loc = $hsp->subject_start."-<br>".$hsp->subject_stop if (length $subject_loc > 8);
    my $subject_length = $hsp->subject_stop > $hsp->subject_start ? (($hsp->subject_stop) - ($hsp->subject_start) + 1) : (($hsp->subject_start) - ($hsp->subject_stop) + 1);
    my $subject_loc = $hsp->subject_stop > $hsp->subject_start ? $hsp->subject_start."-".$hsp->subject_stop : $hsp->subject_stop."-".$hsp->subject_start;
    my $subject_mismatch = $subject_length - $hsp->match;
    
    my $subject_name = "<pre>".$hsp->subject_name."</pre>";
    $subject_name = wrap('','',$subject_name);
    $subject_name =~ s/\n/<br>/g;
    
    my @table1 = ({HSP_EVAL_QUERY=>$hsp->pval,
		   HSP_PID_QUERY=>$hsp->percent_id,
		   HSP_PSIM_QUERY=>$hsp->percent_sim,
		   HSP_GAP_QUERY=>$hsp->query_gaps,
		   HSP_EVAL_SUB=>$hsp->pval,
		   HSP_PID_SUB=>$hsp->percent_id,
		   HSP_PSIM_SUB=>$hsp->percent_sim,
		   HSP_GAP_SUB=>$hsp->subject_gaps,
		   HSP_SCORE_SUB=>$hsp->score,
		   HSP_SCORE_QUERY=>$hsp->score,
		  });
    my @table2 = ({HSP_STRAND_QUERY=>$hsp->strand,
		   HSP_LENGTH_QUERY=>$hsp->length,
		   HSP_STRAND_SUB=>$hsp->strand,
		   HSP_LENGTH_SUB=>$hsp->length,
		   HSP_MATCH_QUERY=>$hsp->match,
		   HSP_POSITION_QUERY=>$query_loc,
		   HSP_POSITION_SUB=>$subject_loc,
		   HSP_MISMATCH_QUERY=>$query_mismatch,
		   HSP_MATCH_SUB=>$hsp->match,
		   HSP_MISMATCH_SUB=>$subject_mismatch,
    		  });
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
     $template->param(HSP_IF=>1);
     $template->param(HSP_NUM=>$hsp->number);
     $template->param(HSP_QUERY=>\@table1);
     $template->param(HSP_SUB=>\@table2);
     
     my $query_seq = $hsp->query_alignment;
     $query_seq = wrap('','',$query_seq);
     my @query = split(/\n/,$query_seq);
     $query_seq =~ s/[^atgc]//ig;
     $query_seq = wrap('','',$query_seq);
     $query_seq =~ s/\n/<br>/g;
     $query_seq =~ tr/atgc/ATGC/;
     $query_seq = qq{<pre>$query_seq</pre>};
     
     my $sub_seq = $hsp->subject_alignment;
     $sub_seq = wrap('','',$sub_seq);
     my @sub = split(/\n/, $sub_seq);
     $sub_seq =~ s/[^atgc]//ig;
     $sub_seq = wrap('','',$sub_seq);
     $sub_seq =~ s/\n/<br>/g;
     $sub_seq =~ tr/atgc/ATGC/;
     $sub_seq = qq{<pre>$sub_seq</pre>};
     
     my ($sub_chr) = $hsp->subject_name =~ /chromosome: (\d+)/;
     my ($sub_dsid) = $hsp->subject_name =~ /id: (\d+)/;
    
     my $alignment = $hsp->alignment;
     $alignment = wrap('','',$alignment);
     my @align = split(/\n/,$alignment);
     my $align_str = "";
     for(my $i=0;$i<scalar(@sub);$i++)
     {
       $align_str .= $query[$i]."<br>".$align[$i]."<br>".$sub[$i]."<br>";
     }
     $align_str =~ s/<br>$//;
     $align_str = "<pre>$align_str</pre>";
     
     $template->param(QUERY_SEQ=>qq{<a href="#" onclick="show_seq('$query_seq','$query_name',1,'seqObj','seqObj','}.$hsp->query_start."','".$hsp->query_stop.qq{')">Click for Query Sequence</a>});
     
     $template->param(SUB_SEQ=>qq{<a href="#" onclick="show_seq('$sub_seq','$subject_name',2,'$sub_dsid','$sub_chr','}.$hsp->subject_start."','".$hsp->subject_stop.qq{')">Click for Subject Sequence</a>});
     
     $template->param(ALIGNMENT=>qq{<a href="#" onclick="show_seq('$align_str','N/A',0,0,0,0)">Click for Alignment Sequence</a>});
     
     my $html = $template->output;
     $template->param(HSP_IF=>0);
     my ($query_image, $subject_image) = generate_hit_image(report=>$report, hsp_num=>$hsp_num, hsp=>$hsp);
     my $query_link = "<img src=$query_image border=0>";
     $query_link =~ s/$TEMPDIR/$TEMPURL/;
     my $subject_link = "<img src=$subject_image border=0>";
     $subject_link =~ s/$TEMPDIR/$TEMPURL/;
#    print STDERR $query_link,"\n", $query_image,"\n";
     return $html, $query_link, $subject_link;
  }
  
sub generate_overview_image
  {
     my %opts = @_;
     my $basename = $opts{basename};
     my $type = $opts{type};
     my $image_width = $opts{width};
     my @set = split/\n/, `ls $TEMPDIR/$basename*.blast`;
     my @reports;
     my $count = 1;
     $BASEFILE = $TEMPDIR."/".$basename;
     $BASEFILENAME = $basename;
     foreach my $blast (@set){
       my $report = new CoGe::Accessory::blast_report({file=>$blast});
       my ($org_name) = $report->hsps->[$count-1]->subject_name =~ /^\s*(.*?)\s*\(/;
       push @reports,{report=>$report,organism=>$org_name,link=>$TEMPURL."/".$basename."-".$count.".blast",};
       $count++;
     }
     #print STDERR Dumper \@reports;
     my $image_filename = $BASEFILE."_".$type;
     my ($chromosome_data, $chromosome_data_large) = generate_chromosome_images(results=>\@reports,hsp_type=>$type,large_width=>$image_width,filename=>$image_filename);
     my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeBlast.tmpl');
     $template->param(CHROMOSOMES_IF=>1);
     #print STDERR Dumper $chromosome_data;
     $template->param(CHROMOSOME_LOOP=>$chromosome_data);
     my $chromosome_element = $template->output;
     $template->param(CHROMOSOMES_IF=>0);
     my $chr_large_element;
     if (scalar(@$chromosome_data_large) > 0)
     {
       $template->param(CHROMOSOMES_LARGE_IF=>1);
       $template->param(CHROMOSOME_LOOP_LARGE=>$chromosome_data_large);
       $chr_large_element = $template->output;
       $template->param(CHROMOSOMES_LARGE_IF=>0);
     }
     my $html = $chromosome_element.$chr_large_element;
     #print STDERR $html,"\n";
     return $html;
  }
     

sub generate_hit_image
  {
    my %opts = @_;
    my $report = $opts{report};
    my $hsp_num = $opts{hsp_num};
    my $file = $opts{file};
    my $hsp = $opts{hsp};
    my $width = $opts{width} || 400;
    return unless ($report || -r $file);
    $report = new CoGe::Accessory::blast_report({file=>$file}) unless $report;
    unless ($hsp)
      {
	foreach my $item (@{$report->hsps})
	  {
	    next unless $item->number eq $hsp_num;
	    $hsp = $item;
	    last;
	  }
      }
#    print STDERR Dumper $hsp, $report;
    #generate_query_image
    my $c = new CoGe::Graphics::Chromosome ();
    $c->chr_length($hsp->query_length);

    $c->iw($width);
    $c->draw_chromosome(1);
    $c->draw_ruler(1);
    $c->draw_chr_end(0);
    $c->mag(0);
    $c->mag_off(1);
    $c->minor_tick_labels(1);
    $c->draw_hi_qual(0);
    $c->set_region(start=>1, stop=>$c->chr_length);
    my $strand = $hsp->strand =~ /-/ ? "-1" : 1;
    my $feat = CoGe::Graphics::Feature::HSP->new({start=>$hsp->qstart, stop=>$hsp->qstop, strand=>$strand});
    $feat->color([255,0,0]);
    $c->add_feature($feat);
    my $query_file = $TEMPDIR."/".$BASEFILENAME.".q.".$hsp->number.".png";
    $c->generate_png(file=>$query_file);
    $c = new CoGe::Graphics::Chromosome ();
    my $graphic = new CoGe::Graphics;
    my ($dsid) = $hsp->subject_name =~ /id: (\d+)/;
    my ($chr) = $hsp->subject_name =~ /chromosome: (\d+)/;
    my $len = $hsp->subject_stop - $hsp->subject_start+1;
    my $start = $hsp->subject_start-2*$len;
    $start = 1 if $start < 1;
    my $stop = $hsp->subject_stop+2*$len;    
    $graphic->initialize_c (
			   ds=>$dsid,
			   chr=>$chr,
			   c=>$c,
			   iw=>$width,
			   start=> $start,
			   stop => $stop,
			   draw_chr=>1,
			   draw_ruler=>1,
			   draw_chr_end=>0,
			   #			    chr_start_height=>$ih,
			   #			    chr_mag_height=>5,
			   #			    feature_start_height=>$fh,
			   mag=>0,
			   mag_off=>1,
			   chr_length => $len,
			   fill_labels=>1,
			   forcefit=>1,
			   minor_tick_labels=>1,
			   #			    overlap_adjustment=>$overlap_adjustment,
			   #			    feature_labels=>$feature_labels,
			   draw_hi_qual=>0,
			   #			    padding=>$padding,
			  );
    my $db = new CoGe::Genome;
#    $graphic->process_features(c=>$c, layers=>{all=>1}, db=>$db);
    my $sub_file = $TEMPDIR."/".$BASEFILENAME.".s.".$hsp->number.".png";
    $c->generate_png(file=>$sub_file);
    return $query_file, $sub_file;
  }
