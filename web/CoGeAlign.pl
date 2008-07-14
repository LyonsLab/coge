#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Ajax;
use CoGeX;
use CoGeX::Feature;
use Data::Dumper;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use HTML::Template;
use URI::Escape;
use Spreadsheet::WriteExcel;
use LWP::Simple;
use LWP::Simple::Post qw(post post_xml);


$ENV{PATH} = "/opt/apache/CoGe/";
use vars qw( $TEMPDIR $TEMPURL $USER $DATE $CLUSTAL $BASEFILE $coge $cogeweb $FORM);

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = "/opt/apache/CoGe/tmp/";
$TEMPURL = "/CoGe/tmp/CoGeAlign";
$FORM = new CGI;

$coge = CoGeX->dbconnect();
$CLUSTAL = "/usr/bin/clustalw2";


my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       refresh_seq=>\&refresh_seq,
		       run=>\&run,
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();

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
    $template->param(TITLE=>'ClustalW2 Alignments');
    #$template->param(HELP=>'BLAST');
   # print STDERR "user is: ",$USER,"\n";
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
    $template->param(USER=>$name);
    $template->param(ADJUST_BOX=>1);

    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    $template->param(DATE=>$DATE);
    $template->param(BOX_NAME=>'CoGe: ClustalW2');
    $template->param(BODY=>$body);
    my $prebox = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeAlign.tmpl');
	$prebox->param(RESULTS_DIV=>1);
	$template->param(PREBOX=>$prebox->output);
    $html .= $template->output;
   }
 }
 
 sub gen_body
  {
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/CoGeAlign.tmpl');
    my $form = $FORM;
    my $no_values;
    $BASEFILE = $form->param('basename');
    my $sort_by_type = $form->param('sort_type');
    my $sort_by_location = $form->param('sort_loc');
    my $feat_list = [];
    $feat_list = read_file() if $BASEFILE;#: $opts{feature_list};
    foreach my $item ($form->param('fid'))
      {
	push @$feat_list, $item if $item =~ /^\d+$/;
      }
      
    #print STDERR Dumper \@$feat_list;
    my $dsid = $form->param('dsid') if $form->param('dsid');
    my $chr = $form->param('chr') if $form->param('chr');
    push @$feat_list, @{get_fids_from_dataset(dsid=>$dsid, chr=>$chr)} if $dsid;
    
    my $seqs = generate_sequence(feature_list=>$feat_list);
    my $fid_string;
    foreach my $fid (@$feat_list)
    {
    	$fid_string .= $fid.":";
    }
    $fid_string =~ s/:$//;

    $seqs =~ s/(\sGenomic.+-\d+)?//g;
    if ($seqs)
      {
      	my $num_seqs = $seqs =~ tr/>/>/;
    $template->param(MAIN=>1);  	
	$template->param(SEQUENCE=>$seqs);
	$template->param(FIDS=>$fid_string);
	return $template->output;
      }
    else
      {
	return "No feature ids were specified.";
      }
  }
  
  sub generate_sequence
  {
  	my %opts = @_;
    my $feat_list = $opts{feature_list};
    return unless @$feat_list;
    $feat_list = [map {$coge->resultset("Feature")->find($_)} @$feat_list];
    my $seqs;
  	foreach my $feat(@$feat_list)
    {
      unless ($feat)
	{
#	  warn "feature id $featid failed to return a valid feature object\n";
	  next;
	}
	$seqs .= $feat->fasta(col=>0,prot=>0, name_only=>1);
    }
    return $seqs;
  }
  
  sub refresh_seq
  {
	my %opts = @_;
	my $fids = $opts{fids};
	my $seq_type = $opts{seq_type};
	my $protein = $seq_type =~ /dna/i ? 0 : 1;
	my $seqs;
	
	foreach my $fid (split /:/,$fids)
	{
		my ($feat) = $coge->resultset("Feature")->find($fid);
		$seqs .= $feat->fasta(col=>0,prot=>$protein, name_only=>1);
	}
	 $seqs =~ s/(\sGenomic.+-\d+)?//g;
	 return $seqs;
  }
  
    sub read_file
  {
    my $file = "$TEMPDIR/$BASEFILE.featlist";
    my @featlist;
    unless (-r $file)
      {
	warn "unable to read file $file for feature ids\n";
	return \@featlist;
      }
    open (IN, $file) || die "can't open $file for reading: $!";
    while (<IN>)
      {
	chomp;
	push @featlist,$_;
      }
    close IN;
    return \@featlist;
  }
  
  sub run
  {
  	my %opts = @_;
  	my $seqs = $opts{seq};
  	
  	my $seq_type = $opts{seq_type};
  	my $matrix = $opts{matrix};
  	my $gap_open = $opts{gap_open};
  	my $gap_ext = $opts{gap_ext};
  	my $iteration = $opts{iteration};
  	my $format = $opts{format};
  	#my $num_interations = $opts{num_iterations};
  	
  	my $num_seqs = $seqs =~ tr/>/>/;

  	my %file_format = (NEXUS=>'nxs',PHYLIP=>'phy',GDE=>'gde',PIR=>'pir');
  	#print STDERR Dumper \%file_format;
  	
  	$cogeweb = initialize_basefile(prog=>"CoGeAlign");
  	
  	my $seq_file = $cogeweb->basefile."_clustalw.infile";
  	#print STDERR $seq_file,"\n";
  	
  	open(NEW,"> $seq_file");
    print  NEW $seqs;
    close NEW;
  	#print STDERR $format,"\n";
  	my $suffix = $format =~/(coge|clustal)/ ? 'aln' : $file_format{$format};
  	#print STDERR $suffix," is your suffix!\n";
  	my $outfile = $cogeweb->basefile."_clustalw.".$suffix;
  	my $tree_out = $TEMPURL."/".$cogeweb->basefilename."_clustalw.dnd";
  	my $phylip_file = $TEMPURL."/".$cogeweb->basefilename."_clustalw.ph";
  	
  	my $pre_command = "-INFILE=$seq_file";
  	
  	$pre_command .= " -TYPE=$seq_type";
  	if($format && $format !~/coge/i && $format !~ /clustal/i)
  	{
  	$pre_command .= " -OUTPUT=$format";		
  	}
  	
  	$pre_command .= $seq_type=~/dna/i ? " -DNAMATRIX=$matrix" : " -MATRIX=$matrix";
  	$pre_command .= " -GAPOPEN=$gap_open";
  	$pre_command .= " -GAPEXT=$gap_ext";
  	$pre_command .= " -ITERATION=$iteration";
  	$pre_command .= " -OUTPUTTREE=phylip";
  	#$pre_command .= " -NUMITER=$num_interations";
 
  	
  	 my $x;
    ($x, $pre_command) = check_taint($pre_command);
  	
  	my $command = "$CLUSTAL $pre_command";
  	
  	print STDERR $command, "\n";
  	
  	`$command`;
  	
  	$command .= " -tree";
  	`$command`;
  	
  	my $output;
  	
  	open (IN, $outfile) || die "$!";
    while (<IN>)
      {
	$output .= $_;
      }
    close IN;
    $outfile =~ s/$TEMPDIR/\/CoGe\/tmp\//;
    
     my ($header_html,$seq_html);
     my $no_color = $format =~ /color/ ? 1 : 0;
     ($header_html,$seq_html) = parse_results(clustal=>$output,num_seqs=>$num_seqs, no_color=>$no_color) if $format =~ /coge/i;
#    my $box_template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
#     $box_template->param(BOX_NAME=>"ClustalW Alignment Results");
#     $box_template->param(BODY=>$html);
#     my $outhtml = $box_template->output;
#print STDERR $html,"\n";
	my $html = qq{<font align=center style="font-size:24px;">CLUSTALW Alignment Results</font>};
	if($format =~/coge/)
	{
		$html .= qq{<div id=max_height style="max-height:200px;overflow:auto;"><table class="resultborder"><tr><td valign="top"><div id=alignment_header style="color:black;font-size:18px;">$header_html</div></td><td valign="top"><div id=alignment_seq style="color:black;font-size:18px;overflow:auto;width:100%;">$seq_html</div></table></div>};
  	
	}
	else{
		$output =~s/\n/<br\/>/g;
		$html .= qq{<div align=left class=resultborder style="overflow:auto;width:700px;max-height:700px;"><br/>$output</div>};
		
	}
	$html .= qq{<table align=center><tr>};
  	$html .= qq{<td valign=top><table style="font-size:12px;"><tr><td>Log File:</tr><tr><td>TODO</tr><tr><td>ClustalW Raw Alignment File:</tr><tr><td><a href="$outfile" target="_blank">Download</a></tr><tr><td>ClustalW Raw Tree File:</tr><tr><td><a href="$tree_out" target="_blank">Download</a></td></table></td>};
  	my $select_menu = qq{<select id=file_types><option value=GCG>GCG<option value=GDE>GDE<option value=PHYLIP>Phylip<option value=PIR>PIR<option value=NEXUS>NEXUS</select>};
  	$html .= qq{<td valign=top><table style="font-size:16px;"><tr><td>View NJ Tree of Results:</tr><tr><td><a href="$phylip_file" target="_blank">NJ Tree</a></table></td>};
  	#</tr><tr><td>Select Alt. Output File for Download:</tr><tr><td>$select_menu<input type=button value="Go">
  	$html .= "</table>";
  	
  	return $html;
  }
  
  sub parse_results
  {
  	my %opts = @_;
  	my $clustal = $opts{clustal};
  	my $num_seqs = $opts{num_seqs};
  	my $no_color = $opts{no_color};
  	my @lines;
  	my @headers;
  	my @seqs;
  	
  	my $clustal_title;
  	
  	push @lines, split(/\n/,$clustal);
  	
  	for(my $i=0;$i< scalar @lines;$i++)
  	{
  		if($lines[$i]=~/clustal/i)
  		{
  			$clustal_title=$lines[$i];
  			next;
  		}
  		next unless $lines[$i] =~ /\w/i;
  		#print $lines[$i],"<br>";
  		for(my $j=0;$j<$num_seqs;$j++)
  		{
  			my ($name,$seq) = $lines[$i+$j] =~/(\w+)?\s+([NATGC-]+)/i;
  			#print "i is $i, j is $j <br>$lines[$i+$j]<br>name: $name, seq: $seq<br>";
  			
  			$headers[$j] = $name;
  			$seqs[$j] .= $seq;
  		}
  		$i += $num_seqs-1;
  	}
  	
  	#print STDERR Dumper \@seqs;
  	my ($header_output,$seq_output);
  	for(my $j=0;$j<$num_seqs;$j++)
  		{
			#print STDERR $seqs[$j],"\n";
			unless($no_color)
			{
			$seqs[$j] = uc($seqs[$j]);
			$seqs[$j] =~ s/T/<font style='background-color:DeepSkyBlue'>T<\/font>/g;
			$seqs[$j] =~ s/A/<font style='background-color:red'>A<\/font>/g;
			$seqs[$j] =~ s/G/<font style='background-color:yellow'>G<\/font>/g;
			$seqs[$j] =~ s/C/<font style='background-color:green'>C<\/font>/g;
			}
			#print STDERR $seqs[$j],"\n";
			$header_output .= $headers[$j]."<br/>";
			$seq_output .= $seqs[$j]."<br/>";
			
  		}
  	return $header_output,$seq_output;
  }