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
use CoGe::Accessory::genetic_code;

$ENV{PATH} = "/opt/apache/CoGe/";
$ENV{THREADS} =8;
use vars qw( $TEMPDIR $TEMPURL $USER $DATE $CLUSTAL $BASEFILE $coge $cogeweb $FORM);

$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$TEMPDIR = "/opt/apache/CoGe/tmp/";
$TEMPURL = "/CoGe/tmp/CoGeAlign";
$FORM = new CGI;

$coge = CoGeX->dbconnect();
$CLUSTAL = "/usr/bin/clustalw2";
#$CLUSTAL = "/usr/bin/clustalw-mtv";


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
	$template->param(LOGO_PNG=>"CoGeAlign-logo.png");
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
	$template->param(codon_align=>1);
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
	$seqs .= $feat->fasta(col=>0,prot=>0, name_only=>1, add_fid=>1);
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
		$seqs .= $feat->fasta(col=>0,prot=>$protein, name_only=>1, add_fid=>1);
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
    my $inseqs = $opts{seq};
    
    my $seq_type = $opts{seq_type};
    my $matrix = $opts{matrix};
    my $gap_open = $opts{gap_open};
    my $gap_ext = $opts{gap_ext};
    my $iteration = $opts{iteration};
    my $format = $opts{format};
    my $gen_matrix = $opts{gen_matrix};
    my $codon = 0;
    if ($seq_type eq "prot")
      {$seq_type = "PROTEIN";}
    elsif ($seq_type eq "codon")
      {
	$codon =1;
	$seq_type = "PROTEIN";
      }
    else {$seq_type = "DNA";}
    #my $num_interations = $opts{num_iterations};
    
    my $num_seqs = $inseqs =~ tr/>/>/;
    
    my %file_format = (NEXUS=>'nxs',PHYLIP=>'phy',GDE=>'gde',PIR=>'pir');
    #print STDERR Dumper \%file_format;
    
    $cogeweb = initialize_basefile(prog=>"CoGeAlign");
    
    my $seq_file = $cogeweb->basefile."_clustalw.infile";
    #print STDERR $seq_file,"\n";
    
    open(NEW,"> $seq_file");
    print  NEW $inseqs;
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
    
    my $color = $format =~ /color/ ? 1 : 0;
    my ($header_html,$seq_html, $seqs, $codon_alignment) = parse_results(clustal=>$output,num_seqs=>$num_seqs, color=>$color, codon_align=>$codon) if $format =~ /coge/i;
#    print STDERR Dumper $codon_alignment;
    #    my $box_template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/box.tmpl');
    #     $box_template->param(BOX_NAME=>"ClustalW Alignment Results");
    #     $box_template->param(BODY=>$html);
    #     my $outhtml = $box_template->output;
    #print STDERR $html,"\n";
    my $html = qq{<font style="font-size:24px;">CLUSTALW Alignment Results</font>};
    if($format =~/coge/)
      {
	$html .= qq{<div id=max_height style="max-height:700px;overflow:auto;"><table class="resultborder"><tr><td valign="top"><div id=alignment_header style="color:black;"><pre>$header_html</pre></div></td><td valign="top"><div id=alignment_seq style="color:black;;overflow:auto;"><pre>$seq_html</pre></div></table></div>};
  	
      }
    else{
      $output =~s/\n/<br\/>/g;
      $html .= qq{<div align=left class=resultborder style="overflow:auto;width:700px;max-height:700px;"><br/>$output</div>};
      
    }
    $html .= qq{<table align=left><tr>};
    $html .= qq{<td valign=top><table style="font-size:12px;"><tr><td>Log File:</tr><tr><td>TODO</tr><tr><td>ClustalW Raw Alignment File:</tr><tr><td><a href="$outfile" target="_blank">Download</a></tr><tr><td>ClustalW Raw Tree File:</tr><tr><td><a href="$tree_out" target="_blank">Download</a></td></table></td>};
    my $select_menu = qq{<select id=file_types><option value=GCG>GCG<option value=GDE>GDE<option value=PHYLIP>Phylip<option value=PIR>PIR<option value=NEXUS>NEXUS</select>};
    $html .= qq{<td valign=top><table style="font-size:16px;"><tr><td>View NJ Tree of Results:</tr><tr><td><a href="$phylip_file" target="_blank">NJ Tree</a></table></td>};
    #</tr><tr><td>Select Alt. Output File for Download:</tr><tr><td>$select_menu<input type=button value="Go">
    $html .= "</table>";
    if ($gen_matrix)
      {
	if ($codon)
	  {
	    my $matrix = gen_matrix(seqs=>$seqs, type=>'p');
	    $html .= "<br>".$matrix if $matrix;
	    $matrix = gen_matrix(seqs=>$codon_alignment, type=>'c');
	    $html .= "<br>".$matrix if $matrix;
	  }
	else
	  {
	    my $matrix = gen_matrix(seqs=>$seqs, type=>$seq_type);
	    $html .= "<br>".$matrix if $matrix;
	  }
      }
    return $html;
  }
  
sub parse_results
  {
    my %opts = @_;
    my $clustal = $opts{clustal};
    my $num_seqs = $opts{num_seqs};
    my $color = $opts{color};
    my $codon_align = $opts{codon_align};
    my $matrix = $opts{gen_matrix};
    my @lines;
    my @headers;
    my @seqs;
    
    my $clustal_title;
    my %seqs;
    my $pad=0;
    foreach my $line (split /\n/,$clustal)
      {
	next unless $line;
	next if $line =~ /^CLUSTAL/;
	if ($line =~ /^\s/)
	  {
	    $line = substr($line,$pad);
	    $seqs{alignment_stuff}.= $line;
	    next;
	  }
	unless ($pad)
	  {
	    my ($tmp) = $line=~ /(^\S+\s+)/;
	    $pad = length $tmp;
	  }
	my ($name, $seq) = split /\s+/,$line,2;
	push @headers, $name unless $seqs{$name};
	$seqs{$name}.= $seq;
      }
    my ($header_output,$seq_output);
    my %codon_alignment;
    foreach my $name (@headers)
      {
	my $seq = uc($seqs{$name});
	if ($codon_align)
	  {
	    my ($fid) = $name =~/fid_(\d+)/;
	    if ($fid)
	      {
		my $feat = $coge->resultset('Feature')->find($fid);
		my $tseq = $feat->genomic_sequence;
		my $pos = 0;
		my $dna;
		my $prot;
		foreach my $chr (split //,$seqs{$name})
		  {
		    if ($chr eq "-")
		      {
			$dna .= "-"x3 ." ";
		      }
		    else
		      {
			$dna .= substr($tseq,$pos,3)." ";
			$pos+=3;
		      }
		    $prot .= "  $chr ";
		  }
		$seq = $prot."<br>".$dna;
		$dna =~ s/\s//g;
		$codon_alignment{$name}=$dna;
		$name .= "<br>";
	      }
	    
	  }
	if($color)
	  {
	    $seq =~ s/(T+)/<span style='background-color:deepskyblue'>$1<\/span>/g;
	    $seq =~ s/(A+)/<span style='background-color:red'>$1<\/span>/g;
	    $seq =~ s/(G+)/<span style='background-color:yellow'>$1<\/span>/g;
	    $seq =~ s/(C+)/<span style='background-color:green'>$1<\/span>/g;
	  }
	$header_output .= $name."<br/>";
	$seq_output .= $seq."<br/>";
	
      }
    $header_output .= "<br/>";
    $codon_alignment{alignment_stuff} = "  ".join ("  ", (split//,$seqs{alignment_stuff})) if $codon_align;
    $seq_output .= $codon_align ? "<pre>"."  ".join ("  ", map{$_." "}(split//,$seqs{alignment_stuff}))."</pre><br/>" : "<pre>".$seqs{alignment_stuff}."</pre><br/>";
    return $header_output,$seq_output, \%seqs, \%codon_alignment;
  }

sub gen_matrix
  {
    my %opts = @_;
    my $seqs = $opts{seqs};
    my $type = $opts{type};
    my %sub_count;
    my $aln = $seqs->{alignment_stuff};
    delete $seqs->{alignment_stuff};
    my $length = length $aln;
    my %freq;    
    my %total_chrs;
    #log-odds score based on Henikoff and Henikoff 1992 paper:  Amino acid substitution matrices from protein blocks
    my $add = $type =~ /^c/ ? 3 : 1;
    for (my $i=0; $i < $length; $i+=$add)
      {
	my %chrs;
	if ($type =~ /^c/)
	  {
	    map{$chrs{$_}++} map {substr($_,$i,3)} values %$seqs;
	  }
	else
	  {
	    map{$chrs{$_}++} map {substr($_,$i,1)} values %$seqs;
	  }
	my @chrs = keys %chrs;
	foreach my $c1 (@chrs)
	  {
	    next if $c1 eq "-" || $c1 eq "---"; #skip gaps
	    #skip non ATCG DNA characters in codons/DNA
	    next if ($type =~ /^c/i || $type =~ /^d/i) && $c1 =~ /[^ATCG]/;
	    
	    #calculate pair count for identical match and add to total
	    my $x = $chrs{$c1}-1;
	    $freq{$c1}{$c1}+= $x*($x+1)/2;
	    #add count to total character count
	    $total_chrs{$c1} += $chrs{$c1};
	    foreach my $c2 (@chrs)
	      {
		next if $c2 eq "-";
		next if $c1 eq $c2; #already done
		#calculate pair count for non-indentical matches
		$freq{$c1}{$c2}+= $chrs{$c2}*$chrs{$c1} unless $freq{$c2}{$c1};
		
	      }
	  }
      }
    my $total_pairs=0;
    foreach my $c1 (keys %total_chrs)
      {
	foreach my $c2 (keys %total_chrs)
	  {
	    $total_pairs += $freq{$c1}{$c2} if $freq{$c1}{$c2};
	    
	  }
      }
    my %q;
    foreach my $c1 (keys %total_chrs)
      {
	foreach my $c2 (keys %total_chrs)
	  {
	    $q{$c1}{$c2} = $freq{$c1}{$c2}/$total_pairs if $freq{$c1}{$c2};
	  }
      }
    my %p;
    foreach my $c1 (keys %total_chrs)
      {
	$p{$c1} = $q{$c1}{$c1};
	foreach my $c2 (keys %total_chrs)
	  {
	    next if $c1 eq $c2;
	    $p{$c1} += $q{$c1}{$c2} if $q{$c1}{$c2};
	  }
      }
#    return "<pre>".Dumper (\%q)."</pre>";
    my %matrix;
    foreach my $c1 (keys %total_chrs)
      {
	foreach my $c2 (keys %total_chrs)
	  {
	    my $q = $q{$c1}{$c2} ? $q{$c1}{$c2} : $q{$c2}{$c1};
	    next unless $q;
#	    print STDERR $c1,"::",$c2," ",$q{$c1}{$c2}," ",$q{$c2}{$c1}," ",$q/($p{$c1}*$p{$c2}),"\n";
	    my $val = sprintf("%.0f",2*log($q/($p{$c1}*$p{$c2})));
	    $val = 0 if $val eq "-0";
	    $matrix{uc($c1)}{uc($c2)} = $val;
	  }
      }
    $seqs->{alignment_stuff} = $aln;
    my $html = gen_matrix_output(matrix=>\%matrix, type=>$type);
#    my $html = "<pre>".Dumper (\%matrix)."</pre>";
    return $html;
  }

sub gen_matrix_output
  {
    my %opts = @_;
    my $data = $opts{matrix};
    my $type = $opts{type};
    my $code = CoGe::Accessory::genetic_code->code;
    $code = $code->{code};
    my $html ="<table>";
    if ($type =~ /^p/i) #proteins
      {
	my $aa_sort = CoGe::Accessory::genetic_code->sort_aa_by_gc();
	foreach my $c1 (keys %$aa_sort)
	  {
	    foreach my $c2 (keys %$aa_sort)
	      {
		$data->{$c1}{$c2}="-50" unless defined $data->{$c1}{$c2};
	      }
	  }
	$html .= "<tr><th><th>".join ("<th>", sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort);
	$html.="<th>Total:";
	$html .= "<tr>";
	foreach my $aa1 (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
	  {	
	    $html .= "<th>$aa1";
	    #	print STDERR join ("\t",  sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %aa_sort),"\n";
	    my %vals;
	    map {$vals{$_}++}  map {$data->{$aa1}{$_}} keys %$aa_sort;
	    my ($max) = sort {$b<=>$a} keys %vals;
	    my ($min1, $min2) = sort {$a<=>$b} keys %vals;
	    $min2 = $min1 unless $min2;
	    my $total = 0;
	    foreach my $aa2 (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
	      {	
		my $val = $data->{$aa1}{$aa2};
		$total+=$val if $val =~ /^\d+$/;
		my $color;
		if ($val >0)
		  {
		    $color = color_by_usage($max, $val);
		    $html .= "<td style=\"background-color: rgb($color,255,$color)\">".$val;
		  }
		else
		  {
		    $color = color_by_usage(abs($min2), $val*-1);
		    $html .= "<td style=\"background-color: rgb(255,$color,$color)\">".$val;
		  }
	      }
	    $html .= "<td>$total<tr>";
	  }
      }
    elsif ($type =~ /^c/i) #codons
      {
	my @dna = qw(A T C G);
	my @codons;
	foreach my $c1 (@dna)
	  {
	    foreach my $c2 (@dna)
	      {
		foreach my $c3 (@dna)
		  {
		    push @codons, $c1.$c2.$c3;
		  }
	      }
	  }
	foreach my $c1 (@codons)
	  {
	    foreach my $c2 (@codons)
	      {
		$data->{$c1}{$c2}="-50" unless defined $data->{$c1}{$c2};
	      }
	  }
	$html .= "<tr><th><th>".join ("<th>", map{$_."(".$code->{$_}.")"}sort { sort_nt1(substr($a, 0, 1)) <=> sort_nt1(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt3(substr($a,2,1)) <=> sort_nt3(substr($b,2,1)) }@codons);
	$html .= "<th>Total";
	$html .= "<tr>";
	foreach my $aa1 (sort { sort_nt1(substr($a, 0, 1)) <=> sort_nt1(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt3(substr($a,2,1)) <=> sort_nt3(substr($b,2,1)) } @codons)
	  {
	    $html .= "<th>$aa1"."(".$code->{$aa1}.")";
	    my %vals;
	    map {$vals{$_}++}  map {$data->{$aa1}{$_}} keys %{$data->{$aa1}};
	    my ($max) = sort {$b<=>$a} keys %vals;
	    my ($min1, $min2) = sort {$a<=>$b} keys %vals;#	print STDERR join ("\t",  sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %aa_sort),"\n";
	    $min2 = $min1 unless $min2;
	    my $total =0;
	    foreach my $aa2 (sort { sort_nt1(substr($a, 0, 1)) <=> sort_nt1(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt3(substr($a,2,1)) <=> sort_nt3(substr($b,2,1)) } @codons)#	    foreach my $aa2 (sort keys %$data)
	      {	
		my $val = $data->{$aa1}{$aa2};
		$total+=$val;
		my $color;
		if ($val == $min1)
		  {
		    $html .= "<td style=\"background-color: rgb(125,125,125)\">".$val." ".$code->{$aa1}."-".$code->{$aa2};
		  }
		elsif ($val >0)
		  {
		    $color = color_by_usage($max, $val);
		    $html .= "<td style=\"background-color: rgb($color,255,$color)\">".$val." ".$code->{$aa1}."-".$code->{$aa2};
		  }
		else
		  {
		    $color = color_by_usage(abs($min2), $val*-1);
		    $html .= "<td style=\"background-color: rgb(255,$color,$color)\">".$val." ".$code->{$aa1}."-".$code->{$aa2};
		  }
	      }
	    $html .= "<td>$total<tr>";
	  }
      }
    else #dna
      {
	my @dna = qw(A T C G);
	$html .= "<tr><th><th>".join ("<th>", @dna);
	$html .= "<th>Total";
	$html .= "<tr>";
	foreach my $c1 (@dna)
	  {
	    $html .= "<th>$c1";
	    my ($max) = sort {$b<=>$a} map {$data->{$c1}{$_}} keys %{$data->{$c1}};
	    my ($min)= sort {$a<=>$b} map {$data->{$c1}{$_}} keys %{$data->{$c1}}; 
	    my $total =0;
	    foreach my $c2 (@dna)
	      {
		my $val = $data->{$c1}{$c2};
		$total+=$val;
		my $color;
		if ($val >0)
		  {
		    $color = color_by_usage($max, $val);
		    $html .= "<td style=\"background-color: rgb($color,255,$color)\">".$val;
		  }
		else
		  {
		    $color = color_by_usage(abs($min), $val*-1);
		    $html .= "<td style=\"background-color: rgb(255,$color,$color)\">".$val;
		  }
	      }
	    $html .= "<td>$total<tr>";
	  }
      }
    $html .= "</table>";
    return $html;
  }

sub color_by_usage
      {
	my ($max,$value, $opt) = @_;
	$opt = 255 unless $opt;
	return $opt unless $max;
	my $g = $opt*(($max - $value) / $max);
	return int($g + .5);
      }

sub sort_nt1
  {
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "G")
      {
	$val = 1;
      }
    elsif ($chr eq "A")
      {
	$val = 2;
      }
    elsif ($chr eq "T")
      {
	$val = 3;
      }
    return $val;
  }

sub sort_nt2
  {
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "G")
      {
	$val = 1;
      }
    elsif ($chr eq "A")
      {
	$val = 2;
      }
    elsif ($chr eq "T")
      {
	$val = 3;
      }
    return $val;
  }

sub sort_nt3
  {
    my $chr = uc(shift);

    $chr = substr($chr, -1,1) if length($chr)>1;
    my $val = 0;
    if ($chr eq "G")
      {
	$val = 1;
      }
    elsif ($chr eq "T")
      {
	$val = 2;
      }
    elsif ($chr eq "C")
      {
	$val = 3;
      }
    return $val;
  }
