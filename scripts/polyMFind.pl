#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long;
use CoGeX;
use POSIX;
use LWP::Simple;
use CoGe::Accessory::Utils qw( commify );

use vars qw($DEBUG $file $coge $INSTRAIN $homo_threshold);

GetOptions(
	   "file|f=s"=>\$file,
	   "debug"=>\$DEBUG,
	   "instrain=s"=>\$INSTRAIN, #name of specific strain for which polymorphisms are tracked
	   "homopolyer_theshold|ht=i"=>\$homo_threshold, #the minimum lenght of a homopolymer before it is scored.
	   );
$homo_threshold = 5 unless defined $homo_threshold;

###COMMENTING OUT AND REPLACING###
#dsgids == sequence_order in alignment file => coge database ids
#my @dsgids = (
#	      {id=>7111,		  name=>'NCM4139'},
#	      {id=>7112,		  name=>'NCM4287'},
#	      {id=>7113,		  name=>'NCM4299'},
#	      {id=>7114, 		 name=>'NCM4300'},
#	      {id=>7115,		 name=>'NCM4370'},
#	      {id=>7116,		 name=>'NCM4384'},
#	      {id=>7117,		  name=>'NCM4401'},
#	      {id=>7118,	          name=>'NCM4781'},
#	      {id=>4242,	          name=>'MG1655'}, #MG1655 REFERENCE STRAIN
#	     );

unless (-r $file)
  {
    print qq{
Welcome to $0!
  };
    exit;
  }
my $coge = CoGeX->dbconnect();
###COMMENTING OUT AND REPLACING###
#foreach my $val (@dsgids)
#  {
#    $val->{dsg} = $coge->resultset('DatasetGroup')->find($val->{id}); #reference strain dsg
#    my ($ds) = $val->{dsg}->datasets;
#    $val->{ds} = $ds;
#  }

my %SYNMORPHS; #global var, I'm lazy
my %ROW_TYPE_COUNT; #global var, I'm lazy
my %HOMOPOLYMERS; #global var, I'm lazy
my %FALSE_POSITIVE_SCORES;
my %NAMES_COUNT; #used to track whether a gene has multiple polymorphisms in it.  This is indiciative of sequencing or assembly problems, and the polymorphisms may be false positives.

my ($dsgs, $output) = process_file($file);
#sorting the dsgs makes life a lot easier

my @dsgs = map {$dsgs->{$_}} sort { $a <=> $b } keys %$dsgs;
print_header();
print "<table id=feat_table>";
print "<THEAD>";
print "<tr bgcolor=#DDFFDD><TH>";
#strain lables
my @strain_lables;
foreach my $item ( @dsgs )
  {
    my $name = $item->{name};
    push @strain_lables, join ("<br>", split//, $name);
  }

print join ("<TH>", qw(POSITION POSITION<br>SORT GENE FALSE<br>POSITIVE<br>SCORE TYPE SUBTYPE SIZE SIZE<br>SORT POLYMORPH<br>SUMMARY), @strain_lables, qw( 10_NT_ALIGN ANNOTATION (STRAND)CODON GEVO_LINK COGE_ALIGN_LINK HOMOPOLYMER<br>FALSE<br>POSITIVE<br>SCORE MULTI-HIT<br>FALSE<br>POSITIVE<br>SCORE )),"\n";
print "</THEAD>";
print "<TBODY id=feat_table_body>\n";

###COMMENTING OUT AND REPLACING###
#process_file($file);
print $output;

print "</TBODY>\n";
print "</table>\n";

print "<br><h3>Polymorphism Summary</h3>";
print "<table id=mut_types><THEAD><TR bgcolor=#DDFFDD><TH>";
print join ("<TH>", qw(TYPE COUNT)),"</THEAD>\n";
print "<TBODY id=mut_types_body>\n";
#need to rewrite using %SYNMORPHS
foreach my $type (sort keys %SYNMORPHS)
  {
    my $sum = 0;
    map {$sum+=$_}values %{$SYNMORPHS{$type}};
    print "<tr><Td>$type<td>".$sum,"\n";
  }
print "</TBODY>";
print "</table>";

print "<br><h3>Polymorphism counts by strain</h3>";
print "<table id=common_muts><THEAD><TR bgcolor=#DDFFDD><TH>Polymorphism Type<TH>";
#print join ("<TH>", qw(COMMON_MUTATIONS COUNT)),"</THEAD>\n";
print join ("<th>", map {$_->{name}} @dsgs),"\n";
print "<TBODY id=common_muts_body>\n";
#print STDERR Dumper \%SYNMORPHS;
#print Dumper \%NAMES_COUNT;
foreach my $mut_type (sort keys %SYNMORPHS)
  {
    print "<tr><td>$mut_type";
    foreach my $strain (map {$_->{name}} @dsgs)
      {
	my $count = $SYNMORPHS{$mut_type}{$strain} ? $SYNMORPHS{$mut_type}{$strain} : 0;
	print "<td>",$count;
      }
  }
print "</TBODY>";
print "</table>";

print "<br><h3>Polymorphism Summary Row Count</h3>";
print "This is a count of the rows in the table with a type of polymorphism<br>";
print "<table id=row_types><THEAD><TR bgcolor=#DDFFDD><TH>";
print join ("<TH>", qw(TYPE COUNT)),"</THEAD>\n";
print "<TBODY id=row_types_body>\n";

foreach my $type (sort keys %ROW_TYPE_COUNT)
  {
    my $val = $ROW_TYPE_COUNT{$type};
    print "<tr><Td>$type<td>".$val,"\n";
  }
print "</TBODY>";
print "</table>";

print "<br><h3>Putative Homopolymer Sequencing Errors:</h3>";
print "<table id=homopolymer_tally><THEAD><TR bgcolor=#DDFFDD><TH>";
print join ("<TH>", qw(HOMOPOLYMER_LENGTH TOTAL_ERRORS MISSING_NT EXTRA_NT)),"</THEAD>\n";
print "<TBODY id=homopolymer_tally_body>\n";
foreach my $len (sort keys %{$HOMOPOLYMERS{frequency}})
  {
    my $missing = $HOMOPOLYMERS{frequency}{$len}{missing};
    my $extra = $HOMOPOLYMERS{frequency}{$len}{extra};
    my $total = $missing+$extra;
    print "<tr><td>", join ("<td>", $len, $total, $missing, $extra),"\n";
  }
print "</TBODY>";
print "</table>";

print "<table id=homopolymer_strains><THEAD><TR bgcolor=#DDFFDD><TH>";
print join ("<TH>", qw(STRAIN_NAME HOMOPOLYMER_ERROR_COUNT)),"</THEAD>\n";
print "<TBODY id=homopolymer_strains_body>\n";
#foreach my $strain (map {$_->{name}} @dsgids)
foreach my $strain (map {$_->{name}} @dsgs)
  {
    my $error_count = $HOMOPOLYMERS{strain_errors}{$strain};
    print "<tr><td>$strain<td>$error_count\n";
  }
print "</TBODY>";
print "</table>";
print "<br><h3>False Positives Table:</h3>";
print "Contig joins are not counted in this table.<br>";
print "<table id=fp_Scores><THEAD><TR bgcolor=#DDFFDD><TH>";
print join ("<TH>", qw(FALSE_POSITIVE_SCORE COUNT)),"</THEAD>\n";
print "<TBODY id=homopolymer_strains_body>\n";
foreach my $fps (sort {$a<=>$b} keys %FALSE_POSITIVE_SCORES)
  {
    print "<tr><td>$fps<td>".$FALSE_POSITIVE_SCORES{$fps},"\n";
  }
print "</TBODY>";
print "</table>";

print "<br><h3>False positive score calculations</h3>";
print qq{
<p><b>Total false positive score</b> is the sum of "homopolymer false positive score" + "multiple hit false positive score".</p>

<p><b>Homopolmer false positive score</b> = length of homopolymer.  Note: only calculated for homopolymers with a length greater than or equal to $homo_threshold.
<p> 454 sequencing chemistry has a known problem with sequencing homopolymers (http://www.broadinstitute.org/crd/wiki/index.php/Homopolymer).  When a single nucleotide indel polymorphism is detected, this metric helps determine whether it may be due to this type of sequencing problem.  Be aware that if a strand-slippage mutagen is used (e.g. an intercalator such as Ethidium Bromide), such polymorphisms may be real.</p>

<p><b>Multiple hit false positive score</b> = number of unique polymorphisms occurring within the same gene in a strain.  If multiple strains have different number of polymorphisms in the same gene, the greatest value is used.
<p>One source of error in genome assembly is when multiple short repeats are tandemly arranged.  An example of such genome structure are tandemly arrayed tRNAs.  Assembly software often has problems correctly determining the correct read-order for these regions and makes errors in assembly.  Since bacterial gene density is quite high, multiple polymorphisms in a gene is used as a proxy for detecting such assembly errors.</p>
};

print "</html>";

sub process_file
  {
    my $file = shift;

    open (IN, $file)|| die "Can't open $file for reading: $!\n";
    $/ = "\n=\n";
    my $count =1;
    my $dsgs={}; #store hash of dsg information
    my $output;
    while (<IN>)
      {
	$output .= process_alignment_block(block=>$_, dsgs=>$dsgs, count=>$count);
#	last; #for testing purposes
	$count++;
      }
    close IN;
    return $dsgs, $output;
  }

sub process_alignment_block
  {
    my %opts = @_;
    my $block= $opts{block}; #alignment block for processing
    my $dsgs = $opts{dsgs}; #hash ref of dsg information
    my $count = $opts{count}; #count of block that is being processed
    print STDERR "Processing alignment block $count\n";
    $block =~ s/#.*?\n//g;
    my @seqs;
    foreach my $item (split /\n>/, $block)
      {
	$item =~ s/>//g;
	my ($name, $seq) = split /\n/, $item, 2;
	my ($genome, $start, $stop, $dsgid) = $name =~ /(\d+):(\d+)-(\d+).*?\/?(\d+)\.faa/;
	unless ($dsgs->{$genome})
	  {
	    my $dsg = $coge->resultset('DatasetGroup')->find($dsgid);
	    print STDERR "Creating genome object for dataset_group $dsgid: ".$dsg->organism->name,"\n";

	    my ($ds) = $dsg->datasets;
	    $dsgs->{$genome} = {
				id=>$dsgid,
				dsg=>$dsg,
				name=>$dsg->organism->name,
				genome_order=>$genome,
				ds=>$ds,
			       };
	  }
	$seq =~ s/\n|\r//g;
	$seq = uc($seq);
	my @seq = split//, $seq;
	push @seqs, {genome=>($genome),
		     start=>$start,
		     stop=>$stop,
		     seq=>\@seq};
      }
    my $output = find_alignment_differences(data=>\@seqs, dsgs=>$dsgs);
    return $output;
  }

sub find_alignment_differences
  {
    print STDERR "Finding alignment differences\n";
    my %opts = @_;
    my $data = $opts{data};
    my $dsgs = $opts{dsgs};

    my $output;  #place to store the good for output
    my %rw_position; #real world genomic positions
    map { $rw_position{$_->{genome}} = ($_->{start}-1)} @$data;  #initialize real world genomic positions, start at position just before the aligned block
    my $gap_flag =0; #are we in a gap?
    my $n_flag =0; #are we in a joined contig region?
    my @prior_char;
    my $prior_str;

    my @output;  #storage for output;
    my @n_trace; #which sequences have n?
    for (my $i=0; $i < scalar(@{$data->[0]{seq}}); $i++)
      {
	print STDERR "Processing position $i\n" unless $i % 100000;
#	last if $i > 50000;
	my ($type, $subtype, $codon, $size, $info, $gevolink, $alignlink) = qw(- - - - - - - -);
	my ($names, $annos) = qw(- -);
	my $homo_false_pos_score = 0;
	my $multi_false_pos_score = 0;
	my @char = map {$_->{seq}[$i]} @$data;
	my @strain_info;
	#increment real world positions
	my $count = 1;
	foreach my $char (@char)
	  {
	    $rw_position{$count}++ unless $char eq "-";
	    $count++;
	  }
	my $str = join "", @char;
#	print join ("\t", $str, length($str), scalar keys %$dsgs),"\n";
	next unless length $str == scalar keys %$dsgs;

	if ($gap_flag)
	  {
	    my $gevosize = $gap_flag*2;
	    $gevosize = 100 if $gevosize < 100;
	    $gevolink = gen_gevolink(positions=>\%rw_position, data=>$data, size=>$gevosize, dsgs=>$dsgs);
	    my $gap_count = $prior_str =~ tr/-/-/;
	    if ($gap_count == (length($str)-1))
	      {
		$type= "insertion";
	      }
	    elsif ($gap_count == 1)
	      {
		$type="deletion";
	      }
	    else
	      {
		$type = "indel";
	      }
	    $size = $gap_flag;
	    @strain_info = @prior_char;
	    if ($str =~/-/)
	      {
		if ($str =~ /n/i)
		  {
		    trace_n($str, \@n_trace);
		    $n_flag++;
		  }
		else
		  {
		    $gap_flag++;
		    @prior_char = @char;
		    $prior_str = $str;
		    next;
		  }
	      }
	    else
	      {
		($codon, $names, $annos, $subtype) = get_position_info(pos=>\%rw_position, range=>$size, dsgs=>$dsgs);
#		$names_count{$names}++ if $names;#
		$subtype = "Frameshift" if $names;
		if ($annos eq "Intergenic")
		  {
		    $alignlink = "-";
		  }
		elsif ($size > 100)
		  {
		    $subtype = "Large deletion/insertion";
		  }
		else
		  {
		    $alignlink = gen_alignlink(pos=>\%rw_position, dsgs=>$dsgs);
		  }
		$info = $prior_str;
		$gap_flag=0;
	      }
	  }
	if ($n_flag)
	  {
	    my $gevosize = $n_flag*2;
	    $gevosize = 100 if $gevosize < 100;
	    $gevolink = gen_gevolink(positions=>\%rw_position, data=>$data, size=>$gevosize, dsgs=>$dsgs);
	    $type="contig join";
	    $#strain_info = 0;
	    $#strain_info = $#char;
	    if ($str=~/n/i || $str=~/-/i)
	      {
		trace_n($str, \@n_trace);
		$n_flag++;
		$gap_flag=0;
		@prior_char = @char;
		$prior_str = $str;
		next;
	      }
	    else
	      {
		$n_flag=0;
	      }
	  }
	my $seqview_links = gen_SeqViewlink(positions=>\%rw_position, data=>$data, char=>\@strain_info, adjust=>-1, dsgs=>$dsgs);

	if ($gevolink ne "-")
	  {
	    my $align_length = 12;
	    $align_length = 300 if $type =~ /contig/;
	    my ($align_str, $tmp_score) = gen_alignment(data=>$data, pos=>$i-1, length => $align_length, dsgs=>$dsgs);
	    $homo_false_pos_score = $tmp_score if $tmp_score;
	    #need to modify false_pos_score and @$seqview_links if we are in a contig join gap
	    if ($type eq "contig join") #need different scoring scheme
	      {
		$homo_false_pos_score=0;
		$multi_false_pos_score=0;

#		$false_pos_score =0;
		my $seqs_with_gaps;
		my $count =0;
		for (my $i=0; $i<scalar @char; $i++)
		  {
		    my $item = $n_trace[$i];
		    $item = 0 unless $item;
		    if ($item>0)
		      {$seqview_links->[$i] = "<span align='center' class='seqview'>N</span>";$count++;}
		    else {$seqview_links->[$i]= " ";}
		    $seqs_with_gaps++ if $item;
		  }
		$subtype = "Contig join in all" if $count == scalar @char;
#		$false_pos_score = 10 if ($seqs_with_gaps == scalar @n_trace || $seqs_with_gaps == 1);
	      }

#	    push @output, [commify($rw_position{$data->[-1]{genome}}), $rw_position{$data->[-1]{genome}}, $names, $homo_false_pos_score, $multi_false_pos_score, $type, $subtype, commify($size), $size, "<span style=\"font-family: 'Courier New', monospace;\">$info</span>", @$seqview_links, $align_str, $annos, $codon, $gevolink, $alignlink];

	    my $in_strain = process_synmorphs(data=>$info, type=>$type, subtype=>$subtype, name=>$names, n_trace=>\@n_trace, dsgs=>$dsgs);
	    undef @n_trace;
	    next if ($INSTRAIN && !$in_strain);  #specific strain specified, no match to it
	    push @output, [commify($rw_position{$data->[-1]{genome}}), $rw_position{$data->[-1]{genome}}, $names, 0, $type, $subtype, commify($size), $size, "<span style=\"font-family: 'Courier New', monospace;\">$info</span>", @$seqview_links, $align_str, $annos, $codon, $gevolink, $alignlink, $homo_false_pos_score, $multi_false_pos_score];

	  }

	#are all the characters the same
	my %seen;
	map {$seen{$_}++} @char;
	next if keys %seen == 1 && $str !~ /N/; #skip if the sequence is all the same, but keep contig breaks
	@prior_char = @char;
	$prior_str = $str;
	if($str=~ /N/)
	  {
#	    print $str,"\n";
	    trace_n($str, \@n_trace);
	    $n_flag++; #set n flag
	    next;
	  }
	elsif ($str =~ /-/)
	  {
	    $gap_flag++; #set gap flag
	    next;
	  }
	else
	  {
	    $type = "SNP";
	    $size = 1;
	    ($codon, $names, $annos, $subtype) = get_position_info(pos=>\%rw_position, dsgs=>$dsgs);
#	    $names_count{$names}++if $names;
	    $info = join ("",@char);
	    $gevolink = gen_gevolink(positions=>\%rw_position, data=>$data, size=>100, dsgs=>$dsgs);
	    $alignlink = $annos eq "INTERGENIC"?"-": gen_alignlink(\%rw_position);

	  }
	my $seqview_links = gen_SeqViewlink(positions=>\%rw_position, data=>$data, char=>\@char);
	if ($gevolink ne "-")
	  {
	    my $align_length = 12;
	    $align_length = 300 if $type =~ /contig/;
	    my ($align_str, $tmp_score) = gen_alignment(data=>$data, pos=>$i, length=>$align_length);
	    $homo_false_pos_score = $tmp_score if $tmp_score;
#	    push @output, [commify($rw_position{$data->[-1]{genome}}), $rw_position{$data->[-1]{genome}}, $names, $homo_false_pos_score, $multi_false_pos_score, $type, $subtype, commify($size), $size, "<span style=\"font-family: 'Courier New', monospace;\">$info</span>", @$seqview_links, $align_str, $annos, $codon, $gevolink, $alignlink];
	    my $in_strain = process_synmorphs(data=>$info, type=>$type, subtype=>$subtype, name=>$names, n_trace=>\@n_trace, dsgs=>$dsgs);
	    @n_trace=undef;
	    next if ($INSTRAIN && !$in_strain);  #specific strain specified, no match to it
	    push @output, [commify($rw_position{$data->[-1]{genome}}), $rw_position{$data->[-1]{genome}}, $names, 0, $type, $subtype, commify($size), $size, "<span style=\"font-family: 'Courier New', monospace;\">$info</span>", @$seqview_links, $align_str, $annos, $codon, $gevolink, $alignlink, $homo_false_pos_score, $multi_false_pos_score];

	  }
      }
    foreach my $item (@output)
      {
#	print "here\n";
	#increase false positive score for genes that have been hit multiple times
#	$item->[4]+=$names_count{$item->[2]} if $names_count{$item->[2]} > 1;
	#find the strain with the most number of mutations in this gene (if applicable)
#	print "!!\/".$item->[2],"\n!!\n";
	if (my $name_hash = $NAMES_COUNT{$item->[2]})
	  {
#	    print $item->[2],"\n",Dumper $name_hash;
	    my ($max) =  sort {$b <=> $a} values %$name_hash;
#	    print "##$max##\n";
	    $item->[-1] =$max if $max> 1;
	  }
	$item->[3] = $item->[-1]+$item->[-2];

	$FALSE_POSITIVE_SCORES{$item->[3]}++ unless $item->[4] =~ /contig/;
	$output .=  "<tr>";
	foreach my $thing (@$item)
	  {
	    $output .=  $thing =~ /seqview/i ? "<td align=center>" : "<td>";
	    $output .=  $thing;
	  }
	$output .=  "\n";

	$ROW_TYPE_COUNT{$item->[4]}++;
	my $st = $item->[5];
	$st = "Nonsynonymous" if $st =~ /nonsynon/i;
	$ROW_TYPE_COUNT{"subtype: ".$st}++ if $item->[5] && $item->[5] ne "-" && $item->[5] ne " ";
#	print "<tr><td>", join ("<td>", @$item),"\n";
      }
    return $output;
  }

sub trace_n
    {
      my $seq = shift;
      my $ncount = shift;
      my $i =0;
      foreach my $chr (split//, $seq)
	{
	  $ncount->[$i] = 0 unless $ncount->[$i];
	  $ncount->[$i]++ if $chr eq "N";
	  $i++;
	}
    }

sub gen_gevolink
  {
    my %opts = @_;
    my $positions = $opts{positions};
    my $data = $opts{data};
    my $size = $opts{size};
    my $dsgs = $opts{dsgs};
    $size = 500 unless $size;
    my $gevolink= "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?";
    my $count =1;
    foreach my $item (@$data)
      {
	$gevolink.= "dsgid".$count."=".$dsgs->{$item->{genome}}{id}.";x".$count."=".($positions->{$item->{genome}}).";chr".$count."=1;";
	$count++;
      }
    $gevolink.= "color_hsp=1;num_seqs=".scalar(@$data);
    $gevolink.= ";autogo=1;apply_all=";
#    $gevolink = get_tiny_url(url=>$gevolink);
    return "<span class = link onclick=window.open('".$gevolink.$size."') >GEvo Link</span>"."<br><span class = link onclick=window.open('".$gevolink."20000')>GEvo_Link_20k</span>";
#    return "=hyperlink(\"$gevolink\")";
  }

sub gen_alignlink
  {
    my %opts = @_;
    my $positions = $opts{pos};
    my $dsgs = $opts{dsgs};
    my $link = "http://synteny.cnr.berkeley.edu/CoGe/CoGeAlign.pl?fid=";
    my @fids;
    foreach my $key (sort keys %$positions)
      {
	foreach my $feat ($coge->get_features_in_region(dataset_id=>$dsgs->{$key}{ds}->id, chr=>1, start=>$positions->{$key}))
	  {
	    push @fids, $feat->id if $feat->type->name eq "CDS";
	  }
      }
    $link .= join (";fid=", @fids);
    $link .= ";align=codon;autogo=1";
    if (@fids)
      {
#	$link = get_tiny_url(url=>$link);
#	return "=hyperlink(\"$link\")";
	return "<a href=\"".$link."\" target=_new>CoGeAlign Link</a>";
      }
    else
      {
	return;
      }
  }

sub gen_SeqViewlink
  {
    my %opts = @_;
    my $positions = $opts{positions};
    my $data = $opts{data};
    my $chr = $opts{char};
    my $adjust=>$opts{adjust}; #amount to offest positions
    my $dsgs = $opts{dsgs};
    my $adjust = 0 unless $adjust;
    my $link= "http://synteny.cnr.berkeley.edu/CoGe/SeqView.pl?";
    my $count =1;
    my @links;
    foreach my $item (@$data)
      {
	push @links, "<span align='center' class = link onclick=window.open('".$link."dsgid=".$dsgs->{$item->{genome}}{id}.";start=".($positions->{$item->{genome}}+$adjust-10).";stop=".($positions->{$item->{genome}}+$adjust+10).";chr=1;') >".$chr->[$count-1]."</span>";
	$count++;
      }
    return \@links;
  }

sub gen_alignment
    {
      my %opts = @_;
      my $data = $opts{data};
      my $pos = $opts{pos};
      my $length = $opts{length};
      my $dsgs = $opts{dsgs};
      $length = 12 unless $length;
      my $score = 0; #for tracking runs of a nucleotide.  Runs of a single nucleotide are a source of sequencing errors.  Want to increse this score as runs containing a polymorphism get longer, and decrease this score if the same polymorphism is seen in more than one organism
      my $aln = "&nbsp;"x($length)."!<br>";
      my @alns = map {join ("", @{$_->{seq}}[($pos-$length)..($pos+$length)])} @$data;
#      print " "x$length."!\n";
      my %types;
      #454 generates sequencing errors when processing homopolymers.  This results in an extra or missing polymer base being called.  Let's track this stuff and get a tabulation of such sequencing errors
      my @homopolymers; #store length of homopolymer runs
      my $skip_homopolymer_check=0;  #flag for skipping homopolymer check.  Needed when "Ns" are present
      foreach my $seq (@alns)
	{
	  my $chr = substr($seq, $length, 1);
	  next unless $chr;
	  $types{$chr}++;
	  my $pre = substr($seq, 0, $length);
	  my $post = substr($seq, $length+1, $length);
	  my ($pre_match) = $pre =~ /($chr+)$/;
	  my ($post_match) = $post =~ /^($chr+)/;
	  my $run = $pre_match.$chr.$post_match;
	  my $len = $chr ne "-" && $chr ne "N" ? length($run) : 0; #need to make sure that contig joins, unsequenced regions, large insertions and large deletions are skipped
	  push @homopolymers, $len;
	  $skip_homopolymer_check = 1 if $chr eq "N"; #in a contig join or unsequenced region.  Can't do this check
	  $skip_homopolymer_check =1 if $chr eq "-" && length ($run) > 1; #in an insertion or deletion. Can't do this check

	  next if $chr eq "-"; #this sequence has an insert polymorphism
	  next unless $len >= $homo_threshold; #only sore if homopolymer is $homo_threshold nt or longer.  $homo_threshold is a global var
	  $score =$len unless $score && $score > $len;
	}
      $skip_homopolymer_check =1 unless $types{"-"}; #need at least one gap, otherwise just a SNP
      $score=0 unless $types{"-"};  #not a homopolymer sequencing error, it is a SNP
      #track homopolymers.  Common source of sequence errors with 454's technology
      my ($max_hp_length) = sort {$b <=> $a} @homopolymers; #get longest homopolymer
      if (!$skip_homopolymer_check && $max_hp_length > 1) #we have a homopolymer.  Otherwise, just a SNP
	{ #need to determine of the SNP is an extra nt or missing a nt.
	  #debugging
#	  for (my $i=0; $i<@homopolymers; $i++)
#	    {
#	      print $homopolymers[$i],"\t", $alns[$i],"\n";
#	    }
	  my %lens;
	  map {$lens{$_}++} @homopolymers;
	  my ($correct_length) = sort {$lens{$b} <=> $lens{$a}} keys %lens; #assume correct length is seen most frequently
	  my $real_correct_length = $correct_length ? $correct_length : $max_hp_length-1; # $correct_length will == 0 if the other positions are gaps or Ns

	  $HOMOPOLYMERS{frequency}{$real_correct_length}{'missing'} = 0 unless defined $HOMOPOLYMERS{frequency}{$real_correct_length}{'missing'};
	  $HOMOPOLYMERS{frequency}{$real_correct_length}{'extra'} = 0 unless defined $HOMOPOLYMERS{frequency}{$real_correct_length}{'extra'};

	  #tally which strains had the putative sequencing error
	  for (my $i=0; $i<@homopolymers; $i++)
	    {
	      my $len = $homopolymers[$i];
	      if ($correct_length != $len)
		{
		  #tally up for the stain
		  $HOMOPOLYMERS{strain_errors}{$dsgs->{$i+1}{name}}++;
		  #tally if this is an extra or missing nt
		  if ($max_hp_length == $correct_length) #we have an underreporting
		    {
		      $HOMOPOLYMERS{frequency}{$real_correct_length}{'missing'}++;
		    }
		  else #we have an extra nt called
		    {
		      $HOMOPOLYMERS{frequency}{$real_correct_length}{'extra'}++;
		    }
		}
	    }
	}
      #modify score if polymorphism occurrence is in more than one sequence
#      my ($min_type) = sort {$a<=>$b}values %types; #which is the least represented character?
#      if ($min_type > 1 && $min_type < scalar (@alns)-1)
#	{
#	  my $min = $min_type < scalar(@alns)/2 ? $min_type : scalar(@alns)-$min_type;
#	  $score /=$min; #
#	}
      $aln .= join ("<br>", @alns);
      $aln = "<span class='a$pos' style=\"display:none; font-family: 'Courier New', monospace;\"><pre>$aln</pre></span>";
#      $aln .= "<pre class='a$pos' style=\"display:none; font-family: 'Courier New', monospace;\">$aln</pre>";
      $aln .= qq{<span class="link a$pos" onclick="\$('.a$pos').toggle()">show</span>};
      $aln .= qq{<br><span class="link a$pos" style='display:none' onclick="\$('.a$pos').toggle()">hide</span>};
      $score = sprintf("%.2f", $score) if $score =~ /\./;
      return ($aln, $score);
    }

sub get_position_info
  {
    my %opts = @_;
    my $pos = $opts{pos};
    my $range = $opts{range};
    my $dsgs = $opts{dsgs};
    $range = 0 unless $range;
    my @feats;
    my %codons;
    my %amino_acids;
    my $codon_position;
    my $strand;
    my $id; #need a unique id number for html stuff
    my $genetic_code;
    foreach my $i (keys %$pos)
      {
	my $start = $pos->{$i};
	$id = $start unless $id;
	my $ds = $dsgs->{$i}{ds};
	unless ($ds)
	  {
	    print STDERR "Unable to get a $ds for genome $i\n";
	    exit;
	  }
	foreach my $feat ($coge->get_features_in_region(dataset_id=>$ds->id, chr=>1, start=>$start-$range, stop=>$start))
	  {
	    push @feats, $feat unless $feat->type->name eq "chromosome";
	    if ($feat->type->name eq "CDS")
	      {
		($genetic_code) = $feat->genetic_code unless $genetic_code;
		$strand = $feat->strand unless $strand;
 		my $genomic_seq = $feat->genomic_sequence;
 		#need to determine the codon that has been changed
		$codon_position = $feat->strand =~ /-/ ? $feat->stop-$start+1 : $start-$feat->start+1;
 		my $codon_start = ceil($codon_position/3)*3;
 		$codon_start-=2;
 		my $codons = substr($genomic_seq, $codon_start-1, 15); #grab 5 codons, with the one of interest on the left
		$codons =~ s/(\S{3})/$1&nbsp;/g; #add space between codons.
		my ($codon) = $codons =~ /(^.{3})/; #grab the first codon;
		$amino_acids{$genetic_code->{uc($codon)}}=1;
		push @{$codons{$i}}, "(".$feat->strand."):".$codons;

	      }
	  }
      }
    my @genes;
    my $codon;
    my $type;
    if (@feats && @feats == 1)
      {
	@genes = ($feats[0]);
	$type = $genes[0]->type->name;

      }
    elsif (@feats)
      {
	foreach my $feat (@feats)
	  {
	    if ($feat->type->name eq "CDS")
	      {
		push @genes, $feat;
	      }
	  }
	if (!@genes)
	  {
	    @genes = @feats;
	    $type = $feats[-1]->type->name;
	  }
      }
    my $annos;
    my $names;
    my @names; #store non-Estrains names;
    if (@genes)
      {
	my %names;
	foreach my $item (@genes)
	  {
	    map {$names{$_}=1} $item->names;
	  }
	foreach my $item (keys %names)
	  {
	    push @names, $item if $item !~ /E\d{4}_/;
	  }
	@names = keys %names unless @names;
	$names = join (", ",  map{qq{<span class=link onclick=window.open('/CoGe/FeatView.pl?accn=$_')>$_</span>}}sort @names);
	my %annos;
	foreach my $item (@genes)#$genes[0]->annotations)
	  {
	    foreach my $anno ($item->annotations)
	      {
		$annos{$anno->annotation}=1 if $anno->type->name =~ /note/i || $anno->type->name =~ /annotation/i|| $anno->type->name =~ /product/i;
	      }
	  }
	$annos = join (", ", sort keys %annos);
      }
    else
      {
	$annos = "INTERGENIC";
      }
    if (%codons)
      {
	$codon = "<span class=\"c$id\" style=\"font-family: 'Courier New', monospace; display:none\">";
	my $rel_codon_position = $codon_position%3;
	$rel_codon_position = 3 unless $rel_codon_position;
	$codon .= "&nbsp;"x($rel_codon_position+length($strand)+2)."!<br>";
	$codon .= join ("<br>", map{@{$codons{$_}}} sort keys %codons)."</span>";
	$codon .= qq{<span class="link c$id" onclick="\$('.c$id').toggle()">show</span>};
	$codon .= qq{<br><span class="link c$id" style='display:none' onclick="\$('.c$id').toggle()">hide</span>};
      }
    if (%amino_acids)
      {
	if (keys %amino_acids == 1)
	  {
	    $type = "Synonymous";
	  }
	else
	  {
	    $type = "Nonsynonymous<br>Amino acids: ";
	    $type .= join (" ", sort keys %amino_acids);
	  }
      }

    return $codon, $names, $annos, $type;
  }

sub get_tiny_url
  {
    my %opts = @_;
    my $url = $opts{url};
    my $html;
    my $tiny = get("http://tinyurl.com/create.php?url=$url");
    unless ($tiny)
      {
        return "Unable to produce tiny url from <a href=tinyurl.com>tinyurl.com</a>";
      }
    ($tiny) = $tiny =~ /<b>(http:\/\/tinyurl.com\/\w+)<\/b>/;

    return $tiny;
  }

sub process_synmorphs
  {
    my %opts = @_;
    my $data = $opts{data}; #string of text corresonding to characters at position
    my $type = $opts{type}; #type of mutation;
    my $subtype = $opts{subtype};
    my $name = $opts{name};
    my $n_trace = $opts{n_trace}; #for tracking Ns in contig joins
    my $dsgs = $opts{dsgs};
    $subtype = 0 if $subtype eq "-";
    my @chars = split //,$data; #break up characters;
    my $i = 1;
    my %counts;
    if ($data eq "-" && $n_trace && scalar @$n_trace) #dealing with a contig join
      {
	foreach my $item (@$n_trace)
	  {
	    push @{$counts{"N"}}, $i if $item;
	    $i++;
	  }
#	$subtype = "Contig join in all" if scalar @{$counts{"N"}} eq scalar @dsgids;
      }
    else
      {
	foreach my $char (@chars)
	  {
	    push @{$counts{$char}}, $i;
	    $i++;
	  }
      }
    #assuming we only have two states per position.
    my ($k) =  sort {scalar@{$counts{$a}} <=> scalar @{$counts{$b}}} keys %counts;
    if ($INSTRAIN) #strain specific run
      {
	next unless scalar @{$counts{$k}} == 1; #only going to process if a single strain is present
      }
    my $in_strain=0;
    foreach my $strain (map {$dsgs->{$_}{id}} @{$counts{$k}})
      {
	if ($INSTRAIN) #is this a strain specific run?
	  {
	    next unless $strain eq $INSTRAIN;
	    $in_strain=1;
	  }
	$SYNMORPHS{$type}{$strain}++;
	my $tmp_type = $subtype;
	$tmp_type = "Nonsynonymous" if $subtype =~ /nonsynon/i;
	$SYNMORPHS{"subtype: ".$tmp_type}{$strain}++ if $tmp_type;
	$NAMES_COUNT{$name}{$strain}++ if $name && $name ne "-"; #increment counter of mutations in this gene for this strain
      }
    return $in_strain;
#    $SYNMORPHS{$str}++;
  }

sub print_header
    {
print qq{<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery-1.3.2.min.js"></SCRIPT>
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery.tablesorter.2.0.3.js"></SCRIPT>
<SCRIPT language="JavaScript" type="text/javascript" src="/CoGe/js/jquery.tablesorter.pager.js"></SCRIPT>
<link rel="stylesheet" type="text/css" href="/CoGe/css/jquery-ui-1.7.2.custom.css" />
<link rel="stylesheet" type="text/css" href="/CoGe/css/coge.css" />

<Script language="JavaScript">
\$(document).ready(function(){
                \$('#feat_table').tablesorter({
                cssAsc: 'headerSortUp',         // Class name for ascending sorting action to header
                cssDesc: 'headerSortDown',      // Class name for descending sorting action to header
                widgets: ['zebra'],
                headers: {
                          1: {sorter:'digit'},
                          3: {sorter:'digit'},
                          7: {sorter:'digit'}
                         }
                });
                \$('#mut_types').tablesorter({
                cssAsc: 'headerSortUp',         // Class name for ascending sorting action to header
                cssDesc: 'headerSortDown',      // Class name for descending sorting action to header
                widgets: ['zebra']
                });
                \$('#row_types').tablesorter({
                cssAsc: 'headerSortUp',         // Class name for ascending sorting action to header
                cssDesc: 'headerSortDown',      // Class name for descending sorting action to header
                widgets: ['zebra']
                });
                \$('#common_muts').tablesorter({
                cssAsc: 'headerSortUp',         // Class name for ascending sorting action to header
                cssDesc: 'headerSortDown',      // Class name for descending sorting action to header
                widgets: ['zebra']
                });
                \$('#pairs').tablesorter({
                cssAsc: 'headerSortUp',         // Class name for ascending sorting action to header
                cssDesc: 'headerSortDown',      // Class name for descending sorting action to header
                widgets: ['zebra']
                });
 });
</script>
};

    }
