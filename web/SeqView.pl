#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use CoGe::Genome;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;

$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $DB );

$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;
$FORM = new CGI;

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seq=>\&get_seq,
		       get_extend_box=>\&get_extend_box,
		       gen_title=>\&gen_title,
		       new_foot=>\&new_foot,
		       find_feats=>\&find_feats,
			);
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header;
#print gen_html();
sub gen_html
  {
    my $form = $FORM;
    #print "When joining multiple locations for feature, (a) remove 'join' from location in fasta header, and (b) probably indicate where discontinuous sequences are being joined. Also, need to open Get Features in new Window, and place highlighting buttons under each annotation description for the features listed. Finally, need to have GeLo send SeqView start/stop values correctly.";
    my $feat_id = $form->param('featid');
    my ($body) = gen_body();
    my $rc = $form->param('rc');
    my $pro;
    my $foot = gen_foot();
    my ($title) = gen_title(protein=>$pro, rc=>$rc);
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
#    unless ($feat_id) 
#     {
    $template->param(TITLE=>'Sequence Viewer');
    #}
    $template->param(HELP=>'SeqView');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"SeqView-logo.png");
    $template->param(BOX_NAME=>qq{<DIV id="box_name">$title</DIV>});
    $template->param(BODY=>$body);
    $template->param(POSTBOX=>qq{<DIV id = buttons>$foot</DIV});
#    $template->param(CLOSE=>1);
    $template->param(HEAD=>qq{<script src="js/kaj.stable.js"></script>});
    #print STDERR gen_foot()."\n";
    my $html;
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $form = $FORM;
    my $feat_id = $form->param('featid') || 0;
    my $chr = $form->param('chr');
    my $dsid = $form->param('dsid');
    my $feat_name = $form->param('featname');
    my $rc = $form->param('rc');
    my $pro = $form->param('pro');   
    my $upstream = $form->param('upstream');
    my $downstream = $form->param('downstream');
    my $start = $form->param('start');
    my $stop = $form->param('stop');
    my $seq;
    #$seq = find_feats(dsid=>$dsid, chr=>$chr, start=>$start, stop=>$stop);
    unless ($feat_id)
    {
      my $strand = $form->param('strand') || 1;
      $seq = get_seq(pro=>$pro,
		      rc=>$rc,
		      chr=>$chr,
		      dsid=>$dsid,
		      strand=>$strand, 
		      upstream=>$upstream, 
		      downstream=>$downstream,
		      start=>$start,
		      stop=>$stop,
		      );
    }
    else
    {
    my $strand = get_strand($feat_id);
    $strand = check_strand(strand=>$strand, rc=>$rc);
    $seq = get_seq(featid=>$feat_id, 
		      pro=>$pro,
		      rc=>$rc,
		      chr=>$chr,
		      dsid=>$dsid,
		      featname=>$feat_name,
		      strand=>$strand, 
		      upstream=>$upstream, 
		      downstream=>$downstream,
		      );
    }
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    #$template->param(TEXT=>1);
    $template->param(SEQ=>$seq);
    #$seq = qq{<TABLE style="width: 628px; height: 300px; overflow: auto;"><TR><TD align=left valign="top">$seq</TD></TR></TABLE>};
#    return qq{<DIV id="seq" style="width: 623px; height: 300px; overflow: auto;">$seq</div>};
    return qq{<DIV id="seq" style="height: 300px; overflow: auto;">$seq</div>};
#    return qq{<DIV id="seq">$seq</div>};
  }
 
sub check_strand
{
    my %opts = @_;
    my $strand = $opts{'strand'};
    my $rc = $opts{'rc'};
    #print STDERR Dumper \%opts;
    if ($rc==1)
      {
        if ($strand =~ /-/)
          {
            $strand = "1";
          }
        else
          {
            $strand = "-1";
          }
      }
    return $strand;
}

sub get_seq
  {
    my %opts = @_;
    my $add_to_seq = $opts{'add'};
    my $feat_id = $opts{'featid'} || 0;
    my $pro = $opts{'pro'};
    my $rc = $opts{'rc'};
    my $chr = $opts{'chr'};
    my $dsid = $opts{'dsid'};
    my $feat_name = $opts{'featname'};
    my $upstream = $opts{'upstream'};
    my $downstream = $opts{'downstream'};
    my $start = $opts{'start'};
    my $stop = $opts{'stop'};
    my $change_strand = $opts{'changestrand'} || 0;
    if($add_to_seq){
      $start = $upstream;
      $stop = $downstream;
    }
    my $ds = $DB->get_dataset_obj->retrieve($dsid);
    #print $rc;
    my $strand;
    my $seq;
    my $fasta;
    my $fasta_no_html;
   # print STDERR Dumper \%opts;
    unless ($feat_id)
    {
      $strand = $opts{'strand'};
      $fasta = ">".$ds->org->name.", Location: ".$start."-".$stop.", Chromosome: ".$chr.", Strand: ".$strand."\n";
      $fasta_no_html = ">".$ds->org->name.", Location: ".$start."-".$stop.", Chromosome: ".$chr;
      $fasta = qq{<FONT class="main"><i>$fasta</i></FONT>};
#      $columns = 80;
#      $fasta = join ("\n", wrap('','',$fasta));
    }
    else
    {
    $strand = get_strand($feat_id);
    $strand = check_strand(strand=>$strand, rc=>$rc);
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    $fasta = ">".$ds->org->name."(v.".$feat->version.") ".", Name: ".$feat_name.", Type: ".$feat->type->name.", Location: ".$feat->genbank_location_string.", Chromosome: ".$chr.", Strand: ".$strand."\n";
    $fasta = qq{<FONT class="main"><i>$fasta</i></FONT>};
#    $columns = 80;
#    $fasta = join ("\n", wrap('','',$fasta));
    }
    unless ($pro)
    {
      $seq .= get_dna_seq_for_feat (featid=>$feat_id,
      				    dsid=>$dsid, 
      				    rc=>$rc, 
      				    upstream=>$upstream,
      				    downstream=>$downstream, 
      				    start=>$start,
      				    stop=>$stop,
      				    chr=>$chr,
      				    fasta=>$fasta_no_html);
    }
    else
    {
      $seq .= get_prot_seq_for_feat($feat_id);
    }
    #print length($seq);
     my $up;
     my $down;
     my $main;
     unless($pro)
     {
      if($feat_id)
      {
       unless ($rc)
       {
        $up = upstream_color(seq=>$seq, upstream=>$upstream);
        $down = downstream_color(seq=>$seq, downstream=>$downstream);
        $main = main_color(seq=>$seq, upstream=>$upstream, downstream=>$downstream);
       }
       else
       {
        $down = upstream_color(seq=>$seq, upstream=>$upstream);
        $up = downstream_color(seq=>$seq, downstream=>$downstream);
        $main = main_color(seq=>$seq, upstream=>$downstream, downstream=>$upstream);
       }
       $seq = join("", $up, $main, $down);
      }
     }
     else{
      $seq = qq{<FONT class="main">$seq</FONT>};
     }
    unless ($rc==2)
     {$seq = ($fasta. $seq);}
    return "<pre>".$seq."</pre>";
  }
  
sub gen_foot
  {
    my $form = $FORM;
    my $feat_id = $form->param('featid');
    my $chr = $form->param('chr');
    my $dsid = $form->param('dsid');
    my $feat_name = $form->param('featname');
    my $rc = $form->param('rc');
    my $pro = $form->param('pro');
    my $upstream = $form->param('upstream');
    my $downstream = $form->param('downstream');
    my $start = $form->param('start');
    my $stop = $form->param('stop');
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    my $strand;
    unless ($feat_id)
    {$strand = $form->param('strand');}
    else
    {$strand = $feat->strand;}
#   my $button = qq{<input type="submit" value="$name">};
    #my $url = $form->url(-full=>1, -query=>1);
    #$url =~ s/rc=.//g;
    #$url =~ s/pro=.//g;
    #my $url2 = $url;
    #my $url3 = $url;
    #$url .= "&rc=0";
    #$url2 .="&rc=1";
    #$url3 .= "&pro=1";
    my $output = new_foot(featid=>$feat_id,
			  chr=>$chr,
			  dsid=>$dsid,
			  featname=>$feat_name,
			  rc=>$rc,
			  pro=>$pro,
			  upstream=>$upstream,
			  downstream=>$downstream,
			  start=>$start,
			  stop=>$stop,
			  strand=>$strand);
    #$button .= qq{<DIV id="Proseq"><input type="button" value="$name2" onClick="window.location='$url2'"> };
    #return qq{<FORM ACTION="SeqView.pl?featid=$feat_id&dsid=$dsid&chr=$chr&name=$feat_name&rc=$rc">$button</FORM>};
   # return qq{<FORM METHOD = "GET", ACTION="$url">$button</FORM>};
   return qq{<DIV id = buttons>$output</DIV};
  }
    
sub get_dna_seq_for_feat
  {
    my %opts = @_;
    my $featid = $opts{featid};
    my $dsid = $opts{dsid};
    my $rc = $opts{rc};
    my $upstream = $opts{upstream};
    my $downstream = $opts{downstream};
    my $start = $opts{start};
    my $stop = $opts{stop};
    my $chr = $opts{chr};
    my $fasta = $opts{fasta};
    my $seq;
    #print STDERR Dumper \%opts;
#    print STDERR "dsid;$dsid\n";
    unless ($featid)
    {
    $seq = $DB->get_genomic_sequence(start=>$start,
    					stop=>$stop,
    					chr=>$chr,
    					dataset_id=>$dsid);
    }
    else
    {
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    unless (ref($feat) =~ /Feature/i)
      {
        return "Unable to retrieve Feature object for id: $featid";
      }
    $seq = $feat->genomic_sequence(upstream=>$upstream, downstream=>$downstream);
    }
#    print STDERR "Done\n";
    if ($rc==1)
      {$seq = reverse_complement($seq);}
    elsif ($rc==2)
      {$seq = sixframe(seq=>$seq, fasta=>$fasta);}
    $columns = 80;
    unless ($rc==2)
     {$seq = join ("\n", wrap('','',$seq));}
    return $seq;
  }
  
sub get_prot_seq_for_feat
  {
    my $featid = shift;
    my ($feat) = $DB->get_feat_obj->retrieve($featid);
    my ($seq) = $DB->get_protein_sequence_for_feature($feat);
    $seq = "No sequence available" unless $seq;
    #print $seq;
    $columns = 60;
    $seq = join ("\n", wrap('','',$seq));
    return $seq;
  }
  
sub reverse_complement
  {
    my $seq = shift;
    my $rcseq = reverse($seq);
    $rcseq =~ tr/ATCG/TAGC/; 
    return $rcseq;
  }
  
sub upstream_color
    {
      my %opts = @_;
      my $seq = $opts{seq};
      my $upstream = $opts{upstream};
      my $rc = $opts{rc};
      my $up = substr($seq, 0, $upstream);
      unless ($rc)
       {return qq{<FONT class="up">$up</FONT>};}
      else
       {return qq{<FONT class="down">$up</FONT>};}
    }

sub downstream_color
    {
      my %opts = @_;
      my $seq = $opts{seq};
      my $downstream = $opts{downstream};
      my $rc = $opts{rc};
      my $down = substr($seq, ((length $seq)-($downstream)), length $seq);
      unless ($rc)
       {return qq{<FONT class="down">$down</FONT>};}
      else
       {return qq{<FONT class="up">$down</FONT>};}
     }
     
sub main_color
    {
      my %opts = @_;
      my $seq = $opts{seq};
      my $upstream = $opts{upstream};
      my $downstream = $opts{downstream};
      my $main = substr($seq, $upstream, ((length $seq) - ($downstream)));
      return qq{<FONT class="main">$main</FONT>};
    }
    
sub gen_title
    {
      my %opts = @_;
      my $rc = $opts{rc};
      my $pro = $opts{pro};
      my $title;
      unless ($pro)
      {
       if ($rc == 2)
        {$title = "Six Frame Translation";}
       elsif ($rc == 1)
        {$title = "Reverse Complement";}
       else
        {$title = "DNA Sequence";}
      }
      else
      {
       $title = "Protein Sequence";
      }
      return $title;
    }

sub sixframe
	{
	  my %opts = @_;
	  my $seq = $opts{seq};
	  my $fasta = $opts{fasta};
	  my $key;
	  my $sixframe;
     	  my $sequence = $DB->get_feat_obj->frame6_trans(seq=>$seq);
          #print STDERR Dumper ($sequence);
          foreach $key (sort {abs($a) <=> abs($b)} keys %$sequence)
           {
      	     $seq = join ("\n", wrap('','',$sequence->{$key}));
      	     $sixframe = join("\n",$sixframe, qq/$fasta Frame $key\n$seq/);
           }
          return "<pre>".$sixframe."</pre>";
        }
	
sub new_foot
{
    my %opts = @_;
    my $feat_id = $opts{'featid'};
    my $chr = $opts{'chr'};
    my $dsid = $opts{'dsid'};
    my $feat_name = $opts{'featname'};
    my $rc = $opts{'rc'};
    my $pro = $opts{'pro'};
    my $upstream = $opts{'upstream'} || 0;
    my $downstream = $opts{'downstream'} || 0;
    my $start = $opts{'start'};
    my $stop = $opts{'stop'};
    my $strand = $opts{'strand'};
    my $change_strand = $opts{'changestrand'} || 0;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    my $title;
    my $DNAButton;
    my $RCButton;
    my $PROButton;
    my @button_loop;
    
    $strand = check_strand(strand=>$strand, rc=>$rc, changestrand=>$change_strand, pro=>$pro);
    #print STDERR Dumper \%opts;
    unless($feat_id){
                      $title = "RANGE";
                      $DNAButton = {BUTTON_NUM=>1, 
                                    BUTTON_NAME=>'DNA Sequence', 
                                    START=>$start,
    			            PRO=>0,
    				    RC=>0,
    			     	    CHR=>$chr,
    				    DSID=>$dsid,
    				    STOP=>$stop,
   				    STRAND=>1,
   				    RC=>0,
    				    PRO=>0,
   				    };
                       $RCButton = {BUTTON_NUM=>2, 
                            	    BUTTON_NAME=>'Reverse Complement', 
                            	    START=>$start,
   				    PRO=>0,
    				    RC=>1,
    			     	    CHR=>$chr,
    				    DSID=>$dsid,	
    				    STOP=>$stop,
   				    STRAND=>-1,
   				    RC=>1,
   				    PRO=>0};
		      $PROButton = {BUTTON_NUM=>3, 
                            	    BUTTON_NAME=>'Six Frame Translation', 
                            	    START=>$start,
   				    PRO=>0,
    				    RC=>2,
    			     	    CHR=>$chr,
    				    DSID=>$dsid,	
    				    STOP=>$stop,
   				    RC=>2,
   				    PRO=>0};
    }
    else{
    $title = "FEATURE"; 
    $DNAButton = {BUTTON_NUM=>1, 
                        BUTTON_NAME=>'DNA Sequence', 
                        FEATID=>$feat_id,
    			PRO=>0,
    			RC=>0,
    			CHR=>$chr,
    			DSID=>$dsid,
    			FEATNAME=>$feat_name,
   			RC=>0,
    			PRO=>0};
     $RCButton = {BUTTON_NUM=>2, 
                       BUTTON_NAME=>'Reverse Complement', 
                       FEATID=>$feat_id,
   			PRO=>0,
    			RC=>1,
    			CHR=>$chr,
    			DSID=>$dsid,	
    			FEATNAME=>$feat_name,
   			RC=>1,
   			PRO=>0};
     $PROButton = {BUTTON_NUM=>3, 
                        BUTTON_NAME=>'Protein Sequence',
                        FEATID=>$feat_id,
   			PRO=>1,
    			CHR=>$chr,
    			DSID=>$dsid,	
    			FEATNAME=>$feat_name,
   			RC=>0,
   			PRO=>1};
      }
    $template->param(BOTTOM_BUTTONS=>1);
    $template->param($title=>1);
        if ($rc == 1)
	{
           push (@button_loop, ($DNAButton, $PROButton));
	}
    	elsif ($pro || $rc == 2)
    	{
	  push (@button_loop, ($DNAButton, $RCButton));
	}
    	else{
    	  push (@button_loop, ($RCButton, $PROButton));
    	}			 
    	$template->param(BUTTON_LOOP=>\@button_loop);
   
   $template->param(START=>$start);
   $template->param(STOP=>$stop);
   $template->param(CHR=>$chr);
   $template->param(DSID=>$dsid);
   unless($feat_id)
    {
      $template->param(EXTEND=>"Sequence Range");
      $template->param(UPSTREAM=>"START: ");
      $template->param(UPVALUE=>$start);
      $template->param(DOWNSTREAM=>"STOP: ");
      $template->param(DOWNVALUE=>$stop);
    }
    else {
      $template->param(FEATID=>$feat_id);
      $template->param(EXTEND=>"Extend Sequence");
      $template->param(UPSTREAM=>"UPSTREAM: ");
      $template->param(UPVALUE=>$upstream);
      $template->param(DOWNSTREAM=>"DOWNSTREAM: ");
      $template->param(DOWNVALUE=>$downstream);
    }
    $template->param(PRO=>$pro);
    $template->param(RC=>$rc);
    $template->param(CHR=>$chr);
    $template->param(DSID=>$dsid);
    $template->param(FEATNAME=>$feat_name);
    $template->param(STRAND=>$strand);
    return $template->output;
  }
  
sub find_feats
{
	#print STDERR "Here";
	my %opts = @_;
	my $start = $opts{'start'};
	my $stop = $opts{'stop'};
	my $chr = $opts{'chr'};
	my $dsid = $opts{'dsid'};
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
	my $html = `$ENV{PATH}/FeatAnno.pl start=$start stop=$stop ds=$dsid chr=$chr`;
	$html = substr($html, (44), length($html));
        $template->param(FEATUREBOX=>1);
        $template->param(LISTFEATURES=>$html);
        $html = $template->output;
        return qq{<DIV id = "" style = "overflow: auto; height: 300px;">$html</DIV>};

# 	my $feat;
# 	my $anno;
# 	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
# 	my @feats = $DB->get_feat_obj->get_features_in_region(start=>$start, stop=>$stop, chr=>$chr, dataset_id=>$dsid);
# 	unless (@feats)
# 	 {$feat = "There are no features in the selected sequence.";}
# 	else
# 	{
# 	 foreach (@feats)
# 	 {
# 	  $feat .= $_->annotation_pretty_print_html()."\n<BR><HR><BR>\n";
# 	 }
# 	}
# 	#"\n<BR>\n".$template->param(HIGHLIGHT=>1)
# 	#print $pretty;
# 	return $feat;
}

sub get_strand
  {
    my $feat_id = shift;
    my ($feat) = $DB->get_feat_obj->retrieve($feat_id);
    my $strand = $feat->strand;
    return $strand;
  }