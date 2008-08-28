#! /usr/bin/perl -w

use strict;
use CGI;
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use CoGe::Accessory::Web;
use CoGeX;
use CoGeX::Feature;
use CoGeX::Dataset;
use HTML::Template;
use Text::Wrap qw($columns &wrap);
use Data::Dumper;
use POSIX;

$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $TEMPDIR $TEMPURL $FORM $USER $DATE $coge);

$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;

$coge = CoGeX->dbconnect();

my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		       get_seq=>\&get_seq,
		       gen_title=>\&gen_title,
		       find_feats=>\&find_feats,
		       parse_url=>\&parse_url,
		       generate_feat_info=>\&generate_feat_info,
		       generate_gc_info=>\&generate_gc_info,
			);
$pj->js_encode_function('escape');
print $pj->build_html($FORM, \&gen_html);
#print $FORM->header, gen_html();

sub gen_html
  {
    my $html;
    unless ($USER)
      {
		$html = login();
      }
    else
     {
    my $form = $FORM;
    my $featid = $form->param('featid');
    my $rc = $form->param('rc');
    my $pro;
    my ($title) = gen_title(protein=>$pro, rc=>$rc);
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
    $template->param(TITLE=>'Sequence Viewer');
    $template->param(HELP=>'SeqView');
    my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"SeqView-logo.png");
    $template->param(BOX_NAME=>qq{<DIV id="box_name">$title</DIV>});
    $template->param(BODY=>gen_body());
    $template->param(POSTBOX=>gen_foot());
    $template->param(LOGON=>1) unless $USER->user_name eq "public";
    #if($featid)
     #{$template->param(CLOSE=>1);}
    #print STDERR gen_foot()."\n";
    $html .= $template->output;
    }
    return $html;
  }

sub gen_body
  {
    my $form = $FORM;
    my $featid = $form->param('featid') || 0;
    my $chr = $form->param('chr');
    my $dsid = $form->param('dsid');
    my $feat_name = $form->param('featname');
    my $rc = $form->param('rc');
    my $pro = $form->param('pro');   
    my $upstream = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $start = $form->param('start');
    my $stop = $form->param('stop');
    $stop = $start unless $stop;
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    $template->param(RC=>$rc);
    $template->param(JS=>1);
    $template->param(SEQ_BOX=>1);
    if ($featid)
    {
      my ($feat) = $coge->resultset('Feature')->find($featid);
      $dsid = $feat->dataset_id;
      $chr = $feat->chromosome;

      $template->param(FEATID=>$featid);
      $template->param(FEATNAME=>$feat_name);
      $template->param(FEAT_INFO=>qq{<td valign=top><input type=button value="Get Feature Info" onClick="generate_feat_info(['args__$featid'],[display_feat_info])">});
      $template->param(DSID=>$dsid); #to make JS happy
    $template->param(CHR=>$chr); #to make JS happy
    }
    else
    {
        $template->param(DSID=>$dsid);
    	$template->param(CHR=>$chr);
    	$template->param(FEATID=>0); #to make JS happy
    	$template->param(FEATNAME=>'null'); #to make JS happy
        #generate_gc_info(chr=>$chr,stop=>$stop,start=>$start,dsid=>$dsid);
    }
#
    $template->param(GC_INFO=>qq{<td valign=top><input type=button value="Calculate GC Content" onClick="generate_gc_info(['seq_text','args__'+myObj.pro],[display_gc_info],'POST')">});
    my $html = $template->output;
    return $html;
  }
 
sub check_strand
{
    my %opts = @_;
    my $strand = $opts{'strand'} || 1;
    my $rc = $opts{'rc'} || 0;
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
     elsif ($strand =~ /-/)
     {
       $strand =~ s/^\-$/-1/;
     }
     else 
     {
       $strand =~ s/^\+$/1/;
     }
    return $strand;
}

sub get_seq
  {
    my %opts = @_;
    my $add_to_seq = $opts{'add'};
    my $featid = $opts{'featid'} || 0;
    $featid = 0 if $featid eq "undefined"; #javascript funkiness
    my $pro = $opts{'pro'};
    #my $pro = 1;
    my $rc = $opts{'rc'} || 0;
    my $chr = $opts{'chr'};
    my $dsid = $opts{'dsid'};
    my $feat_name = $opts{'featname'};
    my $upstream = $opts{'upstream'};
    my $downstream = $opts{'downstream'};
    my $start = $opts{'start'};
    my $stop = $opts{'stop'};
    my $nowrap = $opts{'nowrap'} || 0;
    $nowrap = 0 if $nowrap =~ /undefined/;
    #print STDERR Dumper \%opts;
    #my $change_strand = $opts{'changestrand'} || 0;
    if($add_to_seq){
      $start = $upstream if $upstream;
      $stop = $downstream if $downstream;
    }
    else
      {
	$start-=$upstream;
	$stop+=$downstream;
      }
    #print $rc;
    my $strand;
    my $seq;
    my $fasta;
    #print STDERR Dumper \%opts;
    if ($featid)
      {
	my $feat = $coge->resultset('Feature')->find($featid);
	($fasta,$seq) = ref($feat) =~ /Feature/i ?
	  $feat->fasta(
		       prot=>$pro,
		       rc=>$rc,
		       upstream=>$upstream,
		       downstream=>$downstream,
		       col=>80,
		       sep=>1,
		      )
	    :
	      ">Unable to retrieve Feature object for id: $featid\n";
	$seq = $rc ? color(seq=>$seq, upstream=>$downstream, downstream=>$upstream) : color(seq=>$seq, upstream=>$upstream, downstream=>$downstream);
		unless ($nowrap !~ /undefined/)
		{
    		$columns = 80;
        	$seq = join ("\n", wrap('','',$seq));
        	$fasta = ($fasta. "\n".$seq);  
		}
#	print STDERR join ("\n\n", $feat->genomic_sequence),"\n";
      }
    else
      {
	
	my $col= $nowrap > 0 ? 0 : 80;
	my $ds = $coge->resultset('Dataset')->find($dsid);
	$fasta = ref ($ds) =~ /dataset/i ? 
	  
	  $ds->fasta
	    (
	     start=>$start,
	     stop=>$stop,
	     chr=>$chr,
	     prot=>$pro,
	     rc=>$rc,
	     col=>$col,
	    )
	      :
		">Unable to retrieve dataset object for id: $dsid";
      }
    return $fasta;
  }
  
sub gen_foot
  {
    my $form = $FORM;
    my $featid = $form->param('featid');
    my $chr = $form->param('chr');
    my $dsid = $form->param('dsid');
    my $feat_name = $form->param('featname');
    my $rc = $form->param('rc');
    my $pro = $form->param('pro');
    my $upstream = $form->param('upstream') || 0;
    my $downstream = $form->param('downstream') || 0;
    my $start = $form->param('start');
    my $stop = $form->param('stop');
    $stop = $start unless $stop;
    my $feat = $coge->resultset('Feature')->find($featid);
    my $strand;
    my $DNAButton;
    my $RCButton;
    my $PROButton;
    my @button_loop;
    $strand = $featid ? $feat->strand : $form->param('strand');
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
    $template->param(BOTTOM_BUTTONS=>1);
    #print STDERR $featid if $featid;
   # print STDERR "nuthin" unless $featid;
    $template->param(ADDITION=>1);
    if ($rc)
      {
	$start -= $downstream;
	$stop += $upstream;
      }
    else
      {
	$start -= $upstream;
	$stop += $downstream;
      }

    if ($featid){
      my ($feat) = $coge->resultset('Feature')->find($featid);
      $dsid = $feat->dataset_id;
      $chr = $feat->chromosome;

      $template->param(PROTEIN=>'Protein Sequence');
      $template->param(SIXFRAME=>0);
      $template->param(FIND_FEATS=>1);
      $template->param(EXTEND=>"Extend Sequence");
      $template->param(UPSTREAM=>"UPSTREAM (5'): ");
      $template->param(UPVALUE=>$upstream);
      $template->param(DOWNSTREAM=>"DOWNSTREAM (3'): ");
      $template->param(DOWNVALUE=>$downstream);
      $template->param(FEATURE=>1);
      $start = $feat->start;
      $stop = $feat->stop;
    }
    else{
      $template->param(PROTEIN=>'Six Frame Translation');
      $template->param(SIXFRAME=>1);
      $template->param(FIND_FEATS=>1);
      $template->param(RANGE=>1);
      $template->param(EXTEND=>"Sequence Range");
      $template->param(UPSTREAM=>"START: ");
      $template->param(UPVALUE=>$start);
      $template->param(DOWNSTREAM=>"STOP: ");
      $template->param(DOWNVALUE=>$stop);
      $template->param(ADD_EXTRA=>1);
      $template->param(RANGE=>1);
#      $template->param(ADDUP=>$upstream);
#      $template->param(ADDDOWN=>$downstream);
      }
    if ($rc)
      {
	$start -= $downstream;
	$stop += $upstream;
      }
    else
      {
	$start -= $upstream;
	$stop += $downstream;
      }
    my ($link, $types) = find_feats(dsid=>$dsid, start=>$start, stop=>$stop, chr=>$chr);
    $template->param(FEATLISTLINK=>$link);
    $template->param(FEAT_TYPE_LIST=>$types);
   $template->param(REST=>1);
   #print STDERR $template->output."\n";
   my $html = $template->output;
   return $html;
  }
    
sub color
    {
      my %opts = @_;
      my $seq = $opts{'seq'};
#       my $rc = $opts{'rc'};
      my $upstream = $opts{'upstream'};
      my $downstream = $opts{'downstream'};
      $upstream = 0 if $upstream < 0;
      $downstream = 0 if $downstream < 0;
      my $up;
      my $down;
      my $main;
      my $nl1;
      $nl1 = 0;
      $up = substr($seq, 0, $upstream);
      while ($up=~/\n/g)
	{$nl1++;}
      my $check = substr($seq, $upstream, $nl1);
      
      $nl1++ if $check =~ /\n/;
      $upstream += $nl1;
      $up = substr($seq, 0, $upstream);
      my $nl2 = 0;
      $down = substr($seq, ((length $seq)-($downstream)), length $seq);
      while ($down=~/\n/g)
	{$nl2++;}
      $check = substr($seq, ((length $seq)-($downstream+$nl2)), $nl2);
      
      $nl2++ if $check =~ /\n/;
      $downstream += $nl2;
      $down = substr($seq, ((length $seq)-($downstream)), $downstream);
      $up = lc($up);
      $down = lc($down);
      $main = substr($seq, $upstream, (((length $seq)) - ($downstream+$upstream)));
      $main = uc($main);
      $seq = join("", $up, $main, $down);
      return $seq;
    }
    
sub gen_title
    {
      my %opts = @_;
      #print STDERR Dumper \@_;
      my $rc = $opts{'rc'} || 0;
      my $pro = $opts{'pro'};
      my $sixframe = $opts{sixframe};
      my $title;
      if ($pro)
      {
        $title = $sixframe ? "Six Frame Translation" : "Protein Sequence";
      }
      else
      {
        $title = $rc ? "Reverse Complement" : "DNA Sequence";
      }
      #print STDERR $title, "\n";
      return $title;
    }
	
sub find_feats
{
	#print STDERR "Here";
	my %opts = @_;
	my $start = $opts{'start'};
	my $stop = $opts{'stop'};
	my $chr = $opts{'chr'};
	my $dsid = $opts{'dsid'};
#	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/SeqView.tmpl');
	my $link = "FeatList.pl?";
	my %type;
	$link .="start=$start;stop=$stop;chr=$chr;dsid=$dsid";
	foreach my $ft ($coge->resultset('FeatureType')->search(
								{"features.dataset_id"=>$dsid,
								 "features.chromosome"=>$chr},
								{join=>"features",
								 select=>[{"distinct"=>"me.feature_type_id"},"name"],
								 as=>["feature_type_id","name"],
								}
							   ))
	  {
	    $type{$ft->name}=$ft->id;
	  }
	$type{All}=0;
	my $type = qq{<SELECT ID="feature_type">};
	$type .= join ("\n", map {"<OPTION value=".$type{$_}.">".$_."</option>"} sort keys %type)."\n";
	$type .= "</select>";
        return $link,$type;
}

sub generate_feat_info
  {
    my $featid = shift;
    my ($feat) = $coge->resultset("Feature")->find($featid);
    unless (ref($feat) =~ /Feature/i)
    {
      return "Unable to retrieve Feature object for id: $featid";
    }
    my $html = qq{<a href="#" onClick="\$('#feature_info').slideToggle(pageObj.speed);" style="float: right;"><img src='/CoGe/picts/delete.png' width='16' height='16' border='0'></a>};
    $html .= $feat->annotation_pretty_print_html();
    return $html;
  }
  
sub generate_gc_info
  {
    my $seq = shift;
    my $seq_type = shift;
    return "Cannot Calculate GC content of Protein Sequence" if $seq_type;
    $seq =~ s/>.*?\n//;
    $seq =~ s/\n//g;
    my $length = length($seq);
    return "No sequence" unless $length;
    my $gc = $seq =~ tr/GCgc/GCgc/;
    my $at = $seq =~ tr/ATat/ATat/;
    my $pgc = sprintf("%.2f",$gc/$length*100);
    my $pat = sprintf("%.2f",$at/$length*100);
    my $total_content = "GC: $gc (".$pgc."%)  AT: $at (".$pat."%)  total length: $length";
    return $total_content;    
  }
