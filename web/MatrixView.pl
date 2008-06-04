#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;
use CGI::Ajax;
#use CoGeX;
use Benchmark;
use CoGe::Accessory::Web;
use CoGe::Accessory::genetic_code;
use Statistics::Basic::Mean;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $MATRIXDIR $USER $FORM $coge $connstr);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache2/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$MATRIXDIR = "/opt/apache2/CoGe/data/blast/matrix/aa/";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
($USER) = CoGe::Accessory::LogUser->get_user();

my $pj = new CGI::Ajax(
		       gen_data=>\&gen_data,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n";print gen_html($FORM);

sub gen_html
  {
    my $html;
    unless ($USER && $USER->user_name  !~/public/i)
      {
	$html = login();
      }
    else
      {
	my ($body) = gen_body();
	
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	
	$template->param(TITLE=>'Sequence Alignment Matrix View');
	$template->param(HEAD=>qq{});
	my $name = $USER->user_name;
        $name = $USER->first_name if $USER->first_name;
        $name .= " ".$USER->last_name if $USER->first_name && $USER->last_name;
        $template->param(USER=>$name);

	$template->param(LOGON=>1) unless $USER->user_name eq "public";
	$template->param(DATE=>$DATE);
	$template->param(LOGO_PNG=>"MatrixView-logo.png");
	$template->param(BODY=>$body);
	$html .= $template->output;
      }
    return $html;
  }
  
sub gen_body
    {
      my $form = shift || $FORM;
      my $matrix = $form->param('matrix');
      my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/MatrixView.tmpl');
      my $matrices = process_dir(matrix=>$matrix);
      $template->param(MATRICES=>$matrices);
      my $data = gen_data(file=>$matrix);
      $template->param(DATA=>$data);
      return $template->output;
      
    }

sub gen_data
  {
    my %opts = @_;
    my $form = $opts{form} || $FORM;
    my $code_layout = $opts{code_layout} || 0;
    my $bin_size = $opts{bin_size} || 5;
    my $file = $opts{file};
    return unless $file;
    my $code = CoGe::Accessory::genetic_code->code;
    $code = $code->{code};
    my $html ="<table>";
    my ($data) = process_file(file=>$file);

    if ($data->{P})
      {
	my $aa_sort = CoGe::Accessory::genetic_code->sort_aa_by_gc();
	$html .= "<tr><th><th>".join ("<th>", sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort);
	$html .= "<tr>";
	foreach my $aa1 (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
	  {	
	    $html .= "<th>$aa1";
	    #	print STDERR join ("\t",  sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %aa_sort),"\n";
	    my ($max) = sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %$aa_sort;
	    my ($min1, $min2) = sort {$a<=>$b} map {$data->{$aa1}{$_}} keys %$aa_sort; #min one will be for the stop aa, we'll skip this since it is always rediculous 
	    foreach my $aa2 (sort {$aa_sort->{$b} <=> $aa_sort->{$a} || $a cmp $b}keys %$aa_sort)
	      {	
		my $val = $data->{$aa1}{$aa2};
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
#		my $color = color_by_usage($max+abs($min2), $val+abs($min2));
#		$html .= "<td style=\"background-color: rgb($color,255,$color)\">".$val;
	      }
	    $html .= "<tr>";
	  }
      }
    else
      {
	$html .= "<tr><th><th>".join ("<th>", sort { sort_nt1(substr($a, 0, 1)) <=> sort_nt1(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt3(substr($a,2,1)) <=> sort_nt3(substr($b,2,1)) }keys  %$data);
	$html .= "<tr>";
	foreach my $aa1 (sort { sort_nt1(substr($a, 0, 1)) <=> sort_nt1(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt3(substr($a,2,1)) <=> sort_nt3(substr($b,2,1)) }keys  %$data)
	  {
	    $html .= "<th>$aa1";
	    #	print STDERR join ("\t",  sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %aa_sort),"\n";
	    my ($max) = sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %$data;
	    my ($min1, $min2, $min3, $min4) = sort {$a<=>$b} map {$data->{$aa1}{$_}} keys %$data; #min one will be for the stop aa, we'll skip this since it is always rediculous 
	foreach my $aa2 (sort { sort_nt1(substr($a, 0, 1)) <=> sort_nt1(substr($b,0, 1)) || sort_nt2(substr($a,1,1)) <=> sort_nt2(substr($b,1,1)) || sort_nt3(substr($a,2,1)) <=> sort_nt3(substr($b,2,1)) } keys  %$data)#	    foreach my $aa2 (sort keys %$data)
	      {	
		my $val = $data->{$aa1}{$aa2};
		my $color;
		if ($val >0)
		  {
		    $color = color_by_usage($max, $val);
		    $html .= "<td style=\"background-color: rgb($color,255,$color)\">".$val;
		  }
		else
		  {
		    $color = color_by_usage(abs($min4), $val*-1);
		    $html .= "<td style=\"background-color: rgb(255,$color,$color)\">".$val;
		  }
	      }
	    $html .= "<tr>";
	  }
      }
    $html .= "</table>";
    return $html;
  }

sub process_dir
  {
    my %opts = @_;
    my $matrix = $opts{matrix};
    my %matrices;
    opendir (DIR, $MATRIXDIR);
    while (my $file = readdir(DIR))
      {
	next if $file =~ /^\./;
	$matrices{"$file"}=uc($file) unless $matrices{uc($file)} || $matrices{lc($file)};
      }
    closedir DIR;
    my $html = "<select id=matrix>";
    $html .= join ("\n", map {"<option value=\"".$_."\">".$matrices{$_}."</option>"} sort {$matrices{$a} cmp $matrices{$b}} keys %matrices )."\n";
    $html =~ s/(value=\"$matrix\")/selected $1/i if $matrix;
    $html .= "</select>";
    return $html;
  }

sub process_file
  {
    my %opts = @_;
    my $bin_size = $opts{bin_size} || 5;
    $bin_size = 1 if $bin_size < 1;
    $bin_size = 100 if $bin_size > 100;
    my $skip = $opts{skip} || [];#[qw(mitochond chloroplast virus phage)];
    my $keep = $opts{keep} || [];#[qw(mitochondr)];
    my $file = $opts{file};
    $file = $MATRIXDIR."/".$file unless $file =~ /$MATRIXDIR/i;
    my %data;
    my @head;
    my $max = 0;
    my $min = 10000;
    print STDERR $file,"\n";
    open (IN, $file);
    while (<IN>)
      {
	chomp;
	next if /^#/;
	my @line = split /\s+/;
	unless ($line[0])
	  {
	    shift @line;
	    @head = @line;
	    next;
	  }
	my $aa = shift @line;
	for (my $i = 0; $i <= $#line; $i++)
	  {
	    $data{$aa}{$head[$i]} = $line[$i];
	    $max = $line[$i] if $line[$i] > $max;
	    next if $aa eq "*";
	    next if $head[$i] eq "*";
	    $min = $line[$i] if $line[$i] < $min;
	  }
      }
    close IN;
    return \%data, $max, $min;
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

1;
