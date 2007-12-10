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
#$connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
#$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

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
    unless ($USER && $USER  !~/public/i)
      {
	$html = login();
      }
    else
      {
	my ($body) = gen_body();
	
	my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/generic_page.tmpl');
	
	$template->param(TITLE=>'CoGe: Secret Projects: Matrix View');
	$template->param(HEAD=>qq{});
	$template->param(USER=>$USER);
	$template->param(DATE=>$DATE);
	$template->param(LOGO_PNG=>"CoGe-logo.png");
	$template->param(BODY=>$body);
	$html .= $template->output;
      }
    return $html;
  }
  
sub gen_body
    {
      my $form = shift || $FORM;
      my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/MatrixView.tmpl');
      my $matrices = process_dir();
#      my ($data) = gen_data();
      $template->param(MATRICES=>$matrices);
#      $template->param(DATA=>$data);
#      my $groups = join ("<br>", map {qq{<INPUT TYPE=checkbox NAME=groups id=groups value=$_ checked>$_}} keys %$types)."<br>";
      #      $template->param(GROUPS=>$groups);
      return $template->output;
      
    }

sub gen_data
  {
    my %opts = @_;
    my $form = $opts{form} || $FORM;
    my $code_layout = $opts{code_layout} || 0;
    my $bin_size = $opts{bin_size} || 5;
    my $file = $opts{file};
    my $code = CoGe::Accessory::genetic_code->code;
    $code = $code->{code};
    
    my %aa2codon;
    foreach my $item (keys %$code)
      {
	push @{$aa2codon{$code->{$item}}},$item;
      }
    my $html;
    my %aa_sort = map {$_,{}} keys %aa2codon;
    foreach my $codon (keys %$code)
      {
	my $gc = $codon =~ tr/GCgc/GCgc/;
	my $at = $codon =~ tr/ATat/ATat/;
	$aa_sort{$code->{$codon}}{AT}+=$at;
	$aa_sort{$code->{$codon}}{GC}+=$gc;
      }
    %aa_sort = map {$_,($aa_sort{$_}{GC}/($aa_sort{$_}{AT}+$aa_sort{$_}{GC}))} keys %aa_sort;

    my ($data, $max_val, $min_val) = process_file(file=>$file);
    $html .= "<table>";
    $html .= "<tr><th>".join ("<th>", "GC% (org count)", sort {$aa_sort{$b} <=> $aa_sort{$a} || $a cmp $b}keys %aa_sort);
    $html .= "<tr>";
    foreach my $aa1 (sort {$aa_sort{$b} <=> $aa_sort{$a} || $a cmp $b}keys %aa_sort)
      {	
	$html .= "<th>$aa1";
#	print STDERR join ("\t",  sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %aa_sort),"\n";
	my ($max) = sort {$b<=>$a} map {$data->{$aa1}{$_}} keys %aa_sort;
	my ($min1, $min2) = sort {$a<=>$b} map {$data->{$aa1}{$_}} keys %aa_sort; #min one will be for the stop aa, we'll skip this since it is always rediculous 
	foreach my $aa2 (sort {$aa_sort{$b} <=> $aa_sort{$a} || $a cmp $b}keys %aa_sort)
	  {	
	    my $val = $data->{$aa1}{$aa2};
	    my $color = color_by_usage($max+abs($min2), $val+abs($min2));
	    $html .= "<td style=\"background-color: rgb($color,255,$color)\">".$val;
	  }
	$html .= "<tr>";
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
	$matrices{"$MATRIXDIR/$file"}=uc($file) unless $matrices{"$MATRIXDIR/".uc($file)} || $matrices{"$MATRIXDIR/".lc($file)};
      }
    closedir DIR;
    my $html = "<select id=matrix>";
    $html .= join ("\n", map {"<option value=\"".$_."\">".$matrices{$_}."</option>"} sort {$matrices{$a} cmp $matrices{$b}} keys %matrices )."\n";
    $html .= "</selected>";
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
	@line = sort {$b <=> $a} @line;

	  

      }
    close IN;
    return \%data, $max, $min;
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

sub commify
    {
      my $text = reverse $_[0];
      $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
      return scalar reverse $text;
    }
sub color_by_usage
      {
	my ($max,$value, $opt) = @_;
	$opt = 255 unless $opt;
	return $opt unless $max;
	my $g = $opt*(($max - $value) / $max);
	return int($g + .5);
      }


1;
