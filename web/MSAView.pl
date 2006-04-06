#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use CGI::Ajax;
use CoGe::Accessory::LogUser;
use HTML::Template;
use Data::Dumper;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$USER = CoGe::Accessory::LogUser->get_user();
my $pj = new CGI::Ajax(
		       gen_html=>\&gen_html,
		      );
$pj->JSDEBUG(0);
$pj->DEBUG(0);
print $pj->build_html($FORM, \&gen_html);
#print "Content-Type: text/html\n\n". gen_html();

sub gen_html
  {
    my $body = gen_body();
    my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/MSAView.tmpl');
    $template->param(TITLE=>'Multiple Sequence Alignment Viewer');
    $template->param(USER=>$USER);
    $template->param(DATE=>$DATE);
    $template->param(LOGO_PNG=>"MSAView-logo.png");
    $template->param(BODY=>$body);
    my $html;# =  "Content-Type: text/html\n\n";
    $html .= $template->output;
    return $html;
  }

sub gen_body
  {
    my $html;
    if ($FORM->param('file'))
      {
	my $fh = $FORM->upload('file');
	my $content;
	while (<$fh>) {$content .= $_};
	$html .= gen_alignment_view($content);
      }
    else
      {
	$html .= qq{
<form method="post" enctype="multipart/form-data">
<DIV>Please select an alignment file in fasta format to upload:  
<input type = "file" name= "file"
</DIV>
<input type = "submit" label="GO">
</FORM>
}
      }
    return $html
  }

sub gen_alignment_view
  {
    my $content = shift;
    my ($seqs, $order, $stats) = parse_fasta($content);
    my $html;

    $html .= "<table>";
    $html .= "<tr>";
    $html .= "<td>";
    $html .= "<pre>\n";
    $html .= join "\n", @$order;
    $html .= "</pre>\n";
    $html .= "<td>";
    $html .= "<pre>";
    foreach my $name (@$order)
      {
	$html .= $seqs->{$name}."\n";
      }
    chomp $html;
    $html .= "</pre>";
    $html .= "</tr>";
#    $html .= "<TR>\n";
#    $html .= "<td><hr><td><hr>\n";
#    $html .= "</TR>\n";
#    $html .= "<TR>\n";
#    $html .= "<td><pre>Consensus</pre>";
#    $html .= "<td><pre>".$seqs->{consensus}."</pre>";
#    $html .= "</TR>\n";
    $html .= "</table>";
    $html .= "<hr>";
    $html .= "Stats:";
    $html .= "<table>";
    $html .= "<tr>";
    $html .= "<td><div class=small>";
    $html .= "Number of sequences";
    $html .= "</div>";
    $html .= "<td><div class=small>";
    $html .= $stats->{num_seqs};
    $html .= "</div>";

    $html .= "<tr>";
    $html .= "<td><div class=small>";
    $html .= "Alignment length";
    $html .= "</div>";
    $html .= "<td><div class=small>";
    $html .= $stats->{length}. " characters";
    $html .= "</div>";
    foreach my $num (qw(100 75 50 25 5))
      {
	$html .= "<tr>";
	$html .= "<td><div class=small>";
	$html .= "Characters with ".$num."% identity";
	$html .= "</div>";
	$html .= "<td><div class=small>";
	my $val = sprintf("%.2f", $stats->{percent}{$num}/$stats->{length})*100;
	$html .= $stats->{percent}{$num}."(".$val."%)";
	$html .= "</div>";	
      }
    $html .= "</table>";
    return $html;
  }

sub parse_fasta
  {
    my $text = shift;
    my %seqs;
    my @order;
    foreach my $ent (split /\n>/, $text)
      {
	$ent =~ s/>//g;
	my ($name, $seq) = split /\n/, $ent, 2;
	$seq =~ s/\n//g;
	$seqs{$name}.=$seq;
	push @order, $name;
      }
    my ($cons, $stats) = generate_consensus_seq([values %seqs]);
    $seqs{"consensus"} = $cons;
    push @order, "consensus";
    return (\%seqs, \@order, $stats);
  }

sub generate_consensus_seq
  {
    my $seqs = shift;
    my $num_seqs = scalar @$seqs;
    return unless $num_seqs;
    my $seq_len = length $seqs->[0];
    my $con;
    my %stats;
    $stats{length} = $seq_len;
    $stats{num_seqs} = $num_seqs;
    for (my $i = 0; $i < $seq_len; $i++)
      {
	my $let;
	my %res;
	for (my $j = 0; $j < $num_seqs; $j++)
	  {
	    my $chr = substr $seqs->[$j], $i, 1;
	    $res{$chr}++ unless $chr eq "-";
	    $let = $chr unless $let;
	    $let = " " unless $chr eq $let;
	  }

	my ($max_chr) = sort {$b <=> $a} values %res;
	$stats{percent}{100}++ unless $let eq " ";
	my $pid = $max_chr/$num_seqs;
	$stats{percent}{75}++ if $pid >= 0.75;
	$stats{percent}{50}++ if $pid >= 0.5;
	$stats{percent}{25}++ if $pid >= 0.25;
	$stats{percent}{5}++ if $pid >= 0.05;
	$let = "*" if $pid >= .75 && $let eq " ";
	$let = ":" if $pid >= .5 && $let eq " ";
	$let = "." if $pid >= .25 && $let eq " ";
	$con .= $let;
      }
    return $con, \%stats;
  }
