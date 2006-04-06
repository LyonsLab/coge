#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use HTML::Template;
use Data::Dumper;
use CoGe::Accessory::LogUser;

my $email = CoGe::Accessory::LogUser->get_user();
my $template = HTML::Template->new(filename=>'/opt/apache/CoGe/tmpl/index.tmpl');
my $DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$template->param(TITLE=>'Comparative Genomics Homepage');
$template->param(USER=>$email);
$template->param(DATE=>$DATE);

my @actions = (
	       {
		ACTION=>qq{<a href="./B2V.pl">Comparative BL2SEQ Viewer</a>},
		DESC => qq{The comparative bl2seq viewer uses the bl2seq algorithm to find regions of similarity between two sequences and displays the results graphically.  This program has several unqiue options to aid in the finding of evolutionary conserved non-coding sequences.},
	       },
	       {

		ACTION=>qq{<a href="./FeatView.pl">Genomic Location Feature Viewer</a>},
		DESC => qq{Allows users to find and display information about a genomic feature.  A <font class=oblique>feature</font> in terms of sequence or genomic data denotes some genomic region that has something known about it.  For example, if some genomic region is known to be a gene, and that gene makes an mRNA transcript, and that mRNA is translated into a functional protein, and that protein has several functional domains, it will be likely that there will be several genomic features in that region:  One to denote the gene, one or more to denote the mRNA, and one or more to denote the protein, and one or more to denote the functional domains.  },
	       },
	       {
		ACTION => qq{<a href="./MSAView.pl">Multiple Sequence Alignment Viewer</a>},
		DESC   => qq{Allows users to submit a multiple sequence alignment in FASTA format (if people would like additional formats, please request via e-mail).  This program also generates a consensus sequence from the alignment and displays some basic statistics about the alignment.},
		},
	      );
$template->param(ACTIONS=>[@actions]);



my $html =  "Content-Type: text/html\n\n";
$html .= $template->output;
print $html;
