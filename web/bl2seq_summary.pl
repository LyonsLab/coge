#! /usr/bin/perl -w

## cns/bl2seq_summary.pl
#  This script presents a summary of a bl2seq run and can focus on
#  individual HSPs if required
#
#  It needs the web form param "blast_report" and the web button submit
#  value set to "GO".  Optionally, it can take param "hsp" to show the
#  summary for a specific hsp, and it can also take the params "begin"
#  and "end" and "other", to provide the genome sequence index for the
#  begin, end of the sequences, and any other params that need to come
#  through

use strict;
use CGI;
#use CGI::Carp 'fatalsToBrowser';
#use CNS::MyDB;
use GS::bl2seqReport;
use Data::Dumper;

# for security purposes
$ENV{PATH} = "/opt/apache2/cns/";
delete @ENV{ 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

use vars qw( $DATE $DEBUG );

# set this to 1 to print verbose messages to logs
$DEBUG = 1;

$| = 1; # turn off buffering
my $form = new CGI;
$DATE = sprintf( "%04d%02d%02d%02d%02d%02d",
        sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

# if form with data... start process loop
if ( $form->param() ) {
	if ( ( $form->param("submit") eq "GO" ) and 
			 ( $form->param("blast_report") ne "" ) ) { 
			 		# got an blast_report and the "GO" signal
		if ( not -f $form->param("blast_report") ) {
			&Signal_Error( $form, "FILE" );
		} else {
			&Show_Summary( $form );
		}
	} else {
		# redraw the blank form because there was a problem with user data
		&Show_Form( $form );
	}
} else {
	# draw blank form for first time
	&Show_Form( $form );
}

sub Show_Summary {
	my $form = shift;
	# set up output page
	print $form->header();
	print $form->start_html(-bgcolor=>"#CCCCCC", -title=>"At Pair HSP Summary");
	my($rc,$data) = (0,0);
	($rc,$data) = &parse_bl2seq_x( $form->param("blast_report"),
							$form->param("accnq"), $form->param("accns"));
		# data is an arrayref to the HSPs ($data->[0] is most signif hit)
	if ( not $rc ) {
		&Signal_Error($form, "PARSE");
	} else {
		if ( defined $form->param("hsp") and $form->param("hsp") ne "" ) {
			my $hsp = $form->param("hsp")-1; #to account for array[0]
			if ( exists $data->[ $hsp ] ) {
				my $html = build_table_hsp( $form->param("blast_report"),
						$form->param("accnq"), $form->param("accns"),
						$data->[ $hsp ],
						$form->param("qbegin"), $form->param("qend"),
						$form->param("sbegin"), $form->param("send"),
						);
						#$form->param("other"));
				print "$html\n";
			} else {
				&Signal_Error($form, "NOHSP");
			}
		}
	}
	print $form->end_html();
}

sub build_table_hsp {
	my $blast_report = shift;
	my $accnq = shift;
	my $accns = shift;
	my $hsp = shift;
	my $qbegin = shift;
	my $qend = shift;
	my $sbegin = shift;
	my $send = shift;
	my $other = shift;
	my $html = qq!<TABLE width="70%" border="0" cellspacing="0" cellpadding="2">\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#998899">BLAST REPORT</td>!;
	$html .= qq!<td bgcolor="#998899" colspan=3>!;
	$html .= $blast_report;
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">QUERY ACCN</td>!;
	$html .= qq!<td bgcolor="#999999">!;
	$html .= $accnq;
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">SBJCT HSP NO.</td>!;
	$html .= qq!<td bgcolor="#999999">!;
	$html .= $hsp->{number};
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;

	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">QUERY BEGIN</td>!;
	$html .= qq!<td bgcolor="#999999">!;
	$html .= $hsp->{qb};
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">QUERY END</td>!;
	$html .= qq!<td bgcolor="#999999">!;
	$html .= $hsp->{qe};
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">GENOME BEGIN Query</td>!;
	$html .= qq!<td bgcolor="#999999">!;
	$html .= $qbegin;
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">GENOME END Query</td>!;
	$html .= qq!<td bgcolor="#999999">!;
	$html .= $qend;
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;


	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">SBJCT ACCN</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $accns;
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">SBJCT HSP NO.</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{number};
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">SBJCT BEGIN</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{sb};
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">SBJCT END</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{se};
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">GENOME BEGIN Sbjct</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $sbegin;
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">GENOME END Sbjct</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $send;
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;
	
	$html .= qq!<tr>\n!;
	$html .= qq!<td colspan=4 width="30%" bgcolor="#998888">HSP INFORMATION</td>!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">HSP MATCH</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{hspmatch} . " / " . $hsp->{length};
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">HSP ORIENTATION</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{orientation};
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">HSP IDENTITY</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{identity};
	$html .= qq!</td>\n!;
	$html .= qq!<td width="30%" bgcolor="#999988">HSP EVAL</td>!;
	$html .= qq!<td bgcolor="#999988">!;
	$html .= $hsp->{eval};
	$html .= qq!</td>\n!;
	$html .= qq!</tr>\n!;

	$html .= qq!<tr>\n!;
	$html .= qq!<td colspan=4 width="30%" bgcolor="#998888">HSP SEQUENCE</td>!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">QUERY HSP</td>!;
	$html .= qq!<td colspan=3 bgcolor="#999999"><PRE>!;
	$html .= uc(seqwrap($hsp->{qmatchseq}));
	$html .= qq!</PRE></td>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!<tr>\n!;
	$html .= qq!<td width="30%" bgcolor="#999999">SBJCT HSP</td>!;
	$html .= qq!<td colspan=3 bgcolor="#999999"><PRE>!;
	$html .= uc(seqwrap($hsp->{smatchseq}));
	$html .= qq!</PRE></td>\n!;
	$html .= qq!</tr>\n!;



	$html .= qq!<tr>\n!;
	$html .= qq!</tr>\n!;
	$html .= qq!</TABLE>!;


	return($html);
}

sub Show_Form {
	# this isn't ussually called directly, so it is pretty sparse
	# Most of the time, this script is call from another page as an
	# embedded URL, like
	# http://toxic.berkeley.edu/cns/bl2seq_summary.pl?blast_report="XX"&submit="GO"
	#
	my $form = shift;
	print $form->header();
	print $form->start_html(-title=>'BL2SEQ Summary Form');
	print $form->center($form->h1('BL2SEQ Summary'));
	print $form->start_multipart_form(-method=>'get',
			-action=> $form->url( -relative ) );
	print "\n";

	# accn entry textfield
	print "BL2SEQ Report name: " .
		$form->textfield(-name=>'blast_report', -size=>100 );
	print "\n";
	print $form->p();
	print "<TABLE border=0><TR>\n";
	print "<TD>Query Accn (optional)</TD><TD>" . 
			$form->textfield(-name=>'accnq', -size=>18 );
	print "</TD></TR><TR>";
	print "<TD>Sbjct Accn (optional)</TD><TD>" .
			$form->textfield(-name=>'accns', -size=>18 );
	print "</TD></TR><TR>";
	print "<TD>HSP Number (optional)</TD><TD>" . $form->textfield(-name=>'hsp', -size=>8 );
	print "</TD></TR><TR>";
	print "<TD>Genome Begin Query(optional)</TD><TD>" .
				$form->textfield(-name=>'qbegin', -size=>12 );
	print "</TD></TR><TR>";
	print "<TD>Genome End Query(optional)</TD><TD>" .
				$form->textfield(-name=>'qend', -size=>12 );
	print "</TD></TR><TR>";
	print "<TD>Genome Begin Sbjct (optional)</TD><TD>" .
				$form->textfield(-name=>'sbegin', -size=>12 );
	print "</TD></TR><TR>";
	print "<TD>Genome End Sbjct (optional)</TD><TD>" .
				$form->textfield(-name=>'send', -size=>12 );
	print "</TD></TR>";
	print "<TR>";
#	print "<TD>Other Info (optional)</TD><TD>" .
#					$form->textfield(-name=>'other', -size=>72 );
#	print "</TD></TR>";
	print "</TABLE>\n";
	print "\n";

	# submit and reset buttons
	print $form->submit(-name=>"submit",-value=>"GO"),
		$form->reset(-name=>"Reset", -value=>"Reset");

	print $form->end_form();
	print $form->end_html();
}

sub Signal_Error {
	my $form = shift;
	my $error = shift;

	if ( $error eq "FILE" ) {
		print $form->header();
		print $form->start_html(-title=>'bl2seq summary');
		print $form->center($form->h1('bl2seq summary problem'));
		print $form->p();
		print $form->center($form->h3("FILE PROBLEM"));
		print $form->p();
		print $form->center($form->h3('Please make sure the file exists'));
		print $form->p();
		print "File was: ", $form->param("blast_report");
		print $form->p();
		print $form->end_html();
	} elsif ( $error eq "PARSE" ) {
		print $form->header();
		print $form->start_html(-title=>'bl2seq summary');
		print $form->center($form->h1('bl2seq summary problem'));
		print $form->p();
		print $form->center($form->h3("FILE PARSING PROBLEM"));
		print $form->p();
		print $form->center($form->h3('There was an error in parsing'));
		print $form->p();
		print "File was: ", $form->param("blast_report");
		print $form->p();
		print $form->end_html();
	} elsif ( $error eq "NOHSP" ) {
		print $form->header();
		print $form->start_html(-title=>'bl2seq summary');
		print $form->center($form->h1('bl2seq summary problem'));
		print $form->p();
		print $form->center($form->h3("MISSING HSP"));
		print $form->p();
		print $form->center($form->h3('There is no such HSP'));
		print $form->p();
		print "File was: ", $form->param("blast_report");
		print $form->p();
		print $form->end_html();
	}
}

sub parse_bl2seq_x {
	# takes a bl2seq report, processes it with bl2seqReport module, then
	# builds an array of hash references.  The hashes contain
	# information about each hsp.  Success/Fail and a reference to the
	# final array is returned.  The final array will be ranked by best
	# hsp first, etc.
	# Also sent into this subroutine is the optional parameters from the
	# form about filtering on the length (T/F) and the actual length to filter
	# on (integer)
	my $blast_report = shift;
	my $accn_q = shift;
	my $accn_s = shift;
	#my $bl2seq_params = shift;
	my(@hsplist) = ();

	my $report = "";
	$report = new bl2seqReport( $blast_report );
	my $count = 1;
	while ( my $hsp = $report->nextHSP() ) {
		my %match = ();
		$match{"blast_report"} = $blast_report;
		$match{"number"} = $count;
		$match{"accn_q"} = $accn_q;
		$match{"accn_s"} = $accn_s;
		$match{"qb"} = $hsp->qb;
		$match{"qe"} = $hsp->qe;
		$match{"sb"} = $hsp->sb;
		$match{"se"} = $hsp->se;
		$match{"hspmatch"} = $hsp->match;
		$match{"length"} = $hsp->length;
		$match{"identity"} = $hsp->percent;
		$match{"eval"} = $hsp->P;
		$match{"orientation"} = $hsp->strand;
		$match{"qmatchseq"} = $hsp->queryAlignment;
		$match{"smatchseq"} = $hsp->sbjctAlignment;
		push(@hsplist, \%match);
		$count++;
	}
	if ( @hsplist > 0 ) {
		return( 1, \@hsplist);
	} else {
		return( 0,0 );
	}
}

sub seqwrap {
	my $seq = shift;
	my $wrap = shift || 50;
	my @chars = split( //, $seq );
	my $newseq = "";
	my $counter = 0;
	foreach my $char ( @chars ) {
	if ( $counter < $wrap ) {
		$newseq .= $char;
	} else {
		$newseq .= "\n" . $char;
		$counter = 0;
	}
		$counter++;
	}
	return($newseq);
}
