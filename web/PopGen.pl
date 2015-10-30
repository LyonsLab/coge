#! /usr/bin/perl -w
use strict;
use CoGe::Accessory::Web;
use CGI;
use Data::Dumper;
use HTML::Template;
use JSON::XS;
use POSIX;
no warnings 'redefine';

use vars qw($P $USER $FORM $ACCN $FID $db $PAGE_NAME $PAGE_TITLE $LINK);

$PAGE_TITLE = 'PopGen';
$PAGE_NAME  = "$PAGE_TITLE.pl";

$| = 1;    # turn off buffering

$FORM = new CGI;

my $export = $FORM->Vars->{'export'};
if ($export) {
	print "Content-Disposition: attachment; filename=PopGen.txt\n\n";
	my @columns = (index($export, 'name') != -1, index($export, 'pi') != -1, index($export, 'theta') != -1, index($export, 'td') != -1);
	my $type = $FORM->Vars->{'type'};
	my $done = 0;
	open my $fh, '<', '/home/sdavey/test.txt';
	while (my $row = <$fh>) {
		chomp $row;
		my @tokens = split '\t', $row;
		if ($type eq $tokens[0]) {
			my $line;
			for (my $i=0; $i<4; $i++) {
				if ($columns[$i]) {
					$line .= "\t" if $line;
					$line .= $tokens[$i];
				}
			}
			print $line . "\n";
			$done = 1;
		}
		elsif ($done) {
			return;
		}
	}
}

( $db, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi => $FORM,
    page_title => $PAGE_TITLE
);

my %FUNCTION = (
);
CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub gen_data {
	my $html;
	my $type = '';
	my @types;
	my $set;
	open my $fh, '<', '/home/sdavey/test.txt';
	while (my $row = <$fh>) {
		chomp $row;
		my @tokens = split '\t', $row;
		if ($type ne $tokens[0]) {
			if ($type) {
				$html .= ',' if $html;
				$html .= '\'';
				$html .= $type;
				$html .= '\':[';
				$html .= $set;
				$html .= ']';
				$set = '';
			}
			$type = $tokens[0];
			push @types, $type;
		}
		$set .= ',' if $set;
		shift @tokens;
		$set .= encode_json(\@tokens);
	}
	$html .= ',' if $html;
	$html .= '\'';
	$html .= $type;
	$html .= '\':[';
	$html .= $set;
	$html .= ']};var dt;$(document).ready(function(){dt=$(\'#main\').DataTable({data:data[\'';
	$html .= $types[1];
	$html .= '\'],lengthChange:false,order:[],ordering:false,pageLength:500,searching:false});});';
	my $select = '<select id="types" onchange="var t=this.options[this.selectedIndex].text;dt.clear().rows.add(data[t]).draw();">';
	for (@types) {
		$select .= '<option>';
		$select .= $_;
		$select .= '</option>';
	}
	$select .= '</select>';
	return 'var data={' . $html, $select;
}

sub gen_html {
    my $template = HTML::Template->new( filename => $P->{TMPLDIR} . 'PopGen.tmpl' );
    $template->param( PAGE_TITLE => 'PopGen',
		              TITLE      => 'PopGen',
    				  PAGE_LINK  => $LINK,
    				  HEAD       => qq{},
				      HOME       => $P->{SERVER},
                      HELP       => 'PopGen',
                      WIKI_URL   => $P->{WIKI_URL} || '',
    				  ADMIN_ONLY => $USER->is_admin,
                      USER       => $USER->display_name || '',
                      CAS_URL    => $P->{CAS_URL} || ''
    );
    $template->param( LOGON    => 1 ) unless ($USER->user_name eq "public");
    my ($script, $select) = gen_data();
    $template->param( script => $script );
    $template->param( select => $select );

    return $template->output;
}
