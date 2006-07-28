#! /usr/bin/perl -w
use strict;

use CGI;
use CGI::Carp 'fatalsToBrowser';
use Data::Dumper;
use CoGe::Genome;
use CoGe::Accessory::LogUser;

$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $FID $DI $CHR $LOC);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;
$FID = $FORM->param('fid');
$DI = $FORM->param('di');# || 61;
$CHR = $FORM->param('chr');# || 7;
$LOC = $FORM->param('loc') || $FORM->param('pos') || $FORM->param('x');# || 6049802;
$LOC = 0 unless $LOC;
$DB = new CoGe::Genome;
print "Content-Type: text/html\n\n";
my $rhtml = gen_html(featid=>$FID, loc=>$LOC, chr=>$CHR, di=>$DI) if $LOC > 0;
print "<font class=title3>Position:</font> <font class=data>$LOC</font><br><hr>";
$rhtml = "No annotations" unless $rhtml;
print $rhtml;

sub gen_html
  {
    my %args = @_;
    my $featid = $args{featid};
    my $loc = $args{loc};
    my $chr = $args{chr};
    my $di = $args{di};
    my ($dio) = $DB->get_data_info_obj->search({data_information_id=>$di});
    my @feats;
    foreach my $tdio ($dio->get_associated_data_infos)
      {
	push @feats, $DB->get_feat_obj->get_features_in_region(info_id => $tdio->id, 
							       chr => $chr,
							       start => $loc,
							       stop => $loc,
							      ) if ($chr && $loc);
      }
    push @feats, $DB->get_feat_obj->search(feature_id=>$featid) if $featid;
    my $html;
    foreach my $feat (sort {$a->type->name cmp $b->type->name} @feats)
      {
	next if $feat->type->name =~ /^source$/;
	$html .= $feat->annotation_pretty_print_html();
	unless ($FORM->param('no_org'))
	  {
	    $html .= qq{<font class="title4">Organism: </font>};
	    $html .= qq{<font class="data">}.$feat->data_info->org->name."</font>\n";
	    $html .= qq{<br>};
	  }
	$html .= qq{<font class="title4">Type: </font>};
	$html .= qq{<font class="data">}.$feat->feat_type->name."</font>\n";
	$html .= qq{<br>};
	

	$html .= qq{<HR>};
      }
    return $html;
  }
