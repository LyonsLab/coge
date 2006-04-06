#! /usr/bin/perl -w
use strict;
use CGI;
use CGI::Carp 'fatalsToBrowser';
use Data::Dumper;
use CoGe::Genome;
use CoGe::Accessory::LogUser;

$ENV{PATH} = "/opt/apache2/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $FID $DI $CHR $LOC);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
$DATE = sprintf( "%04d-%02d-%02d %02d:%02d:%02d",
		 sub { ($_[5]+1900, $_[4]+1, $_[3]),$_[2],$_[1],$_[0] }->(localtime));

$FORM = new CGI;
$FID = $FORM->param('fid');
$DI = $FORM->param('di');# || 61;
$CHR = $FORM->param('chr');# || 7;
$LOC = $FORM->param('loc') || $form->param('pos') || $form->param('x');# || 6049802;
$USER = CoGe::Accessory::LogUser->get_user();
$DB = new CoGe::Genome;
gen_html(featid=>$FID, loc=>$LOC, chr=>$CHR, di=>$DI);

sub gen_html
  {
    my %args = @_;
    my $featid = $args{featid};
    my $loc = $args{loc};
    my $chr = $args{chr};
    my $di = $args{di};
    my @feats;
    push @feats, $DB->get_feat_obj->get_features_in_region(info_id => $di, 
							   chr => $chr,
							   start => $loc,
							   stop => $loc,
							  ) if ($di && $chr && $loc);
    push @feats, $DB->get_feat_obj->search(feature_id=>$featid) if $featid;
    foreach my $feat (@feats)
      {
	print qq{<font class="annotation">Organism</font>};
	print qq{<li>};
	print $feat->data_info->org->name,"\n";
	print qq{<br>};
	print qq{<font class="annotation">Type</font>};
	print qq{<li>};
	print $feat->feat_type->name,"\n";
	print qq{<br>};
	
	print $feat->annotation_pretty_print_html();
	print qq{<HR>};
      }
  }
