#! /usr/bin/perl -w
use strict;

use CGI;
use CGI::Carp 'fatalsToBrowser';
use Data::Dumper;
use CoGeX;
use DBIxProfiler;
use CoGe::Accessory::LogUser;
#use CoGeX;
$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $FID $DS $CHR $LOC $ORG $VERSION $START $STOP $NAME_ONLY $coge $GSTID);

# set this to 1 to print verbose messages to logs
$DEBUG = 0;
$TEMPDIR = "/opt/apache/CoGe/tmp";
$TEMPURL = "/CoGe/tmp";
$| = 1; # turn off buffering
($USER) = CoGe::Accessory::LogUser->get_user();
$FORM = new CGI;
$FID = $FORM->param('fid');
$DS = $FORM->param('ds');# || 61;
$CHR = $FORM->param('chr');# || 7;
$LOC = $FORM->param('loc') || $FORM->param('pos') || $FORM->param('x');# || 6049802;
$LOC = 0 unless $LOC;
$START = $FORM->param('start');
$START = $LOC unless defined $START;
$STOP = $FORM->param('stop');
$STOP = $START unless defined $STOP;
$ORG = $FORM->param('org') || $FORM->param('organism');
$VERSION = $FORM->param('version') || $FORM->param('ver');
$NAME_ONLY = $FORM->param('name_only') || 0;
$GSTID = $FORM->param('gstid'); #genomic_sequence_type_id

$coge = CoGeX->dbconnect();
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

print "Content-Type: text/html\n\n";
my $rhtml = gen_html(featid=>$FID, start=>$START, stop=>$STOP, chr=>$CHR, ds=>$DS, org=>$ORG, version=>$VERSION, name_only=>$NAME_ONLY, gstid=>$GSTID) if $START > 0 || $FID;
if ($START && $STOP)
  {
    if ($START == $STOP)
      {
	$rhtml =  "<font class=title3>Position:</font> <font class=data>$START</font><br><hr>".$rhtml;
      }
    else
      {
	$rhtml =  "<font class=title3>Position:</font> <font class=data>$START-$STOP</font><br><hr>".$rhtml;
      }
  }
$rhtml = "No annotations" unless $rhtml;
print $rhtml;

sub gen_html
  {
    my %args = @_;
    my $featid = $args{featid};
    my $start = $args{start};
    my $stop = $args{stop};
    $stop = $start unless $stop;
    my $chr = $args{chr};
    my $ds = $args{ds};
    my $version = $args{version};
    my $name_only = $args{name_only};
    my $gstid = $args{gstid};
    ($ds) = $coge->resultset('Dataset')->resolve($ds) if $ds;
    if ($ds) #there can be additional information about a chromosome for a particular version of the organism that is not in the same data_information.  Let's go find the organism_id and version for the specified data information item.
      {
	$version = $ds->version unless $version;
	($chr) = $ds->get_chromosomes() unless $chr;
      }
    my @feats;
    push @feats, $coge->get_features_in_region(dataset_id => $ds->id, 
					       chr => $chr,
					       start => $start,
					       stop => $stop,
					      ) if ($chr && $start && $stop);

    push @feats, $coge->resultset('Feature')->find($featid) if $featid;
    return " " unless @feats;
    my $html;
    if($name_only)
      {
	$html = "<table>";
	foreach my $feat (sort {$a->type->name cmp $b->type->name} @feats)
	  {
	    next if $feat->type->name =~ /^source$/;
	    $html .= "<tr><td valign='top'>".$feat->type->name."</td><td valign='top'>".(join (", ",$feat->names))."</td></tr>";
	  }
	$html .="</table>";
	return $html;
      }
    foreach my $feat (sort {$a->type->name cmp $b->type->name} @feats)
      {
	next if $feat->type->name eq "chromosome";
	next if $feat->type->name =~ /^source$/;
	my $color;
	if ($feat->type->name eq "CDS")
	  {
	    $color = "#DDFFDD"
	  }
	elsif ($feat->type->name eq "mRNA")
	  {
	    $color = "#DDDDFF";
	  }
	elsif ($feat->type->name =~ "rna")
	  {
	    $color = "#DDDDDD";
	  }
	elsif ($feat->type->name =~ "gene")
	  {
	    $color = "#FFDDDD";
	  }
	else
	  {
	    $color = "#FFDDBB";
	  }
	$html .= "<table bgcolor=$color width=100%><tr><td>".$feat->annotation_pretty_print_html(loc_link=>"SeqView.pl", gstid=>$gstid);
	$html .= qq{<font class="title4">Type: </font>};
	$html .= qq{<font class="data">}.$feat->type->name."</font>\n";
	$html .= qq{</table>};
	$html .= qq{<HR>};
      }
    return $html;
  }
