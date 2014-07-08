#! /usr/bin/perl -w
use strict;

use CGI;
use CGI::Carp 'fatalsToBrowser';
use Data::Dumper;
use CoGeX;
use CoGe::Accessory::LogUser;

$ENV{PATH} = "/opt/apache/bpederse/cogex/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $DB $FORM $FID $DS $CHR $LOC $ORG $VERSION $START $STOP);

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

my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $DB = CoGeX->connect($connstr, 'USER', 'PASSWORD' );

print "Content-Type: text/html\n\n";
my $rhtml = gen_html(featid=>$FID, start=>$START, stop=>$STOP, chr=>$CHR, ds=>$DS, org=>$ORG, version=>$VERSION) if $START > 0;
print "<font class=title3>Position:</font> <font class=data>$START-$STOP</font><br><hr>";
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
    my $org = $args{org};
    my $version = $args{version};
##working here, taken frmo Graphics.pm
    $org = $DB->search('Organism')->resolve($org) if $org;
    $ds =  $DB->search('Dataset')->search(
            { 'dataset_id' => $ds}
            ,{'prefetch' => ['organism']});
                        if $ds;
    if ($ds) #there can be additional information about a chromosome for a particular version of the organism that is not in the same data_information.  Let's go find the organism_id and version for the specified data information item.
      {
	$org = $ds->organism->name unless $org;
	$version = $ds->version unless $version;
	#($chr) = $ds->get_chromosomes() unless $chr;
      }

    if ($org && !$ds)
      {
	$version =
    $DB->get_dataset_obj->get_current_version_for_organism(org=>$org) unless $version;
	($ds) = $DB->get_dataset_obj->search({version=>$version, organism_id=>$org->id});
	($chr) = $ds->get_chromosomes() unless $chr;
      }
###end working section
    my ($dso) = $ds;#$DB->get_dataset_obj->retrieve($ds);
    my @feats;
    foreach my $tdso ($dso->get_associated_datasets)
      {
	push @feats, $DB->get_feat_obj->get_features_in_region(dataset => $tdso->id,
							       chr => $chr,
							       start => $start,
							       stop => $stop,
							      ) if ($chr && $start && $stop);
      }
    push @feats, $DB->get_feat_obj->retrieve($featid) if $featid;
    my $html;
    foreach my $feat (sort {$a->type->name cmp $b->type->name} @feats)
      {
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
	$html .= "<table bgcolor=$color width=100%><tr><td>".$feat->annotation_pretty_print_html(loc_link=>"SeqView.pl");
	unless ($FORM->param('no_org'))
	  {
	    $html .= qq{<font class="title4">Organism: </font>};
	    $html .= qq{<font class="data">}.$feat->dataset->org->name."</font>\n";
	    $html .= qq{<br>};
	  }
	$html .= qq{<font class="title4">Type: </font>};
	$html .= qq{<font class="data">}.$feat->feat_type->name."</font>\n";
	$html .= qq{</table>};
	$html .= qq{<HR>};
      }
    return $html;
  }
