#! /usr/bin/perl -w
use strict;

use CGI;
use CGI::Carp 'fatalsToBrowser';
use Data::Dumper;
#use CoGe::Genome;
use CoGeX;
use DBIxProfiler;
use CoGe::Accessory::LogUser;
#use CoGeX;
$ENV{PATH} = "/opt/apache/CoGe/";

use vars qw( $DATE $DEBUG $TEMPDIR $TEMPURL $USER $FORM $FID $DS $CHR $LOC $ORG $VERSION $START $STOP $NAME_ONLY $coge);

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

#$DB = new CoGe::Genome;
my $connstr = 'dbi:mysql:dbname=genomes;host=biocon;port=3306';
$coge = CoGeX->connect($connstr, 'cnssys', 'CnS' );
#$coge->storage->debugobj(new DBIxProfiler());
#$coge->storage->debug(1);

print "Content-Type: text/html\n\n";
my $rhtml = gen_html(featid=>$FID, start=>$START, stop=>$STOP, chr=>$CHR, ds=>$DS, org=>$ORG, version=>$VERSION, name_only=>$NAME_ONLY) if $START > 0;
if ($START == $STOP)
{
  print "<font class=title3>Position:</font> <font class=data>$START</font><br><hr>";
}
else
{
  print "<font class=title3>Position:</font> <font class=data>$START-$STOP</font><br><hr>";
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
    my $org = $args{org};
    my $version = $args{version};
    my $name_only = $args{name_only};
#    print STDERR Dumper \%args;
    ($org) = $coge->resultset('Organism')->resolve($org) if $org;
    ($ds) = $coge->resultset('Dataset')->resolve($ds) if $ds;
    if ($ds) #there can be additional information about a chromosome for a particular version of the organism that is not in the same data_information.  Let's go find the organism_id and version for the specified data information item.
      {
	($org) = $ds->organism unless $org;
	$version = $ds->version unless $version;
	($chr) = $ds->get_chromosomes() unless $chr;
      }

#    if ($org && !$ds)
#      {
#	$version = $DB->get_dataset_obj->get_current_version_for_organism(org=>$org) unless $version;
#	($ds) = $DB->get_dataset_obj->search({version=>$version, organism_id=>$org->id});
#	($chr) = $ds->get_chromosomes() unless $chr;
#      }
    my @feats;
    #for finding associated datasets for additional annotation
#    foreach my $tmpds ($org->datasets)
#      {
#	next unless $tmpds->version eq $ds->version;
#	my $chrpass=0;
#	foreach my $tmpchr ($tmpds->get_chromosomes)
#	  {
#	    $chrpass = 1 if $tmpchr eq $chr;
#	  }
#	next unless $chrpass;
#	push @feats, $coge->get_features_in_region(dataset_id => $tmpds->id, 
#						   chr => $chr,
#						   start => $start,
#						   stop => $stop,
#						  ) if ($chr && $start && $stop);
#      }
    push @feats, $coge->get_features_in_region(dataset_id => $ds->id, 
					       chr => $chr,
					       start => $start,
					       stop => $stop,
					      ) if ($chr && $start && $stop);

    push @feats, $coge->resultset('Feature')->find($featid) if $featid;
    my $html;
    if($name_only)
      {
	$html = "<table>";
	foreach my $feat (sort {$a->type->name cmp $b->type->name} @feats)
	  {
	    next if $feat->type->name =~ /^source$/;
	    $html .= "<tr><td valign='top'>".$feat->type->name."</td><td valign='top'>".(join (", ", map {$_->name} $feat->names))."</td></tr>";
	  }
	$html .="</table>";
	return $html;
      }
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
	    $html .= qq{<font class="data">}.$feat->dataset->organism->name."</font>\n";
	    $html .= qq{<br>};
	  }
	$html .= qq{<font class="title4">Type: </font>};
	$html .= qq{<font class="data">}.$feat->type->name."</font>\n";
	$html .= qq{</table>};
	$html .= qq{<HR>};
      }
    return $html;
  }
