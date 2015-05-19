#!/usr/bin/perl -w

use strict;
use CoGeX;
use DBIxProfiler;
use Data::Dumper;
my $orgid = shift;
#die "Usage: $0 <org_id>\n" unless $orgid;

$| =1;
my $connstr = 'dbi:mysql:dbname=DB;host=HOST;port=PORT';
my $coge = CoGeX->connect($connstr, 'USER', 'PASSWORD' );
$coge->storage->debugobj(new DBIxProfiler());
$coge->storage->debug(1);

if ($orgid)
  {codon_usage_for_org($orgid);}
else
  {
    for (my $i = 1; $i<50; $i++)
      {
	codon_usage_for_org($i);
      }
  }
sub codon_usage_for_org
  {
    my $orgid = shift;
    my %codon;
    my ($org) = $coge->resultset('Organism')->resolve($orgid);
    return unless $org;
    foreach my $ds ($org->current_datasets())
      {
	print STDERR "DATASET: ". $ds->name,": ", $ds->description," (v",$ds->version,")\n";
#	next;
	foreach my $feat ($ds->features({"feature_type.name"=>"CDS"},{join=>"feature_type"}))
	  {
#	    next unless $feat->type->name eq "CDS";
#	    print STDERR join (", ", sort $feat->names)," . . . ";
	    my ($codon) = $feat->codon_frequency(counts=>1);
	    print Dumper $codon;
#	    my $seq = $feat->genomic_sequence;
#	    my $x=0;
#	    while ($x<length($seq))
#	      {
#		$codon{uc(substr($seq, $x, 3))}++;
#		$x+=3;
#	      }
#	    print STDERR $x/3, " codons processed\n";
	  }
      }

    foreach my $codon (sort keys %codon)
      {
	print $codon,"\t",$codon{$codon},"\n";;
      }
  }

sub get_codon_usage_for_chromosome
  {
    my %args = @_;
    my $dsid = $args{dsid};
    my $chr = $args{chr};
    return unless $dsid;
    my $ds = $coge->resultset('Dataset')->find($dsid);
    my %codons;
    my ($code, $code_type);
    my $codon_total = 0;
    my %aa;
    my $aa_total=0;
    foreach my $feat ($ds->features({"feature_type.name"=>"CDS"},{join=>"feature_type"}))
      {
	($code, $code_type) = $feat->genetic_code() unless $code;
	my ($codon) = $feat->codon_frequency(counts=>1);
	grep {$codon_total+=$_} values %$codon;
	grep {$codons{$_}+=$codon->{$_}} keys %$codon;
	foreach (keys %$codon)
	  {
	    next if length ($_) eq 3;
	    print STDERR join ("\t", $_, $feat->id, $feat->names),"\n"
	  }
	foreach my $tri (keys %$code)
	  {
	    $aa{$code->{$tri}}+=$codon->{$tri};
	    $aa_total+=$codon->{$tri};
	  }
      }
#    foreach my $tri (keys %code)
#      {
#	$codons{$tri} = 0 unless $codons{$tri};
#      }
    my $count = 0;
    my $html = "Codon Usage: $code_type<table>";
    foreach (sort { substr($a, 0, 2) cmp substr ($b, 0, 2) || sort_nt($a) <=> sort_nt($b) } keys %$code)
      {
	my $str = "<tr><td>".$_."(".$code->{$_}.")<td>".$codons{$_}."<td>(".sprintf("%.2f",100*$codons{$_}/$codon_total)."%)";
	delete $codons{$_};
	$html .= "</table><tr><td><table>" unless $count;
	if ($count%4)
	  {
	  }
	else
	  {
	    $html .= "</table><td nospan><table>";
	  }

	$html .= $str;
	$count++;
	$count = 0 if $count == 16;
      }
     $count = 0;
     foreach (sort keys %codons)
       {
 	my $str = "<tr><td>".$_."<td>".$codons{$_}."<td>(".sprintf("%.2f",100*$codons{$_}/$codon_total)."%)";
	$html .= "</table><tr><td><table>" unless $count;
 	if ($count%4)
 	  {
 	  }
 	else
 	  {
 	    $html .= "</table><td nospan><table>";
 	  }
 	$html .= $str;
 	$count++;
 	$count = 0 if $count == 16;
      }
    $html .="</table></table>";
    $html .= "Predicted amino acid usage using $code_type<table>";
    $html .= join ("<tr>",map  {"<td>$_<td>".$aa{$_}."<td>(".sprintf("%.2f",100*$aa{$_}/$aa_total)."%)"} sort keys %aa);
    $html .= "</table>";
    return $html;
#     foreach (map {$_."(".$code->{$_}.") ".$codons{$_}} sort { substr($a, 0, 2) cmp substr ($b, 0, 2) || sort_nt($a) <=> sort_nt($b) }keys %codons)
#       {
# 	$html .= "<tr>" unless $count;
# 	$html .= "<td nospan>" unless $count%4;
# 	$html .= $_."<br>";
# 	$count++;
# 	$count = 0 if $count == 16;
#       }
#     return $html;
  }
