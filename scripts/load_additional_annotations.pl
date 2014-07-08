#!/usr/bin/perl -w

use strict;
use CoGe::Genome;
use Data::Dumper;

my $gdb = CoGe::Genome->new;
my $version = shift;
my $anno_type = $gdb->get_annotation_type_obj->find_or_create({name=>"Description"});

while (<>)
  {
    chomp;
    my @line = split /\t/;
    my $alias = process_alias($line[1], $line[-1]);
    foreach my $feat ($gdb->get_feat_by_name($line[0]))
      {
	if ($version)
	  {
	    next unless $feat->data_info->version == $version;
	  }
	my $feat_names = get_feat_names($feat);
#	print join ("\n", keys %$feat_names, "---",@$alias),"\n\n";
	foreach my $name (@$alias)
	  {
	    next if $feat_names->{$name};
	    print "Creating feature name $name for $line[0]\n";
#	    $gdb->get_feature_name_obj->create({name=>$name, feature_id=>$feat->id});
	  }
#	next if $line[2] eq "hypothetical protein" || $line[2] eq "expressed protein";
	$line[2] =~ s/^\s+//g;
#	print $line[2] if $line[2] =~ /^\w+\s\w+$/;
	my $annos = get_annos($feat);
	next if $annos->{$line[2]};
	print "Creating annotation:  $line[0]: $line[2]\n";
	$gdb->get_annotation_obj->create({
					  annotation=> $line[2],
					  feature_id=> $feat->id,
					  annotation_type_id=>$anno_type->id(),
					  });
      }
  }

sub get_annos
  {
    my $feat = shift;
    my %annos;
    foreach my $anno ($feat->annos())
      {
	$annos{$anno->anno} = 1;
      }
    return \%annos;
  }

sub get_feat_names
  {
    my $feat = shift;
    my %names;
    foreach my $name ($feat->names())
      {
	$names{$name->name} = 1;
      }
    return \%names;
  }

sub process_alias
  {
    my ($sym, $ali) = @_;
    my %alias;
    my ($tmp) = $sym =~ /Symbol:\s+(\S*)/;
    $alias{$tmp}=1 unless $tmp =~ /None/i;
    ($tmp) = $ali=~/Aliases:\s+(.*)/;
    print "Error with alias: $ali\n" unless $tmp;
    foreach (split /,\s*/, $tmp)
      {
	$alias{$_}=1 unless /None/i;
      }
    return [keys %alias];
  }
