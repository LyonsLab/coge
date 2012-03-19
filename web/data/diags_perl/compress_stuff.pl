#!/usr/bin/perl -w

use strict;

my $dir = shift || ".";

process_dir($dir);

sub process_dir
  {
    my $dir = shift;
    my @items;
    opendir (DIR, $dir);
    while (my $item = readdir(DIR))
      {
	next if $item =~ /^\.\.?$/;
	push @items, $item;
      }
    foreach my $item (@items)
      {
	next if $item =~ /^\.\.?$/;
	next if $item =~ /^perl$/;
	if (-d "$dir/$item")
	  {
	    process_dir("$dir/$item");
	  }
	elsif (-r "$dir/$item")
	  {
	    if ($item !~ /gz$/)
	      {
		process_file ("$dir/$item");
	      }
	  }
      }
  }

sub process_file
    {
      my $file = shift;
      return if $file =~ /html/;
      return if $file =~ /sqlite/;
      print $file,"\n";
      `gzip -f $file`;
    }
