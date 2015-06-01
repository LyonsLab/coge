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
	if (-d "$dir/$item")
	  {
	    process_dir("$dir/$item");
	  }
	elsif (-r "$dir/$item")
	  {
	    if ($item =~ /condensed/)
	      {
		process_file ("$dir/$item");
	      }
	  }
      }
  }

sub process_file
    {
      my $file = shift;
      print "deleting $file\n\n";
      `rm -f $file`;
#      `rm -f $dir/$file`;
#      `rm -f $dir/$tmp`;

    }
