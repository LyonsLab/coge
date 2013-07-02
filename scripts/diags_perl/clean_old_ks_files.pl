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
	    if ($item =~ /html/)
	      {
#		print "deleting directory: $dir/$item\n";
#		`rm -rf $dir/$item`;
	      }
	    else
	      {
		process_dir("$dir/$item");
	      }
	  }
	elsif (-r "$dir/$item")
	  {
	    process_file (dir=>$dir, file=>$item);
	  }
      }
  }

sub process_file
    {
      my %opts = @_;
      my $dir = $opts{dir};
      my $file = $opts{file};
      if ($file =~ /\.ks$/)
	{
	  print "deleting $dir/$file\n";
	  `rm -f $dir/$file`;
	}
    }
