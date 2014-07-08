#!/usr/bin/perl -w

use strict;
use Getopt::Long;

use vars qw($dir $DEBUG $GO);

GetOptions ( "debug=s" => \$DEBUG,
             "go=s"    => \$GO,
	     "d|dir=s" =>\$dir
	   );

$dir =  "." unless defined $dir;

print "GO flag is not set:  nothing will be deleted\n" unless $GO;

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
	    if ($item =~ /running$/)
	      {
		process_file (dir=>$dir, file=>$item);
	      }
	  }
      }
  }

sub process_file
    {
      my %opts = @_;
      my $dir = $opts{dir};
      my $file = $opts{file};
      my $tmp = $file;
      $tmp =~ s/\.running$//;
      print "deleting $dir/$file\n";
      print "deleting $dir/$tmp\n\n";
      `rm -f $dir/$file` if $GO;
      `rm -f $dir/$tmp` if $GO;

    }
