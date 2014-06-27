#!/usr/bin/perl -w

use strict;

my $dir = shift;
process_dir($dir);

sub process_dir
  {
    my $dir = shift;
    my @dirs;
    opendir (DIR, $dir);
    while (my $item = readdir(DIR))
      {
	next if $item =~ /^\.\.?$/;
	if ($item =~ /faa$/)
	  {
	    print "processing $dir/$item\n";
	    process_file (dir=>$dir, file=>$item);
	  }
	elsif (-d "$dir/$item")
	  {
	    push @dirs, "$dir/$item";
	  }
      }
    closedir DIR;

    foreach my $item (@dirs)
      {
	process_dir("$item");
      }
  }

sub process_file
    {
      my %opts = @_;
      my $dir = $opts{dir};
      my $file = $opts{file};
      mkdir "$dir/chr" unless -d "$dir/chr";
      $/ = "\n>";
      open (IN, "$dir/$file");
      while (<IN>)
	{
	  s/>//g;
	  my ($name, $seq) = split/\n/,$_,2;
	  $name =~ s/^.*?\|//;
	  $seq =~ s/\n//g;
	  print $name," ",length($seq),"\n";
	  open (OUT, ">$dir/chr/$name");
	  print OUT $seq;
	  close OUT;
	}
      close IN;
      $/="\n";
    }
