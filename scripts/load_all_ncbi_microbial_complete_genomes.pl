#!/usr/bin/perl -w

use strict;
use LWP::Simple;
use Data::Dumper;
use Getopt::Long;

$| =1;

my ($DEBUG, $forks, $pause, $autoupdate);

GetOptions (
            "debug"=>\$DEBUG,
            "forks|f=i"=>\$forks,
            "pause|p"=>\$pause,
            "autoupdate"=>\$autoupdate,
            );
my @skipped; #array to hold commands that had skipped reloading.  Will be printed out at end of program for user to manually go through;

my $data = get_accns_from_ncbi();

if (@skipped)
  {
    print "Commands that had automatially skipped datasets:\n";
    print join ("\n", @skipped),"\n";
  }

sub get_accns_from_ncbi
  {
    my $ftp = "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/lproks_0.txt";
#    my $page = get($ftp);
#    foreach my $item ($page=~/(lproks.*?.txt)/g)
#      {
#	print $item,"\n";
	my $list = get($ftp);
	#@[20] genbank accn
	#@[21] refseq accn
	foreach my $row (split/\n/,$list)
	  {
	    next if $row =~ /^#/;
	    my @line = split /\t/, $row;
#	    print join ("\t", $line[21], $line[22]),"\n";
	    my $accns = $line[22] ? $line[22] : $line[21];
	    next if $accns eq "-";
	    print "!!",$row,"\n" unless $accns;
	    my @accns = split/,/,$accns;
	    if (scalar @accns == 2 && ($accns[0] =~ /_\w\w\d/ || $accns[0] =~ /^\w\w\d/) && $accns[1] =~ /\D\D\D\D\d/)
	      {
		pop @accns;
	      }
#	    print join ("\t", @accns),"\n";
	    print "\n",$line[3],"\n";
	    run_genbank_genome_loader(accns=>\@accns);
	  }
#      }
  }

sub run_genbank_genome_loader
  {
    my %opts = @_;
    my $accns = $opts{accns};
    my $prog = '/home/elyons/projects/CoGeX/scripts/genbank_genome_loader.pl';
    return unless @$accns;
    $prog .= " -autoupdate" if $autoupdate;
    my $run = $prog." -accn ".join (" -accn ", @$accns);

    $run .= " -td '/tmp/gb/'";
    $run .= " -go";
    $run .= " -autoskip";
    $run .= " -delete_src_file";
    print "\n",$run,"\n";
    my $previously_loaded = 0;
    my $skipped = 0;
    #       next;
    open (IN, $run." |");
    while (<IN>)
      {
	$previously_loaded =1 if /previously loaded/;
	$skipped=1 if /skipping/;
	print $_;
      }
    close IN;
    push @skipped, $run if $skipped;
    if ($pause && !$previously_loaded)
      {
	print "Paused:  press enter to continue\n";
	<STDIN>;
      }
  }
