#!/usr/bin/perl

die "DIE $prog: Statistics::Descriptive module is missing! (try: perl -MCPAN -e 'install Statistics::Descriptive')" unless(eval{require Statistics::Descriptive});
die "DIE $prog: Statistics::R module is missing! (try: perl -MCPAN -e 'install Statistics::R')" unless(eval{require Statistics::R});
die "DIE $prog: Parallel::ForkManager module is missing!" unless(eval{require Parallel::ForkManager});
die "DIE $prog: Bio::DB::Sam module is missing!" unless(eval{require Bio::DB::Sam});

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Parallel::ForkManager;
use Hash::Merge;
use Hash::Merge::Simple;
use Storable;
use Cwd;
use Bio::DB::Sam;
use Statistics::Descriptive;
use Statistics::R;

my $prog = basename($0);
my $VERSION = '2.0';
my $lastmodif = '2016-04-04';
my $help;
my ($outside, $inside) = (100,100);
my $gffFile;
my $bamFile;
my $rmInclude;
my $outOfRange;
my $matrice;
my $windowSize = 10;
my $metType;
my $outputPrefix;
my $FEAT_TYPE;
my $CPUS = 1; # default to single processor
my $QUIET = 0;

&GetOptions(
    "mat"        => \$matrice,
    "quiet"      => \$QUIET,
    "outRange"   => \$outOfRange,
    "met:s"      => \$metType,
    "include"    => \$rmInclude,
    "outside:i"  => \$outside,
    "inside:i"   => \$inside,
    "w:i"        => \$windowSize,
    "gff:s"      => \$gffFile,
    "bam:s"      => \$bamFile,
    "featType:s" => \$FEAT_TYPE,
    "cpu:i"      => \$CPUS,
    "o:s"        => \$outputPrefix,
    "h|help"     => \$help
);

$| = 1;
$help and &help;
#@ARGV or &help;

&main(\@ARGV);

sub main {
	my $self = {};
	bless $self;
	my $datestring = localtime();
	print STDERR "INFO $prog: Running script $datestring\n";
	$self->setOptions();
	$self->{totalFeatures} = $self->readGffFile($gffFile);
#	$self->{inputFiles} = shift;
#	foreach (@{$self->{inputFiles}}) {
#		$self->{_currentFile} = $_;
        $self->{_currentFile} = $bamFile;
		unless (-e $self->{_currentFile}){
			print STDERR "INFO $prog: Could not find file: " . $self->{_currentFile} . " -> next file\n"; next;
		}
		$self->makeWindows();
	    $self->getCoverageFromBam2();
		$self->calcMeanAndCi();
#		delete($self->{windows});
#	}
	$self->print();
	$self->R(); 
	$datestring = localtime();
	print STDERR "INFO $prog: Finished succesfully $datestring\n";
	exit(0);
}

sub setOptions {
	my $self = shift;
    if (defined ($gffFile)){
        unless (-e $gffFile){ print STDERR "DIE $prog: Could not find file: " . $gffFile . "\n" and die; }
        $self->{option}->{gffFile} = $gffFile;
    }
	else{
		print STDERR "DIE $prog: Could not find -s and -e and -chr or -gff option.\n" and &help;
	}
	if(defined $rmInclude){$rmInclude = "YES";}
	else{$rmInclude = "NO";}
	my $remainder = ($inside+$outside) % $windowSize;
	my $quotient = ($inside+$outside) / $windowSize;
	if($remainder != 0){
		print STDERR "DIE $prog: Number of windows beetwen outside and inside equal to $quotient. Should be a whole number.\n" and &help;
	}
	if(not defined $outputPrefix){
		print STDERR "DIE $prog: option -o not defined. Need output file\n" and &help;
	}
	else{
		if($outputPrefix =~/.*/.*/){
			$self->{outputMatrice} = $outputPrefix.".tab";
			$self->{outputFigure} = $outputPrefix.".png"; #$outputPrefix.".pdf";
		}
		else{
			my $dir = getcwd;
			$self->{outputMatrice} = $dir."/".$outputPrefix.".tab";
			$self->{outputFigure}  = $dir."/".$outputPrefix.".png"; #$outputPrefix.".pdf";
		}
	}
}

sub readGffFile { # load a CoGe GFF file
    my $self = shift;
    my $file = shift;
    
    $FEAT_TYPE = lc($FEAT_TYPE) if $FEAT_TYPE;
    
    open(my $fh, $file) or 
        die("Error: cannot open file '$file'\n");
    
    my $count = 0;
    while (<$fh>) {
        next if /^#/;
        chomp;
        my ($chr, undef, $type, $start, $stop, undef, $strand, undef, $attr) = split /\t/;
        next unless (defined $chr && defined $type && defined $start && defined $stop && defined $attr);
        
        $type = lc($type);
        next if (defined $FEAT_TYPE && $type ne $FEAT_TYPE);
        
#        $chr = formatChr($chr);
        
        # Parse ID from attributes field
        my ($id) = $attr =~ /ID=([^;]+);*/i;
        print STDERR "readGffFile: $id $type\n" unless $QUIET;
        
        ( $start, $stop ) = ( $stop, $start ) if ( $start > $stop ); # reorder coordinates if needed
        $start += 1;    # put start in base 1     # mdb: what about stop?
        
        $strand = '+' unless ( defined $strand and ( $strand eq '+' or $strand eq '-' ) );
        
        $self->{features}->{$chr}->{$id} = {
            start  => $start, 
            stop   => $stop, 
            strand => $strand 
        };
        
        $count++;
    }
    close($fh);
    
   print STDERR "readGffFile: $count features loaded\n";
   return $count;
}

sub getCoverageFromBam {
	my $self = shift;
	my $count = 0;
	
	my $sam = Bio::DB::Sam->new(-bam => $self->{_currentFile});
	foreach my $chr (sort keys %{$self->{features}} ){
		foreach my $featureId (sort {$self->{features}->{$chr}->{$a}->{start} <=> $self->{features}->{$chr}->{$b}->{start} } keys %{$self->{features}->{$chr}} ) {
			print STDERR "\r$featureId ", sprintf("%.1f", 100*(++$count / $self->{totalFeatures})), '%                     ' unless $QUIET;
			
			my $feature = $self->{features}->{$chr}->{$featureId};
			
			# Check if feature is not included in area of interest
			next if ($rmInclude eq "YES" && ($feature->{start} + $inside) > $feature->{stop});
#			if($rmInclude eq "YES"){ # check if feature is not inclute in area of interest
#                if($self->{features}->{$chr}->{$featureId}->{start}+$inside>$self->{features}->{$chr}->{$featureId}->{stop}){ next ;}
#            }
			
			my %coverage;
			foreach my $typeCoord (sort ('start', 'stop')) {
				# Extract coverage from bam file
				my ($tmpStart, $tmpStop);
				if ($typeCoord eq "start") {
					$tmpStart = $feature->{start} - $outside;
					$tmpStop = $feature->{start} + $inside;
				}
				elsif ($typeCoord eq "stop") {
					$tmpStart = $feature->{stop} - $inside;
					$tmpStop = $feature->{stop} + $outside;
				}
				
				$sam->pileup($chr.':'.$tmpStart.'-'.$tmpStop,
					sub {
				         my ($seqid,$pos,$pileup) = @_;
#				         $self->{tmp}->{$pos} = @$pileup;
				         $coverage{$pos} = $pileup;
				    }
				 );

				# Add 0 and filtered out unwanted position
				my $k = 1;
				for(my $i = $tmpStart; $i < $tmpStop; $i++ ) {
					my ($stOrSt, $pos);
					if ($feature->{strand} eq "+"){
						if($typeCoord eq "start"){
							$stOrSt = "st";
							$pos = ($i-$feature->{start})+$outside+1;
						}
						elsif($typeCoord eq "stop"){
							$stOrSt = "sp";
							$pos = ($i-$feature->{stop})+$inside+1;
						}
					}
					elsif ($feature->{strand} eq "-"){
						if ($typeCoord eq "start"){
							$stOrSt = "sp";
							$pos = $inside+$outside-(($i-$feature->{start})+$outside+1)+1;
						}
						elsif ($typeCoord eq "stop"){
							$stOrSt = "st";
							$pos = $inside+$outside-(($i-$feature->{stop})+$inside);
						}
					}

					if (not defined $coverage{$i}) { #$self->{tmp}->{$i}) { # when cov is equal to 0
						if (defined $matrice) { $self->{mat}->{$stOrSt}->{$k++}->{$featureId} = 0; }
						next if (defined $outOfRange and $typeCoord eq "start" and $i > $feature->{stop} );
						next if (defined $outOfRange and $typeCoord eq "stop" and $i < $feature->{start} );
						push(@{$self->{windows}->{$self->{nc}->{$pos}}->{cov}->{$stOrSt}}, 0);
					}
					else {
						if (defined $matrice) { $self->{mat}->{$stOrSt}->{$k++}->{$featureId} = $coverage{$i}; }#$self->{tmp}->{$i};}
						next if (defined $outOfRange and $typeCoord eq "start" and $i > $feature->{stop} );
						next if (defined $outOfRange and $typeCoord eq "stop" and $i < $feature->{start} );
						push(@{$self->{windows}->{$self->{nc}->{$pos}}->{cov}->{$stOrSt}}, $coverage{$i});#$self->{tmp}->{$i});
					}
				}
#				delete($self->{tmp});
			}
		}
		print STDERR "\n";
		# print Dumper $self->{windows}; exit;
	}
	
	if (defined($matrice)) {
		print join("\t", keys %{$self->{mat}->{up}->{"1"}}), "\n";
		foreach my $startOrStop (sort keys %{$self->{cov}}) {
			foreach my $pos (sort {$a <=> $b} keys %{$self->{mat}->{$startOrStop}}) {
				my @print;
				push (@print, $pos);
				foreach my $featureId (sort keys %{$self->{mat}->{$startOrStop}->{$pos}}) {
					push (@print, $self->{mat}->{$startOrStop}->{$pos}->{$featureId});
				}
				print join("\t", @print), "\n";
			}
		}
	}
}

sub getCoverageFromBam2 {
    my $self = shift;

    print STDERR "Using $CPUS cpus\n" unless $QUIET;
    my $pm = new Parallel::ForkManager($CPUS-1);
    $pm->run_on_finish(
        sub { 
            my ($pid, $exit_code, $chr) = @_;
            
            my $tmpFile = $outputPrefix.'_'.$chr.'.tmp';
            next unless -e $tmpFile;
            
            print STDERR "Merging $chr\n" unless $QUIET;
            my $h = retrieve($tmpFile);
            Hash::Merge::set_behavior('RIGHT_PRECEDENT');
            $self->{windows} = Hash::Merge::merge($self->{windows}, $h);
        }
    );
    
    my @sortedChr = sort keys %{$self->{features}};
    
    foreach my $chr (@sortedChr) {
        print STDERR "Starting $chr\n" unless $QUIET;
        $pm->start($chr) and next; # fork
        
        my $sam = Bio::DB::Sam->new(-bam => $self->{_currentFile});
        
        foreach my $featureId (sort {$self->{features}->{$chr}->{$a}->{start} <=> $self->{features}->{$chr}->{$b}->{start} } keys %{$self->{features}->{$chr}} ) {
            #print STDERR "\r$chr $featureId" unless $QUIET;
            
            my $start  = $self->{features}->{$chr}->{$featureId}->{start};
            my $stop   = $self->{features}->{$chr}->{$featureId}->{stop};
            my $strand = $self->{features}->{$chr}->{$featureId}->{strand};
            
            next if ($rmInclude eq 'YES' && ($start+$inside)>$stop); # check if feature is not inclute in area of interest
            
            foreach my $typeCoord (sort ('start', 'stop') ) {
                # extract coverage from bam file
                my ($tmpStart, $tmpStop);
                if($typeCoord eq 'start'){
                    $tmpStart = $start - $outside;
                    $tmpStop = $start + $inside;
                }
                elsif($typeCoord eq 'stop'){
                    $tmpStart = $stop - $inside;
                    $tmpStop = $stop + $outside;
                }
                
                my %coverage;
                $sam->fast_pileup($chr.':'.$tmpStart.'-'.$tmpStop,
                    sub {
                         my (undef,$pos,$pileup) = @_;
                        $coverage{$pos} = @$pileup;
                    }
                );

                # add 0 and filtered out unwanted position
                my $k = 1 ;
                for(my $i = $tmpStart ; $i < $tmpStop ; $i++ ){
                    my ($stOrSt, $pos);
                    if ($strand eq '+') {
                        if($typeCoord eq 'start') {
                            $stOrSt = 'st' ;
                            $pos = ($i-$start)+$outside+1 ;
                        }
                        elsif($typeCoord eq 'stop') {
                            $stOrSt = 'sp' ;
                            $pos = ($i-$stop)+$inside+1 ;
                        }
                    }
                    elsif ($strand eq '-') {
                        if ($typeCoord eq 'start') {
                            $stOrSt = 'sp' ;
                            $pos = $inside+$outside-(($i-$start)+$outside+1)+1 ;
                        }
                        elsif ($typeCoord eq 'stop') {
                            $stOrSt = 'st' ;
                            $pos = $inside+$outside-(($i-$stop)+$inside) ;
                        }
                    }

                    if (not defined $coverage{$i}) { # when cov is equal to 0
                        if(defined $matrice){ $self->{mat}->{$stOrSt}->{$k++}->{$featureId} = 0 ; }
                        next if (defined $outOfRange and $typeCoord eq 'start' and $i>$stop );
                        next if (defined $outOfRange and $typeCoord eq 'stop' and $i<$start );
                        push(@{$self->{windows}->{$self->{nc}->{$pos}}->{cov}->{$stOrSt}}, 0) ;
                    }
                    else {
                        if(defined $matrice){ $self->{mat}->{$stOrSt}->{$k++}->{$featureId} = $coverage{$i} ;}
                        next if (defined $outOfRange and $typeCoord eq 'start' and $i>$stop );
                        next if (defined $outOfRange and $typeCoord eq 'stop' and $i<$start );
                        push(@{$self->{windows}->{$self->{nc}->{$pos}}->{cov}->{$stOrSt}}, $coverage{$i});
                    }
                }
            }
        }
        print STDERR "\n" unless $QUIET;
        
        print STDERR 'Finished ', $chr, ' ', scalar(keys %{$self->{windows}}), "\n" unless $QUIET;
        store $self->{windows}, $outputPrefix.'_'.$chr.'.tmp';
        $pm->finish; # exit the child process
        
        # print Dumper $self->{windows}; exit;
    }
    $pm->wait_all_children;
    
    if (defined($matrice)) {
        print join("\t", keys %{$self->{mat}->{up}->{"1"}}), "\n" ;
        foreach my $startOrStop (sort keys %{$self->{cov}}) {
            foreach my $pos (sort {$a <=> $b} keys %{$self->{mat}->{$startOrStop}}) {
                my @print ;
                push (@print, $pos) ;
                foreach my $featureId (sort keys %{$self->{mat}->{$startOrStop}->{$pos}}) {
                    push (@print, $self->{mat}->{$startOrStop}->{$pos}->{$featureId}) ;
                }
                print join("\t", @print), "\n" ;
            }
        }
    }
}

sub makeWindows {
	my $self = shift;
	my $k = 0;
	my $nc = 1;
	for(my $i = 1; $i<=$outside+$inside-$windowSize+1; $i=$i+$windowSize){
		$self->{windows}->{++$k}->{start} = $i;
		$self->{windows}->{$k}->{stop} = $i+$windowSize-1;
		for (my $j = 1; $j <= $windowSize; $j++ ){
			$self->{nc}->{$nc++} = $k;
		}
	}
	# print Dumper $self->{windows}; exit;
}

sub calcMeanAndCi {
	my $self = shift;
	foreach my $windK (sort {$a<=>$b} keys %{$self->{windows}}){
		foreach my $stOrSp (sort keys %{$self->{windows}->{$windK}->{cov}}){
			my $statObj = Statistics::Descriptive::Full->new();
		 	$statObj->add_data(@{$self->{windows}->{$windK}->{cov}->{$stOrSp}});
			my $count    = sprintf "%.8f", $statObj->count;
			my $sum      = sprintf "%.8f", $statObj->sum;
			my $median   = sprintf "%.8f", $statObj->median;
			my $stddev   = sprintf "%.8f", $statObj->standard_deviation;
			my $min      = sprintf "%.8f", $statObj->min;
			my $max      = sprintf "%.8f", $statObj->max;
			$self->{print}->{$stOrSp}->{$windK}->{basename($self->{_currentFile})}->{mean} = sprintf "%.8f", $statObj->mean;
			$self->{print}->{$stOrSp}->{$windK}->{basename($self->{_currentFile})}->{ci} = sprintf "%.8f", (1.960*($stddev/sqrt($count)));
			# print join("\t", ($startOrStop, $self->{windows}->{$upOrdown}->{$windK}->{start}, $self->{windows}->{$upOrdown}->{$windK}->{stop}, $count, $avg, $stddev, $ci)), "\n";
		}
	}
}

sub print {
	my $self = shift;
	open(my $out, '>', $self->{outputMatrice}) or die("Could not open file $self->{outputFile}.\n");

	my @header; 
	push(@header, "5or3", "pos", "start", "stop");
	my @fileName = sort keys $self->{print}->{st}->{'1'};
	foreach(@fileName){
		push(@header,$_);
		push(@header,$_);	
	}
	# print join ("\t", @header), "\n";
	foreach my $startOrStop (sort keys %{$self->{print}}){
		foreach my $windK (sort {$a <=> $b} keys %{$self->{print}->{$startOrStop}}){
			my $x;
			if($startOrStop eq "st"){
				$x = $windK-($outside/$windowSize)-1;
			}
			elsif($startOrStop eq "sp"){
				$x = $windK-($inside/$windowSize);
			}
			# my @print;
			# push(@print, $startOrStop, $x, $self->{windows}->{$windK}->{start}, $self->{windows}->{$windK}->{stop});
			foreach my $file (sort keys %{$self->{print}->{$startOrStop}->{$windK}}){
				print $out join("\t", $startOrStop, $file, $x, $self->{print}->{$startOrStop}->{$windK}->{$file}->{mean}, $self->{print}->{$startOrStop}->{$windK}->{$file}->{ci}),"\n";
				# push(@print, $self->{print}->{$startOrStop}->{$windK}->{$file}->{mean}, $self->{print}->{$startOrStop}->{$windK}->{$file}->{ci});
			}
			# print join("\t", @print), "\n";
		}
	}
}

sub R {
	my $self = shift;
	$self->{R} = Statistics::R->new();
	$self->{R}->startR();

	# install R lib if necessary
	$self->{R}->run(q`if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}`); # install package Deriv if necessary
	$self->{R}->run(q`library("RColorBrewer")`);
	$self->{R}->run(q`if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}`); # install package Deriv if necessary
	$self->{R}->run(q`library("ggplot2")`);
	$self->{R}->run(q`if("useful" %in% rownames(installed.packages()) == FALSE) {install.packages("useful")}`); # install package Deriv if necessary
	$self->{R}->run(q`library("useful")`);
	$self->{R}->run(q`if("grid" %in% rownames(installed.packages()) == FALSE) {install.packages("grid")}`); # install package Deriv if necessary
	$self->{R}->run(q`library("grid")`);
	$self->{R}->run(q`if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}`); # install package Deriv if necessary
	$self->{R}->run(q`library("gridExtra")`);

	$self->{R}->run(qq`mat = read.table("$self->{outputMatrice}", h=T, sep = "\t")`);
	$self->{R}->run(q`colnames(mat) = c("storsp", "g", "bin", "avg", "ci")`);
	$self->{R}->run(q`min = min(mat[,4]-mat[,5])`);
	$self->{R}->run(q`max = max(mat[,4]+mat[,5])`);

	$self->{R}->run(q`subMat5 = mat[which(mat[,1]=="st"),-1]`);
	$self->{R}->run(q`col <- brewer.pal(n = length(unique(subMat5[,1])), name = "Set1")`);
	$self->{R}->run(q`end5 = ggplot(subMat5, aes(x=bin, y=avg, fill=g, color=g)) + geom_ribbon(aes(ymin = avg - ci, ymax = avg + ci),alpha = 0.4) + geom_line() + scale_fill_manual(values=col) + scale_colour_manual(values = col ) + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y =element_blank()) + scale_y_continuous(limits = c(min, max)) + geom_vline(xintercept=c(0), linetype=2, colour="black")`);

	$self->{R}->run(q`subMat3 = mat[which(mat[,1]=="sp"),-1]`);
	$self->{R}->run(q`col <- brewer.pal(n = length(unique(subMat3[,1])), name = "Set1")`);
	$self->{R}->run(q`end3 = ggplot(subMat3, aes(x=bin, y=avg, fill=g, color=g)) + geom_ribbon(aes(ymin = avg - ci, ymax = avg + ci),alpha = 0.4) + geom_line() + scale_fill_manual(values=col) + scale_colour_manual(values = col ) + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y =element_blank()) + scale_y_continuous(limits = c(min, max)) + geom_vline(xintercept=c(0), linetype=2, colour="black")`);

#	$self->{R}->run(qq`pdf("$self->{outputFigure}")`);
	$self->{R}->run(qq`png(file="$self->{outputFigure}", width = 800, height = 800, units = "px", type = "cairo")`);
	$self->{R}->run(q`grid.newpage()`);
	$self->{R}->run(q`pushViewport(viewport(layout = grid.layout(1, 2)))`);
	$self->{R}->run(q`print(end5, vp = vplayout(1, 1))`);
	$self->{R}->run(q`print(end3, vp = vplayout(1, 2))`);
	$self->{R}->run(q`dev.off()`);
}

sub help {
my $prog = basename($0);
print STDERR <<EOF;
#### $prog ####
#
# CREATED:    2016-04-01
# LAST MODIF: $lastmodif
# AUTHOR:     Josquin Daron (Ohio State University, Slotkin Lab)
# VERSION:    1.1
#
# This script is used to make metaplot.
# $prog -bed <bed file> <option> file.bam or file.txt

USAGE:

       ### OPTIONS ###

       -h|--help             print this help
       -o         <file>     output file (without extension)
       -bed       <file>     bed file of the list of positions
       -outside   <int>      number of nucleotide outside of the start/stop (defaut 100)
       -inside    <int>      number of nucleotide inside of the start/stop (defaut 100)
       -w         <int>      windows size (defaut 10)
       -met       <string>   methylation context (CG, CHG or CHH)
       -outRange             if features length < <int> inside length, don't report coverage 
                             out of the range of a feature
       -include              remove feature include in the window of interest

EOF
exit(1);
}








