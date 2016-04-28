#!/usr/bin/perl

die "DIE $prog: Statistics::Descriptive module is missing! (try: perl -MCPAN -e 'install Statistics::Descriptive')" unless(eval{require Statistics::Descriptive});
die "DIE $prog: Tie::File module is missing! (try: perl -MCPAN -e 'install Tie::File')" unless(eval{require Tie::File});
die "DIE $prog: Statistics::R module is missing! (try: perl -MCPAN -e 'install Statistics::R')" unless(eval{require Statistics::R});
die "DIE $prog: Bio::DB::Sam module is missing!" unless(eval{require Bio::DB::Sam});

use strict ;
use warnings ;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Cwd;
use Bio::DB::Sam;
use Statistics::Descriptive;
use Tie::File;
use Statistics::R ;


my $prog = basename($0) ;
my $VERSION = '2.0';
my $lastmodif = '2016-04-04';
my $help;
my ($outside, $inside) = (100,100) ;
my $stdin ;
my $fastaId ;
my $bedFile ;
my $gffFile;
my $statCalc ;
my $prefix = '';
my $rmInclude ;
my $outOfRange ;
my $matrice ;
my $windowSize = 10;
my $metType ;
my $outputPrefix ;

&GetOptions(
			"mat"                   => \$matrice,
			"outRange"              => \$outOfRange,
			"met:s"                 => \$metType,
			"include"               => \$rmInclude,
			"outside:i"             => \$outside,
			"inside:i"              => \$inside,
			"w:i"                   => \$windowSize,
			"bed:s"                 => \$bedFile,
			"gff:s"                 => \$gffFile,
			"o:s"                   => \$outputPrefix,
			"h|help"                => \$help);


$help and &help;
@ARGV or &help;

&main(\@ARGV);

sub main {
	my $self = {};
	bless $self;
	my $datestring = localtime();
	print STDERR "INFO $prog: Running script $datestring\n" ;
	$self->setOptions();
	#$self->readBedFile();
	$self->readGffFile($gffFile);
	$self->{inputFiles} = shift;
	foreach (@{$self->{inputFiles}}) {
		$self->{_currentFile} = $_;
		unless (-e $self->{_currentFile}){
			print STDERR "INFO $prog: Could not find file: " . $self->{_currentFile} . " -> next file\n"; next;
		}
		$self->makeWindows();
		if(defined $metType){
			$self->getCoverageFromMet();
			$self->calcMetPercentage();
		}
		else{
			$self->getCoverageFromBam();
		}
		$self->calcMeanAndCi();
		delete($self->{windows});
	} # END OF foreach
	$self->print();
	$self->R() ; 
	$datestring = localtime();
	print STDERR "INFO $prog: Finishied succesfully $datestring\n" ;
	exit(0);
} # END OF main

sub setOptions {
	my $self = shift;
	if (defined ($bedFile)){
		unless (-e $bedFile){ print STDERR "DIE $prog: Could not find file: " . $bedFile . "\n" and die ; }
		$self->{option}->{bedFile} = $bedFile ;
	}
    elsif (defined ($gffFile)){
        unless (-e $gffFile){ print STDERR "DIE $prog: Could not find file: " . $gffFile . "\n" and die ; }
        $self->{option}->{gffFile} = $gffFile ;
    }	
	else{
		print STDERR "DIE $prog: Could not find -s and -e and -chr or -bed option.\n" and &help ;
	}
	if(defined $rmInclude){$rmInclude = "YES" ;}
	else{$rmInclude = "NO" ;}
	my $remainder = ($inside+$outside) % $windowSize ;
	my $quotient = ($inside+$outside) / $windowSize ;
	if($remainder != 0){
		print STDERR "DIE $prog: Number of windows beetwen outside and inside equal to $quotient. Should be a whole number.\n" and &help ;
	}
	if(not defined $outputPrefix){
		print STDERR "DIE $prog: option -o not defined. Need output file\n" and &help ;
	}
	else{
		if($outputPrefix =~/.*/.*/){
			$self->{outputMatrice} = $outputPrefix.".tab" ;
			$self->{outputFigure} = $outputPrefix.".pdf" ;
		}
		else{
			my $dir = getcwd;
			$self->{outputMatrice} = $dir."/".$outputPrefix.".tab" ;
			$self->{outputFigure}  = $dir."/".$outputPrefix.".pdf" ;
		}
	}
}

sub readBedFile{
	my $self = shift;
	open("BED",$self->{option}->{bedFile}) or print STDERR "DIE $prog: Could not open file: " . $self->{option}->{bedFile} . " -> skip\n" and die ;
	while(<BED>){
		chomp $_ ;
		my @line = split("\t", $_) ;
		if(scalar(@line)<4){ print STDERR "ERROR $prog: Bed file must have 4 column in line: $_ \n" ; die ; }
		if($line[1]!~/^\d+$/ or $line[2]!~/^\d+$/){ print STDERR "ERROR $prog: Warning in bed file: start or stop should be a number in line: $_ \n" ; next ; }
		if(defined $line[5] and ($line[5] eq "+" or $line[5] eq "-")){
			$self->{bedFile}->{$line[0]}->{$line[3]}->{coord}->{start} = $line[1]+1 ; # put start in base 1
			$self->{bedFile}->{$line[0]}->{$line[3]}->{coord}->{stop} = $line[2] ;
			$self->{bedFile}->{$line[0]}->{$line[3]}->{strand} = $line[5] ;
		}
		else{
			$self->{bedFile}->{$line[0]}->{$line[3]}->{coord}->{start} = $line[1]+1 ; # put start in base 1
			$self->{bedFile}->{$line[0]}->{$line[3]}->{coord}->{stop} = $line[2] ;
			$self->{bedFile}->{$line[0]}->{$line[3]}->{strand} = "+" ;
		}
	}
	if (defined($self->{option}->{tabFile})){print STDERR "INFO $prog: read tab File $self->{option}->{tabFile} done !\n";}
}

sub readGffFile { # load a CoGe GFF file
    my $self = shift;
    my $file = shift;
    
    open(my $fh, $file) or 
        die("Error: cannot open file '$file'\n");
    
    my $count = 0;
    while (<$fh>) {
        next if /^#/;
        chomp;
        my ($chr, undef, $type, $start, $stop, undef, $strand, undef, $attr) = split /\t/;
        next unless (defined $chr && defined $type && defined $start && defined $stop && defined $attr);
        $type = lc($type);
        ($start, $stop) = ($stop, $start) if ($stop < $start); # reorder coordinates if needed

#        $chr = formatChr($chr);
        
        # Parse ID from attributes field
        my ($id) = $attr =~ /(ID)=([^;]+);*/i;
        
        ( $start, $stop ) = ( $stop, $start ) if ( $start > $stop );
        $start += 1;    # put start in base 1 # mdb: what about stop?
        
        $strand = '+' unless ( defined $strand and ( $strand eq '+' or $strand eq '-' ) );
        
        $self->{bedFile}->{$chr}->{$id}->{coord} = { start => $start, stop => $stop };
        $self->{bedFile}->{$chr}->{$id}->{strand} = $strand;
        
        $count++;
    }
    close($fh);
    
   print STDERR "readGffFile: $count features loaded\n";
}

sub getCoverageFromBam{
	my $self = shift;
	my $sam = Bio::DB::Sam->new(-bam => $self->{_currentFile});
	foreach my $chr (sort keys %{$self->{bedFile}} ){
		foreach my $featureId (sort {$self->{bedFile}->{$chr}->{$a}->{coord}->{start} <=> $self->{bedFile}->{$chr}->{$b}->{coord}->{start} } keys %{$self->{bedFile}->{$chr}} ){
			print STDERR "\r$featureId" ;
			if($rmInclude eq "YES"){ # check if feature is not inclute in area of interest
				if($self->{bedFile}->{$chr}->{$featureId}->{coord}->{start}+$inside>$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop}){ next ;}
			}
			foreach my $typeCoord (sort keys %{$self->{bedFile}->{$chr}->{$featureId}->{coord}} ){

				# extract coverage from bam file
				my $tmpStart ;
				my $tmpStop ;
				if($typeCoord eq "start"){
					$tmpStart = $self->{bedFile}->{$chr}->{$featureId}->{coord}->{$typeCoord}-$outside ;
					$tmpStop = $self->{bedFile}->{$chr}->{$featureId}->{coord}->{$typeCoord}+$inside ;
				}
				elsif($typeCoord eq "stop"){
					$tmpStart = $self->{bedFile}->{$chr}->{$featureId}->{coord}->{$typeCoord}-$inside ;
					$tmpStop = $self->{bedFile}->{$chr}->{$featureId}->{coord}->{$typeCoord}+$outside ;
				}
				$sam->pileup($chr.":".$tmpStart."-".$tmpStop,
					sub {
				         my ($seqid,$pos,$pileup) = @_;
				         $self->{tmp}->{$pos} = @$pileup ;
				 });

				# add 0 and filtered out unwanted position
				my $k = 1 ;
				for(my $i = $tmpStart ; $i < $tmpStop ; $i++ ){
					my $stOrSt ;
					my $pos ;
					if($self->{bedFile}->{$chr}->{$featureId}->{strand} eq "+"){
						if($typeCoord eq "start"){
							$stOrSt = "st" ;
							$pos = ($i-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start})+$outside+1 ;
						}
						elsif($typeCoord eq "stop"){
							$stOrSt = "sp" ;
							$pos = ($i-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop})+$inside+1 ;
						}
					}
					elsif($self->{bedFile}->{$chr}->{$featureId}->{strand} eq "-"){
						if($typeCoord eq "start"){
							$stOrSt = "sp" ;
							$pos = $inside+$outside-(($i-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start})+$outside+1)+1 ;
						}
						elsif($typeCoord eq "stop"){
							$stOrSt = "st" ;
							$pos = $inside+$outside-(($i-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop})+$inside) ;
						}
					}

					if(not defined $self->{tmp}->{$i}){ # when cov is equal to 0
						if(defined $matrice){ $self->{mat}->{$stOrSt}->{$k++}->{$featureId} = 0 ; }
						if(defined $outOfRange and $typeCoord eq "start" and $i>$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop} ){ next ;}
						if(defined $outOfRange and $typeCoord eq "stop" and $i<$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start} ){ next ;}
						push(@{$self->{windows}->{$self->{nc}->{$pos}}->{cov}->{$stOrSt}}, 0) ;
					}
					else{
						if(defined $matrice){ $self->{mat}->{$stOrSt}->{$k++}->{$featureId} = $self->{tmp}->{$i} ;}
						if(defined $outOfRange and $typeCoord eq "start" and $i>$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop} ){ next ;}
						if(defined $outOfRange and $typeCoord eq "stop" and $i<$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start} ){ next ;}
						push(@{$self->{windows}->{$self->{nc}->{$pos}}->{cov}->{$stOrSt}}, $self->{tmp}->{$i}) ;
					}
				}
				delete($self->{tmp}) ;
			}
		}
		print STDERR "\n" ;
		# print Dumper $self->{windows} ; exit ;
	}
	if(defined($matrice)){
		print join("\t", keys %{$self->{mat}->{up}->{"1"}}), "\n" ;
		foreach my $startOrStop (sort keys %{$self->{cov}}){
			foreach my $pos (sort {$a <=> $b} keys %{$self->{mat}->{$startOrStop}}){
				my @print ;
				push (@print, $pos) ;
				foreach my $featureId (sort keys %{$self->{mat}->{$startOrStop}->{$pos}}){
					push (@print, $self->{mat}->{$startOrStop}->{$pos}->{$featureId}) ;
				}
				print join("\t", @print), "\n" ;
			}
		}
	}
}

sub getCoverageFromMet{
	my $self = shift;
	tie my @array, 'Tie::File', $self->{_currentFile} or die("Unable to open file \"$self->{_currentFile}\": $!\n");
	print STDERR "INFO $prog: Reading file $self->{_currentFile} ... done !\n" ;
	my $startFile = 1 ;
	my $endFile = scalar(@array)-1 ;
	# my $endFile = 100011 ;
	my $rememberStartLine = 1;
	foreach my $chr (sort keys %{$self->{bedFile}}){
		my ($number) = $chr =~/Chr([0-9])/ ;
		foreach my $featureId (sort {$self->{bedFile}->{$chr}->{$a}->{coord}->{start} <=> $self->{bedFile}->{$chr}->{$b}->{coord}->{start} } keys %{$self->{bedFile}->{$chr}} ){
			if($self->{bedFile}->{$chr}->{$featureId}->{coord}->{start}-$outside<=0){ next ;} # if start - outside is negative
			# check if feature is not inclute in area
			if($rmInclude eq "YES"){
				if($self->{bedFile}->{$chr}->{$featureId}->{coord}->{start}+$inside>$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop}){ next ;}
			}
			my $test = "TRUE" ;
			for (my $i = $startFile ; $i < $endFile ; $i++){ # going through the txt file
				print STDERR "\r$featureId" ;
				my @line = split("\t", $array[$i]) ;
				if($line[0] < $number){ next ;}
				elsif($line[0] == $number and $line[1]>$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop}+$outside){
					$startFile = $rememberStartLine ;
					last ;
				}
				elsif($line[0] > $number){
					$startFile = $rememberStartLine ;
					last ;}
				#### surrounding start + and stop -
				if($line[1] >= $self->{bedFile}->{$chr}->{$featureId}->{coord}->{start}-$outside and $line[1] < $self->{bedFile}->{$chr}->{$featureId}->{coord}->{start}+$inside and $line[7] eq $metType){
					if($test eq "TRUE"){
						$rememberStartLine = $i ;
						$test = "FALSE" ;
					}
					my $pos ;
					my $stOrSp ;
					if(defined $outOfRange and $line[1]>$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop} ){ goto(HERE) ;} # need to use a goto because area around 5' and 3' could be overlapping
					if ($self->{bedFile}->{$chr}->{$featureId}->{strand} eq "+"){
						$pos = ($line[1]-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start})+$outside+1 ;
						$stOrSp = "st" ;
					}
					elsif ($self->{bedFile}->{$chr}->{$featureId}->{strand} eq "-"){
						$pos = $inside+$outside-(($line[1]-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start})+$outside+1)+1 ;
						$stOrSp = "sp" ;
					}
					push (@{$self->{windows}->{$self->{nc}->{$pos}}->{met}->{$stOrSp}->{$featureId}->{nbCm}}, $line[4]) ;
					push (@{$self->{windows}->{$self->{nc}->{$pos}}->{met}->{$stOrSp}->{$featureId}->{nbR}}, $line[5]) ;
				}
				
				#### surrounding stop + and start -
				HERE:if($line[1] >= $self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop}-$inside and $line[1] < $self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop}+$outside and $line[7] eq $metType){
					my $pos ;
					my $stOrSp ;
					if(defined $outOfRange and $line[1]<$self->{bedFile}->{$chr}->{$featureId}->{coord}->{start} ){ next ;}
					if ($self->{bedFile}->{$chr}->{$featureId}->{strand} eq "+"){
						$pos = ($line[1]-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop})+$inside+1 ;
						$stOrSp = "sp" ;
					}
					if ($self->{bedFile}->{$chr}->{$featureId}->{strand} eq "-"){
						$pos = $inside+$outside-(($line[1]-$self->{bedFile}->{$chr}->{$featureId}->{coord}->{stop})+$inside) ;
						$stOrSp = "st" ;
					}
					push (@{$self->{windows}->{$self->{nc}->{$pos}}->{met}->{$stOrSp}->{$featureId}->{nbCm}}, $line[4]) ;
					push (@{$self->{windows}->{$self->{nc}->{$pos}}->{met}->{$stOrSp}->{$featureId}->{nbR}}, $line[5]) ;
				}
			} # end of For line in map file
		}
	}
	print STDERR "\n" ;
	# print Dumper $self->{windows} ; exit ;
}

sub calcMetPercentage {
	my $self = shift;
	foreach my $windK (sort {$a <=> $b} keys %{$self->{windows}} ){
		if(not defined $self->{windows}->{$windK}->{met}){
			print STDERR "ERROR $prog: Not Cytosine in a given window: $prog can't calculate the % of cytosine metylated foreach windows (increase the window size).\n" and die ;
		}
		foreach my $storsp (sort keys $self->{windows}->{$windK}->{met}){
			foreach my $featureId (sort keys $self->{windows}->{$windK}->{met}->{$storsp}){
				my $statNbCm = Statistics::Descriptive::Full->new();
			 	$statNbCm->add_data(@{$self->{windows}->{$windK}->{met}->{$storsp}->{$featureId}->{nbCm}});
				my $statNbR = Statistics::Descriptive::Full->new();
			 	$statNbR->add_data(@{$self->{windows}->{$windK}->{met}->{$storsp}->{$featureId}->{nbR}});
			 	my $percMet = ($statNbCm->sum/$statNbR->sum)*100 ;
			 	push(@{$self->{windows}->{$windK}->{cov}->{$storsp}}, $percMet) ;
			}
		}
	}
	# print Dumper $self->{windows} ;
}

sub makeWindows{
	my $self = shift;
	my $k = 0 ;
	my $nc = 1 ;
	for(my $i = 1 ; $i<=$outside+$inside-$windowSize+1 ; $i=$i+$windowSize){
		$self->{windows}->{++$k}->{start} = $i ;
		$self->{windows}->{$k}->{stop} = $i+$windowSize-1 ;
		for (my $j = 1 ; $j <= $windowSize ; $j++ ){
			$self->{nc}->{$nc++} = $k ;
		}
	}
	# print Dumper $self->{windows} ; exit ;
}

sub calcMeanAndCi{
	my $self = shift;
	foreach my $windK (sort {$a<=>$b} keys %{$self->{windows}}){
		foreach my $stOrSp (sort keys %{$self->{windows}->{$windK}->{cov}}){
			my $statObj = Statistics::Descriptive::Full->new();
		 	$statObj->add_data(@{$self->{windows}->{$windK}->{cov}->{$stOrSp}});
			my $count    = sprintf "%.8f", $statObj->count ;
			my $sum      = sprintf "%.8f", $statObj->sum ;
			my $median   = sprintf "%.8f", $statObj->median ;
			my $stddev   = sprintf "%.8f", $statObj->standard_deviation ;
			my $min      = sprintf "%.8f", $statObj->min ;
			my $max      = sprintf "%.8f", $statObj->max ;
			$self->{print}->{$stOrSp}->{$windK}->{basename($self->{_currentFile})}->{mean} = sprintf "%.8f", $statObj->mean ;
			$self->{print}->{$stOrSp}->{$windK}->{basename($self->{_currentFile})}->{ci} = sprintf "%.8f", (1.960*($stddev/sqrt($count))) ;
			# print join("\t", ($startOrStop, $self->{windows}->{$upOrdown}->{$windK}->{start}, $self->{windows}->{$upOrdown}->{$windK}->{stop}, $count, $avg, $stddev, $ci)), "\n";
		}
	}
}

sub print {
	my $self = shift;
	open(my $out, '>', $self->{outputMatrice}) or die("Could not open file $self->{outputFile}.\n");

	my @header ; 
	push(@header, "5or3", "pos", "start", "stop") ;
	my @fileName = sort keys $self->{print}->{st}->{'1'} ;
	foreach(@fileName){
		push(@header,$_) ;
		push(@header,$_) ;	
	}
	# print join ("\t", @header), "\n" ;
	foreach my $startOrStop (sort keys %{$self->{print}}){
		foreach my $windK (sort {$a <=> $b} keys %{$self->{print}->{$startOrStop}}){
			my $x ;
			if($startOrStop eq "st"){
				$x = $windK-($outside/$windowSize)-1 ;
			}
			elsif($startOrStop eq "sp"){
				$x = $windK-($inside/$windowSize) ;
			}
			# my @print ;
			# push(@print, $startOrStop, $x, $self->{windows}->{$windK}->{start}, $self->{windows}->{$windK}->{stop}) ;
			foreach my $file (sort keys %{$self->{print}->{$startOrStop}->{$windK}}){
				print $out join("\t", $startOrStop, $file, $x, $self->{print}->{$startOrStop}->{$windK}->{$file}->{mean}, $self->{print}->{$startOrStop}->{$windK}->{$file}->{ci}),"\n" ;
				# push(@print, $self->{print}->{$startOrStop}->{$windK}->{$file}->{mean}, $self->{print}->{$startOrStop}->{$windK}->{$file}->{ci}) ;
			}
			# print join("\t", @print), "\n" ;
		}
	}
}

sub R{
	my $self = shift;
	$self->{R} = Statistics::R->new() ;
	$self->{R}->startR() ;

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

	$self->{R}->run(qq`pdf("$self->{outputFigure}")`);
	$self->{R}->run(q`grid.newpage()`);
	$self->{R}->run(q`pushViewport(viewport(layout = grid.layout(1, 2)))`);
	$self->{R}->run(q`print(end5, vp = vplayout(1, 1))`);
	$self->{R}->run(q`print(end3, vp = vplayout(1, 2))`);
	$self->{R}->run(q`dev.off()`);
}

sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
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
exit(1) ;
}








