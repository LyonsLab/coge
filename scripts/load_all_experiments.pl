#!/usr/bin/perl -w

use strict;
use CoGeX;
use CoGe::Accessory::Web qw(get_defaults);
use File::Path qw(make_path);
use Data::Dumper;

# Parameters
my $DEBUG = 0;
my $user          = 'mbomhoff';
my $source        = 'Brian Gregory';
my $version       = '1';
my $coge_dir      = '/home/mbomhoff/dev/coge';
my $install_dir   = '/storage/coge/data/experiments';
# maize:
#my $gid           = 8062;
#my $data_dir      = '/home/mbomhoff/tmp/fanli/coge_data/GSE39232';
#my $staging_dir   = '/home/mbomhoff/tmp/fanli/staging/GSE39232'; #"/home/elyons/projects/EPIC-coge/Jacobsen_2013/stagging";
#my $metadata_file = '/home/mbomhoff/tmp/fanli/coge_data/GSE39232/metadata_GSE39232.txt'; #"/home/elyons/projects/EPIC-coge/Jacobsen_2013/GSE39901_metadata.txt";
#my $citation_link = "http://www.ncbi.nlm.nih.gov/pubmed/23739895";
# soybean:
my $gid           = 1611;
my $data_dir      = '/home/mbomhoff/data/fanli/coge_data/GSE41753';
my $staging_dir   = '/home/mbomhoff/data/fanli/staging/GSE41753';
my $metadata_file = '/home/mbomhoff/data/fanli/coge_data/GSE41753/metadata_GSE41753.txt';
my $citation_link = "http://www.ncbi.nlm.nih.gov/pubmed/23739894";

# Load config file
my $config = "$coge_dir/web/coge.conf";
my $P      = CoGe::Accessory::Web::get_defaults($config);

# Connect to database
my $connstr = "dbi:mysql:dbname=".$P->{DBNAME}.";host=".$P->{DBHOST}.";port=".$P->{DBPORT}.";";
my $coge = CoGeX->connect( $connstr, $P->{DBUSER}, $P->{DBPASS} );

# Load metadata file
my $metadata = get_metadata($metadata_file);
unless ($metadata) {
	print STDERR "load_all_experiments: no metadata\n";
	exit(-1);
}
#print STDERR Dumper $metadata;

# Load each experiment file
my $exp_count = 0;
process_dir( dir => $data_dir, metadata => $metadata );

# Yay!
print STDERR "All done!\n";
exit;

#-------------------------------------------------------------------------------
sub get_metadata {
    my $file = shift;
    open( IN, $file );
    my %data;

	#Name   Mark    Genome  Protocol        Platform ID     Platform name   Citation        Processing      Raw files       Processed files Series  Samples
    while (<IN>) {
        chomp;
        next unless $_;
        next if /^#/;
        my @line = split(/\t/);
        my $processed_file = $line[9];
        $data{ $processed_file } = \@line;
    }
    close IN;
    return \%data;
}

sub process_dir {
    my %opts = @_;
    my $dir  = $opts{dir};
    opendir( DIR, $dir );
    while ( my $item = readdir(DIR) ) {
        if ( $item =~ /\.csv/ && -r "$dir/$item" ) {
            unless ( $metadata->{$item} ) {
                print "WARNING: no metadata for $item\n";
                next;
            }
            process_file(
                file     => "$dir/$item",
                metadata => $metadata->{$item}
            );
        }
    }
    closedir DIR;
}

sub process_file {
    my %opts = @_;
    my $file = $opts{file};
    my $md   = $opts{metadata};
    
    my (
        $name,        $mark,            $genome,   $protocol,
        $platform_id, $platform_name,   $citation, $processing,
        $raw_files,   $processed_files, $series,   $samples
    ) = @$md;

	$exp_count++;
	my $staging = "$staging_dir/$exp_count";
    if (-e $staging) {
	    print "Clearing staging directory: $staging\n";
	    `rm $staging/*`;
    }
    else {
    	make_path($staging);
    }

    my $cmd = "$coge_dir/scripts/load_experiment.pl " .
    	"-config $config -user_name $user -name '$name' -desc '$mark' " .
    	"-version '$version' -gid $gid -source_name '$source' " .
    	"-staging_dir $staging -install_dir $install_dir -data_file $file";
    print "Running: " . $cmd, "\n";
	return if ($DEBUG);
    
    my $output = qx{ $cmd };
    if ( $? != 0 ) {
        print STDERR "load_experiment.pl failed with rc=$?\n";
        return $?;
    }
    
    my ($eid) = $output =~ /experiment id: (\d+)/;
    if (!$eid) {
    	print STDERR "Unable to retrieve experiment id\n";
    	exit(-1);
    }
    else {
        print "Captured experimemt id $eid\n";
        my $exp = $coge->resultset('Experiment')->find($eid);
        my $exp_type =
          $coge->resultset('ExperimentType')
          ->find_or_create( { name => $protocol } );
        $coge->resultset('ExperimentTypeConnector')->create(
            {
                experiment_type_id => $exp_type->id,
                experiment_id      => $exp->id
            }
        );
        $exp_type =
          $coge->resultset('ExperimentType')
          ->find_or_create( { name => $platform_name } );
        $coge->resultset('ExperimentTypeConnector')->create(
            {
                experiment_type_id => $exp_type->id,
                experiment_id      => $exp->id
            }
        );
        my $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Citation" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $citation,
                link               => $citation_link
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Mark" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $mark,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Platform" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $platform_name,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Platform ID" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $platform_id,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Processed files" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $processed_files,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Processing" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $processing,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Protocol" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $protocol,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Series" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $series,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Samples" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $samples,
            }
        );
        $etype =
          $coge->resultset('AnnotationType')
          ->find_or_create( { name => "Gene mutant" } );
        $exp->add_to_experiment_annotations(
            {
                annotation_type_id => $etype->id,
                annotation         => $genome,
            }
        );
    }
}
