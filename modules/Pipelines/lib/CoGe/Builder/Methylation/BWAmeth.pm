package CoGe::Builder::Methylation::BWAmeth;

use Moose;
extends 'CoGe::Builder::Methylation::Analyzer';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

my ($PICARD, $PILEOMETH);
sub BUILD { # called immediately after constructor
    my $self = shift;
    $PICARD = $self->conf->{PICARD};
    unless ($PICARD) {
        CoGe::Exception::Generic->throw(message => 'Missing PICARD in config file');
    }

    $PILEOMETH = $self->conf->{PILEOMETH} || 'PileOMeth';
}

sub build {
    my $self = shift;
    my %opts = @_;
    my $bam_file = shift @{$opts{data_files}};

    # Validate inputs not already checked in Request
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }

    my $genome = $self->request->genome;

    # Set metadata for the pipeline being used
    my $annotations = generate_additional_metadata();
    my @annotations2 = CoGe::Core::Metadata::to_annotations($self->params->{additional_metadata});
    push @$annotations, @annotations2;

    #
    # Build the workflow
    #
    if ($self->params->{methylation_params}->{'picard-deduplicate'}) {
        $self->add_task(
            $self->picard_deduplicate($bam_file)
        );
        $bam_file = $self->previous_output;
    }

# mdb removed 1/27/15 -- resulting .svg files are not used, reinstate when there is a way to show/download them
#    push @tasks, create_pileometh_plot_job(
#        bam_file => $input_file,
#        gid => $genome->id,
#        staging_dir => $staging_dir
#    );

    $self->add_task(
        $self->pileometh_extraction(
            bam_file => $bam_file,
            gid => $genome->id
        )
    );
    
    my @outputs = @{$self->previous_outputs};
    foreach my $file (@outputs) {
        my ($name) = $file =~ /(CHG|CHH|CpG)/;

        $self->add_task(
            $self->create_pileometh_import_job(
                input_file => $file,
                name => $name
            )
        );

        my $md = clone($metadata);
        $md->{name} .= " ($name methylation)";

        $self->add_task(
            $self->load_experiment(
                metadata => $md,
                gid => $genome->id,
                input_file => $self->previous_output,
                name => $name,
                annotations => $annotations
            )
        );
    }
}

sub picard_deduplicate {
    my $self = shift;
    my $bam_file = shift;

    my $cmd = 'java -jar ' . $PICARD;
    
    my $output_file = $bam_file . '-deduplicated.bam';
    
    return {
        cmd => "$cmd MarkDuplicates REMOVE_DUPLICATES=true INPUT=$bam_file METRICS_FILE=$bam_file.metrics OUTPUT=$output_file.tmp ; mv $output_file.tmp $output_file",
        args => [
#            ['MarkDuplicates', '', 0],
#            ['REMOVE_DUPLICATES=true', '', 0],
#            ["INPUT=$bam_file", '', 0],
#            ["METRICS_FILE=$bam_file.metrics", '', 0],
#            ["OUTPUT=$output_file", '', 0],
        ],
        inputs => [
            $bam_file
        ],
        outputs => [
            $output_file
        ],
        description => "Deduplicating PCR artifacts using Picard"
    };
}

sub pileometh_plot {
    my $self = shift;
    my %opts = @_;
    my $bam_file = $opts{bam_file};
    my $gid = $opts{gid};

    my $BWAMETH_CACHE_FILE = catfile(get_genome_cache_path($gid), 'bwameth_index', 'genome.faa.reheader.faa');
    my $output_prefix = 'pileometh';
    
    return {
        cmd => $PILEOMETH,
        args => [
            ['mbias', '', 0],
            ['--CHG', '', 0],
            ['--CHH', '', 0],
            ['', $BWAMETH_CACHE_FILE, 0],
            ['', $bam_file, 0],
            [$output_prefix, '', 0]
        ],
        inputs => [
            $bam_file
        ],
        outputs => [
            catfile($self->staging_dir, $output_prefix . '_OB.svg'),
            catfile($self->staging_dir, $output_prefix . '_OT.svg')
        ],
        description => "Plotting methylation bias with PileOMeth"
    };
}

sub pileometh_extraction {
    my $self = shift;
    my %opts = @_;
    my $bam_file = $opts{bam_file};
    my $gid = $opts{gid};

    my $params = $self->params->{methylation_params};
    my $q  = $params->{'pileometh-min_converage'} // 10;
    my $ot = $params->{'--OT'} // '0,0,0,0';
    my $ob = $params->{'--OB'} // '0,0,0,0';
    
    my $BWAMETH_CACHE_FILE = catfile(get_genome_cache_path($gid), 'bwameth_index', 'genome.faa.reheader.faa');
    my $output_prefix = to_filename_without_extension($bam_file);
    
    return {
        cmd => $PILEOMETH,
        script => undef,
        args => [
            ['extract', '', 0],
            ['--methylKit', '', 0],
            ['--CHG', '', 0],
            ['--CHH', '', 0],
            ['-q', $q, 0],
            ['--OT', $ot, 0],
            ['--OB', $ob, 0],
            ['', $BWAMETH_CACHE_FILE, 0],
            ['', $bam_file, 1]
        ],
        inputs => [
            $bam_file
        ],
        outputs => [
            catfile($self->staging_dir, $output_prefix . '_CpG.methylKit'),
            catfile($self->staging_dir, $output_prefix . '_CHH.methylKit'),
            catfile($self->staging_dir, $output_prefix . '_CHG.methylKit')
        ],
        description => "Extracting methylation calls with PileOMeth"
    };
}

sub pileometh_import {
    my $self = shift;
    my %opts = @_;
    my $input_file = $opts{input_file};
    my $name = $opts{name};

    my $params = $self->params->{methylation_params};
    my $c = $params->{'pileometh-min_converage'} // 10;
    
    my $output_file = $input_file . '.filtered.coge.csv';
    
    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'methylation', 'coge-import_pileometh.py'),
        script => undef,
        args => [
            ['-u', 'f', 0],
            ['-c', $c, 0],
            ['', $input_file, 1]
        ],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file
        ],
        description => "Converting $name"
    };
}

sub generate_additional_metadata {
    my $self = shift;
    my $methylation_params = $self->params->{methylation_params};

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php/Methylation_Analysis_Pipeline||note|Generated by CoGe's Methylation Analysis Pipeline};

    push @annotations, 'note|picard MarkDuplicates REMOVE_DUPLICATES=true ' if ($methylation_params->{'bismark-deduplicate'});
    #push @annotations, 'note|PileOMeth mbias --CHG --CHH' if ($methylation_params->{'bismark-deduplicate'});
    push @annotations, 'note|PileOMeth extract --CHG --CHH ' . join(' ', map { $_.' '.$methylation_params->{$_} } ('-q', '--OT', '--OB'));

    return \@annotations;
}

1;
