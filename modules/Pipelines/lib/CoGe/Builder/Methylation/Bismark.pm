package CoGe::Builder::Methylation::Bismark;

use Moose;
extends 'CoGe::Builder::Methylation::Analyzer';

use Data::Dumper;
use File::Spec::Functions qw(catdir catfile);

use CoGe::Accessory::Web;
use CoGe::Accessory::Utils;
use CoGe::Core::Storage;
use CoGe::Core::Metadata;
use CoGe::Exception::Generic;

sub build {
    my $self = shift;
    my %opts = @_;
    my $bam_file = shift @{$opts{data_files}}; # IMPORTANT: this should be the unsorted version, see COGE-706 and http://seqanswers.com/forums/showthread.php?t=45192
    unless ($bam_file) {
        CoGe::Exception::Generic->throw(message => 'Missing bam');
    }

    # Validate inputs not already checked in Request
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }

    my $gid = $self->request->genome->id;

    # Set metadata for the pipeline being used
    my $annotations = generate_additional_metadata($self->params->{read_params}, $self->params->{methylation_params});
    my @annotations2 = CoGe::Core::Metadata::to_annotations($self->params->{additional_metadata});
    push @$annotations, @annotations2;

    #
    # Build the workflow
    #
    my (@tasks, @done_files);

    if ($self->params->{methylation_params}->{'bismark-deduplicate'}) {
        $self->add_task(
            $self->bismark_deduplicate($bam_file)
        );
        $bam_file = $self->previous_output;
    }

    $self->add_task(
         $self->extract_methylation($bam_file)
    );
    
    my @outputs = $self->previous_outputs;
    while (@outputs) {
        my $file1 = shift @outputs;
        my $file2 = shift @outputs;
        
        my ($name) = $file1 =~ /(CHG|CHH|CpG)/;

        $self->add_task(
            $self->bismark_import(
                ob_input_file => $file1,
                ot_input_file => $file2,
                name => $name
            )
        );

        my $md = clone($metadata);
        $md->{name} .= " ($name methylation)";

        $self->add_task(
            $self->load_experiment(
                metadata    => $md,
                gid         => $gid,
                input_file  => $self->previous_output,
                name        => $name,
                annotations => $annotations
            )
        );
    }
}

sub bismark_deduplicate {
    my $self = shift;
    my $bam_file = shift;

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    my $cmd = ($self->conf->{BISMARK_DIR} ? catfile($self->conf->{BISMARK_DIR}, 'deduplicate_bismark') : 'deduplicate_bismark');
    $cmd = 'nice ' . $cmd;
    
    my $args;
    if ($read_type eq 'paired') {
        push @$args, ['-p', '', 0];
    }
    else { # single-ended
        push @$args, ['-s', '', 0];
    }
    
    push @$args, ['--bam', $bam_file, 1];
    
    my $output_file = qq[$bam_file.deduplicated.bam];
    
    return {
        cmd => $cmd,
        args => $args,
        inputs => [ $bam_file ],
        outputs => [ $output_file ],
        description => "Deduplicating PCR artifacts using Bismark"
    };
}

sub extract_methylation {
    my $self = shift;
    my $bam_file = shift;

    my $read_params = $self->params->{read_params};
    my $read_type   = $read_params->{read_type} // 'single';

    my $methylation_params = $self->params->{methylation_params};
    my $ignore             = $methylation_params->{'--ignore'} // 0;
    my $ignore_r2          = $methylation_params->{'--ignore_r2'} // 0;
    my $ignore_3prime      = $methylation_params->{'--ignore_3prime'} // 0;
    my $ignore_3prime_r2   = $methylation_params->{'--ignore_3prime_r2'} // 0;

    my $cmd = $self->conf->{BISMARK_DIR} ? catfile($self->conf->{BISMARK_DIR}, 'bismark_methylation_extractor') : 'bismark_methylation_extractor';
    $cmd = 'nice ' . $cmd;
    
    my $name = to_filename_without_extension($bam_file);
    
    my $args = [
        ['--multicore',     4,                  0],
        ['--output',        $self->staging_dir, 0],
        ['--ignore',        $ignore,            0],
        ['--ignore_3prime', $ignore_3prime,     0]
    ];
    
    if ($read_type eq 'paired') {
        push @$args, ['-p',                 '',                0];
        push @$args, ['--ignore_r2',        $ignore_r2,        0];
        push @$args, ['--ignore_3prime_r2', $ignore_3prime_r2, 0];
    }
    
    push @$args, ['', $bam_file, 0];
    
    my $done_file = catfile($self->staging_dir, 'extract_methylation.done');
    push @$args, ['', '&& touch ' . $done_file, 0]; # kludge to ensure proper sequence since --output dir must be absolute
    
    return {
        cmd => $cmd,
        args => $args,
        inputs => [
            $bam_file
        ],
        outputs => [
            catfile($self->staging_dir, 'CHG_OB_' . $name . '.txt'),
            catfile($self->staging_dir, 'CHG_OT_' . $name . '.txt'),
            catfile($self->staging_dir, 'CHH_OB_' . $name . '.txt'),
            catfile($self->staging_dir, 'CHH_OT_' . $name . '.txt'),
            catfile($self->staging_dir, 'CpG_OB_' . $name . '.txt'),
            catfile($self->staging_dir, 'CpG_OT_' . $name . '.txt')
        ],
        description => "Extracting methylation status"
    };
}

sub bismark_import {
    my $self = shift;
    my %opts = @_;
    my $ot_input_file = $opts{ot_input_file};
    my $ob_input_file = $opts{ob_input_file};
    my $name          = $opts{name};

    my $methlyation_params = $self->params->{methylation_params};
    my $min_coverage  = $methlyation_params->{min_coverage} // 5;
    
    my $output_file = qq[$name.filtered.coge.csv];
    
    return {
        cmd => catfile($self->conf->{SCRIPTDIR}, 'methylation', 'coge-import_bismark.py'),
        script => undef,
        args => [
            ['-u',   'f',            0],
            ['-c',   $min_coverage,  0],
            ['--OT', $ot_input_file, 1],
            ['--OB', $ob_input_file, 1],
            ['-o',   $name,          0]
        ],
        inputs => [
            $ot_input_file,
            $ob_input_file,
            catfile($self->staging_dir, 'extract_methylation.done')
        ],
        outputs => [
            catfile($self->staging_dir, $output_file),
        ],
        description => "Converting $name"
    };
}

sub generate_additional_metadata {
    my $self = shift;
    my $read_params        = $self->params->{read_params};
    my $methylation_params = $self->params->{methylation_params};

    my @annotations;
    push @annotations, qq{https://genomevolution.org/wiki/index.php/Methylation_Analysis_Pipeline||note|Generated by CoGe's Methylation Analysis Pipeline};

    if ($read_params->{'read_type'} eq 'paired') {
        if ($methylation_params->{'bismark-deduplicate'}) {
            push @annotations, 'note|deduplicate_bismark -p';
        }
        push @annotations, 'note|bismark_methylation_extractor ' . join(' ', map { $_.' '.$methylation_params->{$_} } ('--ignore', '--ignore_r2', '--ignore_3prime', '--ignore_3prime_r2'));
    }
    else {
        if ($methylation_params->{'bismark-deduplicate'}) {
            push @annotations, 'note|deduplicate_bismark -s';
        }
        push @annotations, 'note|bismark_methylation_extractor ' . join(' ', map { $_.' '.$methylation_params->{$_} } ('--ignore', '--ignore_3prime'));
    }

    return \@annotations;
}

1;