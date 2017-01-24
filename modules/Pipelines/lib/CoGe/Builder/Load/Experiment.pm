package CoGe::Builder::Load::Experiment;

use Moose;
extends 'CoGe::Builder::Buildable';

use Data::Dumper qw(Dumper);

use CoGe::Accessory::Web qw(url_for);
use CoGe::Core::Experiment qw(detect_data_type);
use CoGe::Builder::Alignment::Aligner;
use CoGe::Builder::Expression::qTeller;
use CoGe::Builder::SNP::SNPFinder;
use CoGe::Builder::Methylation::Analyzer;
use CoGe::Builder::Protein::ChIPseq;
use CoGe::Exception::MissingField;
use CoGe::Exception::Generic;

sub get_name {
    my $self = shift;
    my $metadata = $self->params->{metadata};
    my $info = '"' . $metadata->{name};
    $info .= ": " . $metadata->{description} if $metadata->{description};
    $info .= " (v" . $metadata->{version} . ")";
    $info .= '"';
    return "Load Experiment " . $info;
}

sub get_site_url {
    my $self = shift;
    return url_for('LoadExperiment.pl', wid => $self->workflow->id);
}

sub build {
    my $self = shift;
    my %opts = @_;
    my $input_files = $opts{data_files};
    unless ($input_files && @$input_files) {
        CoGe::Exception::Generic->throw(message => 'Missing inputs');
    }

    # Validate inputs not already checked in Request
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw(message => "Missing metadata");
    }
    # mdb added 2/25/15 - convert from Mojolicious boolean: bless( do{\\(my $o = 1)}, 'Mojo::JSON::_Bool' )
    $metadata->{restricted} = $metadata->{restricted} ? 1 : 0;

    my $genome = $self->request->genome;

    # Determine file type if not set #FIXME this is a little kludgey
    my $data = $self->params->{source_data};
    my $file_type = $data->[0]->{file_type}; # type of first data file
    ($file_type) = detect_data_type($file_type, $data->[0]->{path}) unless $file_type;
    
    #
    # Build workflow
    #

    # Add analytical tasks based on file type
    if ( $file_type eq 'fastq' || $file_type eq 'bam' || $file_type eq 'sra' ) {
        my @bam_files;
        my @raw_bam_files; # mdb added 2/29/16 for Bismark, COGE-706

        # Align fastq file or take existing bam
        if ( $file_type && ( $file_type eq 'fastq' || $file_type eq 'sra' ) ) {
            # Add alignment workflow
            my $aligner = CoGe::Builder::Alignment::Aligner->new($self);
            $aligner->build(data_files => $input_files);
            $self->add($aligner);
            @bam_files     = @{$aligner->bam};
            @raw_bam_files = @{$aligner->raw_bam};
            unless (@bam_files) {
                CoGe::Exception::Generic->throw(message => 'Alignment returned no results');
            }
        }
        elsif ( $file_type && $file_type eq 'bam' ) {
            $self->add(
                $self->load_bam(
                    bam_file => $input_files->[0]
                )
            );
            @bam_files = @raw_bam_files = @$input_files;
        }
        else {
            CoGe::Exception::Generic->throw(message => 'Invalid file type');
        }
        
        # Add expression workflow (if specified)
        if ( $self->params->{expression_params} ) {
            my $expr = CoGe::Builder::Expression::qTeller->new($self);
            $expr->build(data_files => \@bam_files);
            $self->add($expr);
        }
        
        # Add SNP workflow (if specified)
        if ( $self->params->{snp_params} ) {
            my $isBamSorted = ($file_type ne 'bam');
            my $snp = CoGe::Builder::SNP::SNPFinder->new($self);
            $snp->build(data_files => \@bam_files, is_sorted => $isBamSorted);
            $self->add($snp);
        }
        
        # Add methylation workflow (if specified)
        if ( $self->params->{methylation_params} ) {
            my $aligner = CoGe::Builder::Methylation::Analyzer->new($self);
            $aligner->build(data_files => [ $bam_files[0], $raw_bam_files[0] ]);
            $self->add($aligner);
        }
        
        # Add ChIP-seq workflow (if specified)
        if ( $self->params->{chipseq_params} ) {
            my $chipseq = CoGe::Builder::Protein::ChIPseq->new($self);
            $chipseq->build(data_files => \@bam_files);
            $self->add($chipseq);
        }
    }
    # Else, all other file types
    else {
        # Add conversion step for BigWig files
        my $input_file = $input_files->[0];
        if ( $file_type eq 'bw' ) {
            ($input_file) = $self->add(
                $self->bigwig_to_wig($input_file)
            );
        }

        $self->add(
            $self->load_experiment(
                gid => $genome->id,
                input_file => $input_file,
                normalize => $self->params->{normalize} ? $self->params->{normalize_method} : 0
            )
        );
    }
}

sub bigwig_to_wig {
    my $self = shift;
    my $input_file  = shift;

    my $filename    = basename($input_file) . '.wig';
    my $output_file = catfile($self->staging_dir, $filename);
    my $done_file   = $output_file . '.done';

    my $cmd = $self->conf->{BIGWIGTOWIG} || 'bigWigToWig';

    return {
        cmd => "$cmd $input_file $output_file && touch $done_file",
        args => [],
        inputs => [
            $input_file
        ],
        outputs => [
            $output_file,
            $done_file
        ],
        description => 'Converting BigWig to WIG format'
    };
}

__PACKAGE__->meta->make_immutable;

1;
