package CoGe::Builder::Methylation::Analyzer;

use Moose;
extends 'CoGe::Builder::Buildable';

use CoGe::Builder::Methylation::Bismark;
use CoGe::Builder::Methylation::BWAmeth;
use CoGe::Builder::Methylation::Metaplot;
use CoGe::Exception::Generic;
use CoGe::Exception::MissingField;

sub build {
    my $self = shift;
    my %opts = @_;
    my ($bam_file) = @{$opts{data_files}};
    unless ($bam_file) {
        CoGe::Exception::Generic->throw( message => 'Missing bam input' );
    }

    # Validate inputs (that weren't already checked in Request)
    my $metadata = $self->params->{metadata};
    unless ($metadata) {
        CoGe::Exception::MissingField->throw( message => "Missing metadata" );
    }

    #
    # Build workflow
    #

    # Add methylation analysis tasks
    my $analyzer;
    my $method = lc($self->params->{methylation_params}->{method});
    if ($method eq 'bismark') {
        $analyzer = CoGe::Builder::Methylation::Bismark->new($self);
        $analyzer->build(data_files => [$bam_file]);
    }
    elsif ($method eq 'bwameth') {
        $analyzer = CoGe::Builder::Methylation::BWAmeth->new($self);
        $analyzer->build(data_files => [$bam_file]);
    }
    else {
        CoGe::Exception::Generic->throw(message => 'Invalid methylation method');
    }
    $self->add($analyzer);

    # Add metaplot workflow (if specified and genome is annotated)
    if ( $self->params->{methylation_params}->{metaplot_params} ) {
        my $metaplot = CoGe::Builder::Methylation::Metaplot->new($self);
        $metaplot->build(data_files => [$bam_file]);
        $self->add($metaplot);
    }
}

__PACKAGE__->meta->make_immutable;

1;
