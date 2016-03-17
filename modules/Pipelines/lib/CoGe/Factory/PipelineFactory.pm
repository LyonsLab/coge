package CoGe::Factory::PipelineFactory;

use Moose;

use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);
use Data::Dumper;

use CoGe::Core::Storage qw(get_workflow_paths);
use CoGe::Builder::Export::Fasta;
use CoGe::Builder::Export::Gff;
use CoGe::Builder::Export::Experiment;
use CoGe::Builder::Load::Experiment;
use CoGe::Builder::Load::BatchExperiment;
use CoGe::Builder::Load::Genome;
use CoGe::Builder::Load::Annotation;
use CoGe::Builder::SNP::IdentifySNPs;
use CoGe::Builder::Tools::SynMap;
use CoGe::Builder::Expression::MeasureExpression;
use CoGe::Builder::Methylation::MeasureMethylation;
use CoGe::Builder::PopGen::MeasureDiversity;

has 'db' => (
    is => 'ro',
    required => 1
);

has 'conf' => (
    is => 'ro',
    required => 1
);

has 'user' => (
    is  => 'ro',
    required => 1
);

has 'jex' => (
    is => 'ro',
    required => 1
);

sub get {
    my ($self, $message) = @_;

    my $request = {
        params    => $message->{parameters},
        requester => $message->{requester},
        db        => $self->db,
        jex       => $self->jex,
        user      => $self->user,
        conf      => $self->conf
    };

    # Select pipeline builder
    my $builder;
    if ($message->{type} eq "export_gff") {
        $builder = CoGe::Builder::Export::Gff->new($request);
    }
    elsif ($message->{type} eq "export_fasta") {
        $builder = CoGe::Builder::Export::Fasta->new($request);
    }
    elsif ($message->{type} eq "export_experiment") {
        $builder = CoGe::Builder::Export::Experiment->new($request);
    }
    elsif ($message->{type} eq "load_experiment") {
        $builder = CoGe::Builder::Load::Experiment->new($request);
    }
    elsif ($message->{type} eq "load_batch") {
        $builder = CoGe::Builder::Load::BatchExperiment->new($request);
    }
    elsif ($message->{type} eq "load_genome") {
        $builder = CoGe::Builder::Load::Genome->new($request);
    }
    elsif ($message->{type} eq "load_annotation") {
        $builder = CoGe::Builder::Load::Annotation->new($request);
    }
    elsif ($message->{type} eq "analyze_snps") {
        $builder = CoGe::Builder::SNP::IdentifySNPs->new($request);
    }
    elsif ($message->{type} eq "synmap") {
        $builder = CoGe::Builder::Tools::SynMap->new($request);
    }
    elsif ($message->{type} eq "analyze_expression") {
        $builder = CoGe::Builder::Expression::MeasureExpression->new($request);
    }
    elsif ($message->{type} eq "analyze_methylation") {
        $builder = CoGe::Builder::Methylation::MeasureMethylation->new($request);
    }    
    elsif ($message->{type} eq "analyze_diversity") {
        $builder = CoGe::Builder::PopGen::MeasureDiversity->new($request);
    }
    else {
        print STDERR "PipelineFactory::get unknown type\n";
        return;
    }
    
    # Initialize workflow
    $builder->workflow( $self->jex->create_workflow(name => $builder->get_name, init => 1) );
    return unless ($builder->workflow && $builder->workflow->id);
    my ($staging_dir, $result_dir) = get_workflow_paths($self->user->name, $builder->workflow->id);
    $builder->staging_dir($staging_dir);
    $builder->result_dir($result_dir);
    $builder->workflow->logfile(catfile($result_dir, "debug.log"));
    
    # Get a tiny URL to a status page #TODO simplyify this
    my ($page, $link);
    if ($message->{requester}) { # request is from internal web page - external API requests will not have a 'requester' field
        $page = $message->{requester}->{page};
        my $url = $message->{requester}->{url}; #FIXME why page and url separate? merge these ...
        if ($url) {
            $link = CoGe::Accessory::Web::get_tiny_link( url => $self->conf->{SERVER} . $url . "&wid=" . $builder->workflow->id );
        }
        elsif ($page) {
            $link = CoGe::Accessory::Web::get_tiny_link( url => $self->conf->{SERVER} . $page . "?wid=" . $builder->workflow->id );
        }
    }
    else { # otherwise infer status page by job type (for the DE)
        $page = 'API';
        if ($message->{type} eq 'load_genome') {
            $link = CoGe::Accessory::Web::get_tiny_link( url => $self->conf->{SERVER} . 'LoadGenome.pl' . "?wid=" . $builder->workflow->id );
        }
        elsif ($message->{type} eq 'load_experiment') {
            $link = CoGe::Accessory::Web::get_tiny_link( url => $self->conf->{SERVER} . 'LoadExperiment.pl' . "?wid=" . $builder->workflow->id );
        }
    }
    $builder->page($page) if $page;
    $builder->site_url($link) if $link;

    # Construct the workflow
    my $rc = $builder->build;
    unless ($rc) {
        $rc = 'undef' unless defined $rc;
        print STDERR "PipelineFactory::get build failed, rc=$rc\n";
        return;
    }
    
    # Dump raw workflow to file for debugging -- mdb added 2/17/16
    make_path($result_dir);
    open(my $fh, '>', catfile($result_dir, 'workflow.log'));
    print $fh Dumper $builder->workflow, "\n";
    close($fh);
    
    return $builder;
}

1;
