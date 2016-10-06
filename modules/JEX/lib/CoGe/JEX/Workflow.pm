package CoGe::JEX::Workflow;
use strict;
use warnings;
use 5.10.0;

use Moose;
use Data::Dumper;
use Data::Compare;

# Attributes
has 'id' => (
    is       => 'ro',
);

has 'name' => (
    is       => 'ro',
    isa      => 'Str',
    required => 1,
    default  => "",
);

has 'priority' => (
    is       => 'ro',
    isa      => 'Int',
    default  => 0,
);

has 'logfile' => (
    is       => 'rw',
    default  => "",
);

has 'jobs' => (
    is       => 'rw',
    isa      => 'ArrayRef',
    default  => sub { [] },
    required => 0,
);

# Public functions
sub add_job {
    my ( $self, $opts ) = @_;
    #print STDERR 'Workflow::add_job ', Dumper $opts, "\n";
    my $cmd         = $opts->{cmd};
    my $script      = "" unless defined( $opts->{script} );
    my $args        = $opts->{args} || [];
    my $inputs      = $opts->{inputs} || [];
    my $options     = $opts->{options};
    my $outputs     = $opts->{outputs} || [];
    my $description = $opts->{description};
    my $overwrite;

    # Set defaults
    $options //= {}; #/

    # Build job object
    if ( defined( $opts->{overwrite} ) && $opts->{overwrite} > 0 ) {
        $overwrite = 1;
    }
    else {
        $overwrite = 0;
    }
    my $job = {
        cmd         => $cmd,
        options     => $options,
        script      => $script,
        args        => $args,
        description => $description,
        overwrite   => $overwrite,
        inputs      => $inputs,
        outputs     => $outputs,
    };

    # mdb added 1/5/15 - prevent duplicate jobs, JEX doesn't handle them correctly
    foreach (@{$self->jobs}) {
        if (Compare($_, $job)) {
            print STDERR "Workflow::add_job warning: skipping duplicate job named '", $job->{description}, "'\n";
            return 1;
        }
    }
    
    # Add job to workflow
    push @{$self->jobs}, $job;

    return 1;
}

sub add_jobs {
    my ($self, $jobs) = @_;
    #print STDERR Dumper "Workflow::add_jobs ", $jobs, "\n";
    
    my $rc = 1;
    foreach (@$jobs) {
        $rc &= $self->add_job($_);
    }
    
    return $rc;
}

sub log {
	my $self = shift;
    $| = 1;
    my $message = shift;
    $message =~ /(.*)/xs; # why is this done? doesn't seem to remove anything
    $message = $1;
    open( OUT, '>>' . $self->logfile ) || return;
    print OUT $message, "\n";
    close OUT;
}

sub log_section{
	my $self = shift;
	$self->log( "#" x (25) );
	$self->log( shift );
	$self->log( "#" x (25) );
	$self->log( "" );
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
__END__
