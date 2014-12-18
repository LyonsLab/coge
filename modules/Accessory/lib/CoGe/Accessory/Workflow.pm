package CoGe::Accessory::Workflow;
use strict;
use warnings;
use 5.10.0;

use Moose;
use Data::Dumper;

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
    my $size        = $self->jobs;
    my $overwrite;

    # Set default
    $options //= {}; #/

    if ( defined( $opts->{overwrite} ) && $opts->{overwrite} > 0 ) {
        $overwrite = 1;
    }
    else {
        $overwrite = 0;
    }

    push(
        @{ $self->jobs },
        {
            cmd         => $cmd,
            options     => $options,
            script      => $script,
            args        => $args,
            description => $description,
            overwrite   => $overwrite,
            inputs      => $inputs,
            outputs     => $outputs,
        }
    );

    return scalar( $self->jobs() ) > $size;
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

__PACKAGE__->meta->make_immutable;
no Moose;
1;
__END__
