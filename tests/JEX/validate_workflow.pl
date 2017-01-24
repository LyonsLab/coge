#!/usr/bin/env perl

use CoGe::JEX::Jex;
use Data::Dumper;
use Time::Local;

my $START_ID = 40_000;
my $HOST = 'localhost';
my $PORT = 5151;

# Connect to JEX
my $jex = CoGe::JEX::Jex->new( host => $HOST, port => $PORT );
unless ($jex) {
    print STDERR "Couldn't connect to JEX\n";
    exit(-1);
}

# Get workflows
my @workflow_ids;
my $workflow_id = $ARGV[0];
if ($workflow_id) {
    push @workflow_ids, $workflow_id;
}
else {
    my $workflows = $jex->find_workflows(undef, "completed");
    @workflow_ids = map { $_->[0] } @$workflows;
}

# Validate workflows
my $totalValidated = 0;
print STDERR scalar(@workflow_ids), " candidate workflow(s)\n";
foreach my $id (@workflow_ids) {
    next if (defined $START_ID && $id < $START_ID);

    # Get workflow
    my $workflow = $jex->get_job($id);

    # Index tasks by outputs
    my %tasksByOutput;
    foreach my $task (@{$workflow->{jobs}}) {
        foreach my $output (@{$task->{outputs}}) {
            $tasksByOutput{$output} = $task;
        }
    }

    # Validate workflow
    next unless ($workflow->{jobs} && @{$workflow->{jobs}});
    if (lc($workflow->{status}) ne 'completed') {
        print STDERR "Skipping incomplete workflow $id\n";
        next;
    }
    print STDERR "Validating workflow ", $id, "\n";
    validate_sequence($workflow, \%tasksByOutput);

    $totalValidated++;
}

print STDERR "All done ... $totalValidated workflows validated\n";

exit;

#-----------------------------------------------------------------------------

sub validate_sequence {
    my $workflow = shift;
    my $taskIndex = shift;

    foreach my $task (@{$workflow->{jobs}}) {
        my $inputs = $task->{inputs};
        #print STDERR "   task: ", ($inputs ? scalar(@$inputs) : 0), " input(s)\n";
        next unless ($inputs && @$inputs);
        next unless ($task->{started});

        foreach my $input (@$inputs) {
            my $task2 = $taskIndex->{$input};
            next unless ($task2 && $task2->{ended});
            next if (lc($task2->{status}) eq 'skipped');

            my ($y1, $mo1, $d1, $h1, $m1, $s1, $p1) = $task2->{ended} =~ /(\d+)\/(\d+)\/(\d+) at (\d+):(\d+):(\d+)(\w+)/;
            $h1 += 11 if ($p1 eq 'PM');
            my $time1 = timelocal( $s1, $m1-1, $h1, $d1, $mo1, $y1 );

            my ($y2, $mo2, $d2, $h2, $m2, $s2, $p2) = $task->{started} =~ /(\d+)\/(\d+)\/(\d+) at (\d+):(\d+):(\d+)(\w+)/;
            $h2 += 11 if ($p2 eq 'PM');
            my $time2 = timelocal( $s2, $m2-1, $h2, $d2, $mo2, $y2 );

            my $diff = $time2 - $time1;
            if ($diff < 0) {
                print STDERR '   ERROR: Task started too soon in workflow ', $workflow->{id}, "\n";
                print STDERR '   overlap:   ', abs($diff), " sec\n";
                print STDERR '   task:      ', $task->{description}, "\n";
                print STDERR '   input:     ', $input, "\n";
                print STDERR '   cur task:  ', $task->{started}, '   ', $task->{ended}, "\n";
                print STDERR '   prev task: ', $task2->{started}, '   ', $task2->{ended}, "\n";
            }
        }
    }
}
