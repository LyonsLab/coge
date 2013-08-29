package CoGeX::Result::Job;

use v5.10;
use strict;
use warnings;

use DateTime;
use DateTime::Format::HTTP;
use base 'DBIx::Class::Core';

=head1 NAME

My::Schema::Result::Log

=cut

__PACKAGE__->table("job");

=head1 ACCESSORS

=cut

__PACKAGE__->add_columns(
    "job_id",
    {
        data_type     => "INT",
        default_value => undef,
        is_nullable   => 0,
        size          => 11
    },
    "user_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
    "page",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 0,
        size          => 255
    },
    "start_time",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
    "end_time",
    { data_type => "TIMESTAMP", default_value => undef, is_nullable => 0 },
    "link",
    {
        data_type     => "VARCHAR",
        default_value => undef,
        is_nullable   => 1,
        size          => 255
    },
    "status",
    { data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 5 },
    "type",
    { data_type => "TINYINT", default_value => 0, is_nullable => 0, size => 5 },
    "process_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 5 },
    "log_id",
    { data_type => "INT", default_value => 0, is_nullable => 0, size => 11 },
);
__PACKAGE__->set_primary_key("job_id");

__PACKAGE__->belongs_to( 'user' => "CoGeX::Result::User", 'user_id' );
__PACKAGE__->belongs_to( 'log'  => "CoGeX::Result::Log",  'log_id' );

################################################ subroutine header begin ##

=head2 info

 Usage     : $self->info
 Purpose   : generate a string of information about the log entry
 Returns   : a string
 Argument  : None
 Throws    : None
 Comments  :

=cut

################################################## subroutine header end ##

sub info {
    my $self = shift;
    return
        $self->start_time . ' | '
      . $self->page . ' | '
      . ( $self->log ? $self->log->description . ' | ' : '')
      . $self->status_description;
}

sub info_html {
    my $self = shift;
    
    my $color = $self->status_color;
    
    return
        $self->start_time . ' | '
      . $self->page . ' | '
      . ( $self->log ? $self->log->description . ' | ' : '')
      . ($color ? '<span style="padding-bottom:1px;padding-right:5px;padding-left:5px;border-radius:15px;color:white;background-color:'.$color.';">' . $self->status_description . '</span>' : $self->status_description);
}

sub status_description {
    my $self = shift;
    given ( $self->status ) {
        when (0) { return 'Scheduled';  }
        when (1) { return 'Running';	}
        when (2) { return 'Completed';	}
        when (3) { return 'Cancelled';	}
        when (4) { return 'Terminated';	}
        when (5) { return 'Failed';		}
        default  { return 'Unknown';	}
    }
}

sub status_color {
    my $self = shift;
    given ( $self->status ) {
        when (1) { return 'yellowgreen'; }
        when (3) { return 'salmon';	}
        when (4) { return 'salmon';	}
        when (5) { return 'salmon';	}
        default  { return; }
    }
}

sub elapsed_time {
    my $self = shift;
    my $diff;

    my $start_time = DateTime::Format::HTTP->parse_datetime($self->start_time);
    if ($self->end_time) {
        my $end_time = DateTime::Format::HTTP->parse_datetime($self->end_time);
        $diff = $end_time->subtract_datetime($start_time);
    } else {
        $diff = DateTime->now()->subtract_datetime($start_time);
    }

    my $days = $diff->in_units('days');
    my $hours = $diff->hours();

    my $elapsed = '';
    $elapsed .= "$days days " if $days;
    $elapsed .= "$hours hrs " if $hours;
    $elapsed .= $diff->minutes() . " mins " if $diff->minutes() && not $days;
    $elapsed .= $diff->seconds() . " secs" if $diff->seconds() && not $days;

    return $elapsed
}

1;
