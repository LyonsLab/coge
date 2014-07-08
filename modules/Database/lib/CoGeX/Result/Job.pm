package CoGeX::Result::Job;

use v5.10;
use strict;
use warnings;

#use DateTime;
#use DateTime::Format::HTTP;
use Time::Local;
#use DateTime::Duration;
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
      . $self->elapsed_time . ' | '
      . $self->page . ' | '
      . ( $self->log ? $self->log->description . ' | ' : '')
      . $self->status_description;
}

sub info_html {
    my $self = shift;

    my $color = $self->status_color;

    return
        $self->start_time . ' | '
      . $self->elapsed_time . ' | '
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

# mdb removed 10/15/13 -- too slow
#    my $start_time = DateTime::Format::HTTP->parse_datetime($self->start_time);
#    my $diff;
#    if ($self->end_time and $self->status > 2) {
#        my $end_time = DateTime::Format::HTTP->parse_datetime($self->end_time);
#        $diff = $end_time->subtract_datetime($start_time);
#    } else {
#        my $timezone = DateTime::TimeZone->new( name => 'local' );
#        my $dt = DateTime->now(time_zone => $timezone);
#
#        $diff = $dt->subtract_datetime($start_time);
#    }
#
#    my ($days, $hours, $minutes, $seconds) = $diff->in_units('days', 'hours', 'minutes', 'seconds');
#
#    my $elapsed = '';
#    $elapsed .= "${days}d " if $days;
#    $elapsed .= "${hours}h " if $hours;
#    $elapsed .= "${minutes}m " if $minutes && not $days;
#    $elapsed .= "${seconds}s" if $seconds && not $days;
#    return $elapsed;

	my ($y1, $mo1, $d1, $h1, $m1, $s1) = $self->start_time =~ /(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+)/;
	my $time1 = timelocal( $s1, $m1, $h1, $d1, $mo1-1, $y1 );

	my $time2;
	if ($self->end_time and $self->status > 2) {
		my ($y2, $mo2, $d2, $h2, $m2, $s2) = $self->end_time =~ /(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+)/;
		$time2 = timelocal( $s2, $m2, $h2, $d2, $mo2-1, $y2 );
	}
	else {
		$time2 = time;
	}

	my $diff = $time2 - $time1;
	my $d = int($diff / (60*60*24));
	$diff -= $d * (60*60*24);
	my $h = int($diff / (60*60));
	$diff -= $h * (60*60);
	my $m = int($diff / 60);
	$diff -= $m * 60;
	my $s = $diff % 60;

    my $elapsed = '';
    $elapsed .= "${d}d " if $d > 0;
    $elapsed .= "${h}h " if $h > 0;
    $elapsed .= "${m}m " if $m > 0 && $d <= 0;
    $elapsed .= "${s}s" if $s > 0 && $d <= 0;
	return $elapsed;
}

1;
