#! /usr/bin/perl -w
use v5.10;
use strict;

use CGI;
use Digest::MD5 qw(md5_base64);
use Time::Piece;

use HTML::Template;
use JSON::XS;

# CoGe packages
use CoGeX;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Utils qw(format_time_diff);
use CoGe::Accessory::Web;

no warnings 'redefine';

our ( $P, $PAGE_TITLE, $USER, $BASEFILE, $coge, %FUNCTION, $FORM, $JEX,
    $LINK );

$PAGE_TITLE = 'Jobs';

$FORM = new CGI;
( $coge, $USER, $P, $LINK ) = CoGe::Accessory::Web->init(
    cgi        => $FORM,
    page_title => $PAGE_TITLE
);

$JEX =
  CoGe::Accessory::Jex->new( host => $P->{JOBSERVER}, port => $P->{JOBPORT} );

%FUNCTION = (
    cancel_job   => \&cancel_job,
    restart_job  => \&restart_job,
    get_jobs     => \&get_jobs_for_user,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub get_jobs_for_user {
    my @entries;

    if ( $USER->is_admin ) {
        @entries = $coge->resultset('Log')->search(
            #{ description => { 'not like' => 'page access' } },
            {
                type     => { '!='  => 0 },
                workflow_id  => { "!=" => undef}
            },
            { order_by => { -desc => 'time' } },
        );
    }
    else {
        @entries = $coge->resultset('Log')->search(
            {
                user_id => $USER->id,
                workflow_id  => { "!=" => undef},

                #{ description => { 'not like' => 'page access' } }
                type => { '!=' => 0 }
            },
            { order_by => { -desc => 'time' } },
        );
    }

    my %users = map { $_->user_id => $_->name } $coge->resultset('User')->all;
    my @workflows = map { $_->workflow_id } @entries;
    my $workflows = $JEX->find_workflows(($USER->is_admin ? undef : \@workflows));

    my @job_items;
    my %workflow_results;

    foreach (@{$workflows}) {
        my($id, $name, $submitted, $completed, $status) = @{$_};

        my $start_time = localtime($submitted)->strftime('%F %I:%M%P');
        my $end_time = "";
        my $diff;

        if ($completed) {
            $end_time = localtime($completed)->strftime('%F %I:%M%P') if $completed;
            $diff = $completed - $submitted;
        } else {
            $diff = time - $submitted;
        }

        $workflow_results{$id} = {
            status    => $status,
            started   => $start_time,
            completed => $end_time,
            elapsed   => format_time_diff($diff)
        };
    }

    my $index = 1;
    foreach (@entries) {
        my $entry = $workflow_results{$_->workflow_id};

        # A log entry must correspond to a workflow
        next unless $entry;

        push @job_items, {
            id => int($index++),
            workflow_id => $_->workflow_id,
            user  => $users{$_->user_id} || "public",
            tool  => $_->page,
            link  => $_->link,
            %{$entry}
        };
    }

    my @filtered;

    # Filter repeated entries
    foreach (reverse @job_items) {
        my $wid = $_->{workflow_id};
        next if (defined $wid and defined $workflow_results{$wid}{seen});
        $workflow_results{$wid}{seen}++ if (defined $wid);

        unshift @filtered, $_;
    }

    return encode_json({ jobs => \@filtered });
}

sub gen_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER => $name );
    $template->param(
        TITLE      => qq{Jobs},
        PAGE_TITLE => $PAGE_TITLE,
        PAGE_LINK  => $LINK,
        HOME       => $P->{SERVER},
        HELP       => 'Jobs',
        WIKI_URL   => $P->{WIKI_URL} || ''
    );
    $template->param( LOGON => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY => gen_body() );
    $template->param( ADMIN_ONLY => $USER->is_admin );
    $template->param( CAS_URL    => $P->{CAS_URL} || '' );

    $html .= $template->output;
}

sub gen_body {
    unless ( $USER->is_admin ) {
        my $template =
          HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
        $template->param( ADMIN_ONLY => 1 );
        return $template->output;
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( PAGE_NAME  => "$PAGE_TITLE.pl" );
    $template->param( MAIN       => 1 );

    return $template->output;
}

sub cancel_job {
    my $job_id = _check_job_args(@_);

    return encode_json( {} ) unless defined($job_id);

    my $status = $JEX->get_status( $job_id );

    if ( $status =~ /scheduled|running|notfound/i ) {
        return encode_json({ status => $JEX->terminate( $job_id ) });
    } else {
        return encode_json( {} );
    }
}

sub restart_job {
    my $job_id = _check_job_args(@_);

    return encode_json( {} ) unless defined($job_id);

    my $status = $JEX->get_status( $job_id );

    if ( $status =~ /running/i ) {
        return encode_json( {} );
    } else {
        return encode_json({ status => $JEX->restart( $job_id ) });
    }
}

sub cmp_by_start_time {
    my $job1 = shift;
    my $job2 = shift;

    $job1->start_time cmp $job2->start_time;
}

sub _check_job_args {
    my %args   = @_;
    my $job_id = $args{job};

    if ( not defined($job_id) ) {
        say STDERR "Job.pl: a job id was not given to cancel_job.";
    }

    return $job_id;
}
