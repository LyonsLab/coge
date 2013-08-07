#! /usr/bin/perl -w
use v5.10;
use strict;

use CGI;
use Digest::MD5 qw(md5_base64);

use HTML::Template;
use JSON::XS;

# CoGe packages
use CoGeX;
use CoGe::Accessory::Jex;
use CoGe::Accessory::Web;

no warnings 'redefine';

our ( $P, $PAGE_TITLE, $USER, $BASEFILE, $coge, %FUNCTION, $FORM, $YERBA );

$PAGE_TITLE = 'Jobs';
$FORM       = new CGI;

( $coge, $USER, $P ) = CoGe::Accessory::Web->init(
    ticket     => $FORM->param('ticket'),
    url        => $FORM->url,
    page_title => $PAGE_TITLE
);

$YERBA = CoGe::Accessory::Jex->new( host => "localhost", port => 5151 );

%FUNCTION = (
    cancel_job   => \&cancel_job,
    schedule_job => \&schedule_job,
);

CoGe::Accessory::Web->dispatch( $FORM, \%FUNCTION, \&gen_html );

sub get_jobs_for_user {

    #my %opts = @_;
    my @jobs;

    if ( $USER->is_admin ) {
        @jobs =
          $coge->resultset('Job')
          ->search( undef, { order_by => 'job_id DESC' } );
    }
    elsif ( $USER->is_public ) {
        @jobs =
          $coge->resultset('Job')
          ->search( { user_id => 0 }, { order_by => 'job_id ASC', } );
    }
    else {
        @jobs = $USER->jobs->search( undef, { order_by => 'job_id DESC' } );
    }

    my @job_items;

    foreach my $job (@jobs) {
        push @job_items, {
            ID     => $job->job_id,
            LINK   => $job->link,
            PAGE   => $job->page,
            STATUS => get_status_message($job),

            #TODO: Should return the job duration
            RUNTIME => $job->start_time,
        };
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( LIST_STUFF => 1 );
    $template->param( JOB_LOOP   => \@job_items );

    return $template->output;
}

sub gen_html {
    my $html;
    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . 'generic_page.tmpl' );
    $template->param( HELP => "/wiki/index.php?title=$PAGE_TITLE" );
    my $name = $USER->user_name;
    $name = $USER->first_name if $USER->first_name;
    $name .= " " . $USER->last_name if $USER->first_name && $USER->last_name;
    $template->param( USER       => $name );
    $template->param( TITLE      => qq{} );
    $template->param( PAGE_TITLE => $PAGE_TITLE );
    $template->param( LOGO_PNG   => "$PAGE_TITLE-logo.png" );
    $template->param( LOGON      => 1 ) unless $USER->user_name eq "public";
    $template->param( BODY       => gen_body() );

    #	$name .= $name =~ /s$/ ? "'" : "'s";
    #	$template->param( BOX_NAME   => $name . " Data Lists:" );
    $template->param( ADJUST_BOX => 1 );
    $html .= $template->output;
}

sub gen_body {
    if ( $USER->is_public ) {
        my $template =
          HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
        $template->param( PAGE_NAME => "$PAGE_TITLE.pl" );
        $template->param( LOGIN     => 1 );
        return $template->output;
    }

    my $template =
      HTML::Template->new( filename => $P->{TMPLDIR} . "$PAGE_TITLE.tmpl" );
    $template->param( PAGE_NAME  => "$PAGE_TITLE.pl" );
    $template->param( MAIN       => 1 );
    $template->param( LIST_INFO  => get_jobs_for_user() );
    $template->param( ADMIN_AREA => 1 ) if $USER->is_admin;
    return $template->output;
}

sub cancel_job {
    my $job_id = _check_job_args(@_);
    my $job    = _get_validated_job($job_id);

    return return encode_json( {} ) unless defined($job);

    my $status = $YERBA->get_status( $job->id );
    say STDERR $status;

    if ( lc($status) eq 'running' ) {
        $job->update( { status => 3 } );
        return encode_json( $YERBA->terminate( $job->id ) );
    }
    else {
        return encode_json( {} );
    }

}

sub schedule_job {
    my $job_id = _check_job_args(@_);
    my $job    = _get_validated_job($job_id);

    return "fail" unless defined($job);
    return "true";
}

sub get_status_message {
    my $job = shift;

    given ( $job->status ) {
        when (1) { return 'Running'; }
        when (2) { return 'Complete'; }
        when (3) { return 'Cancelled'; }
        when (4) { return 'Terminated'; }
        default  { return 'Running'; }
    }
}

sub cmp_by_start_time {
    my $job1 = shift;
    my $job2 = shift;

    $job1->start_time cmp $job2->start_time;
}

# private functions

sub _get_validated_job {
    my $job_id = shift;
    my $job    = $coge->resultset('Job')->find($job_id);

    if ( ( not defined($job) || $job->user_id == $USER->id )
        && not $USER->is_admin )
    {
        say STDERR "Job.pl: job $job->id expected user id "
          . "$job->user_id but received $USER->id";
        return;
    }

    return $job;
}

sub _check_job_args {
    my %args   = @_;
    my $job_id = $args{job};

    if ( not defined($job_id) ) {
        say STDERR "Job.pl: a job id was not given to cancel_job.";
    }

    return $job_id;
}
