package CoGe::Services::Error;

BEGIN {
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw(Exporter);
    @EXPORT  = qw(
      API_STATUS_UNAUTHORIZED API_STATUS_NOTFOUND API_STATUS_SEARCHTERM
      API_STATUS_MISSING_ID API_STATUS_MISSING_BODY API_STATUS_MISSING_DATA
      API_STATUS_BAD_REQUEST API_STATUS_CUSTOM
    );
}

sub API_STATUS_UNAUTHORIZED {
    return (
        status => 401,
        json   => {
            error => { message => "Access denied" }
        }
    );
}

sub API_STATUS_NOTFOUND {
    return (
        status => 404,
        json   => {
            error => { message => "Resource not found" }
        }
    );
}

sub API_STATUS_SEARCHTERM {
    return (
        status => 400,
        json   => {
            error => { message => 'Search term is shorter than 3 characters' }
        }
    );
}

sub API_STATUS_MISSING_ID {
    return (
        status => 400,
        json   => {
            error => { message => 'Missing ID parameter' }
        }
    );
}

sub API_STATUS_MISSING_BODY {
    return (
        status => 400,
        json   => {
            error => { message => 'No request body specified' }
        }
    );
}

sub API_STATUS_MISSING_DATA {
    return (
        status => 400,
        json   => {
            error => { message => 'No data items specified' }
        }
    );
}

sub API_STATUS_BAD_REQUEST {
    my $message = shift;
    return (
        status => 400,
        json   => {
            error => { message => $message }
        }
    );
}

sub API_STATUS_CUSTOM {
    my ($status_code, $message) = @_;
    return (
        status => $status_code,
        json   => {
            error => { message => $message }
        }
    );
}

1;