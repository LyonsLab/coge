package CoGe::Services::Error;

BEGIN {
    require Exporter;

    $VERSION = 0.1;
    @ISA     = qw(Exporter);
    @EXPORT = qw(
      API_STATUS_UNAUTHORIZED API_STATUS_NOTFOUND
    );
}

sub API_STATUS_UNAUTHORIZED {
    return (
        status => 401,
        json   => {
            error => { Auth => "Access denied" }
        }
    );
}

sub API_STATUS_NOTFOUND {
    return (
        status => 400,
        json   => {
            error => { Error => "Not found" }
        }
    );
}

1;