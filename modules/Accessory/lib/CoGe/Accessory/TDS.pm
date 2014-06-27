package CoGe::Accessory::TDS;

=head1 NAME

CoGe::Accessory::TDS

=head1 SYNOPSIS

Tiny Document Store

=head1 DESCRIPTION

Provide simple serialized access to JSON documents.

=head1 AUTHOR

Matt Bomhoff

=head1 COPYRIGHT

The full text of the license can be found in the
LICENSE file included with this module.

=head1 SEE ALSO

=cut

use strict;
use warnings;

use File::Basename qw(dirname);
use File::Path qw(mkpath);
use JSON::XS;

BEGIN {
    use vars qw ($VERSION $MAX_DOCUMENT_SZ @ISA @EXPORT);
    require Exporter;

    $VERSION = 0.1;
    $MAX_DOCUMENT_SZ = 1024*1024*1024; # 1G
    @ISA     = qw( Exporter );
    @EXPORT  = qw( read write );
}

sub read {
    my $filepath = shift;
    my $pData;

    open(my $fh, $filepath) || return;
    flock($fh, 1); # shared lock
    my $json = <$fh>;
    $pData = decode_json($json);
    close($fh);

    return $pData;
}

sub write {
    my ($filepath, $pData) = @_;
    mkpath(dirname($filepath), 0, 0777);

    my $json = encode_json($pData);
    return if (length($json) > $MAX_DOCUMENT_SZ);

    open(my $fh, ">", $filepath) || return;
    flock($fh, 2); # exclusive lock
    seek($fh, 0, 0); truncate(MYFILE, 0); # clear file
    print $fh $json;
    close($fh);

    return 1;
}

1;
