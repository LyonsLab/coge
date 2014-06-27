package CoGe::Pipelines::Misc::Tbl;

use strict;
use warnings;

use File::Spec::Functions;

use CoGe::Accessory::Utils;
use CoGe::Core::Genome qw(get_download_path);

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );

    $VERSION = 0.01;
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw (generate_tbl);
}

sub generate_tbl {
    my %args = @_;

    # Generate filename
    my $organism = $args{basename};
    my $filename = $organism . "id" . $args{gid} . "_tbl.txt";
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd     => catfile($args{script_dir}, "export_NCBI_TBL.pl"),
        args    => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ["-config", $args{conf}, 0]
        ],
        outputs => [$output_file]
    );
}

1;
