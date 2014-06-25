package CoGe::Pipelines::Misc::Copy;

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
    @EXPORT = qw (copy_and_mask);
}

sub copy_and_mask {
    my %args = (
        mask => 0,
        sequence_only => 0,
        @_
    );

    my $cmd = "/copy_genome/copy_load_mask_genome.pl";

    return (
        cmd   => catfile($args{script_dir}, $cmd),
        args  => [
            ["-conf", $args{conf}, 0],
            ["-gid", $args{gid}, 0],
            ["-uid", $args{uid}, 0],
            ["-mask", $args{mask}, 0],
            ["-staging_dir", $args{staging_dir}, 0],
            ["-result_dir", $args{result_dir}, 0],
            ["-sequence_only", $args{sequence_only}, 0]
        ],
        description => "Copy and masking genome..."
    );
}

1;
