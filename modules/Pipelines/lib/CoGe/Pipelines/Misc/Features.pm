package CoGe::Pipelines::Misc::Features;

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
    @EXPORT = qw (generate_features);
}

sub generate_features {
    my %args = @_;

    my $filename = $args{basename} . "-gid-" . $args{gid};
    $filename .= "-prot" if $args{protein};
    $filename .= ".fasta";
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd    => catfile($args{script_dir}, "export_features_by_type.pl"),
        args   => [
            ["-config", $args{conf}, 0],
            ["-f", $filename, 0],
            ["-gid", $args{gid}, 0],
            ["-ftid", $args{fid}, 0],
            ["-prot", $args{protein}, 0],
        ],
        outputs => [$output_file]
    );
}

1;
