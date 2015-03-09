package CoGe::Pipelines::Misc::Bed;

use strict;
use warnings;

use File::Spec::Functions;

use CoGe::Core::Genome qw(get_download_path);

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );

    $VERSION = 0.01;
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw (generate_bed);
}

sub generate_bed {
    my %args = @_;

    # Check for a genome or dataset id
    return unless $args{gid};

    # Generate file name
    my $basename = $args{basename};
    my $filename = "$basename" . "_id" . $args{gid} . ".bed";
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd  => catfile($args{script_dir}, "coge2bed.pl"),
        args => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $args{conf}, 0],
        ],
        outputs => [$output_file]
    );
}

1;
