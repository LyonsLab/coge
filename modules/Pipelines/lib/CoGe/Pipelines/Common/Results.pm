package CoGe::Pipelines::Common::Results;

use strict;
use warnings;

use File::Basename qw(basename);
use File::Spec::Functions;
use URI::Escape::JavaScript qw(escape);

use CoGe::Accessory::Utils;
use CoGe::Core::Genome qw(get_download_path);

BEGIN {
    our ( @EXPORT, @ISA, $VERSION );

    $VERSION = 0.01;
    require Exporter;

    $VERSION = 0.1;
    @ISA = qw ( Exporter );
    @EXPORT = qw (generate_results link_results);
}

sub link_results {
   my ($input, $output, $result_dir, $conf) = @_;

   return (
        cmd     => catfile($conf->{SCRIPTDIR}, "link_results.pl"),
        args    => [
            ['-input_files', escape($input), 0],
            ['-output_files', escape($output), 0],
            ['-result_dir', $result_dir, 0]
        ],
        inputs  => [$input],
        outputs => [catfile($result_dir, basename($output))],
        description => "Generating results..."
   );
}

sub generate_results {
   my ($input, $type, $result_dir, $conf, $dependency) = @_;

   return (
        cmd     => catfile($conf->{SCRIPTDIR}, "generate_results.pl"),
        args    => [
            ['-input_files', escape($input), 0],
            ['-type', $type, 0],
            ['-result_dir', $result_dir, 0]
        ],
        inputs  => [$dependency],
        outputs => [catfile($result_dir, "1")],
        description => "Generating results..."
   );
}

1;
