package CoGe::Pipelines::Misc::Gff;

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
    @EXPORT = qw (generate_gff);
}


sub generate_gff {
    my %args = (
        annos   => 0,
        id_type => 0,
        cds     => 0,
        nu      => 0,
        upa     => 0,
        @_
    );

    # Check for a genome or dataset id
    return unless $args{gid};

    # Generate the output filename
    my $organism = sanitize_name($args{basename});
    my @attributes = qw(annos cds id_type nu upa);
    my $param_string = join "-", map { $_ . $args{$_} } @attributes;
    my $filename = $organism . "_" . $param_string . ".gff";
    my $path = get_download_path($args{secure_tmp}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd     => catfile($args{script_dir}, "coge_gff.pl"),
        args    => [
            ['-id', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $args{conf}, 0],
            # Parameters
            ['-cds', $args{cds}, 0],
            ['-annos', $args{annos}, 0],
            ['-nu', $args{nu}, 0],
            ['-id_type', $args{id_type}, 0],
            ['-upa', $args{upa}, 0],
        ],
        outputs => [$output_file]
    );
}

1;

