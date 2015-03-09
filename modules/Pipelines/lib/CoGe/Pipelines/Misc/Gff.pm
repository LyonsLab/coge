package CoGe::Pipelines::Misc::Gff;

use strict;
use warnings;

use File::Basename qw(basename);
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
    my ($inputs, $conf) = @_;

    my %args = (
        annos   => 0,
        id_type => 0,
        cds     => 0,
        nu      => 0,
        upa     => 0,
    );

    @args{(keys $inputs)} = (values $inputs);

    # Check for a genome or dataset id
    return unless $args{gid};

    # Set the default basename as the id if the basename is not set
    $args{basename} = $args{gid} unless $args{basename};

    # Generate the output filename
    my $organism = "gff";
    my @attributes = qw(annos cds id_type nu upa);
    my $param_string = join "-", map { $_ . $args{$_} } @attributes;
    my $filename = $args{basename} . "_" . $param_string . ".gff";
    $filename =~ s/\s+/_/g;
    $filename =~ s/\)|\(/_/g;
    my $path = get_download_path($conf->{SECTEMPDIR}, $args{gid});
    my $output_file = catfile($path, $filename);

    return $output_file, (
        cmd     => catfile($conf->{SCRIPTDIR}, "coge_gff.pl"),
        args    => [
            ['-gid', $args{gid}, 0],
            ['-f', $filename, 0],
            ['-config', $conf->{_CONFIG_PATH}, 0],
            # Parameters
            ['-cds', $args{cds}, 0],
            ['-annos', $args{annos}, 0],
            ['-nu', $args{nu}, 0],
            ['-id_type', $args{id_type}, 0],
            ['-upa', $args{upa}, 0],
        ],
        outputs => [$output_file],
        description => "Generating gff..."
    );
}
1;
