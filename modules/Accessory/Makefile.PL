

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Accessory::Tile::Cache',
    VERSION_FROM => '', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons@nature.berkeley.edu)',
    ABSTRACT     => 'CoGe::Accessory::Tile::Cache',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
