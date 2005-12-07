

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'Genome::DB',
    VERSION_FROM => 'lib/CoGe/Genome/DB.pm', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons@nature.berkeley.edu)',
    ABSTRACT     => 'Genome::DB',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
