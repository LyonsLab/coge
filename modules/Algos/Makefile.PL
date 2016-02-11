

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Algos',
    VERSION_FROM => undef, # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons@nature.berkeley.edu)',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
    INSTALLMAN1DIR => 'none', # mdb added 11/19/15
);
