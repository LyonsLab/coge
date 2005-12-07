

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Algos::Pairwise',
    VERSION_FROM => 'lib/CoGe/Algos/Pairwise.pm', # finds \$VERSION
    AUTHOR       => 'Eric Lyons ',
    ABSTRACT     => 'CoGe::Algos::Pairwise',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
