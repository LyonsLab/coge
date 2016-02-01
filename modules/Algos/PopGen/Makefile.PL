

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Algos::PopGen',
    VERSION_FROM => undef, # finds \$VERSION
    AUTHOR       => 'Eric Lyons ',
    ABSTRACT     => 'CoGe::Algos::PopGen',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
