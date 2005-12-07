

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::ECNCS',
    VERSION_FROM => 'lib/CoGe/ECNCS.pm', # finds \$VERSION
    AUTHOR       => 'HASH(0x813d9e0) (a.u.thor@a.galaxy.far.far.away)',
    ABSTRACT     => 'CoGe::ECNCS',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
