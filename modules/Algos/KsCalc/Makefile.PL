use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Algos::KsCalc',
    VERSION_FROM => 'lib/CoGe/Algos/KsCalc.pm', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons.uoa@gmail.com)',
    ABSTRACT     => 'CoGe::Algos::KsCalc',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
    INSTALLMAN1DIR => 'none', # mdb added 11/19/15
);
