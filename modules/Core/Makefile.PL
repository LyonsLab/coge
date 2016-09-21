use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Core',
    VERSION_FROM => '', # finds \$VERSION
    AUTHOR       => '',
    ABSTRACT     => 'CoGe::Core',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
