use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::JEX',
    VERSION_FROM => '', # finds \$VERSION
    AUTHOR       => '',
    ABSTRACT     => 'CoGe::JEX',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
