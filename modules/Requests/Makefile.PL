use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Requests',
    VERSION_FROM => '', # finds \$VERSION
    AUTHOR       => 'Evan Briones <evrick@email.arizona.edu>',
    ABSTRACT     => 'CoGe::Requests',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
