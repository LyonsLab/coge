use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Services',
    VERSION_FROM => '', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons.uoa@gmail.com)',
    ABSTRACT     => 'CoGe::Services',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
