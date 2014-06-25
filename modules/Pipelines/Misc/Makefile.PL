use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Pipelines::Misc',
    VERSION_FROM => 'lib/CoGe/Pipelines/Misc/Bed.pm', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons.uoa@gmail.com)',
    ABSTRACT     => 'CoGe::Pipelines::Misc',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
