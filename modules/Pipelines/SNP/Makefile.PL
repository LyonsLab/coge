use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGe::Pipelines::FindSNPs',
    VERSION_FROM => 'lib/CoGe/Pipelines/FindSNPs.pm', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons.uoa@gmail.com)',
    ABSTRACT     => 'CoGe::Pipelines::FindSNPs',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
