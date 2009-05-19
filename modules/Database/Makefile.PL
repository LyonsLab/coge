

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'CoGeX',
    VERSION_FROM => 'lib/CoGeX.pm', # finds \$VERSION
    AUTHOR       => ' (elyons@nature.berkeley.edu)',
    ABSTRACT     => 'CoGeX',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                      'DBIx::Class'  => 0.08

                    },
);
