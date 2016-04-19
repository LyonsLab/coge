use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
# http://search.cpan.org/~bingos/ExtUtils-MakeMaker-7.10/lib/ExtUtils/MakeMaker.pm
WriteMakefile(
    NAME         => 'CoGe',
    VERSION_FROM => '', # finds \$VERSION
    AUTHOR       => 'Eric Lyons (elyons@email.arizona.edu), et al.',
    INSTALLMAN1DIR => 'none', # mdb added 9/30/15
    INSTALLMAN3DIR => 'none', # mdb added 9/30/15
    #ABSTRACT     => '',
    #PREREQ_PM    => {}
);
