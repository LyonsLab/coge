#cd ./modules; perl Makefile.PL lib=./perl; make install; perl Makefile.PL lib=/usr/local/lib/perl/5.18.2/; make install;
#cd ./modules; perl Makefile.PL lib=./perl; make install; 
cd ./modules; perl Makefile.PL INSTALL_BASE=./perl; make install; 
 
