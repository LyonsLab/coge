
PROG=dagchainer
PYTHON_INCLUDE=/usr/include/python2.6

#CCFLAGS=-Wno-deprecated
#CCFLAGS=-O3 -s -fomit-frame-pointer -msse2 -Wno-deprecated
CCFLAGS=-O3 -s -fomit-frame-pointer -Wno-deprecated
all: python dagchainer

dagchainer: dagchainer.cpp
	g++ -o ${PROG} dagchainer.cpp ${CCFLAGS}

python: 
	cython dagtools/cdagline.pyx
	python setup.py build_ext -i
	rm -r build

pythongcc:
	cd dagtools && cython cdagline.pyx && \
		gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I${PYTHON_INCLUDE} -c cdagline.c -o cdagline.o && \
		gcc -pthread -shared -Wl,-O3 -Wl,-Bsymbolic-functions cdagline.o -o cdagline.so && \
		rm cdagline.o && rm cdagline.c 

