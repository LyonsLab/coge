#uncomment and locally adapt next line to fix the full helpfile path name
HELPFILENAME = -DHELPFILENAME=\"/usr/share/njplot/njplot.help\"

#comment out next line and uncomment next 2 to use the PDFLib Lite library 
NO_PDF = -DNO_PDF
#PDF = ../PDFlibLite/libs/pdflib
#PDFLIB = -L$(PDF) -lpdf

# c compiler and linker
CC = gcc

# Vibrant top directory
VIBRANT = /banques0/ncbiJun04


# X11 include directory
X11INCL = /usr/X11R6/include

# motif library directory
MOTIFLIB = /sw/lib
MOTIFINCL = /sw/include

OBJECTS = njplot-vib.o
OBJUNROOTED = unrooted-vib.o preptree.o


CFLAGS  = -c -DWIN_MOTIF -Dunix -I$(VIBRANT)/include -I$(VIBRANT)/vibrant -I$(VIBRANT)/corelib \
	-I$(VIBRANT) -I$(PDF) -I$(X11INCL) -I$(MOTIFINCL) $(HELPFILENAME) $(NO_PDF) \
	-Wimplicit-function-declaration -g
	
all: njplot  unrooted newicktops newicktotxt


njplot : $(OBJECTS)
	$(CC) -g -o njplot $(OBJECTS) \
        -L$(VIBRANT)/lib \
        -lvibrant -lncbi \
        $(PDFLIB)   \
        -L$(MOTIFLIB) -lXm  \
        -L/usr/X11R6/lib -lXmu -lXt -lX11 -lm
        
unrooted : $(OBJUNROOTED)
	$(CC) -g -o unrooted $(OBJUNROOTED) \
        -L$(VIBRANT)/lib -lvibrant -lncbi \
        -L$(MOTIFLIB) -lXm \
        -L/usr/X11R6/lib -lXmu -lXt -lX11  -lm

newicktops: njplot-vib.c
	$(CC) -DNO_GUI  -DNO_PDF -o $@ njplot-vib.c -lm

newicktotxt: njplot-vib.c
	$(CC) -DTTY  -o $@ njplot-vib.c -lm


.c.o : 
	$(CC) $(CFLAGS)  $?

clean:
	rm -f *.o

distclean: clean
	rm -f njplot unrooted newicktops newicktotxt
