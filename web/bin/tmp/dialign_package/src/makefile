


    ###############################
    #                             #
    #                             #
    #    Makefile for DIALIGN     #
    #                             #
    #                             #
    ###############################


CC = gcc
CFLAGS =  -c -O -I$ -DCONS
#CFLAGS =  -g -c -I$ -DCONS
LIBS = -lm
#

OBJS = dialign.o functions.o  para.o \
	output.o wgt.o anchor.o regex.o \
        input.o frag_chain.o \
	alig_graph_closure.o pratique.o

#
dialign2-2: $(OBJS)  
	$(CC) $(OBJS) $(LIBS) -o dialign2-2
#	$(CC) -g $(OBJS) $(LIBS) -o dialign2-2_db
#
#
# Subroutines
#

dialign.o: dialign.c
	$(CC) $(CFLAGS) dialign.c
functions.o: functions.c 
	$(CC) $(CFLAGS) functions.c
input.o: input.c
	$(CC) $(CFLAGS) input.c
frag_chain.o: frag_chain.c
	$(CC) $(CFLAGS) frag_chain.c
para.o: para.c
	$(CC) $(CFLAGS) para.c
output.o: output.c
	$(CC) $(CFLAGS) output.c
wgt.o: wgt.c
	$(CC) $(CFLAGS) wgt.c
regex.o: regex.c
	$(CC) $(CFLAGS) regex.c
anchor.o: anchor.c
	$(CC) $(CFLAGS) anchor.c

#



