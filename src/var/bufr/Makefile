# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

#=============================================================
#  LAPS environment:
#=============================================================

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

RM = rm -f
FLAG = $(FFLAGS)

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90
.f.o:;      $(FC) $(FLAG) -c $*.f
.c.o:;	    $(CC) $(CPPFLAGS) -c $*.c


LIB      = libbufr.a
DEBUGLIB = libdebugbufr.a

all:  
	tar xvf BUFRLIB09.tar
	preproc.sh $(CC) cpp
	$(FC) $(FLAG) -c *.f
	$(CC) $(CPPFLAGS) -c *.c
	ar crv $(LIB) *.o
	rm -f *.f *.F *.c *.h *.sh *.prm README* *.o

debug: 
	tar xvf BUFRLIB09.tar
	preproc.sh $(CC) cpp
	$(FC) $(DBFLAGS) -c *.f
	$(CC) $(DBCFLAGS) -c *.c
	ar crv $(DEBUGLIB) *.o
	rm -f *.f *.F *.c *.h *.sh *.prm README* *.o

