# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)

DRTUFLAG = -Drtu








.SUFFIXES:
.SUFFIXES: .o .c .f .f90 .a

LIB      = librutil.a
DEBUGLIB = libdebugrutil.a




SRC= \
    kessler.f90 dbzcalc.f       

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#


distclean: clean


.f90.o:;    $(FC) $(DBFLAGS) $(FREE) $(STINC) -c $*.f90
.f.o:;      $(FC) $(DBFLAGS) -c $*.f
