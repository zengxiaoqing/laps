# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)

DRTUFLAG = -Drtu








.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = liblapack.a
DEBUGLIB = libdebuglapack.a




SRC=dgecon.f  dgetri.f  dlamch.f  dlaswp.f  dtrti2.f  lsame.f\
    dgetf2.f  dlabad.f  dlange.f  dlatrs.f  dtrtri.f  xerbla.f\
    dgetrf.f  dlacon.f  dlassq.f  drscl.f   ilaenv.f


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
