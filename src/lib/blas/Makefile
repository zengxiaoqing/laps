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

LIB      = libblas.a
DEBUGLIB = libdebugblas.a




SRC= dasum.f \
daxpy.f \
dcopy.f \
ddot.f \
dgemm.f \
dgemv.f \
dger.f \
dscal.f \
dswap.f \
dtrmm.f \
dtrmv.f \
dtrsm.f \
dtrsv.f \
idamax.f



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
