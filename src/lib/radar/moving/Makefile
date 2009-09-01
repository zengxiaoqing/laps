# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..


include $(SRCROOT)/src/include/makefile.inc

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a
.f90.o:;    $(FC) $(DBFLAGS) $(FREE) -c $*.f90

LIB1      = libairad.a
DEBUGLIB1 = libdebugairad.a

SRC1= wrapper.f90

ARCHIVE1 = $(LIB1)($(addsuffix .o, $(basename $(SRC1))))
DBARCHIVE1 = $(DEBUGLIB1)($(addsuffix .o, $(basename $(SRC1))))

all: lib

lib: $(LIB1)


$(LIB1): $(ARCHIVE1) 
	$(RANLIB) $@

debug: $(DBARCHIVE1)
	$(MAKE) lib LIB='$(DEBUGLIB1)' ARCHIVE='$(DBARCHIVE1)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: clean
