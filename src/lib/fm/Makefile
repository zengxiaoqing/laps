# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)









.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libfm.a
DEBUGLIB = libdebugfm.a




SRC=britgo.f \
    gimrad.f \
    plango.f \
    taugim.f \
    zenith.f



ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

install: lib

$(LIB): $(ARCHIVE) 
	$(RANLIB) $(LIB)

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#


distclean: clean
