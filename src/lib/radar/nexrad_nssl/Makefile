# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..


include $(SRCROOT)/src/include/makefile.inc

.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB1      = libnexrad.a
DEBUGLIB1 = libdebugnexrad.a

SRC1= wideband.c     wsr88dutil.c     getconfig.c     ip_utility.c

LIB2      = libnssl.a
DEBUGLIB2 = libdebugnssl.a

SRC2=  	c_buff.c 	\
	checksum.c	\
	config.c   	\
	initbuff.c	\
	mess.c		\
	messages.c	\
	strtime.c	\
	strtoken.c	\
	strutil.c	


ARCHIVE1 = $(LIB1)($(addsuffix .o, $(basename $(SRC1))))
ARCHIVE2 = $(LIB2)($(addsuffix .o, $(basename $(SRC2))))
DBARCHIVE1 = $(DEBUGLIB1)($(addsuffix .o, $(basename $(SRC1))))
DBARCHIVE2 = $(DEBUGLIB2)($(addsuffix .o, $(basename $(SRC2))))


all: lib

lib: $(LIB1) $(LIB2) 


$(LIB1): $(ARCHIVE1) 
	$(RANLIB) $@

$(LIB2): $(ARCHIVE2) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: clean

FORCE:
