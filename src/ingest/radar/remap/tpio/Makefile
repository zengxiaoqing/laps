#
#	Build libtpio.a
#

CC=	cc
CFLAGS=	-DSUN -DRMT_IO_ALLOWED -DDISK_IO_ALLOWED
AR=	/usr/5bin/ar

OBJ=	newtpcfor.o newtplib.o rtapelib.o

CFILES=	newtpcfor.c newtplib.c rtapelib.c
HFILES=	newtpcfor.h newtpdefs.h newtplib.h rmt.h rmtdefs.h
MFILES=	Makefile_sun Makefile_hp

DIST=	tpio.tar

all:	libtpio.a

libtpio.a:	$(OBJ)
	$(AR) rcv libtpio.a $(OBJ)

dist:
	tar cf $(DIST) $(CFILES) $(HFILES) $(MFILES)

clean:
	/bin/rm -rf $(OBJ) $(DIST) libtpio.a core
