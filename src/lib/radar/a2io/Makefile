# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#	A2IO Makefile --- build liba2tp.a and liba2rt.a
#
#	"@(#)Makefile_sun	1.2	2/21/95"
#

#CC=		acc $(FAST)

SRCROOT=../../../..

include $(SRCROOT)/src/include/makefile.inc

OBJA2TP=	a2io_tp.o precip_tp.o read88da2.o merge.o
OBJA2RT=	a2io_rt.o precip_rt.o read88da2.o merge.o

CFILES=		a2io.c precip.c read88da2.c merge.c
HFILES=		a2io.h defs.h vcp.h merge.h
MFILES=		Makefile_sun Makefile_hp
DFILES=		a2io.doc radarinfo.dat

DIST=		a2io.tar

LIBA2TP=	liba2tp.a
LIBA2RT=	liba2rt.a

ALL=		$(LIBA2RT) $(LIBA2TP)

TP_FLAGS=
RT_FLAGS=	-DRTS

all:		$(ALL)

a2io_tp.o:	a2io.c
		$(CC) $(TP_FLAGS) $(CPPFLAGS) $(CFLAGS) -c a2io.c
		mv a2io.o a2io_tp.o

a2io_rt.o:	a2io.c
		$(CC) $(RT_FLAGS) $(CPPFLAGS) $(CFLAGS) -c a2io.c
		mv a2io.o a2io_rt.o

precip_tp.o:	precip.c
		$(CC) $(TP_FLAGS) $(CPPFLAGS) $(CFLAGS) -c precip.c
		mv precip.o precip_tp.o

precip_rt.o:	precip.c
		$(CC) $(RT_FLAGS) $(CPPFLAGS) $(CFLAGS) -c precip.c
		mv precip.o precip_rt.o

read88da2.o:	read88da2.c
		$(CC) $(CPPFLAGS) $(CFLAGS) -c read88da2.c

merge.o:	merge.c
		$(CC) $(CPPFLAGS) $(CFLAGS) -c merge.c

$(LIBA2RT):	$(OBJA2RT)
		$(AR) rcv $(LIBA2RT) $(OBJA2RT)

$(LIBA2TP):	$(OBJA2TP)
		$(AR) rcv $(LIBA2TP) $(OBJA2TP)

clean:
		/bin/rm -rf $(OBJA2RT) $(OBJA2TP) $(ALL) $(DIST) core

dist:
		tar cf $(DIST) $(CFILES) $(HFILES) $(MFILES) $(DFILES)

depend:
		makedepend $(CPPFLAGS) $(CFILES)

# DO NOT DELETE THIS LINE -- make depend depends on it.

a2io.o: a2io.h defs.h vcp.h
precip.o: a2io.h
read88da2.o: merge.h
merge.o: defs.h merge.h
