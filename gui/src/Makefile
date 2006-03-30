# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

ALLLIBS = $(SRCROOT)/src/lib/liblaps.a /usr/lib/libm.a
DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a

RM=rm -f

.SUFFIXES: .o .c .a .exe

LIB = libcgui.a
DEBUGLIB = libdebugcgui.a
ARCHIVE = $(AR) rv
.c.o:
			$(RM) $@
			$(CC) -c $(CFLAGS) $(BYTEFLAG) $(INC) $*.c
			$(ARCHIVE) $(LIB) $(<:.c=.o)

EXE1 = gen_map_bkgnd.exe
EXE2 = pwrap_ll_xy_convert.exe

SRC = create_bcd_bkgnd_xy.c \
ll_to_xy.c \
location_lines.c \
setsup.c \
xy_to_ll.c 

MAIN1 = gen_map_bkgnd.c
MAIN2 = pwrap_ll_xy_convert.c

OBJS=$(SRC:.c=.o)

all:	$(EXE1) $(EXE2)

$(EXE1):$(OBJS)
	$(CC) -o $(EXE1) $(LDFLAGS) $(MAIN1) -I$(SRCROOT)/src/include $(LIB) $(LAPSLIBS) -lm; 
	rm -f *.o

$(EXE2):$(OBJS)
	$(CC) -o $(EXE2) $(LDFLAGS) $(MAIN2) -I$(SRCROOT)/src/include $(LIB) $(LAPSLIBS) -lm; 
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' LIB='$(DEBUGLIB)' \
	CFLAGS='$(DBCFLAGS)'

install: all
	 $(INSTALL) $(EXE1) $(INSTALLROOT)/bin/$(EXE1);
	 $(INSTALL) $(EXE2) $(INSTALLROOT)/bin/$(EXE2);

clean:
	rm -f $(EXE1) $(EXE2)  *.o *.a *~ *\#

distclean: clean
