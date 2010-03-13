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

LIB      = libmthermo.a
DEBUGLIB = libdebugmthermo.a




SRC= alcl.f     es.f       esrw.f     os.f       ptlcl.f    thm.f      tv.f   \
     ct.f       esat.f     esw.f      ow.f       satlft.f   tlcl.f     tw.f   \
     dewpt.f    esgg.f     heatl.f    pccl.f     ssh.f      tlcl1.f    w.f    \
     dpt.f      esice.f    humidity.f pcon.f     tcon.f     tmlaps.f   wmr.f  \
     dwpt.f     esilo.f    o.f        powt.f     tda.f      tmr.f      wobf.f \
     ept.f      eslo.f     oe.f       precpw.f   te.f       tsa.f      z.f        tdew.f


ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(F77INC) $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#


distclean: clean
