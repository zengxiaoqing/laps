SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

LIB = libgrib.a
DEBUGLIB = libdebuggrib.a

SRC  = 	fi7501.f \
        fi7502.f \
        fi7503.f \
        fi7505.f \
        fi7513.f \
        fi7516.f \
        fi7517.f \
        fi7518.f \
        w3fi01.f \
        w3fi58.F \
        w3fi59.F \
        w3fi68.f \
        w3fi71.f \
        w3fi72.F \
        w3fi73.f \
        w3fi74.f \
        w3fi75.f \
        w3fi76.f \
        w3fi82.f \
        gbyte.f \
        sbyte.f \
        sbytes.f \
        xmovex.f \
        xstore.f \
        io_utils.c

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

all:	lib

lib:	$(LIB)

$(LIB):	$(ARCHIVE)
	$(RANLIB) $@

debug: 	
	$(MAKE) all LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
     	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: 	lib

clean:
	$(RM) $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: 	clean

ifdef USECPP
%.o : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $<
	mv $*.i $*.f
	$(FC) $(FFLAGS) -c $*.f
endif

