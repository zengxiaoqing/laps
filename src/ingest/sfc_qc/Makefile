# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#

LAPSROOT=../../..


include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS=$(LAPSLIB) $(BGDATALIB) $(MTHERMOLIB)
DBLAPSLIBS=$(DBLAPSLIB) $(DBBGDATALIB) $(DBMTHERMOLIB)


EXE = sfc_qc.exe

SRC   = kalman.f    matrix_tools.f   \
        QC_tools.f       model_tools.f    write_qc.c    \
        write_tools.f     lapsvanl_sub.f  


MAIN=QC_main.f

LIB=libqc.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o


debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F













