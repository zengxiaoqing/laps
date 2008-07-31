# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#
SRCROOT=../../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/satellite/libsatellite.a \
           $(SRCROOT)/src/lib/goesinav/libgoesinav.a  \
           $(SRCROOT)/src/lib/util/libutil.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/satellite/libdebugsatellite.a \
             $(SRCROOT)/src/lib/goesinav/libdebuggoesinav.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lsr_driver.exe

SRC =  getsnddata.f \
       gen_gvrsndr_lut_lsr.f \
       rddata_line_elem_sub.f \
       read_sounding_db_cdf.f \
       satsndr2laps.f \
       set_missing_sndr.f \
       count_range.f \
       write_lsr.f

MAIN=lsr_driver.f

LIB=libgenlsr.a

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
















