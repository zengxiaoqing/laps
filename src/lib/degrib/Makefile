# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

# g2lib is needed for rd_grib1.F90 only
INCLUDIR =  $(SRCROOT)/src/lib/g2lib

# CDEGRIBFLAGS are -DLINUX, -DIBM, etc; -DSWAPBYTE; etc; 
# DEGRIBFLAGS are USE_PNG, USE_JPEG2000 only;
CFLAGS += $(CDEGRIBFLAGS) $(DEGRIBFLAGS)
FFLAGS += $(FREE) $(DEGRIBFLAGS) -I$(INCLUDIR)
CPPFLAGS += $(CDEGRIBFLAGS) $(DEGRIBFLAGS) -I$(INCLUDIR)

LIB      = libdegrib.a
DEBUGLIB = libdebugdegrib.a

.SUFFIXES: 
.SUFFIXES: .o .c .f .F90 .a .i .exe

#
# machine_dependent needs to be first
#

SRC	= \
          table.F90 \
          build_hdate.F90 \
          cio.c \
          debug_cio.c \
          file_delete.F90 \
          filelist.F90 \
          gbytesys.F90 \
          geth_idts.F90 \
          geth_newdate.F90 \
          gridinfo.F90 \
          module_datarray.F90 \
          module_debug.F90 \
          module_stringutil.F90 \
          new_storage.F90 \
          output.F90 \
          parse_table.F90 \
          rd_grib1.F90 \
          swap.F90 \
          readgrib.F90 \
          datint.F90 \
          gribcode.F90 \
          rd_grib2.F90 


ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS) $(CDEGRIBFLAGS) $(DEGRIBFLAGS)' \
	FFLAGS='$(DBFLAGS) $(FREE) $(DEGRIBFLAGS) -I$(INCLUDIR)'

install: lib

  
ifdef USECPP 
%.o : 	%.F90
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f90
	$(FC) $(FFLAGS) -c $*.f90
	$(RM) $*.f90
else

%.o : 	%.F90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $< 
endif

cleanup:
	rm -f *.o *.mod *.a

clean:
	$(RM) $(LIB) $(DEBUGLIB) *.mod *.o *~ *# 

distclean: clean
