# Generated automatically from Makefile.in by configure.
#
#

LAPSROOT=../../../..
INSTALLROOT=/usr/nfs/lapb/parallel/laps

# Specify Fortran, C, and loader options for default SUN configuration.

CC=xlc

CFLAGS=-g $(INC)  -O3 -qstrict

FC       = xlf


FFLAGS =  -qfixed $(INC)  -O3 -qstrict

# Define pre-processor command line options for .F files.

CPP = xlc -P

# Define include paths.

INC = -I$(LAPSROOT)/src/include -I/usr/nfs/netcdf/include

# Define local library name, this houses the module's code

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  \
           $(LAPSROOT)/src/lib/radar/remap_ftn/libremap_ftn.a  \
           $(LAPSROOT)/src/lib/util/libutil.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(LAPSROOT)/src/lib/radar/remap_ftn/libdebugremap_ftn.a  \
             $(LAPSROOT)/src/lib/util/libdebugutil.a  


OTHERLIBS = /usr/nfs/netcdf/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = remap_polar_netcdf.exe


# Define other variables 

AR  = ar
RANLIB   = ranlib
INSTALL = $(LAPSROOT)/util/install-sh

# Define all source code to be compiled and placed in local library.

SRC=remap.f netcdfio.f


all: exe

exe: $(EXE)

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	rm -f $(EXE) *.o;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' FFLAGS=' -qfixed $(INC) -g'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean
	rm -f Makefile







