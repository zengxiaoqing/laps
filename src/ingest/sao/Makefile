# Generated automatically from Makefile.in by configure.
#
#

LAPSROOT=../../..
INSTALLROOT=/usr/nfs/lapb/parallel/laps

# Specify Fortran, C, and loader options for default SUN configuration.

CC=xlc

CFLAGS=-g $(INC)  -O3 -qstrict

FC       = xlf


FFLAGS =  -qfixed $(INC)  -O3 -qstrict

# Define pre-processor command line options for .F files.

CPP = xlc -P
CPPFLAGS= -Drs6000 -DFORTRANNOUNDERSCORE -DDYNAMIC

# Define include paths.

INC = -I$(LAPSROOT)/src/include -I/usr/local/apps/netcdf-3.3.1/include

# Define local library name, this houses the module's code

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  



OTHERLIBS = /usr/local/apps/netcdf-3.3.1/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = obs_driver.x


# Define other variables 

AR  = ar
RANLIB   = ranlib
INSTALL=$(LAPSROOT)/util/install-sh

# Define all source code to be compiled and placed in local library.

SRC   = get_metar_obs.f \
	read_metar.f \
        get_local_obs.f \
        read_local.f \
        get_buoy_obs.f \
	read_buoy.f \
	write_surface_obs.f

MAIN=obs_driver.f



LIB=libsao.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

#$(EXE.F): $(MAIN) $(ARCHIVE)
#	rm -f $*.f
#	$(CPP) $(CPPFLAGS) $*.F 
#	mv $*.i $*.f
#	$(FC) $(FFLAGS) $*.f $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
#	rm -f *.o $*.f

debug: clean 
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' FFLAGS=' -qfixed $(INC) -g'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean
	rm -f Makefile


%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F










get_metar_obs.f:  $(LAPSROOT)/src/include/surface_obs.inc

obs_driver.f:  $(LAPSROOT)/src/include/surface_obs.inc

read_metar.f:  /usr/local/apps/netcdf-3.3.1/include/netcdf.inc

