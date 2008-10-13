#dis   
#dis    Open Source License/Disclaimer, Forecast Systems Laboratory
#dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
#dis    
#dis    This software is distributed under the Open Source Definition,
#dis    which may be found at http://www.opensource.org/osd.html.
#dis    
#dis    In particular, redistribution and use in source and binary forms,
#dis    with or without modification, are permitted provided that the
#dis    following conditions are met:
#dis    
#dis    - Redistributions of source code must retain this notice, this
#dis    list of conditions and the following disclaimer.
#dis    
#dis    - Redistributions in binary form must provide access to this
#dis    notice, this list of conditions and the following disclaimer, and
#dis    the underlying source code.
#dis    
#dis    - All modifications to this software must be clearly documented,
#dis    and are solely the responsibility of the agent making the
#dis    modifications.
#dis    
#dis    - If significant modifications or enhancements are made to this
#dis    software, the FSL Software Policy Manager
#dis    (softwaremgr@fsl.noaa.gov) should be notified.
#dis    
#dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
#dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
#dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
#dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
#dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
#dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
#dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
#dis   

SRCROOT=../..
EXTMODS=-L$(SRCROOT)/src/lib/modules -lmodules
DBEXTMODS=-L$(SRCROOT)/src/lib/modules -ldebugmodules
LAPSLIB=-L$(SRCROOT)/src/lib llaps
DBLAPSLIB=-L$(SRCROOT)/src/lib -ldebuglaps
GRIBLIB=-L$(SRCROOT)/src/lib/grib -lgrib
DBGRIBLIB=-L$(SRCROOT)/src/lib/grib -ldebuggrib
GRIB2LIB=-L$(SRCROOT)/src/lib/grib2 -lgrib2
DBGRIB2LIB=-L$(SRCROOT)/src/lib/grib2 -ldebuggrib2
NETCDF=/opt/netcdf/3.6.2-pgi-7.1-3
#OTHERLIBS=-L$(NETCDF)/lib -lnetcdf 

LIBS=$(SRCROOT)/src/lib/modules/libmodules.a $(SRCROOT)/src/lib/liblaps.a $(SRCROOT)/src/lib/grib/libgrib.a $(NETCDF)/lib/libnetcdf.a

include $(SRCROOT)/src/include/makefile.inc

FFLAGS=$(DBFLAGS)

RM=rm -f 

.SUFFIXES:              .f90 .o

.f90.o:
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90

.f90.mod:
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90

.mod.o:
#	
EXE=lfmpost.exe
FSRC=lfmpost.f90 lfmutil.f90 lfmgrid.f90 interp.f90 metfunc.f90 constants.f90  \
     mm5util.f90 wrfutil.f90 nmmutil.f90 st4util.f90  \
     ptfcst.f90 lfmwrite.f90

CSRC = 

FOBJS=$(FSRC:.f90=.o)

COBJS=$(CSRC:.c=.o)

all:			$(EXE)

# Note: if optimization is needed then change the DBFLAGS to FFLAGS in $(EXE)

$(EXE):			$(FOBJS) $(COBJS)
			$(FC) $(FREE) -o $@ $(FFLAGS) $(FOBJS) $(COBJS) \
			$(EXTMODS) $(GRIBLIB) $(GRIB2LIB) $(LAPSLIB) $(HDF5LIB) $(LIBS) $(INC) -lz -lm

debug:	
	 $(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
	"LAPSLIB	=	$(DBLAPSLIB)" 	\
        "GRIBLIB        =       $(DBGRIBLIB)"    \
        "GRIB2LIB       =       $(DBGRIB2LIB)"    \
	"FFLAGS 	=	$(DBFLAGS)"    \
	"EXTMODS        =       $(DBEXTMODS)" 	 \
        "LIBS           =       $(NETCDF)/lib/libnetcdf.a 

lfmpost.o:		lfmgrid.o 

lfmutil.o:		lfmgrid.o constants.o

interp.o:		lfmgrid.o

mm5util.o:		lfmgrid.o 

wrfutil.o:		lfmgrid.o constants.o

nmmutil.o:		lfmgrid.o constants.o

st4util.o:		lfmgrid.o constants.o

metfunc.o:		constants.o

ptfcst.o:		lfmgrid.o

lfmwrite.o:		lfmgrid.o

$(EXE):			$(SRCROOT)/src/lib/liblaps.a $(SRCROOT)/src/lib/modules/libmodules.a

clean:
			$(RM) *.o $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)  
