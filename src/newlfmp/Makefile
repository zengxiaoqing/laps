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
LAPSLIB=-L$(SRCROOT)/src/lib -llaps
DBLAPSLIB=-L$(SRCROOT)/src/lib -ldebuglaps
SYNPLIB=-L$(SRCROOT)/src/lib/radar/synp -lsynp
DBSYNPLIB=-L$(SRCROOT)/src/lib/radar/synp -ldebugsynp
GRIBLIB=-L$(SRCROOT)/src/lib/grib -lgrib
DBGRIBLIB=-L$(SRCROOT)/src/lib/grib -ldebuggrib
GRIB2LIB=-L$(SRCROOT)/src/lib/grib2 -lgrib2
DBGRIB2LIB=-L$(SRCROOT)/src/lib/grib2 -ldebuggrib2
#NETCDF=/opt/netcdf/3.6.2-pgi-7.1-3
OTHERLIBS=-L$(NETCDF)/lib -lnetcdf 

include $(SRCROOT)/src/include/makefile.inc

# This can be uncommented if desired for the regular 'make' to have the 
# optimization turned off. An example is with 'xlf' on the IBM where the
# -O3 optimization can take very long to compile.

#FFLAGS=$(DBFLAGS)

RM=rm -f 

.SUFFIXES:              .f90 .o

%.o:%.mod

.f90.o:
			$(RM) $@
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90

EXE=lfmpost.exe
FSRC=lfmpost.f90 lfmutil.f90 lfmgrid.f90 interp.f90 metfunc.f90 constants.f90  \
     mm5util.f90 wrfutil.f90 nmmutil.f90 st4util.f90  \
     ptfcst.f90 lfmwrite.f90

CSRC = 

FOBJS=$(FSRC:.f90=.o)

COBJS=$(CSRC:.c=.o)


FMOD = constants.o \
       lfmgrid.o \
       mm5util.o \
       nmmutil.o \
       st4util.o \
       wrfutil.o 

all:			$(EXE)

$(EXE):			$(FOBJS) $(COBJS) $(FMOD)
			$(FC) $(FREE) -o $@ $(FFLAGS) $(FOBJS) $(COBJS) \
			$(LAPSLIB) $(HDF5LIB) $(OTHERLIBS) $(SYNPLIB) $(EXTMODS) $(GRIBLIB) $(GRIB2LIB) $(INC) -lm

debug:	
	 ($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
         "FFLAGS        =       $(DBFLAGS)"     \
 	 "LAPSLIB	=	$(DBLAPSLIB)" 	\
 	 "SYNPLIB	=	$(DBSYNPLIB)" 	\
         "GRIBLIB       =       $(DBGRIBLIB)"    \
         "GRIB2LIB      =       $(DBGRIB2LIB)"   \
         "FFLAGS90 	=	$(DBFLAGS90)"    \
         "EXTMODS       =       $(DBEXTMODS)" 	 )

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

$(EXE):			$(SRCROOT)/src/lib/liblaps.a $(SRCROOT)/src/lib/modules/libmodules.a $(SRCROOT)/src/lib/radar/synp/libsynp.a

clean:
			$(RM) *.o $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)  
