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
OTHERLIBS=-L$(NETCDF)/lib -lnetcdf 

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
	   $(SRCROOT)/src/lib/grib/libgrib.a  \
	   $(SRCROOT)/src/lib/grib2/libgrib2.a  \
	   $(SRCROOT)/src/lib/radar/rutil/librutil.a \
	   $(SRCROOT)/src/lib/radar/synp/libsynp.a

DBLAPSLIBS = $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
	     $(SRCROOT)/src/lib/grib/libdebuggrib.a  \
	     $(SRCROOT)/src/lib/grib2/libdebuggrib2.a  \
	     $(SRCROOT)/src/lib/radar/rutil/libdebugrutil.a \
	     $(SRCROOT)/src/lib/radar/synp/libdebugsynp.a
             

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
FSRC=lfmpost.f90 lfmutil.f90 start_up_mix.f90 lfmgrid.f90 interp.f90 metfunc.f90 constants.f90  \
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
			$(LAPSLIBS) $(HDF5LIB) $(OTHERLIBS) $(INC) -lm

debug:	
	 ($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
         "FFLAGS        =       $(DBFLAGS)"     \
 	 "LAPSLIB	=	$(DBLAPSLIB)" 	\
 	 "MTHERMOLIB	=	$(DBMTHERMOLIB)" 	\
 	 "UTILLIB	=	$(DBUTILLIB)" 	\
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

$(EXE):			$(SRCROOT)/src/lib/util/libutil.a $(SRCROOT)/src/lib/liblaps.a $(SRCROOT)/src/lib/modules/libmodules.a $(SRCROOT)/src/lib/radar/synp/libsynp.a                                          

clean:
			$(RM) *.o $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)  
