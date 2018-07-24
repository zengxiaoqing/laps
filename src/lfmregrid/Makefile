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

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/satellite/libsatellite.a  \
	   $(SRCROOT)/src/lib/degrib/libdegrib.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
	   $(SRCROOT)/src/lib/bgdata/libbgdata.a  \
	   $(SRCROOT)/src/lib/radar/rutil/librutil.a \
	   $(SRCROOT)/src/lib/radar/synp/libsynp.a \
           $(SRCROOT)/src/lib/w3lib/libw3.a \
           $(SRCROOT)/src/lib/g2lib/libg2_4.a


DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/satellite/libdebugsatellite.a   \
	     $(SRCROOT)/src/lib/degrib/libdebugdegrib.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
	     $(SRCROOT)/src/lib/bgdata/libdebugbgdata.a  \
	     $(SRCROOT)/src/lib/degrib/libdebugdegrib.a  \
	     $(SRCROOT)/src/lib/radar/rutil/libdebugrutil.a \
	     $(SRCROOT)/src/lib/radar/synp/libdebugsynp.a \
             $(SRCROOT)/src/lib/w3lib/libdebugw3.a \
             $(SRCROOT)/src/lib/g2lib/libdebugg2_4.a

             

# This can be uncommented if desired for the regular 'make' to have the 
# optimization turned off. An example is with 'xlf' on the IBM where the
# -O3 optimization can take very long to compile.

#FFLAGS=$(DBFLAGS)

RM=rm -f 

.SUFFIXES:              .f90 .o .c .f .a .F .i .exe

%.o:%.mod

.f90.o:
			$(RM) $@
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90

EXE=lfmregrid.exe
FSRC=lfmregrid.f90 lfmregrid_sub.f90 lfmwrite_rg.f90

CSRC = 

FOBJS=$(FSRC:.f90=.o)

COBJS=$(CSRC:.c=.o)


FMOD = 

all:			$(EXE)

$(EXE):			$(FOBJS) $(COBJS) $(FMOD)
			$(FC) $(FREE) -o $@ $(FFLAGS) $(FOBJS) $(COBJS) \
			$(BGDATALIB) $(LAPSLIBS) $(HDF5LIB) $(OTHERLIBS) $(DEGRIBLIBS) $(INC) -lm

debug:	
	 ($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
         "FFLAGS        =       $(DBFLAGS)"     \
 	 "LAPSLIBS	=	$(DBLAPSLIBS)" 	\
 	 "MTHERMOLIB	=	$(DBMTHERMOLIB)" 	\
 	 "UTILLIB	=	$(DBUTILLIB)" 	\
 	 "SYNPLIB	=	$(DBSYNPLIB)" 	\
         "BGDATALIB     =       $(DBBGDATALIB)"    \
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
