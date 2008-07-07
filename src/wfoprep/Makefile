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

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/bgdata/libbgdata.a \
           $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/bgdata/libdebugbgdata.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a

RM=rm -f 

.SUFFIXES:              .f90 .o

.f90.o:
			$(RM) $@
			$(FC) -c $(FFLAGS) $(FREE) $*.f90

EXE=wfoprep.exe
FSRC=wfoprep.f90 \
     fill_missing_levs.f90

FOBJS=$(FSRC:.f90=.o)

FMOD = module_wfoprep_setup.o  \
       module_wfo_models.o \
       module_map_utils.o \
       module_wfoprep_mm5.o \
       module_wfoprep_wrf.o

all:			$(EXE)

$(EXE):			$(FOBJS) $(FMOD) $(LAPSLIBS)
			$(FC) -o $@ $(FFLAGS) $(FOBJS) $(FMOD) $(LAPSLIBS) $(OTHERLIBS)

debug:	
	($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) \
	"LAPSLIBS	=	$(DBLAPSLIBS)" \
	"FFLAGS         =	$(DBFLAGS)")

wfoprep.o:		$(FMOD)

module_wfo_models.o:	module_map_utils.o

module_wfoprep_mm5.o:	module_map_utils.o

module_wfoprep_mm5.o:	module_map_utils.o

clean:
			$(RM) $(FOBJS) $(FMOD) $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)
