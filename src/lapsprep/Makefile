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

RM=rm -f 

.SUFFIXES:              .f .o

.f.o:
			$(RM) $@
			$(FC) -c $(FFLAGS90) $*.f

EXE=lapsprep.exe
FSRC=lapsprep.f   

FOBJS=$(FSRC:.f=.o)

FMOD = module_constants.o \
       module_date_pack.o \
       module_laps_static.o \
       module_lapsprep_mm5.o \
       module_lapsprep_rams.o \
       module_lapsprep_wrf.o \
       module_setup.o


all:			$(EXE)

$(EXE):			$(FOBJS) $(FMOD)
			$(FC) -o $@ $(FOBJS) $(FMOD)  $(LAPSLIB) $(OTHERLIBS) 

debug:	
	($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
 	"LAPSLIBS	=	$(DBLAPSLIB)" 	\
	"FFLAGS90	=	$(DBFLAGS90)" 	)

lapsprep.o:		$(FMOD)

module_lapsprep_mm5.o:	module_setup.o module_laps_static.o module_date_pack.o

module_lapsprep_rams.o:	module_setup.o module_laps_static.o module_date_pack.o

module_lapsprep_wrf.o:	module_setup.o module_laps_static.o module_date_pack.o


clean:
			$(RM) $(FOBJS) $(FMOD) $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)  
