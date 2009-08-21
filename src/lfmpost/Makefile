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
GRIBLIB=-L$(SRCROOT)/src/lib/grib -lgrib
DBGRIBLIB=-L$(SRCROOT)/src/lib/grib -ldebuggrib
OTHERLIBS=-L$(NETCDF)/lib -lnetcdf

include $(SRCROOT)/src/include/makefile.inc

RM=rm -f 

.SUFFIXES:              .f90 .o

.f90.o:
			$(RM) $@
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90

EXE=lfmpost.exe
FSRC=lfmpost.f90 bint.f90 clouds.f90 integrated_liquid.f90 reflectivity.f90 \
     met_functions.f90 interp_press_to_z.f90 model_pblhgt.f90 \
     snowfall.f90 wintprec.f90 helicity.f90 capecin.f90 the2t.f90 \
     output_laps_format.f90 smooth.f90 grib_sfc_vars.f90 grib_ua_vars.f90 \
     get_file_unit.f90 tke.f90

CSRC = v5d.c binio.c

FOBJS=$(FSRC:.f90=.o)

COBJS=$(CSRC:.c=.o)

FMOD = module_constants.o \
       module_postproc_lfm.o \
       module_setup.o \
       module_vis5d.o \
       module_fire.o

all:			$(EXE)

$(EXE):			$(FOBJS) $(COBJS) $(FMOD)
			$(FC) $(FREE) -o $@ $(FFLAGS) $(FOBJS) $(COBJS) \
                        $(FMOD) $(LAPSLIB) $(OTHERLIBS) $(EXTMODS) $(GRIBLIB) $(INC)

debug:	
	($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
	"LAPSLIB	=	$(DBLAPSLIB)" 	\
        "GRIBLIB        =       $(DBGRIBLIB)"    \
	"FFLAGS90	=	$(DBFLAGS90)"    \
	"EXTMODS        =       $(DBEXTMODS)" 	)

lfmpost.o:		module_setup.o \
			module_postproc_lfm.o module_constants.o \
       			module_vis5d.o


module_postproc_lfm.o:	module_constants.o module_setup.o module_fire.o tke.f90

met_functions.o:	module_constants.o
 
module_vis5d.o:	module_setup.o

slprs_from_sigma.o:	module_constants.o

wintprec.o:	module_constants.o

snowfall.o:	module_constants.o

helicity.o:	module_constants.o

capecin.o:	module_constants.o

the2t.o:	module_constants.o

clean:
			$(RM) $(FOBJS) $(COBJS) $(FMOD) $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)  
