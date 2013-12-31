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

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

RM=rm -f 

.SUFFIXES:

.SUFFIXES:              .f90 .F90 .o

LIB=libmodules.a
DEBUGLIB=libdebugmodules.a
ARCHIVE = ar rv

.F90.o:
			$(RM) $@
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.F90
			$(ARCHIVE) $(LIB) $(<:.F90=.o)

.f90.o:
			$(RM) $@
			$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90
			$(ARCHIVE) $(LIB) $(<:.f90=.o)


FMOD = mem_namelist.o mem_grid.o mem_static.o mem_wind.o mem_sfcanl.o mem_temp.o mem_vol.o module_auto_mod.o \
       module_cloud_rad.o module_constants_laps.o \
       module_horiz_interp.o  module_grib.o module_map_utils.o \
       module_mm5v3_io.o module_time_utils.o module_vinterp_utils.o \
       module_wrfsi_static.o module_wrf_netcdf.o module_grid_utils.o \
       module_grib2.o prmtrs_stmas_cloud.o ppm.o spline.o

all:			lib

lib:			$(LIB)

module_grib.o:	module_map_utils.o

module_wrf_netcdf.o:	module_grid_utils.o

module_grib2.o:	module_map_utils.o

$(LIB):			$(FMOD)
			$(RANLIB) $(LIB)
			$(RM) *.o

debug:	
	($(RM) *.o *.mod ; $(MAKE) lib	\
 	"LIB	=	$(DEBUGLIB)" 	\
	"FFLAGS =	$(DBFLAGS)" 	)

install: lib


clean:
			$(RM) $(FMOD) *.a core *.mod
