# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)

DRTUFLAG = -Drtu








.SUFFIXES:
.SUFFIXES: .o .c .f .f90 .F .F90 .a

LIB      = libcloud.a
DEBUGLIB = libdebugcloud.a




SRC= \
    cloud_var.F90 compare_analysis_to_radiometer.f crtm_kmatrix.f90 cvr_to_tb8_effective.f get_cloud_rad.f90 \
    get_sky_rgb.f90 hotstart_perturbation.f phase_func_moon.f90 places_routines.f \
    solar_routines.f allsky.f90 interp_1d.f90 skyglow_cyl.f90 vi.f90

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#


distclean: clean


.f90.o:;    $(FC) $(DBFLAGS) $(FREE) $(STINC) -c $*.f90
.F90.o:;    $(FC) $(DBFLAGS) $(FREE) $(STINC) -c $*.F90
.f.o:;      $(FC) $(DBFLAGS) -c $*.f
.F.o:;      $(FC) $(DBFLAGS) -c $*.F
