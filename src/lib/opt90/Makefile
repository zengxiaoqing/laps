# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)


#FFLAGS= -c \
#	-qcheck \
#	-qdbg \
#	-qextchk \
#	-qfloat=nomaf:rndsngl \
#	-qflttrap=ov:und:zero:inv \
#	-qfree=f90 \
#	-qhalt=W \
#	-qlanglvl=90pure \
#	-qsuffix=f=f90


.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f

LIB      = libopt90.a
DEBUGLIB = libdebugopt90.a

SRC= absorber_profile.f90 \
	coefficient_utility.f90 \
	error_handler.f90 \
	file_utility.f90 \
	forward_model.f90 \
	initialize.f90 \
	parameters.f90 \
	predictors.f90 \
	radiance.f90 \
	sensor_planck_routines.f90 \
	spectral_coefficients.f90 \
	transmittance.f90 \
	transmittance_coefficients.f90 \
	type_kinds.f90 \
	optran90_fm.f90 \
        forward_model.f90






ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

all: lib

lib: $(LIB)

install: lib

$(LIB): $(ARCHIVE)
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: clean

.f90.o:
	$(FC) $(FFLAGS) -c $(FREE) $<


type_kinds.o : type_kinds.f90 
transmittance_coefficients.o : transmittance_coefficients.f90 type_kinds.o error_handler.o parameters.o coefficient_utility.o 
transmittance.o : transmittance.f90 type_kinds.o parameters.o transmittance_coefficients.o 
tangent_linear_model.o : tangent_linear_model.f90 type_kinds.o error_handler.o parameters.o spectral_coefficients.o absorber_profile.o predictors.o transmittance.o radiance.o 
spectral_coefficients.o : spectral_coefficients.f90 type_kinds.o file_utility.o error_handler.o parameters.o coefficient_utility.o 
sensor_planck_routines.o : sensor_planck_routines.f90 type_kinds.o parameters.o spectral_coefficients.o 
$(EXE_FILE).o : $(EXE_FILE).f90 type_kinds.o error_handler.o parameters.o initialize.o 
radiance.o : radiance.f90 type_kinds.o parameters.o spectral_coefficients.o sensor_planck_routines.o 
predictors.o : predictors.f90 type_kinds.o parameters.o error_handler.o 
parameters.o : parameters.f90 type_kinds.o 
k_matrix_model.o : k_matrix_model.f90 type_kinds.o error_handler.o parameters.o spectral_coefficients.o absorber_profile.o predictors.o transmittance.o radiance.o forward_model.o 
initialize.o : initialize.f90 error_handler.o transmittance_coefficients.o spectral_coefficients.o 
forward_model.o : forward_model.f90 type_kinds.o error_handler.o parameters.o spectral_coefficients.o absorber_profile.o predictors.o transmittance.o radiance.o 
file_utility.o : file_utility.f90 
error_handler.o : error_handler.f90 file_utility.o 
coefficient_utility.o : coefficient_utility.f90 type_kinds.o file_utility.o error_handler.o 
adjoint_model.o : adjoint_model.f90 type_kinds.o error_handler.o parameters.o spectral_coefficients.o absorber_profile.o predictors.o transmittance.o radiance.o forward_model.o
absorber_profile.o : absorber_profile.f90 type_kinds.o parameters.o transmittance_coefficients.o 
optran90_fm.o : optran90_fm.f90 type_kinds.o file_utility.o error_handler.o parameters.o initialize.o forward_model.o spectral_coefficients.o

