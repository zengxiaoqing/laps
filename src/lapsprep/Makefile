SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc


.SUFFIXES:              .f .o

.f.o:
			$(RM) $@
			$(FC) -c $(FFLAGS90) $*.f

FSRC		=	lapsprep.f   

FOBJS		=	$(FSRC:.f=.o)

FMOD		=	module_constants.o module_date_pack.o module_laps_static.o \
			module_lapsprep_output.o module_setup.o


all:			lapsprep.exe

lapsprep.exe:		$(FOBJS) $(FMOD)
			$(FC) -o $@ $(FOBJS) $(FMOD)  $(LAPSLIB) $(OTHERLIBS) 

debug:	
	($(RM) *.o *.exe *.mod ; $(MAKE) lapsprep.exe 	\
 	"LAPSLIBS	=	$(DBLAPSLIB)" 	\
	"FFLAGS90	=	$(DBFLAGS90)" 	)

lapsprep.o:		$(FMOD)

module_lapsprep_output.o:	module_setup.o module_laps_static.o module_date_pack.o 


clean:
			$(RM) $(FOBJS) $(FMOD) lapsprep.exe core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi
