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
       module_lapsprep_output.o \
       module_setup.o


all:			$(EXE)

$(EXE):			$(FOBJS) $(FMOD)
			$(FC) -o $@ $(FOBJS) $(FMOD)  $(LAPSLIB) $(OTHERLIBS) 

debug:	
	($(RM) *.o *.exe *.mod ; $(MAKE) $(EXE) 	\
 	"LAPSLIBS	=	$(DBLAPSLIB)" 	\
	"FFLAGS90	=	$(DBFLAGS90)" 	)

lapsprep.o:		$(FMOD)

module_lapsprep_output.o:	module_setup.o module_laps_static.o module_date_pack.o 


clean:
			$(RM) $(FOBJS) $(FMOD) $(EXE) core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi

install:	$(EXE) 
		$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)  
