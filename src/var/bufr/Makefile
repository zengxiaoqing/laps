# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  


DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a 


.f.o:;      $(FC) $(FLAG) -c $*.f

LIB=libbufr.a

all: $(LIB)

$(LIB):
	tar xvf BUFRLIB.tar
	$(FC) $(FLAG) -c *.f
	$(AR) -r $(LIB) *.o
	rm -f *.f *.o README_BUFRLIB

clean:
	rm -f $(LIB)


distclean: clean

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F

