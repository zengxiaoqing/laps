# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../


include $(SRCROOT)/src/include/makefile.inc


LIB      = libbufr.a
DEBUGLIB = libdebugbufr.a


#
# machine_dependent needs to be first
#

SRC	= \
	addate.f \
	adn30.f \
	bfrini.f \
	bort.f \
	chekstab.f \
	chrtrn.f \
	chrtrna.f \
	cktaba.f \
	closbf.f \
	closmg.f \
	conwin.f \
	copybf.f \
	copymg.f \
	copysb.f \
	cpbfdx.f \
	cpymem.f \
	cpyupd.f \
	datebf.f \
	digit.f \
	drstpl.f \
	dumpbf.f \
	dxinit.f \
	dxmini.f \
	elemdx.f \
	getwin.f \
	i4dy.f \
	idn30.f \
	ifbget.f \
	ifxy.f \
	inctab.f \
	invcon.f \
	invmrg.f \
	invtag.f \
	invwin.f \
	ipkm.f \
	irderm.f \
	ireaderm.f \
	ireadmg.f \
	irev.f \
	iupb.f \
	iupbs1.f \
	iupm.f \
	jstify.f \
	ljust.f \
	lmsg.f \
	lstjpb.f \
	lstrpc.f \
	lstrps.f \
	makestab.f \
	mesgbf.f \
	msgini.f \
	msgupd.f \
	msgwrt.f \
	mvb.f \
	nemock.f \
	nemtab.f \
	nemtba.f \
	nemtbax.f \
	nemtbb.f \
	nemtbd.f \
	nenuck.f \
	newwin.f \
	nmbyt.f \
	nmsub.f \
	numbck.f \
	numtab.f \
	nvnwin.f \
	nwords.f \
	nxtwin.f \
	openbf.f \
	openbt.f \
	openmb.f \
	openmg.f \
	pad.f \
	parseq.f \
	parusr.f \
	parutg.f \
	pkb.f \
	pkc.f \
	pktdd.f \
	posapn.f \
	posapx.f \
	rcstpl.f \
	rdbfdx.f \
	rdcmps.f \
	rdmemm.f \
	rdmems.f \
	rdtree.f \
	rdusdx.f \
	readdx.f \
	readerm.f \
	readerme.f \
	readft.f \
	readibm.f \
	readmg.f \
	readmm.f \
	readns.f \
	readsb.f \
	readtj.f \
	rjust.f \
	rsvfvm.f \
	seqsdx.f \
	standard.f \
	status.f \
	strcln.f \
	string.f \
	strnum.f \
	strsuc.f \
	tabent.f \
	tabsub.f \
	trybump.f \
	ufbcnt.f \
	ufbcpy.f \
	ufbcup.f \
	ufbdmp.f \
	ufbevn.f \
	ufbget.f \
	ufbint.f \
	ufbmem.f \
	ufbmms.f \
	ufbmns.f \
	ufbovr.f \
	ufbqcd.f \
	ufbqcp.f \
	ufbrep.f \
	ufbrms.f \
	ufbrp.f \
	ufbrw.f \
	ufbseq.f \
	ufbsp.f \
	ufbstp.f \
	ufbtab.f \
	ufbtam.f \
	upb.f \
	upbb.f \
	upc.f \
	upds3.f \
	uptdd.f \
	usrtpl.f \
	valx.f \
	wrdlen.f \
	writdx.f \
	writsa.f \
	writsb.f \
	wrtree.f \
	wtstat.f


ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@
	rm -f *.f *.o README_BUFRLIB

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	$(RM) $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: clean

ifdef USECPP
%.o : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f
	$(FC) $(FFLAGS) -c $*.f
endif


addate.f:
	tar xvf BUFRLIB.tar
