# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov


LAPSROOT=.
include $(LAPSROOT)/src/include/makefile.inc

LINK=$(LAPSROOT)/util/link.pl
CWD = $(shell pwd)

MACHDEP = config.log config.status config.cache src/include/config.h

LIBDIRS = src/lib \
          src/lib/bgdata \
          src/lib/blas \
          src/lib/cloud \
          src/lib/fm \
          src/lib/goesinav \
          src/lib/goeslib \
          src/lib/lapack \
          src/lib/mthermo \
          src/lib/nav \
          src/lib/powell  \
          src/lib/radar/remap_ftn \
          src/lib/satellite \
          src/lib/temp   \
          src/lib/util  

EXEDIRS = src/accum \
          src/background \
          src/balance \
          src/cloud \
          src/deriv \
	  src/dprep \
          src/grid \
          src/humid \
          src/ingest/acars \
          src/ingest/pireps \
          src/ingest/profiler \
          src/ingest/radar/ingest_rrv \
          src/ingest/radar/remap \
          src/ingest/radar/vad \
          src/ingest/radar/wsiremap/vrc \
          src/ingest/radar/wsiremap/vrc/table \
          src/ingest/raob \
          src/ingest/rass \
          src/ingest/rass/blp \
          src/ingest/sao \
          src/ingest/satellite/cloud_drift \
          src/ingest/satellite/lvd \
          src/ingest/satellite/lvd/table \
          src/ingest/satellite/sounding \
          src/sched \
          src/sfc \
          src/sfc/table \
          src/soil \
          src/temp \
          src/wind \
          src/WFO/post

DATADIRS = log/qc \
           static \
           time \
           cdl \
           lapsprd/balance/lt1 \
           lapsprd/balance/lw3 \
           lapsprd/cdw \
           lapsprd/dprep \
           lapsprd/d01 \
           lapsprd/d02 \
           lapsprd/d03 \
           lapsprd/d04 \
           lapsprd/d05 \
           lapsprd/d06 \
           lapsprd/d07 \
           lapsprd/d08 \
           lapsprd/d09 \
           lapsprd/d10 \
           lapsprd/d11 \
           lapsprd/d12 \
           lapsprd/d13 \
           lapsprd/d14 \
           lapsprd/d15 \
           lapsprd/d16 \
           lapsprd/d17 \
           lapsprd/d18 \
           lapsprd/d19 \
           lapsprd/d20 \
           lapsprd/fsf \
           lapsprd/fua \
           lapsprd/grid \
           lapsprd/l1s \
           lapsprd/lc3 \
           lapsprd/lcb \
           lapsprd/lco \
           lapsprd/lcp \
           lapsprd/lct \
           lapsprd/lcv \
           lapsprd/lf1 \
           lapsprd/lga \
           lapsprd/lgb \
           lapsprd/lh3 \
           lapsprd/lh4 \
           lapsprd/lhe \
           lapsprd/lil \
           lapsprd/liw \
           lapsprd/lm1 \
           lapsprd/lm2 \
           lapsprd/lmd \
           lapsprd/lmr \
           lapsprd/lmt \
           lapsprd/lpbl \
           lapsprd/lps \
           lapsprd/lq3 \
           lapsprd/lrp \
           lapsprd/lrs \
           lapsprd/lso \
           lapsprd/lsr \
           lapsprd/lsr/dmsp01 \
           lapsprd/lsr/dmsp02 \
           lapsprd/lsr/goes08 \
           lapsprd/lsr/goes09 \
           lapsprd/lsr/goes10 \
           lapsprd/lsr/tros12 \
           lapsprd/lsr/tros14 \
           lapsprd/lsx \
           lapsprd/lt1 \
           lapsprd/lty \
           lapsprd/lvd/goes08 \
           lapsprd/lvd/goes09 \
           lapsprd/lvd/goes10 \
           lapsprd/lvd/gmssat \
           lapsprd/lw3 \
           lapsprd/lwc \
           lapsprd/lwm \
           lapsprd/msg \
	   lapsprd/model/varfiles \
	   lapsprd/model/output \
	   lapsprd/model/sfc \
           lapsprd/pig \
           lapsprd/pin \
           lapsprd/prg \
           lapsprd/pro \
	   lapsprd/ram \
	   lapsprd/rsf \
           lapsprd/rdr \
           lapsprd/rdr/001 \
           lapsprd/rdr/002 \
           lapsprd/rdr/003 \
           lapsprd/rdr/004 \
           lapsprd/rdr/005 \
           lapsprd/rdr/006 \
           lapsprd/rdr/007 \
           lapsprd/rdr/008 \
           lapsprd/rdr/009 \
           lapsprd/rdr/001/vrc \
           lapsprd/rdr/002/vrc \
           lapsprd/rdr/003/vrc \
           lapsprd/rdr/004/vrc \
           lapsprd/rdr/005/vrc \
           lapsprd/rdr/006/vrc \
           lapsprd/rdr/007/vrc \
           lapsprd/rdr/008/vrc \
           lapsprd/rdr/009/vrc \
           lapsprd/sag \
           lapsprd/snd \
           lapsprd/v01 \
           lapsprd/v02 \
           lapsprd/v03 \
           lapsprd/v04 \
           lapsprd/v05 \
           lapsprd/v06 \
           lapsprd/v07 \
           lapsprd/v08 \
           lapsprd/v09 \
           lapsprd/v10 \
           lapsprd/v11 \
           lapsprd/v12 \
           lapsprd/v13 \
           lapsprd/v14 \
           lapsprd/v15 \
           lapsprd/v16 \
           lapsprd/v17 \
           lapsprd/v18 \
           lapsprd/v19 \
           lapsprd/vdr \
           lapsprd/vrc 


all: exe
debug: debuglib

localize: mkdatadirs
	/usr/nfs/bin/perl $(LAPSROOT)/etc/laps_localization.pl --lapsroot=$(INSTALLROOT) \
               --dataroot=$(DATAROOT) --srcroot=$(LAPSROOT)

exe: lib
	@for dir in $(EXEDIRS) ;\
	  do \
	  echo " ";\
	  echo Making Laps in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) all ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

lib:
	@for dir in $(LIBDIRS) ;\
	  do \
	  echo " ";\
	  echo Making Laps in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) all ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

debuglib:
	@for dir in $(LIBDIRS) ;\
	  do \
	  echo " ";\
	  echo Making Laps in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) debug ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

lapsplot: lib
	cd src/lapsplot; $(MAKE)

install_lapsplot: lib
	cd src/lapsplot; $(MAKE) install


install: mkdirs
	@for dir in $(LIBDIRS) $(EXEDIRS);\
	  do \
	  echo " ";\
	  echo Making Laps in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) install ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

mkdirs: mkdatadirs
	if [ ! -d  $(INSTALLROOT) ] ; then  \
	mkdir -p $(INSTALLROOT) ; fi 
	if [ ! -d $(INSTALLROOT)/util ] ; then  \
	cp -r  $(LAPSROOT)/util $(INSTALLROOT)/util ; fi 
	if [ ! -d $(INSTALLROOT)/etc ] ; then  \
	cp -r $(LAPSROOT)/etc  $(INSTALLROOT) ; fi

mkdatadirs :
	if [ ! -d $(DATAROOT) ] ; then  \
	cp -r $(LAPSROOT)/data $(DATAROOT) ; fi 
	@for dir in $(DATADIRS) ;\
	  do \
	  if [ ! -d $(DATAROOT)/$$dir ] ; then  \
	    echo Creating Laps data directory $(DATAROOT)/$$dir ;\
	    (cd $(DATAROOT); \
	    mkdir -p $$dir ; if [ $$? != 0 ] ; then \
	          echo "Exit status from mkdir was $$?" ; exit 1 ; fi ;) ;\
	    fi ; done ; 




clean:
	$(RM) *~ *# *.o 
	cd src/include; $(RM) *~ *# *.o 
	@for dir in $(EXEDIRS);\
	  do \
	  echo " ";\
	  echo Making LAPS clean in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

realclean: clean
	@for dir in $(LIBDIRS);\
	  do \
	  echo " ";\
	  echo Making LAPS clean in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

distclean: realclean cleandirs
	@for dir in $(LIBDIRS) $(EXEDIRS);\
	  do \
	  echo " ";\
	  echo Making LAPS distclean in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) distclean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	$(RM) $(MACHDEP) \
	      ./util/cronfile ./etc/*.pl ./etc/*.csh ; \
	$(RM) -r ./bin ./log

cleandirs:
	@for dir in $(DATADIRS) ;\
	  do \
	  echo " ";\
	  echo Removing Laps data directory $(LAPSROOT)/$$dir ;\
	  (cd $(LAPSROOT); \
	   $(RM) -r $$dir ; if [ $$? != 0 ] ; then \
	        echo "Exit status from rm was $$?" ; exit 1 ; fi ;) ;\
	  done

links:
	$(LINK) $(CWD) $(LAPSROOT) $(INSTALLROOT)/lib $(MACHDEP)
	@for dir in $(LIBDIRS) $(EXEDIRS);\
	  do \
	  echo " ";\
	  echo Making LAPS links in directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) links ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

