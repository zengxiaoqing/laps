      subroutine get_raob_data(i4time_sys,ilaps_cycle_time,NX_L,NY_L
     1                                                     ,filename
     1                                                     ,istatus)

!     Ken Dritz     28-Jul-1997       Added NX_L, NY_L to dummy argument list.
!     Ken Dritz     28-Jul-1997       Added call to get_r_missing_data.
!     Ken Dritz     28-Jul-1997       Changed LAPS_DOMAIN_FILE to 'nest7grid'.
!     Ken Dritz     28-Jul-1997       Removed include of lapsparms.for.
!     Ken Dritz     28-Jul-1997       Removed comment about "non-automatic
!                                     declarations" (above arrays dimensioned
!                                     by NX_L, NY_L); they are now automatic.

C     FORTRAN TEMPLATE FOR FILE= 9614912000300o                          
      PARAMETER (NVARS=39) !NUMBER OF VARIABLES
      PARAMETER (NREC=   200)   !CHANGE THIS TO GENERALIZE
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 39
      INTEGER*4 RCODE
      INTEGER*4 RECDIM
C     ****VARIABLES FOR THIS NETCDF FILE****
C
      INTEGER*4   wmoStaNum                      (NREC)
      CHARACTER*1 staName                        (   6,NREC)
      REAL*4      staLat                         (NREC)
      REAL*4      staLon                         (NREC)
      REAL*4      staElev                        (NREC)
      REAL*8      synTime                        (NREC)
      INTEGER*4   numMand                        (NREC)
      INTEGER*4   numSigT                        (NREC)
      INTEGER*4   numSigW                        (NREC)
      INTEGER*4   numMwnd                        (NREC)
      INTEGER*4   numTrop                        (NREC)
      REAL*8      relTime                        (NREC)
      INTEGER*4   sondTyp                        (NREC)
      REAL*4      prMan                          (  22,NREC)
      REAL*4      htMan                          (  22,NREC)
      REAL*4      tpMan                          (  22,NREC)
      REAL*4      tdMan                          (  22,NREC)
      REAL*4      wdMan                          (  22,NREC)
      REAL*4      wsMan                          (  22,NREC)
      REAL*4      prSigT                         ( 150,NREC)
      REAL*4      tpSigT                         ( 150,NREC)
      REAL*4      tdSigT                         ( 150,NREC)
      REAL*4      htSigW                         (  75,NREC)
      REAL*4      wdSigW                         (  75,NREC)
      REAL*4      wsSigW                         (  75,NREC)
      REAL*4      prTrop                         (   3,NREC)
      REAL*4      tpTrop                         (   3,NREC)
      REAL*4      tdTrop                         (   3,NREC)
      REAL*4      wdTrop                         (   3,NREC)
      REAL*4      wsTrop                         (   3,NREC)
      REAL*4      prMaxW                         (   4,NREC)
      REAL*4      wdMaxW                         (   4,NREC)
      REAL*4      wsMaxW                         (   4,NREC)
      CHARACTER*1 rawTTAA                        ( 601,NREC)
      CHARACTER*1 rawTTBB                        ( 601,NREC)
      CHARACTER*1 rawTTCC                        ( 601,NREC)
      CHARACTER*1 rawTTDD                        ( 601,NREC)
      CHARACTER*1 rawPPBB                        ( 601,NREC)
      CHARACTER*1 rawPPDD                        ( 601,NREC)
C*************************************
      INTEGER*4 START(10)
      INTEGER*4 COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY

      character*70 filename
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,  
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

      NCID=NCOPN(filename
     +,NCNOWRIT,RCODE)
      CALL NCINQ(NCID,NDIMS,NVARS,NGATTS,RECDIM,RCODE)
      CALL NCDINQ(NCID,RECDIM,DUMMY,NRECS,RCODE)
C     !NRECS! NOW CONTAINS NUM RECORDS FOR THIS FILE
C
C    statements to fill wmoStaNum                      
C
      CALL NCVINQ(NCID, 1,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  10 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  10  CONTINUE
      CALL NCVGT(NCID, 1,START,COUNT,
     +wmoStaNum                      ,RCODE)
C
C    statements to fill staName                        
C
      CALL NCVINQ(NCID, 2,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  20 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  20  CONTINUE
      CALL NCVGTC(NCID, 2,START,COUNT,
     +staName                        ,LENSTR,RCODE)
C
C    statements to fill staLat                         
C
      CALL NCVINQ(NCID, 3,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  30 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  30  CONTINUE
      CALL NCVGT(NCID, 3,START,COUNT,
     +staLat                         ,RCODE)
C
C    statements to fill staLon                         
C
      CALL NCVINQ(NCID, 4,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  40 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  40  CONTINUE
      CALL NCVGT(NCID, 4,START,COUNT,
     +staLon                         ,RCODE)
C
C    statements to fill staElev                        
C
      CALL NCVINQ(NCID, 5,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  50 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  50  CONTINUE
      CALL NCVGT(NCID, 5,START,COUNT,
     +staElev                        ,RCODE)
C
C    statements to fill synTime                        
C
      CALL NCVINQ(NCID, 6,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  60 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  60  CONTINUE
      CALL NCVGT(NCID, 6,START,COUNT,
     +synTime                        ,RCODE)
C
C    statements to fill numMand                        
C
      CALL NCVINQ(NCID, 7,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  70 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  70  CONTINUE
      CALL NCVGT(NCID, 7,START,COUNT,
     +numMand                        ,RCODE)
C
C    statements to fill numSigT                        
C
      CALL NCVINQ(NCID, 8,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  80 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  80  CONTINUE
      CALL NCVGT(NCID, 8,START,COUNT,
     +numSigT                        ,RCODE)
C
C    statements to fill numSigW                        
C
      CALL NCVINQ(NCID, 9,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  90 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  90  CONTINUE
      CALL NCVGT(NCID, 9,START,COUNT,
     +numSigW                        ,RCODE)
C
C    statements to fill numMwnd                        
C
      CALL NCVINQ(NCID,10,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 100 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 100  CONTINUE
      CALL NCVGT(NCID,10,START,COUNT,
     +numMwnd                        ,RCODE)
C
C    statements to fill numTrop                        
C
      CALL NCVINQ(NCID,11,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 110 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 110  CONTINUE
      CALL NCVGT(NCID,11,START,COUNT,
     +numTrop                        ,RCODE)
C
C    statements to fill relTime                        
C
      CALL NCVINQ(NCID,12,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 120 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 120  CONTINUE
      CALL NCVGT(NCID,12,START,COUNT,
     +relTime                        ,RCODE)
C
C    statements to fill sondTyp                        
C
      CALL NCVINQ(NCID,13,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 130 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 130  CONTINUE
      CALL NCVGT(NCID,13,START,COUNT,
     +sondTyp                        ,RCODE)
C
C    statements to fill prMan                          
C
      CALL NCVINQ(NCID,14,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 140 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 140  CONTINUE
      CALL NCVGT(NCID,14,START,COUNT,
     +prMan                          ,RCODE)
C
C    statements to fill htMan                          
C
      CALL NCVINQ(NCID,15,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 150 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 150  CONTINUE
      CALL NCVGT(NCID,15,START,COUNT,
     +htMan                          ,RCODE)
C
C    statements to fill tpMan                          
C
      CALL NCVINQ(NCID,16,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 160 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 160  CONTINUE
      CALL NCVGT(NCID,16,START,COUNT,
     +tpMan                          ,RCODE)
C
C    statements to fill tdMan                          
C
      CALL NCVINQ(NCID,17,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 170 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 170  CONTINUE
      CALL NCVGT(NCID,17,START,COUNT,
     +tdMan                          ,RCODE)
C
C    statements to fill wdMan                          
C
      CALL NCVINQ(NCID,18,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 180 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 180  CONTINUE
      CALL NCVGT(NCID,18,START,COUNT,
     +wdMan                          ,RCODE)
C
C    statements to fill wsMan                          
C
      CALL NCVINQ(NCID,19,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 190 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 190  CONTINUE
      CALL NCVGT(NCID,19,START,COUNT,
     +wsMan                          ,RCODE)
C
C    statements to fill prSigT                         
C
      CALL NCVINQ(NCID,20,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 200 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 200  CONTINUE
      CALL NCVGT(NCID,20,START,COUNT,
     +prSigT                         ,RCODE)
C
C    statements to fill tpSigT                         
C
      CALL NCVINQ(NCID,21,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 210 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 210  CONTINUE
      CALL NCVGT(NCID,21,START,COUNT,
     +tpSigT                         ,RCODE)
C
C    statements to fill tdSigT                         
C
      CALL NCVINQ(NCID,22,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 220 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 220  CONTINUE
      CALL NCVGT(NCID,22,START,COUNT,
     +tdSigT                         ,RCODE)
C
C    statements to fill htSigW                         
C
      CALL NCVINQ(NCID,23,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 230 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 230  CONTINUE
      CALL NCVGT(NCID,23,START,COUNT,
     +htSigW                         ,RCODE)
C
C    statements to fill wdSigW                         
C
      CALL NCVINQ(NCID,24,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 240 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 240  CONTINUE
      CALL NCVGT(NCID,24,START,COUNT,
     +wdSigW                         ,RCODE)
C
C    statements to fill wsSigW                         
C
      CALL NCVINQ(NCID,25,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 250 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 250  CONTINUE
      CALL NCVGT(NCID,25,START,COUNT,
     +wsSigW                         ,RCODE)
C
C    statements to fill prTrop                         
C
      CALL NCVINQ(NCID,26,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 260 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 260  CONTINUE
      CALL NCVGT(NCID,26,START,COUNT,
     +prTrop                         ,RCODE)
C
C    statements to fill tpTrop                         
C
      CALL NCVINQ(NCID,27,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 270 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 270  CONTINUE
      CALL NCVGT(NCID,27,START,COUNT,
     +tpTrop                         ,RCODE)
C
C    statements to fill tdTrop                         
C
      CALL NCVINQ(NCID,28,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 280 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 280  CONTINUE
      CALL NCVGT(NCID,28,START,COUNT,
     +tdTrop                         ,RCODE)
C
C    statements to fill wdTrop                         
C
      CALL NCVINQ(NCID,29,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 290 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 290  CONTINUE
      CALL NCVGT(NCID,29,START,COUNT,
     +wdTrop                         ,RCODE)
C
C    statements to fill wsTrop                         
C
      CALL NCVINQ(NCID,30,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 300 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 300  CONTINUE
      CALL NCVGT(NCID,30,START,COUNT,
     +wsTrop                         ,RCODE)
C
C    statements to fill prMaxW                         
C
      CALL NCVINQ(NCID,31,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 310 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 310  CONTINUE
      CALL NCVGT(NCID,31,START,COUNT,
     +prMaxW                         ,RCODE)
C
C    statements to fill wdMaxW                         
C
      CALL NCVINQ(NCID,32,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 320 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 320  CONTINUE
      CALL NCVGT(NCID,32,START,COUNT,
     +wdMaxW                         ,RCODE)
C
C    statements to fill wsMaxW                         
C
      CALL NCVINQ(NCID,33,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 330 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 330  CONTINUE
      CALL NCVGT(NCID,33,START,COUNT,
     +wsMaxW                         ,RCODE)
C
C    statements to fill rawTTAA                        
C
      CALL NCVINQ(NCID,34,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 340 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 340  CONTINUE
      CALL NCVGTC(NCID,34,START,COUNT,
     +rawTTAA                        ,LENSTR,RCODE)
C
C    statements to fill rawTTBB                        
C
      CALL NCVINQ(NCID,35,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 350 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 350  CONTINUE
      CALL NCVGTC(NCID,35,START,COUNT,
     +rawTTBB                        ,LENSTR,RCODE)
C
C    statements to fill rawTTCC                        
C
      CALL NCVINQ(NCID,36,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 360 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 360  CONTINUE
      CALL NCVGTC(NCID,36,START,COUNT,
     +rawTTCC                        ,LENSTR,RCODE)
C
C    statements to fill rawTTDD                        
C
      CALL NCVINQ(NCID,37,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 370 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 370  CONTINUE
      CALL NCVGTC(NCID,37,START,COUNT,
     +rawTTDD                        ,LENSTR,RCODE)
C
C    statements to fill rawPPBB                        
C
      CALL NCVINQ(NCID,38,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 380 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 380  CONTINUE
      CALL NCVGTC(NCID,38,START,COUNT,
     +rawPPBB                        ,LENSTR,RCODE)
C
C    statements to fill rawPPDD                        
C
      CALL NCVINQ(NCID,39,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 390 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 390  CONTINUE
      CALL NCVGTC(NCID,39,START,COUNT,
     +rawPPDD                        ,LENSTR,RCODE)
C
C     HERE IS WHERE YOU WRITE STATEMENTS TO USE THE DATA
C
C
C
!     Write All Raobs to LAPS SND file

      r_nc_missing_data = 1e20

      n_snd = min(NREC,NRECS)

      do isnd = 1,n_snd

!         QC and write out the sounding

          if(abs(reltime(isnd)) .lt. 1e10)then
              i4time_release = int(sngl(reltime(isnd)))+315619200
          else
              i4time_release = 0
          endif

          if(abs(syntime(isnd)) .lt. 1e10)then
              i4time_syn     = int(sngl(syntime(isnd)))+315619200
          else
              i4time_syn = 0
          endif

          i4time_diff    = i4time_release - i4time_sys

          write(6,*)
          write(6,*)' Raob #',isnd,i4time_sys,i4time_release,i4time_diff       
     1                         ,i4time_syn

          if(stalat(isnd) .ge. r_nc_missing_data)then
              write(6,*)' Missing first latitude',i
              goto 999
          endif
          if(stalon(isnd) .ge. r_nc_missing_data)then
              write(6,*)' Missing first longitude',i
              goto 999
          endif

          if(stalat(isnd) .le. rnorth .and. stalat(isnd) .ge. south 
     1                                .AND.      
     1       stalon(isnd) .ge. west   .and. stalon(isnd) .le. east      
     1                                                            )then       
              write(6,*)' Raob is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' Outside domain lat/lon perimeter - reject'
              goto 999
          endif

          call sort_and_write(
     1                        NREC,isnd,r_missing_data
     1                       ,wmostanum,staname,stalat,stalon,staelev
     1                       ,nummand,htman,prman,tpman,tdman      
     1                       ,wdman,wsman
     1                       ,numsigt,prsigt,tpsigt,tdsigt
     1                       ,numsigw,htsigw,wdsigw,wssigw
     1                       ,istatus)

          go to 999

 998      write(6,*)' Error writing out RAOB'

 999      continue

      enddo ! i

      return
      END


      subroutine sort_and_write(
     1                        NREC,isnd,r_missing_data
     1                       ,wmostanum,staname,stalat,stalon,staelev
     1                       ,nummand,htman,prman,tpman,tdman      
     1                       ,wdman,wsman
     1                       ,numsigt,prsigt,tpsigt,tdsigt
     1                       ,numsigw,htsigw,wdsigw,wssigw
     1                       ,istatus)

      INTEGER*4   wmoStaNum                      (NREC)
      CHARACTER*1 staName                        (   6,NREC)
      REAL*4      staLat                         (NREC)
      REAL*4      staLon                         (NREC)
      REAL*4      staElev                        (NREC)

      INTEGER*4   numMand                        (NREC)
      REAL*4      prMan                          (  22,NREC)
      REAL*4      htMan                          (  22,NREC)
      REAL*4      tpMan                          (  22,NREC)
      REAL*4      tdMan                          (  22,NREC)
      REAL*4      wdMan                          (  22,NREC)
      REAL*4      wsMan                          (  22,NREC)

      INTEGER*4   numsigt                        (NREC)
      REAL*4      prSigT                         ( 150,NREC)
      REAL*4      tpSigT                         ( 150,NREC)
      REAL*4      tdSigT                         ( 150,NREC)

      INTEGER*4   numsigw                        (NREC)
      REAL*4      htSigW                         (  75,NREC)
      REAL*4      wdSigW                         (  75,NREC)
      REAL*4      wsSigW                         (  75,NREC)

      PARAMETER (NLVL_OUT=   150)         
      integer*4   indx(NLVL_OUT)  
      REAL*4      prout                          (NLVL_OUT)
      REAL*4      htout                          (NLVL_OUT)
      REAL*4      tpout                          (NLVL_OUT)
      REAL*4      tdout                          (NLVL_OUT)
      REAL*4      wdout                          (NLVL_OUT)
      REAL*4      wsout                          (NLVL_OUT)

!     Generate info for Sorting/QC, write original mandatory data to log file
      write(6,*)
      n_good_levels = 0

      write(6,*)' Subroutine sort_and_write - initial mandatory data'       
      do ilevel = 1,nummand(isnd)
          if(htman(ilevel,isnd) .lt. 90000.)then
              n_good_levels = n_good_levels + 1
              write(6,*) htman(ilevel,isnd),prman(ilevel,isnd)
     1                  ,tpman(ilevel,isnd),tdman(ilevel,isnd)
     1                  ,wdman(ilevel,isnd),wsman(ilevel,isnd)

              indx(n_good_levels) = n_good_levels
              htout(n_good_levels) = htman(ilevel,isnd)
              prout(n_good_levels) = prman(ilevel,isnd)
              tpout(n_good_levels) = tpman(ilevel,isnd)
              tdout(n_good_levels) = tdman(ilevel,isnd)
              wdout(n_good_levels) = wdman(ilevel,isnd)
              wsout(n_good_levels) = wsman(ilevel,isnd)
          endif
      enddo

      write(6,*)' Subroutine sort_and_write - sig wind data'       
      do ilevel = 1,numsigw(isnd)
          if(htsigw(ilevel,isnd) .lt. 90000. .and.
     1       htsigw(ilevel,isnd) .ne. 0.            )then
              n_good_levels = n_good_levels + 1
              write(6,*) htsigw(ilevel,isnd),r_missing_data
     1                  ,r_missing_data,r_missing_data
     1                  ,wdsigw(ilevel,isnd),wssigw(ilevel,isnd)

              indx(n_good_levels) = n_good_levels
              htout(n_good_levels) = htsigw(ilevel,isnd)
              prout(n_good_levels) = r_missing_data
              tpout(n_good_levels) = r_missing_data
              tdout(n_good_levels) = r_missing_data
              wdout(n_good_levels) = wdsigw(ilevel,isnd)
              wsout(n_good_levels) = wssigw(ilevel,isnd)
          endif
      enddo

!     Bubble sort the levels by height
 400  iswitch = 0
      do i = 2,n_good_levels
          if(htout(indx(i)) .lt. htout(indx(i-1)))then
              izz = indx(i-1)
              indx(i-1) = indx(i)
              indx(i) = izz
              iswitch = 1
          endif
      enddo

      if(iswitch .eq. 1)go to 400

      write(6,*)
      write(6,511,err=998)
     1             wmostanum(isnd),n_good_levels,stalat(isnd)
     1            ,stalon(isnd),staelev(isnd),(staname(ic,isnd),ic=1,5)       
      write(11,511,err=998)
     1             wmostanum(isnd),n_good_levels,stalat(isnd)
     1            ,stalon(isnd),staelev(isnd),(staname(ic,isnd),ic=1,5)       

  511 format(i12,i12,f11.4,f15.4,f15.0,1x,5a1)


!     Write out all sorted data for mandatory + sigw levels. 
!     T and Td are in deg C
      do i = 1,n_good_levels
          ilevel = indx(i)

          if(tpout(ilevel) .eq. 99999. .or.
     1       tpout(ilevel) .eq. r_missing_data     )then
              t_c = r_missing_data
          else
              t_c = tpout(ilevel) - 273.15
          endif

          if(tpout(ilevel) .eq. 99999. .or.
     1       tdout(ilevel) .eq. 99999. .or. t_c .eq. r_missing_data)then       
              td_c = r_missing_data
          else
              td_c = tpout(ilevel) - 273.15 - tdout(ilevel)
          endif

          if(wdout(ilevel) .eq. 99999. .or.
     1       wsout(ilevel) .eq. 99999.)then
              wdout(ilevel) = r_missing_data
              wsout(ilevel) = r_missing_data
          endif

          write(6,*) htout(ilevel),prout(ilevel)
     1              ,t_c
     1              ,td_c
     1              ,wdout(ilevel),wsout(ilevel),ilevel
          write(11,*)htout(ilevel),prout(ilevel)
     1              ,t_c
     1              ,td_c
     1              ,wdout(ilevel),wsout(ilevel) 
      enddo

      go to 999

 998  write(6,*)' Error writing out RAOB'

 999  continue

      return
      end


