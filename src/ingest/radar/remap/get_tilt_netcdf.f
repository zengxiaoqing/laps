
      subroutine get_tilt_netcdf_data(filename
     1                               ,siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,istatus)

!.............................................................................

C     FORTRAN TEMPLATE FOR FILE= test.nc                                 
      PARAMETER (NVARS=29) !NUMBER OF VARIABLES
      PARAMETER (NREC=   367)   !CHANGE THIS TO GENERALIZE
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 29
      INTEGER*4 RCODE
      INTEGER*4 RECDIM
C     ****VARIABLES FOR THIS NETCDF FILE****
C
      LOGICAL*1   Z                              ( 460,NREC)
      LOGICAL*1   V                              ( 920,NREC)
      LOGICAL*1   W                              ( 920,NREC)
      INTEGER*2   elevationNumber                
      REAL*4      elevationAngle                 
      INTEGER*2   numRadials                     
      REAL*4      radialAzim                     (NREC)
      REAL*4      radialElev                     (NREC)
      REAL*8      radialTime                     (NREC)
      CHARACTER*1 siteName                       ( 132)
      CHARACTER*1 radarName                      (   5)
      REAL*4      siteLat                        
      REAL*4      siteLon                        
      REAL*4      siteAlt                        
      INTEGER*2   VCP                            
      REAL*8      esStartTime                    
      REAL*8      esEndTime                      
      REAL*4      unambigRange                   
      REAL*4      firstGateRangeZ                
      REAL*4      firstGateRangeV                
      REAL*4      gateSizeZ                      
      REAL*4      gateSizeV                      
      INTEGER*2   numGatesZ                      
      INTEGER*2   numGatesV                      
      REAL*4      resolutionV                    
      REAL*4      nyquist                        
      REAL*4      calibConst                     
      REAL*4      atmosAttenFactor               
      REAL*4      powDiffThreshold               
C*************************************
      INTEGER*4 START(10)
      INTEGER*4 COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY

!     Non-automatic declarations...............................................
      character*70 filename

      NCID=NCOPN(filename
!............................................................................
     +,NCNOWRIT,RCODE)
      CALL NCINQ(NCID,NDIMS,NVARS,NGATTS,RECDIM,RCODE)
      CALL NCDINQ(NCID,RECDIM,DUMMY,NRECS,RCODE)
C     !NRECS! NOW CONTAINS NUM RECORDS FOR THIS FILE
C
C    statements to fill Z                              
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
     +Z                              ,RCODE)
C
C    statements to fill V                              
C
      CALL NCVINQ(NCID, 2,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  20 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  20  CONTINUE
      CALL NCVGT(NCID, 2,START,COUNT,
     +V                              ,RCODE)
C
C    statements to fill W                              
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
     +W                              ,RCODE)
C
C    statements to fill elevationNumber                
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
     +elevationNumber                ,RCODE)
C
C    statements to fill elevationAngle                 
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
     +elevationAngle                 ,RCODE)
C
C    statements to fill numRadials                     
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
     +numRadials                     ,RCODE)
C
C    statements to fill radialAzim                     
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
     +radialAzim                     ,RCODE)
C
C    statements to fill radialElev                     
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
     +radialElev                     ,RCODE)
C
C    statements to fill radialTime                     
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
     +radialTime                     ,RCODE)
C
C    statements to fill siteName                       
C
      CALL NCVINQ(NCID,10,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 100 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 100  CONTINUE
      CALL NCVGTC(NCID,10,START,COUNT,
     +siteName                       ,LENSTR,RCODE)
C
C    statements to fill radarName                      
C
      CALL NCVINQ(NCID,11,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 110 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 110  CONTINUE
      CALL NCVGTC(NCID,11,START,COUNT,
     +radarName                      ,LENSTR,RCODE)
C
C    statements to fill siteLat                        
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
     +siteLat                        ,RCODE)
C
C    statements to fill siteLon                        
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
     +siteLon                        ,RCODE)
C
C    statements to fill siteAlt                        
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
     +siteAlt                        ,RCODE)
C
C    statements to fill VCP                            
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
     +VCP                            ,RCODE)
C
C    statements to fill esStartTime                    
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
     +esStartTime                    ,RCODE)
C
C    statements to fill esEndTime                      
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
     +esEndTime                      ,RCODE)
C
C    statements to fill unambigRange                   
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
     +unambigRange                   ,RCODE)
C
C    statements to fill firstGateRangeZ                
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
     +firstGateRangeZ                ,RCODE)
C
C    statements to fill firstGateRangeV                
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
     +firstGateRangeV                ,RCODE)
C
C    statements to fill gateSizeZ                      
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
     +gateSizeZ                      ,RCODE)
C
C    statements to fill gateSizeV                      
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
     +gateSizeV                      ,RCODE)
C
C    statements to fill numGatesZ                      
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
     +numGatesZ                      ,RCODE)
C
C    statements to fill numGatesV                      
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
     +numGatesV                      ,RCODE)
C
C    statements to fill resolutionV                    
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
     +resolutionV                    ,RCODE)
C
C    statements to fill nyquist                        
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
     +nyquist                        ,RCODE)
C
C    statements to fill calibConst                     
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
     +calibConst                     ,RCODE)
C
C    statements to fill atmosAttenFactor               
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
     +atmosAttenFactor               ,RCODE)
C
C    statements to fill powDiffThreshold               
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
     +powDiffThreshold               ,RCODE)
C
C     HERE IS WHERE YOU WRITE STATEMENTS TO USE THE DATA
C
C
C
!.............................................................................

      return

!.............................................................................
      END
