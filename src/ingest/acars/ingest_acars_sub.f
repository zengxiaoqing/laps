
      subroutine get_acars_data(i4time_sys,ilaps_cycle_time
     1                                    ,NX_L,NY_L
     1                                    ,filename,istatus)

!     Ken Dritz     15-Jul-1997      Changed LAPS_DOMAIN_FILE to 'nest7grid'
!     Ken Dritz     15-Jul-1997      Removed include of lapsparms.for
!     Steve Albers  22-Jul-1997      Added NX_L, NY_L to dummy argument list

!.............................................................................

C     FORTRAN TEMPLATE FOR FILE= latest_q.cdf                            
      PARAMETER (NVARS=32) !NUMBER OF VARIABLES
      PARAMETER (NREC=   50000)   !CHANGE THIS TO GENERALIZE
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 32
      INTEGER*4 RCODE
      INTEGER*4 RECDIM
C     ****VARIABLES FOR THIS NETCDF FILE****
C
      CHARACTER*1 minDate                        (  30)
      CHARACTER*1 maxDate                        (  30)
      REAL*8      minSecs                        
      REAL*8      maxSecs                        
      REAL*4      latitude                       (NREC)
      REAL*4      longitude                      (NREC)
      REAL*4      altitude                       (NREC)
      REAL*8      timeObs                        (NREC)
      REAL*4      temperature                    (NREC)
      REAL*4      windDir                        (NREC)
      REAL*4      windSpeed                      (NREC)
      REAL*4      heading                        (NREC)
      INTEGER*4   turbIntens                     (NREC)
      REAL*4      waterVaporMR                   (NREC)
      REAL*4      vertAccel                      (NREC)
      CHARACTER*1 tailNumber                     (   6,NREC)
      LOGICAL*1   airline                        (NREC)
      LOGICAL*1   dataDescriptor                 (NREC)
      LOGICAL*1   errorType                      (NREC)
      LOGICAL*1   rollFlag                       (NREC)
      LOGICAL*1   waterVaporQC                   (NREC)
      LOGICAL*1   interpolatedTime               (NREC)
      LOGICAL*1   interpolatedLL                 (NREC)
      LOGICAL*1   tempError                      (NREC)
      LOGICAL*1   windDirError                   (NREC)
      LOGICAL*1   windSpeedError                 (NREC)
      LOGICAL*1   speedError                     (NREC)
      LOGICAL*1   bounceError                    (NREC)
      LOGICAL*1   correctedFlag                  (NREC)
      CHARACTER*1 flight                         (  13,NREC)
      CHARACTER*1 rptStation                     (   4,NREC)
      REAL*8      timeReceived                   (NREC)
C*************************************
      INTEGER*4 START(10)
      INTEGER*4 COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY

!     Non-automatic declarations...............................................

      character*1 c1_dataDescriptor                 (NREC)
      equivalence(c1_dataDescriptor,dataDescriptor)

      character*1 c1_errorType                      (NREC)
      equivalence(c1_errorType,errorType)

      character*70 filename
      character*9 a9_timeObs,a9_recptTime 
      character*7 c7_skycover
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif

      NCID=NCOPN(filename
!............................................................................
     +,NCNOWRIT,RCODE)
      CALL NCINQ(NCID,NDIMS,NVARS,NGATTS,RECDIM,RCODE)
      CALL NCDINQ(NCID,RECDIM,DUMMY,NRECS,RCODE)
C     !NRECS! NOW CONTAINS NUM RECORDS FOR THIS FILE
C
C    statements to fill minDate                        
C
      CALL NCVINQ(NCID, 1,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  10 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  10  CONTINUE
      CALL NCVGTC(NCID, 1,START,COUNT,
     +minDate                        ,LENSTR,RCODE)
C
C    statements to fill maxDate                        
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
     +maxDate                        ,LENSTR,RCODE)
C
C    statements to fill minSecs                        
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
     +minSecs                        ,RCODE)
C
C    statements to fill maxSecs                        
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
     +maxSecs                        ,RCODE)
C
C    statements to fill latitude                       
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
     +latitude                       ,RCODE)
C
C    statements to fill longitude                      
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
     +longitude                      ,RCODE)
C
C    statements to fill altitude                       
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
     +altitude                       ,RCODE)
C
C    statements to fill timeObs                        
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
     +timeObs                        ,RCODE)
C
C    statements to fill temperature                    
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
     +temperature                    ,RCODE)
C
C    statements to fill windDir                        
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
     +windDir                        ,RCODE)
C
C    statements to fill windSpeed                      
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
     +windSpeed                      ,RCODE)
C
C    statements to fill heading                        
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
     +heading                        ,RCODE)
C
C    statements to fill turbIntens                     
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
     +turbIntens                     ,RCODE)
C
C    statements to fill waterVaporMR                   
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
     +waterVaporMR                   ,RCODE)
C
C    statements to fill vertAccel                      
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
     +vertAccel                      ,RCODE)
C
C    statements to fill tailNumber                     
C
      CALL NCVINQ(NCID,16,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 160 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 160  CONTINUE
      CALL NCVGTC(NCID,16,START,COUNT,
     +tailNumber                     ,LENSTR,RCODE)
C
C    statements to fill airline                        
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
     +airline                        ,RCODE)
C
C    statements to fill dataDescriptor                 
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
     +dataDescriptor                 ,RCODE)
C
C    statements to fill errorType                      
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
     +errorType                      ,RCODE)
C
C    statements to fill rollFlag                       
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
     +rollFlag                       ,RCODE)
C
C    statements to fill waterVaporQC                   
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
     +waterVaporQC                   ,RCODE)
C
C    statements to fill interpolatedTime               
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
     +interpolatedTime               ,RCODE)
C
C    statements to fill interpolatedLL                 
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
     +interpolatedLL                 ,RCODE)
C
C    statements to fill tempError                      
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
     +tempError                      ,RCODE)
C
C    statements to fill windDirError                   
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
     +windDirError                   ,RCODE)
C
C    statements to fill windSpeedError                 
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
     +windSpeedError                 ,RCODE)
C
C    statements to fill speedError                     
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
     +speedError                     ,RCODE)
C
C    statements to fill bounceError                    
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
     +bounceError                    ,RCODE)
C
C    statements to fill correctedFlag                  
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
     +correctedFlag                  ,RCODE)
C
C    statements to fill flight                         
C
      CALL NCVINQ(NCID,30,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 300 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 300  CONTINUE
      CALL NCVGTC(NCID,30,START,COUNT,
     +flight                         ,LENSTR,RCODE)
C
C    statements to fill rptStation                     
C
      CALL NCVINQ(NCID,31,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 310 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 310  CONTINUE
      CALL NCVGTC(NCID,31,START,COUNT,
     +rptStation                     ,LENSTR,RCODE)
C
C    statements to fill timeReceived                   
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
     +timeReceived                   ,RCODE)
C
C     HERE IS WHERE YOU WRITE STATEMENTS TO USE THE DATA
C
C
C
!.............................................................................

      if(NRECS .gt. NREC)then
          write(6,*)' Warning: NRECS > NREC ',NRECS,NREC
      else
          write(6,*)' NRECS, NREC = ',NRECS,NREC
      endif

      num_acars = min(NREC,NRECS)

      do i = 1,num_acars

          write(6,*)
          write(6,*)' acars #',i,'  ',c1_dataDescriptor(i)
     1                               ,c1_errorType(i)
!         write(6,*)' location = '
!    1             ,latitude(i),longitude(i),altitude(i)

          if(c1_dataDescriptor(i) .eq. 'X')then
            if(c1_errorType(i) .eq. 'W' .or. 
     1         c1_errorType(i) .eq. 'B'                         )then
              write(6,*)' QC flag is bad - reject ',c1_dataDescriptor(i)   
     1                                             ,c1_errorType(i)
              goto 900
            endif
          endif


          if(latitude(i) .le. rnorth .and. latitude(i)  .ge. south .and.
     1       longitude(i) .ge. west  .and. longitude(i) .le. east      
     1                                                             )then       
              continue
          else ! Outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,latitude(i),longitude(i)
              goto 900
          endif


          if(abs(timeObs(i))      .lt. 3d9       .and.
     1       abs(timereceived(i)) .lt. 3d9              )then
              call c_time2fname(nint(timeObs(i)),a9_timeObs)
              call c_time2fname(nint(timereceived(i)),a9_recptTime)
          else
              write(6,*)' Bad observation time - reject'       
     1                   ,timeObs(i),timereceived(i)
              goto 900
          endif


          call cv_asc_i4time(a9_timeObs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_timeObs,i4_resid,ilaps_cycle_time / 2
              goto 900        
          endif

          write(6,1)a9_timeObs,a9_recptTime 
          write(11,1)a9_timeObs,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)latitude(i),longitude(i),altitude(i)
          write(11,2)latitude(i),longitude(i),altitude(i)
 2        format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  

          write(6,*)' Wind = ',windDir(i),windSpeed(i)
 
          if(abs(windSpeed(i)) .gt. 250.)then
              write(6,*)' Above wind speed is suspect - reject'
              go to 900
          endif

          write(6,3)int(windDir(i)),windSpeed(i)
          write(11,3)int(windDir(i)),windSpeed(i)
 3        format(' Wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')

!        if(string(2:5) .eq. 'Wind')then
!            read(lun,202)idir_deg,ispd_kt
!202         format (' ', i3, ' deg @ ', i3, ' kt.')
!            write(6,202)idir_deg,ispd_kt
!            dd = idir_deg
!            ff = ispd_kt * .518
!            return
!        endif


 900  enddo ! i

      return

!.............................................................................

      END
