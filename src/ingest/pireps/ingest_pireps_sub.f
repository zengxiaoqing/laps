
      subroutine get_pirep_data(i4time_sys,ilaps_cycle_time,filename
     1                                                     ,NX_L,NY_L
     1                                                     ,istatus)

!     Ken Dritz     15-Jul-1997       Added NX_L, NY_L to dummy arg. list
!     Ken Dritz     15-Jul-1997       Changed LAPS_DOMAIN_FILE to 'nest7grid'
!     Ken Dritz     15-Jul-1997       Removed include of lapsparms.for

C     FORTRAN TEMPLATE FOR FILE= test.nc                                 
      PARAMETER (NVARS=33) !NUMBER OF VARIABLES
C     PARAMETER (NREC=    97)   !CHANGE THIS TO GENERALIZE
      PARAMETER (NREC=    97)   !CHANGE THIS TO GENERALIZE
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 33
      INTEGER*4 RCODE
      INTEGER*4 RECDIM
C     ****VARIABLES FOR THIS NETCDF FILE****
C
      CHARACTER*1 collSite                       (   4,NREC)
      CHARACTER*1 reportType                     (   4,NREC)
      REAL*4      lat                            (   4,NREC)
      REAL*4      lon                            (   4,NREC)
      CHARACTER*1 locStr                         (  16,   4,NREC)
      REAL*8      timeObs                        (NREC)
      REAL*8      recptTime                      (NREC)
      REAL*4      fltLvlBottom                   (NREC)
      REAL*4      fltLvlTop                      (NREC)
      CHARACTER*1 fltLvlInd                      (   4,NREC)
      CHARACTER*1 aircraftType                   (   5,NREC)
      REAL*4      skyCvrBottom                   (   3,NREC)
      REAL*4      skyCvrTop                      (   3,NREC)
      CHARACTER*1 skyCvrHtInd                    (   4,   3,NREC)
      CHARACTER*1 skyCvrAmt                      (   8,   3,NREC)
      CHARACTER*1 fltWthr                        (  33,NREC)
      INTEGER*4   vis                            (NREC)
      INTEGER*4   temptr                         (NREC)
      INTEGER*4   windDir                        (NREC)
      INTEGER*4   windSpd                        (NREC)
      REAL*4      turbBottom                     (   2,NREC)
      REAL*4      turbTop                        (   2,NREC)
      CHARACTER*1 turbHtInd                      (   4,   2,NREC)
      INTEGER*4   turbIntens                     (   2,   2,NREC)
      CHARACTER*1 turbFreq                       (   4,   2,   2,NREC)
      CHARACTER*1 turbType                       (   5,   2,NREC)
      REAL*4      icingBottom                    (   2,NREC)
      REAL*4      icingTop                       (   2,NREC)
      CHARACTER*1 icingHtInd                     (   4,   2,NREC)
      INTEGER*4   icingIntens                    (   2,NREC)
      CHARACTER*1 icingType                      (   5,   2,NREC)
      INTEGER*4   lowLvlWndShr                   (NREC)
      CHARACTER*1 origReport                     ( 201,NREC)
C*************************************
      INTEGER*4 START(10)
      INTEGER*4 COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY

!     Non-automatic declarations
      character*70 filename
      character*9 a9_timeObs,a9_recptTime 
      character*7 c7_skycover
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

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
C    statements to fill collSite                       
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
     +collSite                       ,LENSTR,RCODE)
C
C    statements to fill reportType                     
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
     +reportType                     ,LENSTR,RCODE)
C
C    statements to fill lat                            
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
     +lat                            ,RCODE)
C
C    statements to fill lon                            
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
     +lon                            ,RCODE)
C
C    statements to fill locStr                         
C
      CALL NCVINQ(NCID, 5,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  50 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  50  CONTINUE
      CALL NCVGTC(NCID, 5,START,COUNT,
     +locStr                         ,LENSTR,RCODE)
C
C    statements to fill timeObs                        
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
     +timeObs                        ,RCODE)
C
C    statements to fill recptTime                      
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
     +recptTime                      ,RCODE)
C
C    statements to fill fltLvlBottom                   
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
     +fltLvlBottom                   ,RCODE)
C
C    statements to fill fltLvlTop                      
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
     +fltLvlTop                      ,RCODE)
C
C    statements to fill fltLvlInd                      
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
     +fltLvlInd                      ,LENSTR,RCODE)
C
C    statements to fill aircraftType                   
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
     +aircraftType                   ,LENSTR,RCODE)
C
C    statements to fill skyCvrBottom                   
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
     +skyCvrBottom                   ,RCODE)
C
C    statements to fill skyCvrTop                      
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
     +skyCvrTop                      ,RCODE)
C
C    statements to fill skyCvrHtInd                    
C
      CALL NCVINQ(NCID,14,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 140 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 140  CONTINUE
      CALL NCVGTC(NCID,14,START,COUNT,
     +skyCvrHtInd                    ,LENSTR,RCODE)
C
C    statements to fill skyCvrAmt                      
C
      CALL NCVINQ(NCID,15,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 150 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 150  CONTINUE
      CALL NCVGTC(NCID,15,START,COUNT,
     +skyCvrAmt                      ,LENSTR,RCODE)
C
C    statements to fill fltWthr                        
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
     +fltWthr                        ,LENSTR,RCODE)
C
C    statements to fill vis                            
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
     +vis                            ,RCODE)
C
C    statements to fill temptr                         
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
     +temptr                         ,RCODE)
C
C    statements to fill windDir                        
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
     +windDir                        ,RCODE)
C
C    statements to fill windSpd                        
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
     +windSpd                        ,RCODE)
C
C    statements to fill turbBottom                     
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
     +turbBottom                     ,RCODE)
C
C    statements to fill turbTop                        
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
     +turbTop                        ,RCODE)
C
C    statements to fill turbHtInd                      
C
      CALL NCVINQ(NCID,23,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 230 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 230  CONTINUE
      CALL NCVGTC(NCID,23,START,COUNT,
     +turbHtInd                      ,LENSTR,RCODE)
C
C    statements to fill turbIntens                     
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
     +turbIntens                     ,RCODE)
C
C    statements to fill turbFreq                       
C
      CALL NCVINQ(NCID,25,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 250 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 250  CONTINUE
      CALL NCVGTC(NCID,25,START,COUNT,
     +turbFreq                       ,LENSTR,RCODE)
C
C    statements to fill turbType                       
C
      CALL NCVINQ(NCID,26,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 260 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 260  CONTINUE
      CALL NCVGTC(NCID,26,START,COUNT,
     +turbType                       ,LENSTR,RCODE)
C
C    statements to fill icingBottom                    
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
     +icingBottom                    ,RCODE)
C
C    statements to fill icingTop                       
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
     +icingTop                       ,RCODE)
C
C    statements to fill icingHtInd                     
C
      CALL NCVINQ(NCID,29,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 290 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 290  CONTINUE
      CALL NCVGTC(NCID,29,START,COUNT,
     +icingHtInd                     ,LENSTR,RCODE)
C
C    statements to fill icingIntens                    
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
     +icingIntens                    ,RCODE)
C
C    statements to fill icingType                      
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
     +icingType                      ,LENSTR,RCODE)
C
C    statements to fill lowLvlWndShr                   
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
     +lowLvlWndShr                   ,RCODE)
C
C    statements to fill origReport                     
C
      CALL NCVINQ(NCID,33,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 330 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 330  CONTINUE
      CALL NCVGTC(NCID,33,START,COUNT,
     +origReport                     ,LENSTR,RCODE)
C
C     HERE IS WHERE YOU WRITE STATEMENTS TO USE THE DATA
C
C
C
!     Write All Pireps to LAPS PIN file
      altitude = 0.

      r_nc_missing_data = 1e20

      if(NRECS .gt. NREC)then
          write(6,*)' Warning: NRECS > NREC ',NRECS,NREC
      else
          write(6,*)' NRECS, NREC = ',NRECS,NREC
      endif

      num_pireps = min(NREC,NRECS)

      do i = 1,num_pireps

          write(6,*)
          write(6,*)' Pirep #',i

          if(lat(1,i) .ge. r_nc_missing_data)then
              write(6,*)' Missing first latitude',i
              goto 999
          endif
          if(lon(1,i) .ge. r_nc_missing_data)then
              write(6,*)' Missing first longitude',i
              goto 999
          endif
          if(timeObs(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing ob time',i
              goto 999
          endif
          if(recptTime(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing received time',i
              goto 999
          endif

!         Test to see how many lat/lons are present
          n_latitude_present = 0
          do j = 1,4
            if(lat(j,i) .lt. r_nc_missing_data)then
              write(6,*)' lat/lon = ',lat(j,i),lon(j,i)
              n_latitude_present = n_latitude_present + 1
            endif
          enddo ! j

          write(6,*)' Num locations = ',n_latitude_present       

          if(n_latitude_present .gt. 1)then
              write(6,*)' Multiple locations, reject ob',i
              goto 999
          endif

          if(lat(1,i) .le. rnorth .and. lat(1,i) .ge. south .and.
     1       lon(1,i) .ge. west   .and. lon(1,i) .le. east      )then        
              write(6,*)' Pirep is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' Outside domain lat/lon perimeter - reject'
              goto 999
          endif

!         Write Individual Pirep to LAPS PIN file
          call c_time2fname(nint(timeObs(i)),a9_timeObs)
          call c_time2fname(nint(recptTime(i)),a9_recptTime)

          call cv_asc_i4time(a9_timeObs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then
              write(6,*)' Outside time window - reject '
     1                              ,a9_timeObs,i4_resid
              goto 999        
          endif


          write(6,1)a9_timeObs,a9_recptTime 
          write(11,1)a9_timeObs,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)lat(1,i),lon(1,i),altitude
          write(11,2)lat(1,i),lon(1,i),altitude
 2        format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  

          write(6,33)
          write(11,33)
 33       format(' Cloud layer')

!         Write out cloud base/top in feet and cloud amount in eighths
          do ilyr = 1,3

!             Use -1000. for missing value of cloud base/top
              if(skyCvrBottom(ilyr,i) .ge. 1e10)then
                  rbase = -1000.
              else
                  rbase = skyCvrBottom(ilyr,i) * 3.281
              endif
              if(skyCvrTop(ilyr,i) .ge. 1e10)then
                  rtop = -1000.
              else
                  rtop = skyCvrTop(ilyr,i) * 3.281
              endif

              do j = 1,7
                  c7_skycover(j:j) = skyCvrAmt(j,ilyr,i)
              enddo ! j

              call skycover_to_frac(c7_skycover,fraction,istatus)

!             Test for missing or unusable cloud cover
              if(istatus .ne. 1)then 
                  ieighths = -999
              else 
                  ieighths = nint(fraction * 8.0)
              endif                  

              write(6,3)rbase,rtop,ieighths,(skyCvrAmt(j,ilyr,i),j=1,8)
     1                                     ,fraction
 3            format(' Cloud layer',2f8.0,i5,1x,8a1,f5.2)

!             if(istatus .eq. 1)then
!                 write(6,*)' Above layer written to PIN file'
                  write(11,3)rbase,rtop,ieighths
!             endif

          enddo ! ilyr

 999      continue

      enddo ! i

      return
      END


