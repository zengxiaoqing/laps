
      subroutine get_cloud_drift_afwa(i4time_sys,i4_ob_window
     1                                ,NX_L,NY_L
     1                                ,filename,istatus)

      character*(*) filename

!.............................................................................

      character*6 C6_A1ACID
      character*9 a9_timeObs,a9_recptTime 
!     character*7 c7_skycover
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)
      real*4 latitude,longitude

!............................................................................

      lun_in = 21
      lun_out = 11

      open(lun_in,file=filename,status='old')

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif
  
      i = 0

      do while (.true.)

          read(lun_in,101,err=890,end=999) !    NAME          UNITS & FACTOR
     1         I_A1CYCC,                       
     1         I_A1MIN,                ! TIME-REPORT-MINUTES
     1         I_A1JUL,                ! JULIAN HOUR          HR since 673650000
     1         I_A1LAT,                ! LATITUDE             DEG * 100
     1         I_A1LON,                ! LONGITUDE            DEG * -100 
     1         I_A1PRES,
     1         I_A1DUM1,
     1         I_A1DUM2,
     1         I_A1DUM3,
     1         I_A1DUM4,               ! FLIGHT-ALT (true)    Meters MSL 
     1         I_A1DUM5,
     1         I_A1DUM6,
     1         I_A1WD,                 ! WIND-DIRECTION       DEG
     1         I_A1WFLS                ! WIND-SPEED           M/S * 10
 101      format(14(i9,2x),i9)

          i = i + 1

          write(6,*)
          write(6,*)' cdw #',i

          latitude  =  float(I_A1LAT)/100.
          longitude = +float(I_A1LON)/100.
          pres_pa   =  I_A1ALT

          write(6,2)latitude,longitude,pres_pa
 2        format(' Lat, lon, pres_pa'/f8.3,f10.3,f8.0)  

          if(latitude  .le. rnorth .and. latitude  .ge. south .and.
     1       longitude .ge. west   .and. longitude .le. east      
     1                                                             )then       
              continue
          else ! Outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,latitude,longitude
              goto 900
          endif

          if(pres_pa .gt. 1e10)then
              write(6,*)' pres_pa is suspect - reject',pres_pa
              goto 900
          endif

          call afwa_julhr_i4time(I_A1JUL,I_A1MIN,i4time_ob)

          call make_fnam_lp(i4time_ob,a9_timeObs,istatus)
          if(istatus .ne. 1)goto900

          a9_recptTime = '         '

!         call cv_asc_i4time(a9_timeObs,i4time_ob)

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_ob_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_timeObs,i4_resid,i4_ob_window
              goto 900        
          endif

!         Test for bad winds
!         if(char(dataDescriptor) .eq. 'X')then
!           if(char(errorType) .eq. 'W' .or. 
!    1         char(errorType) .eq. 'B'                         )then
!             write(6,*)' QC flag is bad - reject '
!    1                 ,char(dataDescriptor),char(errorType)
!             goto 850
!           endif
!         endif

          windSpeed = float(I_A1WFLS) / 10.
          windDIR = I_A1WD

          if(abs(windSpeed) .gt. 250. .or. windDir .gt. 360.)then
              write(6,*)' wind is suspect - reject',windDir,windSpeed

          else ! write out valid wind
              write(6      ,21)latitude,longitude,pres_pa
     1                        ,winddir,windspd,a9_timeObs       
!             write(lun_out,21)latitude,longitude,pres_pa
!    1                        ,winddir,windspd,a9_timeObs       
 21           format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

          endif

          go to 900

 890      write(6,*)' Warning: read error'

 900  enddo ! read line of AFWA file

!............................................................................

 999  write(6,*)' End of AFWA file detected'

      close(lun_in)
      istatus = 1
      return
      end
