
      subroutine get_acars_data(i4time_sys,i4_acars_window
     1                                    ,NX_L,NY_L
     1                                    ,c8_project
     1                                    ,ext
     1                                    ,l_use_tamdar
     1                                    ,filename,istatus)

      character*(*) filename,ext
      logical l_use_tamdar

!.............................................................................

      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status
      character*8 c8_project

C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename
      endif
C
C  Fill all dimension values
C
C
C Get size of recNum
C
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      call acars_sub(nf_fid, recNum, c8_project, ext, l_use_tamdar,
!.............................................................................
     1              i4time_sys,i4_acars_window,NX_L,NY_L,istatus)
      return
!.............................................................................
      end
C
C
      subroutine acars_sub(nf_fid, recNum, c8_project, ext, 
     1                     l_use_tamdar,
!.............................................................................
     1              i4time_sys,i4_acars_window,NX_L,NY_L,istatus)
!.............................................................................

      include 'netcdf.inc'
      character*8 c8_project
      character*(*) ext
      integer recNum,nf_fid, nf_vid, nf_status
      integer airline(recNum), bounceError(recNum),
     +     correctedFlag(recNum), dataDescriptor(recNum),
     +     dataSource(recNum),
     +     errorType(recNum), interpolatedLL(recNum),
     +     interpolatedTime(recNum), missingInputMinutes,
     +     rollFlag(recNum), speedError(recNum), tempError(recNum),
     +     waterVaporQC(recNum), windDirError(recNum),
     +     windSpeedError(recNum)
      real altitude(recNum), heading(recNum), latitude(recNum),
     +     longitude(recNum), maxTurbulence(recNum),
     +     medTurbulence(recNum), temperature(recNum),
     +     vertAccel(recNum), downlinkedRH(recNum), windDir(recNum),
     +     windSpeed(recNum)
      double precision maxSecs, minSecs, timeObs(recNum),
     +     timeReceived(recNum)
      character*4 rptStation(recNum)
      character*6 destAirport(recNum)
      character*30 minDate
      character*6 origAirport(recNum)
      character*9 tailNumber(recNum)
      character*30 maxDate
      character*13 flight(recNum)

!.............................................................................

      character*9 a9_timeObs,a9_recptTime 
      character*7 c7_skycover
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)
      logical l_use_tamdar, l_debug

!............................................................................

      if (c8_project(1:3) .eq. 'WFO' .or. 
     1    c8_project(1:3) .eq. 'RSA'      ) then     
        call read_acars_netcdf_wfo(nf_fid, recNum, airline, bounceError, 
     +     correctedFlag, dataDescriptor, errorType, interpolatedLL, 
     +     interpolatedTime, missingInputMinutes, rollFlag, 
     +     speedError, tempError, waterVaporQC, windDirError, 
     +     windSpeedError, altitude, heading, latitude, longitude, 
     +     maxTurbulence, medTurbulence, temperature, vertAccel, 
     +     downlinkedRH, windDir, windSpeed, maxSecs, minSecs, 
     +     timeObs, timeReceived, destAirport, flight, maxDate, 
     +     minDate, origAirport, rptStation, tailNumber)

         dataSource = 0

      else
        call read_acars_netcdf(nf_fid, recNum, airline, bounceError, 
     +     correctedFlag, dataDescriptor, dataSource,
     +     errorType, interpolatedLL, 
     +     interpolatedTime, missingInputMinutes, rollFlag, 
     +     speedError, tempError, waterVaporQC, windDirError, 
     +     windSpeedError, altitude, heading, latitude, longitude, 
     +     maxTurbulence, medTurbulence, temperature, vertAccel, 
     +     downlinkedRH, windDir, windSpeed, maxSecs, minSecs, 
     +     timeObs, timeReceived, destAirport, flight, maxDate, 
     +     minDate, origAirport, rptStation, tailNumber)

      endif
C
C The netcdf variables are filled - your code goes here
C
!............................................................................

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_laps_perimeter'
          return
      endif

      num_acars = recNum
      write(6,*)' # of acars = ',recNum

      do i = 1,num_acars

          if(i .le. 1000 .or. i .eq. ((i/10)*10) )then
              l_debug = .true.
          else
              l_debug = .false.
          endif

          if(l_debug)write(6,*)
          if(l_debug)write(6,*)' acars #',i,'  ',char(dataDescriptor(i))       
     1                                          ,char(errorType(i))
     1                                          ,dataSource(i)
!         write(6,*)' location = '
!    1             ,latitude(i),longitude(i),altitude(i)

          if((.not. l_use_tamdar) .and. dataSource(i) .eq. 4)then
              if(l_debug)write(6,*)' TAMDAR observation - reject'
              goto 900
          endif

          if(latitude(i) .le. rnorth .and. latitude(i)  .ge. south .and.
     1       longitude(i) .ge. west  .and. longitude(i) .le. east      
     1                                                             )then       
              continue
          else ! Outside lat/lon perimeter - reject
              if(l_debug)write(6,*)' lat/lon - reject'       
!    1                 ,latitude(i),longitude(i)
              goto 900
          endif

          if(nanf(altitude(i)) .eq. 1)then
              if(l_debug)write(6,*)' Altitude failed Nan test - reject'       
     1                            ,altitude(i)
              goto 900
          endif

          if(altitude(i) .gt. 20000. .or. altitude(i) .lt. -1000.
     1                               .or. altitude(i) .eq.     0.)then
              if(l_debug)write(6,*)' Altitude is suspect - reject'
     1                            ,altitude(i)
              goto 900
          endif

          if(abs(timeObs(i))      .lt. 3d9       .and.
     1       abs(timereceived(i)) .lt. 3d9              )then
              call c_time2fname(nint(timeObs(i)),a9_timeObs)
              call c_time2fname(nint(timereceived(i)),a9_recptTime)
          else
              if(l_debug)write(6,*)' Bad observation time - reject'       
     1                   ,timeObs(i),timereceived(i)
              goto 900
          endif


          call cv_asc_i4time(a9_timeObs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_acars_window)then ! outside time window
              if(l_debug)write(6,*)' time - reject '
     1           ,a9_timeObs,i4_resid,i4_acars_window
              goto 900        
          endif

          lun = 31

          call open_ext(lun,i4time_sys,ext(1:3),istatus)       

          if(l_debug)write(6,1)a9_timeObs,a9_recptTime 
                     write(lun,1)a9_timeObs,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

!         l_geoalt is implicitly false with pressure altitude data
          if(l_debug)write(6,2)latitude(i),longitude(i),altitude(i)
          write(lun,2)          latitude(i),longitude(i),altitude(i)
 2        format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  

!         Test for bad winds
          if(char(dataDescriptor(i)) .eq. 'X')then
            if(char(errorType(i)) .eq. 'W' .or. 
     1         char(errorType(i)) .eq. 'B'                         )then
              if(l_debug)write(6,*)' QC flag is bad - reject wind'
     1                 ,char(dataDescriptor(i)),char(errorType(i))
              goto 850
            endif
          endif

          if(abs(windSpeed(i)) .gt. 250.)then
              if(l_debug)write(6,*)' wind speed is suspect - reject'
     1                              ,windSpeed(i)

          elseif(int(windDir(i)).lt.0 .or. int(windDir(i)).gt.360)then     
              if(l_debug)write(6,*)' wind direction is suspect - reject'       
     1                              ,windDir(i)

          else ! write out valid wind
              if(l_debug)write(6 ,3)int(windDir(i)),windSpeed(i)
                         write(lun,3)int(windDir(i)),windSpeed(i)
 3            format(' Wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')

          endif

 850      continue

!         Test for bad temps
          if(char(dataDescriptor(i)) .eq. 'X')then
            if(char(errorType(i)) .eq. 'T' .or. 
     1         char(errorType(i)) .eq. 'B'                         )then
              if(l_debug)write(6,*)' QC flag is bad - reject temp'
     1                 ,char(dataDescriptor(i)),char(errorType(i))
              goto 860
            endif
          endif

          if(abs(temperature(i)) .lt. 400.)then
              if(l_debug)write(6,13)temperature(i)
                         write(lun,13)temperature(i)
 13           format(' Temp:'/1x,f10.1)
       
          else
              if(l_debug)write(6,*)' Temperature is suspect - reject'
     1                             , temperature(i)

          endif

 860      continue

          if(downlinkedRH(i) .ge. 0.   .and. 
     1       downlinkedRH(i) .le. 1.00 .and.
     1       waterVaporQC(i) .le. 2    .and.
     1       waterVaporQC(i) .ge. 0             )then
              if(l_debug)write(6,23)downlinkedRH(i)
              if(l_debug)write(6,*)' RH QC value = ',waterVaporQC(i)
              write(lun,23)downlinkedRH(i)
 23           format(' RH:'/1x,f10.3)

          else
              if(l_debug)write(6,*)' RH rejected: '
     1                             ,downlinkedRH(i),waterVaporQC(i)

          endif

 900  enddo ! i

!............................................................................

      return
      end
C
C  Subroutine to read the file "ACARS data" 
C
      subroutine read_acars_netcdf(nf_fid, recNum, airline, bounceError, 
     +     correctedFlag, dataDescriptor, dataSource,
     +     errorType, interpolatedLL, 
     +     interpolatedTime, missingInputMinutes, rollFlag, 
     +     speedError, tempError, waterVaporQC, windDirError, 
     +     windSpeedError, altitude, heading, latitude, longitude, 
     +     maxTurbulence, medTurbulence, temperature, vertAccel, 
     +     downlinkedRH, windDir, windSpeed, maxSecs, minSecs, 
     +     timeObs, timeReceived, destAirport, flight, maxDate, 
     +     minDate, origAirport, rptStation, tailNumber)
C
      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status
      integer airline(recNum), bounceError(recNum),
     +     correctedFlag(recNum), dataDescriptor(recNum),
     +     dataSource(recNum),
     +     errorType(recNum), interpolatedLL(recNum),
     +     interpolatedTime(recNum), missingInputMinutes,
     +     rollFlag(recNum), speedError(recNum), tempError(recNum),
     +     waterVaporQC(recNum), windDirError(recNum),
     +     windSpeedError(recNum)
      real altitude(recNum), heading(recNum), latitude(recNum),
     +     longitude(recNum), maxTurbulence(recNum),
     +     medTurbulence(recNum), temperature(recNum),
     +     vertAccel(recNum), downlinkedRH(recNum), windDir(recNum),
     +     windSpeed(recNum)
      double precision maxSecs, minSecs, timeObs(recNum),
     +     timeReceived(recNum)
      character*4 rptStation(recNum)
      character*6 destAirport(recNum)
      character*30 minDate
      character*6 origAirport(recNum)
      character*9 tailNumber(recNum)
      character*30 maxDate
      character*13 flight(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      altitude     
C
        nf_status = NF_INQ_VARID(nf_fid,'altitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,altitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altitude'
      endif
C
C     Variable        NETCDF Long Name
C      heading      "heading of flight path over ground"
C
        nf_status = NF_INQ_VARID(nf_fid,'heading',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var heading'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,heading)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var heading'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     
C
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
C
C     Variable        NETCDF Long Name
C      maxTurbulence"Maximum eddy dissipation rate"
C
        nf_status = NF_INQ_VARID(nf_fid,'maxTurbulence',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxTurbulence'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,maxTurbulence)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxTurbulence'
      endif
C
C     Variable        NETCDF Long Name
C      medTurbulence"Median eddy dissipation rate"
C
        nf_status = NF_INQ_VARID(nf_fid,'medTurbulence',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var medTurbulence'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,medTurbulence)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var medTurbulence'
      endif
C
C     Variable        NETCDF Long Name
C      temperature  
C
        nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      vertAccel    "peak vertical acceleration"
C
        nf_status = NF_INQ_VARID(nf_fid,'vertAccel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vertAccel'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,vertAccel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vertAccel'
      endif
C
C     Variable        NETCDF Long Name
C      downlinkedRH "downlinked relative humidity"
C
        nf_status = NF_INQ_VARID(nf_fid,'downlinkedRH',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var downlinkedRH'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,downlinkedRH)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var downlinkedRH'
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "Wind Direction"
C
        nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "Wind Speed"
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      airline      "Airline"
C
        nf_status = NF_INQ_VARID(nf_fid,'airline',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var airline'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,airline)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var airline'
      endif
C
C     Variable        NETCDF Long Name
C      bounceError  "Aircraft altitude variance error"
C
        nf_status = NF_INQ_VARID(nf_fid,'bounceError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bounceError'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,bounceError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bounceError'
      endif
C
C     Variable        NETCDF Long Name
C      correctedFlag"Corrected data indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'correctedFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var correctedFlag'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,correctedFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var correctedFlag'
      endif
C
C     Variable        NETCDF Long Name
C      dataDescriptor"AWIPS-type data descriptor"
C
        nf_status = NF_INQ_VARID(nf_fid,'dataDescriptor',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataDescriptor'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,dataDescriptor)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataDescriptor'
      endif
C
C     Variable        NETCDF Long Name
C      dataSource"AWIPS-type data source"
C
        nf_status = NF_INQ_VARID(nf_fid,'dataSource',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataSource'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,dataSource)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataSource'
      endif
C
C     Variable        NETCDF Long Name
C      errorType    
C
        nf_status = NF_INQ_VARID(nf_fid,'errorType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var errorType'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,errorType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var errorType'
      endif
C
C     Variable        NETCDF Long Name
C      interpolatedLL"UPS ascent/descent lat&lon interpolation indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'interpolatedLL',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var interpolatedLL'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,interpolatedLL)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var interpolatedLL'
      endif
C
C     Variable        NETCDF Long Name
C      interpolatedTime"UPS ascent/descent time interpolation indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'interpolatedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var interpolatedTime'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,interpolatedTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var interpolatedTime'
      endif
C
C     Variable        NETCDF Long Name
C      missingInputMinutes"missing minutes of input data"
C
        nf_status = NF_INQ_VARID(nf_fid,'missingInputMinutes',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var missingInputMinutes'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,missingInputMinutes)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var missingInputMinutes'
      endif
C
C     Variable        NETCDF Long Name
C      rollFlag     "Aircraft roll angle flag "
C
        nf_status = NF_INQ_VARID(nf_fid,'rollFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rollFlag'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,rollFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rollFlag'
      endif
C
C     Variable        NETCDF Long Name
C      speedError   "Aircraft ground speed error"
C
        nf_status = NF_INQ_VARID(nf_fid,'speedError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var speedError'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,speedError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var speedError'
      endif
C
C     Variable        NETCDF Long Name
C      tempError    
C
        nf_status = NF_INQ_VARID(nf_fid,'tempError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempError'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,tempError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempError'
      endif
C
C     Variable        NETCDF Long Name
C      waterVaporQC "water vapor mixing ratio QC character"
C
        nf_status = NF_INQ_VARID(nf_fid,'waterVaporQC',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVaporQC'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,waterVaporQC)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVaporQC'
      endif
C
C     Variable        NETCDF Long Name
C      windDirError 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDirError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirError'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,windDirError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirError'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedError
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeedError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedError'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedError'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      maxSecs      "maximum observation time"
C
        nf_status = NF_INQ_VARID(nf_fid,'maxSecs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxSecs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,maxSecs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxSecs'
      endif
C
C     Variable        NETCDF Long Name
C      minSecs      "minimum observation time"
C
        nf_status = NF_INQ_VARID(nf_fid,'minSecs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var minSecs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,minSecs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var minSecs'
      endif
C
C     Variable        NETCDF Long Name
C      timeObs      "time of observation"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
C
C     Variable        NETCDF Long Name
C      timeReceived "time data was received at ground station"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeReceived',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeReceived'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeReceived)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeReceived'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      destAirport  "Destination Airport"
C
        nf_status = NF_INQ_VARID(nf_fid,'destAirport',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var destAirport'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,destAirport)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var destAirport'
      endif
C
C     Variable        NETCDF Long Name
C      flight       "Flight number"
C
        nf_status = NF_INQ_VARID(nf_fid,'flight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var flight'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,flight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var flight'
      endif
C
C     Variable        NETCDF Long Name
C      maxDate      "maximum observation date"
C
        nf_status = NF_INQ_VARID(nf_fid,'maxDate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxDate'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,maxDate)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var maxDate'
      endif
C
C     Variable        NETCDF Long Name
C      minDate      "minimum observation date"
C
        nf_status = NF_INQ_VARID(nf_fid,'minDate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var minDate'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,minDate)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var minDate'
      endif
C
C     Variable        NETCDF Long Name
C      origAirport  "Originating Airport"
C
        nf_status = NF_INQ_VARID(nf_fid,'origAirport',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origAirport'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,origAirport)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origAirport'
      endif
C
C     Variable        NETCDF Long Name
C      rptStation   "Station reporting through"
C
        nf_status = NF_INQ_VARID(nf_fid,'rptStation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rptStation'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,rptStation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rptStation'
      endif
C
C     Variable        NETCDF Long Name
C      tailNumber   "tail number"
C
        nf_status = NF_INQ_VARID(nf_fid,'tailNumber',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tailNumber'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,tailNumber)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tailNumber'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
C
C  Subroutine to read WFO format file "ACARS data" 
C
      subroutine read_acars_netcdf_wfo(nf_fid, recNum, airline, 
     +     bounceError, 
     +     correctedFlag, dataDescriptor, errorType, interpolatedLL, 
     +     interpolatedTime, missingInputMinutes, rollFlag, 
     +     speedError, tempError, waterVaporQC, windDirError, 
     +     windSpeedError, altitude, heading, latitude, longitude, 
     +     maxTurbulence, medTurbulence, temperature, vertAccel, 
     +     downlinkedRH, windDir, windSpeed, maxSecs, minSecs, 
     +     timeObs, timeReceived, destAirport, flight, maxDate, 
     +     minDate, origAirport, rptStation, tailNumber)
C
      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status
      integer airline(recNum), bounceError(recNum),
     +     correctedFlag(recNum), dataDescriptor(recNum),
     +     errorType(recNum), interpolatedLL(recNum),
     +     interpolatedTime(recNum), missingInputMinutes,
     +     rollFlag(recNum), speedError(recNum), tempError(recNum),
     +     waterVaporQC(recNum), windDirError(recNum),
     +     windSpeedError(recNum)
      real altitude(recNum), heading(recNum), latitude(recNum),
     +     longitude(recNum), maxTurbulence(recNum),
     +     medTurbulence(recNum), temperature(recNum),
     +     vertAccel(recNum), downlinkedRH(recNum), windDir(recNum),
     +     windSpeed(recNum)
      double precision maxSecs, minSecs, timeObs(recNum),
     +     timeReceived(recNum)
      character*4 rptStation(recNum)
      character*6 destAirport(recNum)
      character*30 minDate
      character*6 origAirport(recNum)
      character*9 tailNumber(recNum)
      character*30 maxDate
      character*13 flight(recNum)
      integer i, wind_err, temp_err
      real ws, temp, alt, max_temp, min_temp, wmax

C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      altitude     
C
        nf_status = NF_INQ_VARID(nf_fid,'indAltitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var indAltitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,altitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altitude'
      endif
C
C     Variable        NETCDF Long Name
C      heading      "heading of flight path over ground"
C      heading not available in WFO file....fill with 99999.0
C
      do i = 1, recNum
        heading(i) = 99999.0
      enddo
C
C     Variable        NETCDF Long Name
C      latitude     
C
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
C
C     Variable        NETCDF Long Name
C      maxTurbulence"Maximum eddy dissipation rate"
C      maxTurbulence not available in WFO file....fill with -9.99 
C
      do i = 1, recNum
        maxTurbulence(i) = -9.99
      enddo
C
C     Variable        NETCDF Long Name
C      medTurbulence"Median eddy dissipation rate"
C      medTurbulence not available in WFO file....fill with -9.99 
C
      do i = 1, recNum
        medTurbulence(i) = -9.99
      enddo
C
C     Variable        NETCDF Long Name
C      temperature  
C
        nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      vertAccel    "peak vertical acceleration"
C      vertAccel not available in WFO file....fill with 0.0  
C
      do i = 1, recNum
        vertAccel(i) = 0.0 
      enddo
C
C     Variable        NETCDF Long Name
C      downlinkedRH "downlinked relative humidity"
C
        nf_status = NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,downlinkedRH)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
C
C reset values in downlinkedRH to 0-1 range from percent
C
      do i = 1, recNum
        downlinkedRH(i) = downlinkedRH(i)/100.0
      enddo
C
C     Variable        NETCDF Long Name
C      windDir      "Wind Direction"
C
        nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "Wind Speed"
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      airline      "Airline"
C      airline not available in WFO file....fill with 0
C
      do i = 1, recNum
        airline(i) = 0
      enddo
C
C     Variable        NETCDF Long Name
C      bounceError  "Aircraft altitude variance error"
C      bounceError not available in WFO file....fill with 45
C
      do i = 1, recNum
        bounceError(i) = 45
      enddo
C
C     Variable        NETCDF Long Name
C      correctedFlag"Corrected data indicator"
C      correctedFlag not available in WFO file....fill with 114
C
      do i = 1, recNum
        correctedFlag(i) = 114
      enddo
C
C     Variable        NETCDF Long Name
C      interpolatedLL "UPS ascent/descent lat&lon interpolation indicator"
C      interpolatedLL not available in WFO file....fill with 114
C
      do i = 1, recNum
        interpolatedLL(i) = 114
      enddo
C
C     Variable        NETCDF Long Name
C      interpolatedTime "UPS ascent/descent time interpolation indicator"
C      interpolatedTime not available in WFO file....fill with 114
C
      do i = 1, recNum
        interpolatedTime(i) = 114
      enddo
C
C     Variable        NETCDF Long Name
C      missingInputMinutes "missing minutes of input data"
C      missingInputMinutes not available in WFO file....fill with 0
C
      missingInputMinutes = 0
C
C     Variable        NETCDF Long Name
C      rollFlag     "Aircraft roll angle flag "
C
        nf_status = NF_INQ_VARID(nf_fid,'rollQuality',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rollQuality'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,rollFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rollFlag'
      endif

      do i = 1, recNum
        if (rollFlag(i) .eq. 0) then
          rollFlag(i) = 71
        elseif (rollFlag(i) .eq. 1) then
          rollFlag(i) = 66
        else
          rollFlag(i) = 78
        endif
      enddo
 
C
C     Variable        NETCDF Long Name
C      speedError   "Aircraft ground speed error"
C      speedError not available in WFO file....fill with 45
C
      do i = 1, recNum
        speedError(i) = 45
      enddo
C
C     Variable        NETCDF Long Name
C      tempError    
C      tempError not available in WFO file....fill with 45
C
      do i = 1, recNum
        tempError(i) = 45
      enddo
C
C     Variable        NETCDF Long Name
C      windDirError 
C      windDirError not available in WFO file....fill with 45
C
      do i = 1, recNum
        windDirError(i) = 45
      enddo
C
C     Variable        NETCDF Long Name
C      windSpeedError
C      windSpeedError not available in WFO file....fill with 45
C
      do i = 1, recNum
        windSpeedError(i) = 45
      enddo
C
C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      maxSecs      "maximum observation time"
C      maxSecs not available in WFO file....fill with 0
C
      maxSecs = 0
C
C     Variable        NETCDF Long Name
C      minSecs      "minimum observation time"
C      minSecs not available in WFO file....fill with 0
C
      minSecs = 0
C
C     Variable        NETCDF Long Name
C      timeObs      "time of observation"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
C
C     Variable        NETCDF Long Name
C      timeReceived "time data was received at ground station"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeReceived',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeReceived'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeReceived)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeReceived'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      destAirport  "Destination Airport"
C      destAirport not available in WFO file....fill with "   " (3 spaces)
C
      do i = 1, recNum
        destAirport(i) = "   "
      enddo
C
C     Variable        NETCDF Long Name
C      flight       "Flight number"
C
        nf_status = NF_INQ_VARID(nf_fid,'flight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var flight'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,flight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var flight'
      endif
C
C     Variable        NETCDF Long Name
C      maxDate      "maximum observation date"
C      maxDate not available in WFO file....fill with 30 blank spaces
C
C     Variable        NETCDF Long Name
C      minDate      "minimum observation date"
C      minDate not available in WFO file....fill with 30 blank spaces
C
      maxDate = "                              "
      minDate = "                              "
C
C     Variable        NETCDF Long Name
C      origAirport  "Originating Airport"
C      origAirport not available in WFO file....fill with "   " (3 spaces)
C
      do i = 1, recNum
        origAirport(i) = "   "
      enddo
C
C     Variable        NETCDF Long Name
C      rptStation   "Station reporting through"
C      rptStation not available in WFO file....fill with "    " (4 spaces)
C
      do i = 1, recNum
        rptStation(i) = "    "
      enddo
C
C     Variable        NETCDF Long Name
C      tailNumber   "tail number"
C
        nf_status = NF_INQ_VARID(nf_fid,'tailNumber',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tailNumber'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,tailNumber)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tailNumber'
      endif

C     used by acars program...not available on wfo
C       use FSL ACARS Quality Controls to fill
C
C     Max temp:  if altitude > 35000ft, T < -20C
C                else   T < 60 - 80 * (altitude / 35000)
C     Min temp:  if altitude < 18000ft, T > -60C
C                if altitude > 35000ft, T > -100C
C                else T > 60 - 40 * (altitude -18000) / 17000
C     Wind dir:  0 <= windDir <= 360
C     Wind Spd:  windSpeed(knots) >= 0
C     Max wind Spd:  bad = 0
C                    if (altitude < 30000.) {
C                      wmax = 70. + 230.*altitude / 30000.;
C                    } else if (altitude < 40000.) {
C                      wmax = 300.;
C                    } else if (altitude < 45000.) {
C                      wmax = 300. - 100 * (altitude - 40000.) / 5000. ;
C                    } else {
C                      wmax = 200.;
C                    }
C                    if (windSpeed > wmax) {
C                      bad = 1;
C                    }
C
C     dataDescriptor: 'R' if temp and wind within ranges above
C                     'X' if temp and wind failed any test    
C     if dataDescriptor = 'X':
C     errorType:      'W' for wind
C                     'T' for temperature
C                     'B' for both
C
      do i = 1, recNum
        wind_err = 0  ! no error
        temp_err = 0  ! no error
        temp = temperature(i) - 273.15  ! convert K to C

        if(abs(windSpeed(i)) .lt. 1000.)then
            ws = windSpeed(i) * 1.9438 ! convert m/s to knots
        endif

        if(abs(altitude(i)) .lt. 1e10)then
            alt = altitude(i) * 3.280839895 ! convert m to ft
        endif

        if (alt .gt. 35000.0) then
          max_temp = -20.0
        else
          max_temp = 80 * (alt / 35000.)
        endif
        
        if (alt .lt. 18000.0) then
          min_temp = -60.0
        elseif (alt .gt. 35000.0) then
          min_temp = -100.0
        else
          min_temp = 40 * (alt - 18000.) / 17000.
        endif
        
        if ((temp .le. min_temp) .or. (temp .ge. max_temp)) 
     1    temp_err = 1    ! bad temp

        if ((windDir(i) .lt. 0.0) .or. (windDir(i) .gt. 360.0))
     1    wind_err = 1    ! bad wind dir

        if (alt .lt. 30000.) then
          wmax = 70. + 230.*alt / 30000.
        elseif (alt .lt. 40000.) then
          wmax = 300.
        elseif (alt .lt. 45000.) then
          wmax = 300. - 100 * (alt - 40000.) / 5000. 
        else 
          wmax = 200.
        endif 

        if ((ws .lt. 0.0) .or. (ws .gt. wmax))
     1    wind_err = 1  ! bad wind speed

        if ((temp_err .eq. 0) .and. (wind_err .eq. 0)) then
          dataDescriptor(i) = 82
          errorType(i) = 46
        elseif ((temp_err .eq. 1) .and. (wind_err .eq. 1)) then
          dataDescriptor(i) = 88
          errorType(i) = 66
        elseif (temp_err .eq. 1) then
          dataDescriptor(i) = 88
          errorType(i) = 84
        else  ! (wind_err .eq. 1)
          dataDescriptor(i) = 88
          errorType(i) = 87
        endif
 
      enddo
C
C     Variable        NETCDF Long Name
C      waterVaporQC "water vapor mixing ratio QC character"
C      waterVaporQC not available in WFO file....fill with 0
C
      do i = 1, recNum
        waterVaporQC(i) = 0
      enddo
C
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
