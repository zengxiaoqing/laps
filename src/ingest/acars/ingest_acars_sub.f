
      subroutine get_acars_data(i4time_sys,i4_acars_window
     1                                    ,NX_L,NY_L
     1                                    ,filename,istatus)

      character*(*) filename

!.............................................................................

      include 'netcdf.inc'
      integer recNum,nf_fid, nf_vid, nf_status
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
      call acars_sub(nf_fid, recNum,
!.............................................................................
     1              i4time_sys,i4_acars_window,NX_L,NY_L,istatus)
      return
!.............................................................................
      end
C
C
      subroutine acars_sub(nf_fid, recNum,
!.............................................................................
     1              i4time_sys,i4_acars_window,NX_L,NY_L,istatus)
!.............................................................................

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
     +     vertAccel(recNum), waterVaporMR(recNum), windDir(recNum),
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
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

!............................................................................

      call read_acars_netcdf(nf_fid, recNum, airline, bounceError, 
     +     correctedFlag, dataDescriptor, errorType, interpolatedLL, 
     +     interpolatedTime, missingInputMinutes, rollFlag, 
     +     speedError, tempError, waterVaporQC, windDirError, 
     +     windSpeedError, altitude, heading, latitude, longitude, 
     +     maxTurbulence, medTurbulence, temperature, vertAccel, 
     +     waterVaporMR, windDir, windSpeed, maxSecs, minSecs, 
     +     timeObs, timeReceived, destAirport, flight, maxDate, 
     +     minDate, origAirport, rptStation, tailNumber)
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

          write(6,*)
          write(6,*)' acars #',i,'  ',char(dataDescriptor(i))
     1                               ,char(errorType(i))
!         write(6,*)' location = '
!    1             ,latitude(i),longitude(i),altitude(i)


          if(latitude(i) .le. rnorth .and. latitude(i)  .ge. south .and.
     1       longitude(i) .ge. west  .and. longitude(i) .le. east      
     1                                                             )then       
              continue
          else ! Outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,latitude(i),longitude(i)
              goto 900
          endif

          if(altitude(i) .gt. 20000.)then
              write(6,*)' Altitude is suspect - reject',altitude(i)
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
          if(i4_resid .gt. i4_acars_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9_timeObs,i4_resid,i4_acars_window
              goto 900        
          endif

          write(6,1)a9_timeObs,a9_recptTime 
          write(11,1)a9_timeObs,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)latitude(i),longitude(i),altitude(i)
          write(11,2)latitude(i),longitude(i),altitude(i)
 2        format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  

!         Test for bad winds
          if(char(dataDescriptor(i)) .eq. 'X')then
            if(char(errorType(i)) .eq. 'W' .or. 
     1         char(errorType(i)) .eq. 'B'                         )then
              write(6,*)' QC flag is bad - reject wind'
     1                 ,char(dataDescriptor(i)),char(errorType(i))
              goto 850
            endif
          endif

          if(abs(windSpeed(i)) .gt. 250.)then
              write(6,*)' wind speed is suspect - reject',windSpeed(i)

          else ! write out valid wind
              write(6,3)int(windDir(i)),windSpeed(i)
              write(11,3)int(windDir(i)),windSpeed(i)
 3            format(' Wind:'/' ', i3, ' deg @ ', f6.1, ' m/s')

          endif

 850      continue

!         Test for bad temps
          if(char(dataDescriptor(i)) .eq. 'X')then
            if(char(errorType(i)) .eq. 'T' .or. 
     1         char(errorType(i)) .eq. 'B'                         )then
              write(6,*)' QC flag is bad - reject temp'
     1                 ,char(dataDescriptor(i)),char(errorType(i))
              goto 860
            endif
          endif

          if(abs(temperature(i)) .lt. 400.)then
              write(6,13)temperature(i)
              write(11,13)temperature(i)
 13           format(' Temp:'/1x,f10.1)
       
          else
              write(6,*)' Temperature is suspect - reject'
     1                , temperature(i)

          endif

 860      continue

          if(waterVaporMR(i) .ge. 0. .and. 
     1       waterVaporMR(i) .le. 100.)then
              write(6,23)waterVaporMR(i)
              write(11,23)waterVaporMR(i)
 23           format(' MixR:'/1x,f10.3)

          else
              write(6,*)' water vapor rejected: ',waterVaporMR(i)

          endif

 900  enddo ! i

!............................................................................

      return
      end
C
C  Subroutine to read the file "ACARS data" 
C
      subroutine read_acars_netcdf(nf_fid, recNum, airline, bounceError, 
     +     correctedFlag, dataDescriptor, errorType, interpolatedLL, 
     +     interpolatedTime, missingInputMinutes, rollFlag, 
     +     speedError, tempError, waterVaporQC, windDirError, 
     +     windSpeedError, altitude, heading, latitude, longitude, 
     +     maxTurbulence, medTurbulence, temperature, vertAccel, 
     +     waterVaporMR, windDir, windSpeed, maxSecs, minSecs, 
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
     +     vertAccel(recNum), waterVaporMR(recNum), windDir(recNum),
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
C      waterVaporMR "water vapor mixing ratio"
C
        nf_status = NF_INQ_VARID(nf_fid,'waterVaporMR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVaporMR'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,waterVaporMR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var waterVaporMR'
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
