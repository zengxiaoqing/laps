c
c
      subroutine read_local(nf_fid , recNum, dataProvider, dewpoint,
     +     elevation, latitude, longitude, reportTime, 
     +     stationId, stationPressChangeTime, stationPressure,
     +     tempChangeTime, temperature, windDir,
     +     windDirChangeTime, windGust, windGustChangeTime,
     +     windSpeed, windSpeedChangeTime, istatus)
c
c======================================================================
c
c     Routine to read the NetCDF local data files created by LDAD. 
c     Code created with 'xgennet.pl' by J. Edwards, NOAA/FSL.
c
c     Original:  P. Stamus, NOAA/FSL  02-05-98
c
c======================================================================
c
      include 'netcdf.inc'
      integer recNum, nf_fid, nf_vid, nf_status

      character*11 dataProvider(recNum)
      character*6 stationId(recNum)
      real dewpoint(recNum), elevation(recNum), latitude(recNum),
     +     longitude(recNum), 
     +     stationPressure(recNum), temperature(recNum),
     +     windDir(recNum), windGust(recNum), windSpeed(recNum)

      double precision 
     +     reportTime(recNum), stationPressChangeTime(recNum),
     +     tempChangeTime(recNum), windDirChangeTime(recNum),
     +     windGustChangeTime(recNum), windSpeedChangeTime(recNum)


	istatus = 0
C
C     Variable        NETCDF Long Name
C      dataProvider "LDAD data provider" 
C
        nf_status = NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ dataProvider '
      endif
C
C     Variable        NETCDF Long Name
C      stationId    "alphanumeric station Id" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationId '
      endif
C
C     Variable        NETCDF Long Name
C      dewpoint     "dew point temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dewpoint)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ dewpoint '
      endif
C
C     Variable        NETCDF Long Name
C      elevation    "elevation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ elevation '
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "latitude" 
C
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ latitude '
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "longitude" 
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ longitude '
      endif
C
C     Variable        NETCDF Long Name
C      stationPressure"station pressure" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationPressure '
      endif
C
C     Variable        NETCDF Long Name
C      temperature  "temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ temperature '
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "wind direction" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windDir '
      endif
C
C     Variable        NETCDF Long Name
C      windGust     "wind gust" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windGust',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windGust)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windGust '
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "wind speed" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windSpeed '
      endif
C
C     Variable        NETCDF Long Name
C      reportTime   "time data was processed by the provider" 
C
        nf_status = NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ reportTime '
      endif
C
C     Variable        NETCDF Long Name
C      stationPressChangeTime"station press time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationPressChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressChangeTime'
      endif
        nf_status = 
     &    NF_GET_VAR_DOUBLE(nf_fid,nf_vid,stationPressChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationPressChangeTime '
      endif
C
C     Variable        NETCDF Long Name
C      tempChangeTime"temperature time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'tempChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,tempChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ tempChangeTime '
      endif
C
C     Variable        NETCDF Long Name
C      windDirChangeTime"wind direction time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDirChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windDirChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windDirChangeTime '
      endif
C
C     Variable        NETCDF Long Name
C      windGustChangeTime"wind gust time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windGustChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGustChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windGustChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windGustChangeTime '
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedChangeTime"wind speed time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeedChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windSpeedChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windSpeedChangeTime '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

	istatus = 1
      return
      end
