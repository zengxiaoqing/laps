c
c
      subroutine read_local_obs(nf_fid , recNum, altimeter, 
     +     dataProvider, solarRadiation, seaSurfaceTemp,    
     +     soilTemperature,
     +     dewpoint, elevation, latitude, longitude, observationTime,
     +     presWeather, relHumidity, rhChangeTime, seaLevelPressure,
     +     stationId, stationPressChangeTime, stationPressure,
     +     stationType, tempChangeTime, temperature, visibility,
     +     windDir, windDirChangeTime, windDirMax, windGust,
     +     windGustChangeTime, windSpeed, windSpeedChangeTime,
     &     badflag, istatus)
c
c**********************************************************************
c
c     Routine to read the LDAD NetCDF mesonet observation files at FSL.
c     Code created with 'xgennet.pl' by J. Edwards, NOAA/FSL.
c     
c     Original:  P. Stamus, NOAA/FSL  28 Aug 1998
c	         09-30-98  P. Stamus
c                     Housekeeping changes.
c                10-18-99  P. Stamus
c                     Add check for missing values.
c
c**********************************************************************
c
      include 'netcdf.inc'
      integer recNum, nf_fid, nf_vid, nf_status !, ifilval

      character*11 dataProvider(recNum)
      character*25 presWeather(recNum)
      character*6 stationId(recNum)
      character*11 stationType(recNum)
      real altimeter(recNum), dewpoint(recNum), elevation(recNum),
     +     latitude(recNum), longitude(recNum), relHumidity(recNum),
     +     seaLevelPressure(recNum), stationPressure(recNum),
     +     temperature(recNum), visibility(recNum), windDir(recNum),
     +     windDirMax(recNum), windGust(recNum), windSpeed(recNum),
     +     solarRadiation(recNum), seaSurfaceTemp(recNum),
     +     soilTemperature(recNum)
      real filval, misval

      double precision observationTime(recNum), rhChangeTime(recNum),
     +     stationPressChangeTime(recNum), tempChangeTime(recNum),
     +     windDirChangeTime(recNum), windGustChangeTime(recNum),
     +     windSpeedChangeTime(recNum), dfilval, dmisval
c
c..... Start here.
c
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
C      presWeather  "present weather" 
C
        nf_status = NF_INQ_VARID(nf_fid,'presWeather',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var presWeather'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,presWeather)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ presWeather '
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
C      stationType  "LDAD station type" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationType '
      endif
C
C     Variable        NETCDF Long Name
C      altimeter    "altimeter setting" 
C
        nf_status = NF_INQ_VARID(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeter'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,altimeter)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ altimeter '
      endif
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var altimeter'
      endif
      call ck_array_real(altimeter, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var altimeter'
      endif
      call ck_array_real(altimeter, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var dewpoint'
      endif
      call ck_array_real(dewpoint, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var dewpoint'
      endif
      call ck_array_real(dewpoint, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      dewpoint     "solar radiation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'solarRadiation',nf_vid)       
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadiation'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,solarRadiation)       
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ solarRadiation '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var solarRadiation'
      endif
      call ck_array_real(solarRadiation, recNum, filval, badflag)       
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var solarRadiation'
      endif
      call ck_array_real(solarRadiation, recNum, misval, badflag)       
C
C     Variable        NETCDF Long Name
C      dewpoint     "sea surface temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'seaSurfaceTemp',nf_vid)       
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaSurfaceTemp'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,seaSurfaceTemp)       
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ seaSurfaceTemp '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaSurfaceTemp'
      endif
      call ck_array_real(seaSurfaceTemp, recNum, filval, badflag)       
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaSurfaceTemp'
      endif
      call ck_array_real(seaSurfaceTemp, recNum, misval, badflag)       
C
C     Variable        NETCDF Long Name
C      "soil temperature" 
C
        nf_status = NF_INQ_VARID(nf_fid,'soilTemperature',nf_vid)       
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilTemperature'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,soilTemperature)       
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ soilTemperature '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var soilTemperature'
      endif
      call ck_array_real(soilTemperature, recNum, filval, badflag)       
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var soilTemperature'
      endif
      call ck_array_real(soilTemperature, recNum, misval, badflag)       
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var elevation'
      endif
      call ck_array_real(elevation, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var elevation'
      endif
      call ck_array_real(elevation, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var latitude'
      endif
      call ck_array_real(latitude, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var latitude'
      endif
      call ck_array_real(latitude, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var longitude'
      endif
      call ck_array_real(longitude, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var longitude'
      endif
      call ck_array_real(longitude, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      relHumidity  "relative humidity" 
C
        nf_status = NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ relHumidity '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var relHumidity'
      endif
      call ck_array_real(relHumidity, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var relHumidity'
      endif
      call ck_array_real(relHumidity, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      seaLevelPressure"sea level pressure" 
C
        nf_status = NF_INQ_VARID(nf_fid,'seaLevelPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ seaLevelPressure '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaLevelPressure'
      endif
      call ck_array_real(seaLevelPressure, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var seaLevelPressure'
      endif
      call ck_array_real(seaLevelPressure, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressure'
      endif
      call ck_array_real(stationPressure, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressure'
      endif
      call ck_array_real(stationPressure, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var temperature'
      endif
      call ck_array_real(temperature, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var temperature'
      endif
      call ck_array_real(temperature, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      visibility   "visibility" 
C
        nf_status = NF_INQ_VARID(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,visibility)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ visibility '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var visibility'
      endif
      call ck_array_real(visibility, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var visibility'
      endif
      call ck_array_real(visibility, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDir'
      endif
      call ck_array_real(windDir, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDir'
      endif
      call ck_array_real(windDir, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      windDirMax   "wind direction at gust" 
C
        nf_status = NF_INQ_VARID(nf_fid,'windDirMax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windDirMax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ windDirMax '
      endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirMax'
      endif
      call ck_array_real(windDirMax, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirMax'
      endif
      call ck_array_real(windDirMax, recNum, misval, badflag)
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
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGust'
      endif
      call ck_array_real(windGust, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGust'
      endif
      call ck_array_real(windGust, recNum, misval, badflag)
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
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',filval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeed'
      endif
      call ck_array_real(windSpeed, recNum, filval, badflag)
       nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',misval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeed'
      endif
      call ck_array_real(windSpeed, recNum, misval, badflag)
C
C     Variable        NETCDF Long Name
C      observationTime"time of observation" 
C
        nf_status = NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ observationTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var observationTime'
      endif
      do i=1,recNum
        if(observationTime(i) .eq. dfilval) observationTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var observationTime'
      endif
      do i=1,recNum
        if(observationTime(i) .eq. dmisval) observationTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      rhChangeTime "relative humidity time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'rhChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,rhChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ rhChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var rhChangeTime'
      endif
      do i=1,recNum
         if(rhChangeTime(i) .eq. dfilval) rhChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var rhChangeTime'
      endif
      do i=1,recNum
         if(rhChangeTime(i) .eq. dmisval) rhChangeTime(i) = badflag
      enddo !i
C
C     Variable        NETCDF Long Name
C      stationPressChangeTime"station press time of last change" 
C
        nf_status = NF_INQ_VARID(nf_fid,'stationPressChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressChangeTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,
     &                                        stationPressChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ stationPressChangeTime '
      endif
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressChangeTime'
      endif
      do i=1,recNum
         if(stationPressChangeTime(i) .eq. dfilval) 
     &                             stationPressChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var stationPressChangeTime'
      endif
      do i=1,recNum
         if(stationPressChangeTime(i) .eq. dmisval) 
     &                             stationPressChangeTime(i) = badflag
      enddo !i
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
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var tempChangeTime'
      endif
      do i=1,recNum
         if(tempChangeTime(i) .eq. dfilval) tempChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var tempChangeTime'
      endif
      do i=1,recNum
         if(tempChangeTime(i) .eq. dmisval) tempChangeTime(i) = badflag
      enddo !i
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
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirChangeTime'
      endif
      do i=1,recNum
         if(windDirChangeTime(i) .eq. dfilval) 
     &                                 windDirChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windDirChangeTime'
      endif
      do i=1,recNum
         if(windDirChangeTime(i) .eq. dmisval) 
     &                                 windDirChangeTime(i) = badflag
      enddo !i
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
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGustChangeTime'
      endif
      do i=1,recNum
         if(windGustChangeTime(i) .eq. dfilval) 
     &                               windGustChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windGustChangeTime'
      endif
      do i=1,recNum
         if(windGustChangeTime(i) .eq. dmisval) 
     &                               windGustChangeTime(i) = badflag
      enddo !i
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
       nf_status = NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dfilval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeedChangeTime'
      endif
      do i=1,recNum
         if(windSpeedChangeTime(i) .eq. dfilval) 
     &                                windSpeedChangeTime(i) = badflag
      enddo !i
      nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value',dmisval)
      if(nf_status .ne. NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *, ' in var windSpeedChangeTime'
      endif
      do i=1,recNum
         if(windSpeedChangeTime(i) .eq. dmisval) 
     &                                windSpeedChangeTime(i) = badflag
      enddo !i
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif
c
      istatus = 1
c
      return
      end
