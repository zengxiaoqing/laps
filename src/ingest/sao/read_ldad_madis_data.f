
C
C  Subroutine to read the file "LDAD automated mesonet data " 
C
      subroutine read_ldadmadis_netcdf(nf_fid, maxSensor, recNum, 
     +     filterSetNum, firstOverflow, globalInventory, nStaticIds, 
     +     numPST, numericWMOid, precipIntensity, precipType, 
     +     pressChangeChar, altimeter, dewpoint, elevation, latitude, 
     +     longitude, meanWeightedTemperature, precipAccum, 
     +     precipRate, pressChange3Hour, rawPrecip, relHumidity, 
     +     seaLevelPressure, soilMoisture, soilTemperature, 
     +     solarRadiation, stationPressure, temperature, visibility, 
     +     windDir, windDirMax, windGust, windSpeed, observationTime, 
     +     receivedTime, reportTime, rhChangeTime, 
     +     stationPressChangeTime, tempChangeTime, windDirChangeTime, 
     +     windGustChangeTime, windSpeedChangeTime, altimeterDD, 
     +     dataProvider, dewpointDD, precipAccumDD, precipRateDD, 
     +     pressChange3HourDD, providerId, relHumidityDD, 
     +     seaLevelPressureDD, staticIds, stationId, stationName, 
     +     stationPressureDD, stationType, temperatureDD, test1, 
     +     visibilityDD, windDirDD, windSpeedDD)
C
      include 'netcdf.inc'
      integer maxSensor, recNum,nf_fid, nf_vid, nf_status
      integer filterSetNum, firstOverflow, globalInventory,
     +     nStaticIds, numPST, numericWMOid(recNum), precipIntensity(
     +     maxSensor, recNum), precipType( maxSensor, recNum),
     +     pressChangeChar(recNum)
      real altimeter(recNum), dewpoint(recNum), elevation(recNum),
     +     latitude(recNum), longitude(recNum),
     +     meanWeightedTemperature(recNum), precipAccum(recNum),
     +     precipRate(recNum), pressChange3Hour(recNum),
     +     rawPrecip(recNum), relHumidity(recNum),
     +     seaLevelPressure(recNum), soilMoisture(recNum),
     +     soilTemperature(recNum), solarRadiation(recNum),
     +     stationPressure(recNum), temperature(recNum),
     +     visibility(recNum), windDir(recNum), windDirMax(recNum),
     +     windGust(recNum), windSpeed(recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum), rhChangeTime(recNum),
     +     stationPressChangeTime(recNum), tempChangeTime(recNum),
     +     windDirChangeTime(recNum), windGustChangeTime(recNum),
     +     windSpeedChangeTime(recNum)
      character windDirDD(recNum)
      character*11 stationType(recNum)
      character*51 test1(recNum)
      character*24 staticIds
      character windSpeedDD(recNum)
      character relHumidityDD(recNum)
      character stationPressureDD(recNum)
      character altimeterDD(recNum)
      character pressChange3HourDD(recNum)
      character precipRateDD(recNum)
      character*11 dataProvider(recNum)
      character dewpointDD(recNum)
      character*6 stationId(recNum)
      character seaLevelPressureDD(recNum)
      character visibilityDD(recNum)
      character precipAccumDD(recNum)
      character*51 stationName(recNum)
      character*12 providerId(recNum)
      character temperatureDD(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      altimeter    "altimeter setting"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeter'
      endif
      call ck_array_real(altimeter,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,altimeter)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeter'
      endif
      call ck_array_real(altimeter,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      dewpoint     "dew point temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
      endif
      call ck_array_real(dewpoint,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dewpoint)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
      endif
      call ck_array_real(dewpoint,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      elevation    "elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
      call ck_array_real(elevation,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
      call ck_array_real(elevation,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      latitude     "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
      call ck_array_real(latitude,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
      call ck_array_real(latitude,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      longitude    "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
      call ck_array_real(longitude,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
      call ck_array_real(longitude,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      meanWeightedTemperature"Mean weighted temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'meanWeightedTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var meanWeightedTemperature'
      endif
      call ck_array_real(meanWeightedTemperature,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,meanWeightedTemperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var meanWeightedTemperature'
      endif
      call ck_array_real(meanWeightedTemperature,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      precipAccum  "precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccum'
      endif
      call ck_array_real(precipAccum,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccum'
      endif
      call ck_array_real(precipAccum,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      precipRate   "precipitation rate"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRate'
      endif
      call ck_array_real(precipRate,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipRate)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRate'
      endif
      call ck_array_real(precipRate,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      pressChange3Hour"3 hour pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3Hour'
      endif
      call ck_array_real(pressChange3Hour,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressChange3Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3Hour'
      endif
      call ck_array_real(pressChange3Hour,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      rawPrecip    "raw precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawPrecip',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawPrecip'
      endif
      call ck_array_real(rawPrecip,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rawPrecip)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawPrecip'
      endif
      call ck_array_real(rawPrecip,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      relHumidity  "relative humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
      call ck_array_real(relHumidity,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
      call ck_array_real(relHumidity,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      seaLevelPressure"sea level pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
      endif
      call ck_array_real(seaLevelPressure,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
      endif
      call ck_array_real(seaLevelPressure,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      soilMoisture "Soil moisture"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilMoisture',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilMoisture'
      endif
      call ck_array_real(soilMoisture,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilMoisture)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilMoisture'
      endif
      call ck_array_real(soilMoisture,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      soilTemperature"Soil temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilTemperature'
      endif
      call ck_array_real(soilTemperature,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilTemperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilTemperature'
      endif
      call ck_array_real(soilTemperature,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      solarRadiation"solar radiation"
C
      nf_status=NF_INQ_VARID(nf_fid,'solarRadiation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadiation'
      endif
      call ck_array_real(solarRadiation,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,solarRadiation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadiation'
      endif
      call ck_array_real(solarRadiation,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      stationPressure"station pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
      endif
      call ck_array_real(stationPressure,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
      endif
      call ck_array_real(stationPressure,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      temperature  "temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
      call ck_array_real(temperature,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
      call ck_array_real(temperature,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      visibility   "visibility"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
      endif
      call ck_array_real(visibility,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,visibility)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
      endif
      call ck_array_real(visibility,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      windDir      "wind direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
      call ck_array_real(windDir,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
      call ck_array_real(windDir,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      windDirMax   "wind direction at gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirMax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
      endif
      call ck_array_real(windDirMax,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDirMax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
      endif
      call ck_array_real(windDirMax,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      windGust     "wind gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGust',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
      endif
      call ck_array_real(windGust,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windGust)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
      endif
      call ck_array_real(windGust,recNum,misval,badflag)
C
C     Variable        NETCDF Long Name
C      windSpeed    "wind speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
      call ck_array_real(windSpeed,recNum,filval,badflag)
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
      call ck_array_real(windSpeed,recNum,misval,badflag)

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      filterSetNum 
C
      nf_status=NF_INQ_VARID(nf_fid,'filterSetNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var filterSetNum'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,filterSetNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var filterSetNum'
      endif
C
C     Variable        NETCDF Long Name
C      firstOverflow
C
      nf_status=NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      globalInventory
C
      nf_status=NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
C
C     Variable        NETCDF Long Name
C      nStaticIds   
C
      nf_status=NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
C
C     Variable        NETCDF Long Name
C      numPST       "Number of entries in Provider-Subprovider Table"
C
      nf_status=NF_INQ_VARID(nf_fid,'numPST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numPST'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numPST)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numPST'
      endif
C
C     Variable        NETCDF Long Name
C      numericWMOid "numeric WMO identification"
C
      nf_status=NF_INQ_VARID(nf_fid,'numericWMOid',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numericWMOid'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numericWMOid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numericWMOid'
      endif
C
C     Variable        NETCDF Long Name
C      precipIntensity"precipitation intensity"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipIntensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipIntensity'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipIntensity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipIntensity'
      endif
C
C     Variable        NETCDF Long Name
C      precipType   "precipitation type"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipType'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipType'
      endif
C
C     Variable        NETCDF Long Name
C      pressChangeChar"character of pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChangeChar',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChangeChar'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChangeChar)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChangeChar'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      observationTime"time of observation"
C
      nf_status=NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
C
C     Variable        NETCDF Long Name
C      receivedTime "time data was received"
C
      nf_status=NF_INQ_VARID(nf_fid,'receivedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receivedTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receivedTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receivedTime'
      endif
C
C     Variable        NETCDF Long Name
C      reportTime   "time data was processed by the provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif
C
C     Variable        NETCDF Long Name
C      rhChangeTime "relative humidity time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'rhChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,rhChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhChangeTime'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressChangeTime"station press time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,stationPressChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressChangeTime'
      endif
C
C     Variable        NETCDF Long Name
C      tempChangeTime"temperature time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'tempChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,tempChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempChangeTime'
      endif
C
C     Variable        NETCDF Long Name
C      windDirChangeTime"wind direction time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windDirChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirChangeTime'
      endif
C
C     Variable        NETCDF Long Name
C      windGustChangeTime"wind gust time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGustChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGustChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windGustChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGustChangeTime'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedChangeTime"wind speed time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windSpeedChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedChangeTime'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      altimeterDD  "altimeter setting QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,altimeterDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterDD'
      endif
C
C     Variable        NETCDF Long Name
C      dataProvider "Local data provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
C
C     Variable        NETCDF Long Name
C      dewpointDD   "dew point QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dewpointDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointDD'
      endif
C
C     Variable        NETCDF Long Name
C      precipAccumDD"precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipAccumDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumDD'
      endif
C
C     Variable        NETCDF Long Name
C      precipRateDD "precip rate QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipRateDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateDD'
      endif
C
C     Variable        NETCDF Long Name
C      pressChange3HourDD"3h pressure change QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,pressChange3HourDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourDD'
      endif
C
C     Variable        NETCDF Long Name
C      providerId   "Data Provider station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityDD"relative humidity QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,relHumidityDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityDD'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressureDD"Sea level pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,seaLevelPressureDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureDD'
      endif
C
C     Variable        NETCDF Long Name
C      staticIds    
C
      nf_status=NF_INQ_VARID(nf_fid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,staticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
C
C     Variable        NETCDF Long Name
C      stationId    "alphanumeric station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      endif
C
C     Variable        NETCDF Long Name
C      stationName  "alphanumeric station name"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressureDD"station pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationPressureDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureDD'
      endif
C
C     Variable        NETCDF Long Name
C      stationType  "LDAD station type"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureDD"temperature QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,temperatureDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureDD'
      endif
C
C     Variable        NETCDF Long Name
C      test1        "User defined parameter - test # 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'test1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var test1'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,test1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var test1'
      endif
C
C     Variable        NETCDF Long Name
C      visibilityDD "visibility QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,visibilityDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityDD'
      endif
C
C     Variable        NETCDF Long Name
C      windDirDD    "wind direction QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,windDirDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirDD'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedDD  "wind speed QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,windSpeedDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedDD'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
