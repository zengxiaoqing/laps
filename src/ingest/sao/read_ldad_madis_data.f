C
C  Subroutine to read the file "LDAD automated mesonet data " 
C
      subroutine read_ldad_madis_netcdf(nf_fid, maxSensor, recNum, 
     +     firstOverflow, globalInventory, nStaticIds, numPST, 
     +     numericWMOid, precipIntensity, precipType, 
     +     pressChangeChar, altimeter, dewpoint, elevation, latitude, 
     +     longitude, meanWeightedTemperature, precipAccum, 
     +     precipRate, pressChange3Hour, rawPrecip, relHumidity, 
     +     seaLevelPressure, seaSurfaceTemp, soilMoisture, 
     +     soilTemperature, solarRadiation, stationPressure, 
     +     temperature, visibility, windDir, windDirMax, windGust, 
     +     windSpeed, altimeterDD, dataProvider, dewpointDD, 
     +     precipAccumDD, precipRateDD, presWeather, 
     +     pressChange3HourDD, providerId, relHumidityDD, 
     +     seaLevelPressureDD, stationId, stationName, 
     +     stationPressureDD, stationType, temperatureDD, 
     +     visibilityDD, windDirDD, windSpeedDD, observationTime, 
     +     receivedTime, reportTime, rhChangeTime, 
     +     stationPressChangeTime, tempChangeTime, windDirChangeTime, 
     +     windGustChangeTime, windSpeedChangeTime,badflag)
C
      include 'netcdf.inc'
      integer maxSensor, recNum,nf_fid, nf_vid, nf_status
      integer firstOverflow, globalInventory, nStaticIds, numPST,
     +     numericWMOid(recNum), precipIntensity( maxSensor, recNum),
     +     precipType( maxSensor, recNum), pressChangeChar(recNum)
      real altimeter(recNum), dewpoint(recNum), elevation(recNum),
     +     latitude(recNum), longitude(recNum),
     +     meanWeightedTemperature(recNum), precipAccum(recNum),
     +     precipRate(recNum), pressChange3Hour(recNum),
     +     rawPrecip(recNum), relHumidity(recNum),
     +     seaLevelPressure(recNum), seaSurfaceTemp(recNum),
     +     soilMoisture(recNum), soilTemperature(recNum),
     +     solarRadiation(recNum), stationPressure(recNum),
     +     temperature(recNum), visibility(recNum), windDir(recNum),
     +     windDirMax(recNum), windGust(recNum), windSpeed(recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum), rhChangeTime(recNum),
     +     stationPressChangeTime(recNum), tempChangeTime(recNum),
     +     windDirChangeTime(recNum), windGustChangeTime(recNum),
     +     windSpeedChangeTime(recNum)
      character*12 providerId(recNum)
      character windSpeedDD(recNum)
      character seaLevelPressureDD(recNum)
      character*6 stationId(recNum)
      character temperatureDD(recNum)
      character visibilityDD(recNum)
      character precipAccumDD(recNum)
      character pressChange3HourDD(recNum)
      character*11 dataProvider(recNum)
      character altimeterDD(recNum)
      character stationPressureDD(recNum)
      character dewpointDD(recNum)
      character*25 presWeather(recNum)
      character precipRateDD(recNum)
      character*11 stationType(recNum)
      character relHumidityDD(recNum)
      character windDirDD(recNum)
      character*51 stationName(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     altimeter     "altimeter setting"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var altimeter'
      else
       call ck_array_real(altimeter,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,altimeter)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeter'
       endif
       call ck_array_real(altimeter,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     dewpoint      "dew point temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var dewpoint'
      else
       call ck_array_real(dewpoint,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dewpoint)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
       endif
       call ck_array_real(dewpoint,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     elevation     "elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var elevation'
      else
       call ck_array_real(elevation,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
       endif
       call ck_array_real(elevation,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     latitude      "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var latitude'
      else
       call ck_array_real(latitude,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
       endif
       call ck_array_real(latitude,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     longitude     "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var longitude'
      else
       call ck_array_real(longitude,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
       endif
       call ck_array_real(longitude,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     meanWeightedTemperature"Mean weighted temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'meanWeightedTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var meanWeightedTemperature'
      else
       call ck_array_real(meanWeightedTemperature,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,meanWeightedTemperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var meanWeightedTemperature'
       endif
       call ck_array_real(meanWeightedTemperature,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     precipAccum   "precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var precipAccum'
      else
       call ck_array_real(precipAccum,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccum)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccum'
       endif
       call ck_array_real(precipAccum,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     precipRate    "precipitation rate"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var precipRate'
      else
       call ck_array_real(precipRate,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipRate)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRate'
       endif
       call ck_array_real(precipRate,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     pressChange3Hour"3 hour pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var pressChange3Hour'
      else
       call ck_array_real(pressChange3Hour,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressChange3Hour)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3Hour'
       endif
       call ck_array_real(pressChange3Hour,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     rawPrecip     "raw precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawPrecip',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var rawPrecip'
      else
       call ck_array_real(rawPrecip,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rawPrecip)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawPrecip'
       endif
       call ck_array_real(rawPrecip,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     relHumidity   "relative humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var relHumidity'
      else
       call ck_array_real(relHumidity,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
       endif
       call ck_array_real(relHumidity,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     seaLevelPressure"sea level pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var seaLevelPressure'
      else
       call ck_array_real(seaLevelPressure,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
       endif
       call ck_array_real(seaLevelPressure,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     seaSurfaceTemp"sea surface temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaSurfaceTemp',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var seaSurfaceTemp'
      else
       call ck_array_real(seaSurfaceTemp,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaSurfaceTemp)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaSurfaceTemp'
       endif
       call ck_array_real(seaSurfaceTemp,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     soilMoisture  "Soil moisture"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilMoisture',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var soilMoisture'
      else
       call ck_array_real(soilMoisture,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilMoisture)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilMoisture'
       endif
       call ck_array_real(soilMoisture,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     soilTemperature"Soil temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var soilTemperature'
      else
       call ck_array_real(soilTemperature,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilTemperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilTemperature'
       endif
       call ck_array_real(soilTemperature,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     solarRadiation"solar radiation"
C
      nf_status=NF_INQ_VARID(nf_fid,'solarRadiation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var solarRadiation'
      else
       call ck_array_real(solarRadiation,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,solarRadiation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadiation'
       endif
       call ck_array_real(solarRadiation,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     stationPressure"station pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var stationPressure'
      else
       call ck_array_real(stationPressure,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
       endif
       call ck_array_real(stationPressure,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     temperature   "temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var temperature'
      else
       call ck_array_real(temperature,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
       endif
       call ck_array_real(temperature,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     visibility    "visibility"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var visibility'
      else
       call ck_array_real(visibility,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,visibility)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
       endif
       call ck_array_real(visibility,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windDir       "wind direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windDir'
      else
       call ck_array_real(windDir,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
       endif
       call ck_array_real(windDir,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windDirMax    "wind direction at gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirMax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windDirMax'
      else
       call ck_array_real(windDirMax,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDirMax)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
       endif
       call ck_array_real(windDirMax,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windGust      "wind gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGust',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windGust'
      else
       call ck_array_real(windGust,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windGust)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
       endif
       call ck_array_real(windGust,recNum,misval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windSpeed     "wind speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windSpeed'
      else
       call ck_array_real(windSpeed,recNum,filval,badflag)
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
       endif
       call ck_array_real(windSpeed,recNum,misval,badflag)
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     firstOverflow 
C
      nf_status=NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var firstOverflow'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     globalInventory
C
      nf_status=NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var globalInventory'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     nStaticIds    
C
      nf_status=NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var nStaticIds'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numPST        "Number of entries in Provider-Subprovider Table"
C
      nf_status=NF_INQ_VARID(nf_fid,'numPST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var numPST'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numPST)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numPST'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numericWMOid  "numeric WMO identification"
C
      nf_status=NF_INQ_VARID(nf_fid,'numericWMOid',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var numericWMOid'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numericWMOid)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var numericWMOid'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipIntensity"precipitation intensity"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipIntensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var precipIntensity'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipIntensity)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipIntensity'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipType    "precipitation type"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var precipType'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipType)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipType'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressChangeChar"character of pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChangeChar',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var pressChangeChar'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChangeChar)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChangeChar'
       endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     observationTime"time of observation"
C
      nf_status=NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var observationTime'
      else
       call ck_array_dble(observationTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
       endif
       call ck_array_dble(observationTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     receivedTime  "time data was received"
C
      nf_status=NF_INQ_VARID(nf_fid,'receivedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var receivedTime'
      else
       call ck_array_dble(receivedTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receivedTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receivedTime'
       endif
       call ck_array_dble(receivedTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     reportTime    "time data was processed by the provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var reportTime'
      else
       call ck_array_dble(reportTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
       endif
       call ck_array_dble(reportTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     rhChangeTime  "relative humidity time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'rhChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var rhChangeTime'
      else
       call ck_array_dble(rhChangeTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,rhChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhChangeTime'
       endif
       call ck_array_dble(rhChangeTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     stationPressChangeTime"station press time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var stationPressChangeTime'
      else
       call ck_array_dble(stationPressChangeTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,stationPressChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressChangeTime'
       endif
       call ck_array_dble(stationPressChangeTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     tempChangeTime"temperature time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'tempChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var tempChangeTime'
      else
       call ck_array_dble(tempChangeTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,tempChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempChangeTime'
       endif
       call ck_array_dble(tempChangeTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windDirChangeTime"wind direction time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windDirChangeTime'
      else
       call ck_array_dble(windDirChangeTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windDirChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirChangeTime'
       endif
       call ck_array_dble(windDirChangeTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windGustChangeTime"wind gust time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGustChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windGustChangeTime'
      else
       call ck_array_dble(windGustChangeTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windGustChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGustChangeTime'
       endif
       call ck_array_dble(windGustChangeTime,recNum,dmisval,badflag)
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedChangeTime"wind speed time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windSpeedChangeTime'
      else
       call ck_array_dble(windSpeedChangeTime,recNum,dfilval,badflag)
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windSpeedChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedChangeTime'
       endif
       call ck_array_dble(windSpeedChangeTime,recNum,dmisval,badflag)
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     altimeterDD   "altimeter setting QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var altimeterDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,altimeterDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dataProvider  "Local data provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var dataProvider'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dewpointDD    "dew point QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var dewpointDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dewpointDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumDD "precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var precipAccumDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipAccumDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipRateDD  "precip rate QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var precipRateDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipRateDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     presWeather   "present weather"
C
      nf_status=NF_INQ_VARID(nf_fid,'presWeather',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var presWeather'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,presWeather)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var presWeather'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressChange3HourDD"3h pressure change QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var pressChange3HourDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,pressChange3HourDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     providerId    "Data Provider station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var providerId'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityDD "relative humidity QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var relHumidityDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,relHumidityDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     seaLevelPressureDD"Sea level pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var seaLevelPressureDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,seaLevelPressureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationId     "alphanumeric station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var stationId'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationName   "alphanumeric station name"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var stationName'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationPressureDD"station pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var stationPressureDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationPressureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationType   "LDAD station type"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var stationType'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperatureDD "temperature QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var temperatureDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,temperatureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     visibilityDD  "visibility QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var visibilityDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,visibilityDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirDD     "wind direction QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windDirDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,windDirDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedDD   "wind speed QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status)
       print *,'in var windSpeedDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,windSpeedDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedDD'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
