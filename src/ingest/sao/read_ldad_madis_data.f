C
C  Subroutine to read the file "LDAD automated mesonet data " 
C
      subroutine read_ldad_madis_netcdf(nf_fid, maxSensor, recNum, 
     +     maxPSTEntries, code1PST, code2PST, namePST,
     +     altimeterQCR, dewpointQCR, firstOverflow, globalInventory, 
     +     nStaticIds, numericWMOid, precipAccumQCR, precipIntensity, 
     +     precipRateQCR, precipType, pressChange3HourQCR, 
     +     pressChangeChar, relHumidityQCR, seaLevelPressureQCR, 
     +     stationPressureQCR, temperatureQCR, visibilityQCR, 
     +     windDirQCR, windSpeedQCR, altimeter, dewpoint, elevation, 
     +     latitude, longitude, meanWeightedTemperature, precipAccum, 
     +     precipRate, pressChange3Hour, relHumidity, 
     +     seaLevelPressure, seaSurfaceTemp, soilMoisturePercent, 
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
      integer code1PST(maxPSTEntries), code2PST(maxPSTEntries), 
     +     code4PST(maxPSTEntries),
     +     altimeterQCR(recNum), dewpointQCR(recNum),
     +     firstOverflow, globalInventory, nStaticIds,
     +     numericWMOid(recNum), precipAccumQCR(recNum),
     +     precipIntensity( maxSensor, recNum),
     +     precipRateQCR(recNum), precipType( maxSensor, recNum),
     +     pressChange3HourQCR(recNum), pressChangeChar(recNum),
     +     relHumidityQCR(recNum), seaLevelPressureQCR(recNum),
     +     stationPressureQCR(recNum), temperatureQCR(recNum),
     +     visibilityQCR(recNum), windDirQCR(recNum),
     +     windSpeedQCR(recNum)
      real altimeter(recNum), dewpoint(recNum), elevation(recNum),
     +     latitude(recNum), longitude(recNum),
     +     meanWeightedTemperature(recNum), precipAccum(recNum),
     +     precipRate(recNum), pressChange3Hour(recNum),
     +     relHumidity(recNum), seaLevelPressure(recNum),
     +     seaSurfaceTemp(recNum), soilMoisturePercent(recNum),
     +     soilTemperature(recNum), solarRadiation(recNum),
     +     stationPressure(recNum), temperature(recNum),
     +     visibility(recNum), windDir(recNum), windDirMax(recNum),
     +     windGust(recNum), windSpeed(recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum), rhChangeTime(recNum),
     +     stationPressChangeTime(recNum), tempChangeTime(recNum),
     +     windDirChangeTime(recNum), windGustChangeTime(recNum),
     +     windSpeedChangeTime(recNum)
      character temperatureDD(recNum)
      character precipAccumDD(recNum)
      character dewpointDD(recNum)
      character*6 stationId(recNum)
      character pressChange3HourDD(recNum)
      character relHumidityDD(recNum)
      character*25 presWeather(recNum)
      character*12 providerId(recNum)
      character visibilityDD(recNum)
      character windDirDD(recNum)
      character*11 dataProvider(recNum)
      character*11 namePST(maxPSTEntries)
      character altimeterDD(recNum)
      character*51 stationName(recNum)
      character precipRateDD(recNum)
      character*11 stationType(recNum)
      character stationPressureDD(recNum)
      character seaLevelPressureDD(recNum)
      character windSpeedDD(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     altimeter     "altimeter setting"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeter',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for altimeter'
       print *,'Set altimeter to badflag'
       altimeter = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,altimeter)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for altimeter'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - altimeter'
       else
        call ck_array_real(altimeter,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - altimeter'
       else
        call ck_array_real(altimeter,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dewpoint      "dew point temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dewpoint'
       print *,'Set dewpoint to badflag'
       dewpoint = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dewpoint)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dewpoint'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - dewpoint'
       else
        call ck_array_real(dewpoint,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - dewpoint'
       else
        call ck_array_real(dewpoint,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     elevation     "elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for elevation'
       print *,'Set elevation to badflag'
       elevation = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for elevation'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - elevation'
       else
        call ck_array_real(elevation,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - elevation'
       else
        call ck_array_real(elevation,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitude      "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitude'
       print *,'Set latitude to badflag'
       latitude = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitude'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - latitude'
       else
        call ck_array_real(latitude,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - latitude'
       else
        call ck_array_real(latitude,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitude     "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitude'
       print *,'Set longitude to badflag'
       longitude = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitude'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - longitude'
       else
        call ck_array_real(longitude,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - longitude'
       else
        call ck_array_real(longitude,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     meanWeightedTemperature"Mean weighted temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'meanWeightedTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for meanWeightedTemperature'
       print *,'Set meanWeightedTemperature to badflag'
       meanWeightedTemperature = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,meanWeightedTemperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for meanWeightedTemperature'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - meanWeightedTemperature'
       else
        call ck_array_real(meanWeightedTemperature,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - meanWeightedTemperature'
       else
        call ck_array_real(meanWeightedTemperature,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccum   "precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccum'
       print *,'Set precipAccum to badflag'
       precipAccum = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccum)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccum'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precipAccum'
       else
        call ck_array_real(precipAccum,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precipAccum'
       else
        call ck_array_real(precipAccum,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipRate    "precipitation rate"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipRate'
       print *,'Set precipRate to badflag'
       precipRate = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipRate)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipRate'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - precipRate'
       else
        call ck_array_real(precipRate,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - precipRate'
       else
        call ck_array_real(precipRate,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressChange3Hour"3 hour pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressChange3Hour'
       print *,'Set pressChange3Hour to badflag'
       pressChange3Hour = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressChange3Hour)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressChange3Hour'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - pressChange3Hour'
       else
        call ck_array_real(pressChange3Hour,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - pressChange3Hour'
       else
        call ck_array_real(pressChange3Hour,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidity   "relative humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidity'
       print *,'Set relHumidity to badflag'
       relHumidity = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidity'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - relHumidity'
       else
        call ck_array_real(relHumidity,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - relHumidity'
       else
        call ck_array_real(relHumidity,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     seaLevelPressure"sea level pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for seaLevelPressure'
       print *,'Set seaLevelPressure to badflag'
       seaLevelPressure = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for seaLevelPressure'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - seaLevelPressure'
       else
        call ck_array_real(seaLevelPressure,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - seaLevelPressure'
       else
        call ck_array_real(seaLevelPressure,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     seaSurfaceTemp"sea surface temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaSurfaceTemp',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for seaSurfaceTemp'
       print *,'Set seaSurfaceTemp to badflag'
       seaSurfaceTemp = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaSurfaceTemp)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for seaSurfaceTemp'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - seaSurfaceTemp'
       else
        call ck_array_real(seaSurfaceTemp,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - seaSurfaceTemp'
       else
        call ck_array_real(seaSurfaceTemp,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     soilMoisturePercent  "Soil moisture"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilMoisturePercent',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for soilMoisturePercent'
       print *,'Set soilMoisturePercent to badflag'
       soilMoisturePercent = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilMoisturePercent)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for soilMoisturePercent'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - soilMoisturePercent'
       else
        call ck_array_real(soilMoisturePercent,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - soilMoisturePercent'
       else
        call ck_array_real(soilMoisturePercent,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     soilTemperature"Soil temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for soilTemperature'
       print *,'Set soilTemperature to badflag'
       soilTemperature = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilTemperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for soilTemperature'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - soilTemperature'
       else
        call ck_array_real(soilTemperature,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - soilTemperature'
       else
        call ck_array_real(soilTemperature,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     solarRadiation"solar radiation"
C
      nf_status=NF_INQ_VARID(nf_fid,'solarRadiation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for solarRadiation'
       print *,'Set solarRadiation to badflag'
       solarRadiation = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,solarRadiation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for solarRadiation'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - solarRadiation'
       else
        call ck_array_real(solarRadiation,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - solarRadiation'
       else
        call ck_array_real(solarRadiation,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationPressure"station pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationPressure'
       print *,'Set stationPressure to badflag'
       stationPressure = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationPressure'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - stationPressure'
       else
        call ck_array_real(stationPressure,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - stationPressure'
       else
        call ck_array_real(stationPressure,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperature   "temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperature'
       print *,'Set temperature to badflag'
       temperature = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperature'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - temperature'
       else
        call ck_array_real(temperature,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - temperature'
       else
        call ck_array_real(temperature,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     visibility    "visibility"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for visibility'
       print *,'Set visibility to badflag'
       visibility = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,visibility)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for visibility'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - visibility'
       else
        call ck_array_real(visibility,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - visibility'
       else
        call ck_array_real(visibility,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDir       "wind direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDir'
       print *,'Set windDir to badflag'
       windDir = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDir'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windDir'
       else
        call ck_array_real(windDir,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windDir'
       else
        call ck_array_real(windDir,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirMax    "wind direction at gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirMax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirMax'
       print *,'Set windDirMax to badflag'
       windDirMax = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDirMax)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirMax'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windDirMax'
       else
        call ck_array_real(windDirMax,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windDirMax'
       else
        call ck_array_real(windDirMax,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windGust      "wind gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGust',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windGust'
       print *,'Set windGust to badflag'
       windGust = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windGust)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windGust'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windGust'
       else
        call ck_array_real(windGust,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windGust'
       else
        call ck_array_real(windGust,recNum,valmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeed     "wind speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeed'
       print *,'Set windSpeed to badflag'
       windSpeed = badflag
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeed'
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',valfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windSpeed'
       else
        call ck_array_real(windSpeed,recNum,valfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_REAL(nf_fid,nf_vid,'missing_value',valmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windSpeed'
       else
        call ck_array_real(windSpeed,recNum,valmis
     1                    ,badflag)
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     altimeterQCR  "altimeter setting QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for altimeterQCR'
       print *,'Set altimeterQCR to -99'
       altimeterQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,altimeterQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for altimeterQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dewpointQCR   "dew point QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dewpointQCR'
       print *,'Set dewpointQCR to -99'
       dewpointQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,dewpointQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dewpointQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     firstOverflow 
C
      nf_status=NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for firstOverflow'
       print *,'Set firstOverflow to -99'
       firstOverflow = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for firstOverflow'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     globalInventory
C
      nf_status=NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for globalInventory'
       print *,'Set globalInventory to -99'
       globalInventory = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for globalInventory'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     nStaticIds    
C
      nf_status=NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for nStaticIds'
       print *,'Set nStaticIds to -99'
       nStaticIds = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for nStaticIds'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numericWMOid  "numeric WMO identification"
C
      nf_status=NF_INQ_VARID(nf_fid,'numericWMOid',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numericWMOid'
       print *,'Set numericWMOid to -99'
       numericWMOid = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numericWMOid)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numericWMOid'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumQCR"precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumQCR'
       print *,'Set precipAccumQCR to -99'
       precipAccumQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipIntensity"precipitation intensity"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipIntensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipIntensity'
       print *,'Set precipIntensity to -99'
       precipIntensity = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipIntensity)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipIntensity'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipRateQCR "precip rate QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipRateQCR'
       print *,'Set precipRateQCR to -99'
       precipRateQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipRateQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipRateQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipType    "precipitation type"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipType'
       print *,'Set precipType to -99'
       precipType = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipType)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipType'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressChange3HourQCR"3h pressure change QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressChange3HourQCR'
       print *,'Set pressChange3HourQCR to -99'
       pressChange3HourQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChange3HourQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressChange3HourQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressChangeChar"character of pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChangeChar',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressChangeChar'
       print *,'Set pressChangeChar to -99'
       pressChangeChar = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChangeChar)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressChangeChar'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityQCR"relative humidity QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidityQCR'
       print *,'Set relHumidityQCR to -99'
       relHumidityQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidityQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     seaLevelPressureQCR"sea level pressure QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for seaLevelPressureQCR'
       print *,'Set seaLevelPressureQCR to -99'
       seaLevelPressureQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,seaLevelPressureQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for seaLevelPressureQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationPressureQCR"station pressure QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationPressureQCR'
       print *,'Set stationPressureQCR to -99'
       stationPressureQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,stationPressureQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationPressureQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperatureQCR"temperature QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperatureQCR'
       print *,'Set temperatureQCR to -99'
       temperatureQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperatureQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     visibilityQCR "visibility QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for visibilityQCR'
       print *,'Set visibilityQCR to -99'
       visibilityQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,visibilityQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for visibilityQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirQCR    "wind direction QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirQCR'
       print *,'Set windDirQCR to -99'
       windDirQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedQCR  "wind speed QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeedQCR'
       print *,'Set windSpeedQCR to -99'
       windSpeedQCR = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeedQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     code1PST      "precip variable definition"
C
      nf_status=NF_INQ_VARID(nf_fid,'code1PST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for code1PST'
       print *,'Set code1PST to -99'
       code1PST = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,code1PST)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for code1PST'
       endif
      endif

C
C     Variable        NETCDF Long Name
C     code2PST      "solarRadiation variable definition"
C
      nf_status=NF_INQ_VARID(nf_fid,'code2PST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for code2PST'
       print *,'Set code2PST to -99'
       code2PST = -99
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,code2PST)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for code2PST'
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
       print *, NF_STRERROR(nf_status),' for observationTime'
       print *,'Set observationTime to badflag'
       observationTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for observationTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - observationTime'
       else
        call ck_array_dble(observationTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - observationTime'
       else
        call ck_array_dble(observationTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     receivedTime  "time data was received"
C
      nf_status=NF_INQ_VARID(nf_fid,'receivedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for receivedTime'
       print *,'Set receivedTime to badflag'
       receivedTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receivedTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for receivedTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - receivedTime'
       else
        call ck_array_dble(receivedTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - receivedTime'
       else
        call ck_array_dble(receivedTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     reportTime    "time data was processed by the provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for reportTime'
       print *,'Set reportTime to badflag'
       reportTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for reportTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - reportTime'
       else
        call ck_array_dble(reportTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - reportTime'
       else
        call ck_array_dble(reportTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rhChangeTime  "relative humidity time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'rhChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rhChangeTime'
       print *,'Set rhChangeTime to badflag'
       rhChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,rhChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rhChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - rhChangeTime'
       else
        call ck_array_dble(rhChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - rhChangeTime'
       else
        call ck_array_dble(rhChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationPressChangeTime"station press time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationPressChangeTime'
       print *,'Set stationPressChangeTime to badflag'
       stationPressChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,stationPressChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationPressChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - stationPressChangeTime'
       else
        call ck_array_dble(stationPressChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - stationPressChangeTime'
       else
        call ck_array_dble(stationPressChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tempChangeTime"temperature time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'tempChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tempChangeTime'
       print *,'Set tempChangeTime to badflag'
       tempChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,tempChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tempChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - tempChangeTime'
       else
        call ck_array_dble(tempChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - tempChangeTime'
       else
        call ck_array_dble(tempChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirChangeTime"wind direction time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirChangeTime'
       print *,'Set windDirChangeTime to badflag'
       windDirChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windDirChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windDirChangeTime'
       else
        call ck_array_dble(windDirChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windDirChangeTime'
       else
        call ck_array_dble(windDirChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windGustChangeTime"wind gust time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGustChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windGustChangeTime'
       print *,'Set windGustChangeTime to badflag'
       windGustChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windGustChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windGustChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windGustChangeTime'
       else
        call ck_array_dble(windGustChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windGustChangeTime'
       else
        call ck_array_dble(windGustChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedChangeTime"wind speed time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeedChangeTime'
       print *,'Set windSpeedChangeTime to badflag'
       windSpeedChangeTime = badflag
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,windSpeedChangeTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeedChangeTime'
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'_FillValue',dvalfil)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for _FillValue - windSpeedChangeTime'
       else
        call ck_array_dble(windSpeedChangeTime,recNum,dvalfil
     1                    ,badflag)
       endif
       nf_status=NF_GET_ATT_DOUBLE(nf_fid,nf_vid,'missing_value'
     1                            ,dvalmis)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
     1         ,' for missing_value - windSpeedChangeTime'
       else
        call ck_array_dble(windSpeedChangeTime,recNum,dvalmis
     1                    ,badflag)
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     altimeterDD   "altimeter setting QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for altimeterDD'
       print *,'Set altimeterDD to " "'
       altimeterDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,altimeterDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for altimeterDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dataProvider  "Local data provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dataProvider'
       print *,'Set dataProvider to " "'
       dataProvider = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dataProvider'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dewpointDD    "dew point QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dewpointDD'
       print *,'Set dewpointDD to " "'
       dewpointDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dewpointDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dewpointDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     namePST       "PST Provider or Subprovider name"
C
      nf_status=NF_INQ_VARID(nf_fid,'namePST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for namePST'
       print *,'Set namePST to " "'
       namePST = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,namePST)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for namePST'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipAccumDD "precip amount QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipAccumDD'
       print *,'Set precipAccumDD to " "'
       precipAccumDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipAccumDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipAccumDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipRateDD  "precip rate QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipRateDD'
       print *,'Set precipRateDD to " "'
       precipRateDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,precipRateDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipRateDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     presWeather   "present weather"
C
      nf_status=NF_INQ_VARID(nf_fid,'presWeather',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for presWeather'
       print *,'Set presWeather to " "'
       presWeather = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,presWeather)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for presWeather'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressChange3HourDD"3h pressure change QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressChange3HourDD'
       print *,'Set pressChange3HourDD to " "'
       pressChange3HourDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,pressChange3HourDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressChange3HourDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     providerId    "Data Provider station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for providerId'
       print *,'Set providerId to " "'
       providerId = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for providerId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityDD "relative humidity QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidityDD'
       print *,'Set relHumidityDD to " "'
       relHumidityDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,relHumidityDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidityDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     seaLevelPressureDD"Sea level pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for seaLevelPressureDD'
       print *,'Set seaLevelPressureDD to " "'
       seaLevelPressureDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,seaLevelPressureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for seaLevelPressureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationId     "alphanumeric station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationId'
       print *,'Set stationId to " "'
       stationId = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationName   "alphanumeric station name"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationName'
       print *,'Set stationName to " "'
       stationName = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationName'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationPressureDD"station pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationPressureDD'
       print *,'Set stationPressureDD to " "'
       stationPressureDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationPressureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationPressureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationType   "LDAD station type"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationType'
       print *,'Set stationType to " "'
       stationType = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationType'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperatureDD "temperature QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperatureDD'
       print *,'Set temperatureDD to " "'
       temperatureDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,temperatureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperatureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     visibilityDD  "visibility QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for visibilityDD'
       print *,'Set visibilityDD to " "'
       visibilityDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,visibilityDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for visibilityDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirDD     "wind direction QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirDD'
       print *,'Set windDirDD to " "'
       windDirDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,windDirDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedDD   "wind speed QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeedDD'
       print *,'Set windSpeedDD to " "'
       windSpeedDD = ' '
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,windSpeedDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeedDD'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
