C
C  Subroutine to read the file "LDAD automated mesonet data " 
C
      subroutine read_ldadmadis_netcdf(nf_fid, ICcheckNum, QCcheckNum, 
     +     maxPSTEntries, maxSensor, maxStaticIds, nInventoryBins, 
     +     recNum, altimeterQCA, altimeterQCR, code1PST, code2PST, 
     +     code3PST, dewpointICA, dewpointICR, dewpointQCA, 
     +     dewpointQCR, filterSetNum, firstInBin, firstOverflow, 
     +     globalInventory, invTime, inventory, isOverflow, 
     +     lastInBin, lastRecord, nStaticIds, numPST, numericWMOid, 
     +     precipAccumICA, precipAccumICR, precipAccumQCA, 
     +     precipAccumQCR, precipIntensity, precipRateQCA, 
     +     precipRateQCR, precipType, pressChange3HourICA, 
     +     pressChange3HourICR, pressChange3HourQCA, 
     +     pressChange3HourQCR, pressChangeChar, prevRecord, 
     +     relHumidityQCA, relHumidityQCR, roadState1, roadState2, 
     +     roadState3, roadState4, seaLevelPressureICA, 
     +     seaLevelPressureICR, seaLevelPressureQCA, 
     +     seaLevelPressureQCR, stationPressureICA, 
     +     stationPressureICR, stationPressureQCA, 
     +     stationPressureQCR, temperatureICA, temperatureICR, 
     +     temperatureQCA, temperatureQCR, visibilityICA, 
     +     visibilityICR, visibilityQCA, visibilityQCR, windDirICA, 
     +     windDirICR, windDirQCA, windDirQCR, windSpeedICA, 
     +     windSpeedICR, windSpeedQCA, windSpeedQCR, altimeter, 
     +     altimeterQCD, capPi, dewpoint, dewpointQCD, 
     +     drySignalDelay, elevation, formalError, fuelMoisture, 
     +     fuelTemperature, latitude, longitude, 
     +     meanWeightedTemperature, precipAccum, precipAccumQCD, 
     +     precipRate, precipRateQCD, pressChange3Hour, 
     +     pressChange3HourQCD, rawPrecip, relHumidity, 
     +     relHumidityQCD, roadLiquidChemFactor1, 
     +     roadLiquidChemFactor2, roadLiquidChemFactor3, 
     +     roadLiquidChemFactor4, roadLiquidChemPercent1, 
     +     roadLiquidChemPercent2, roadLiquidChemPercent3, 
     +     roadLiquidChemPercent4, roadLiquidDepth1, 
     +     roadLiquidDepth2, roadLiquidDepth3, roadLiquidDepth4, 
     +     roadLiquidFreezeTemp1, roadLiquidFreezeTemp2, 
     +     roadLiquidFreezeTemp3, roadLiquidFreezeTemp4, 
     +     roadLiquidIcePercent1, roadLiquidIcePercent2, 
     +     roadLiquidIcePercent3, roadLiquidIcePercent4, 
     +     roadSubsurfaceTemp1, roadSubsurfaceTemp2, 
     +     roadSubsurfaceTemp3, roadSubsurfaceTemp4, 
     +     roadTemperature1, roadTemperature2, roadTemperature3, 
     +     roadTemperature4, seaLevelPressure, seaLevelPressureQCD, 
     +     soilMoisture, soilTemperature, solarRadiation, 
     +     stationPressure, stationPressureQCD, temperature, 
     +     temperatureQCD, totalColumnPWV, totalSignalDelay, 
     +     visibility, visibilityQCD, wetSignalDelay, windDir, 
     +     windDirMax, windDirQCD, windGust, windSpeed, windSpeedQCD, 
     +     fuelMoistChangeTime, fuelTempChangeTime, observationTime, 
     +     receivedTime, reportTime, rhChangeTime, 
     +     solarRadChangeTime, stationPressChangeTime, 
     +     tempChangeTime, timeSinceLastPcp, windDirChangeTime, 
     +     windGustChangeTime, windSpeedChangeTime, ICT, QCT, 
     +     altimeterDD, dataProvider, dewpointDD, handbook5Id, 
     +     homeWFO, namePST, precipAccumDD, precipRateDD, 
     +     pressChange3HourDD, providerId, rawMessage, relHumidityDD, 
     +     seaLevelPressureDD, staticIds, stationId, stationName, 
     +     stationPressureDD, stationType, temperatureDD, test1, 
     +     typePST, visibilityDD, windDirDD, windSpeedDD)
C
      include 'netcdf.inc'
      integer ICcheckNum, QCcheckNum, maxPSTEntries, maxSensor, 
     +     maxStaticIds, nInventoryBins, recNum,nf_fid, nf_vid, 
     +     nf_status
      integer altimeterQCA(recNum), altimeterQCR(recNum),
     +     code1PST(maxPSTEntries), code2PST(maxPSTEntries),
     +     code3PST(maxPSTEntries), dewpointICA(recNum),
     +     dewpointICR(recNum), dewpointQCA(recNum),
     +     dewpointQCR(recNum), filterSetNum,
     +     firstInBin(nInventoryBins), firstOverflow,
     +     globalInventory, invTime(recNum), inventory(maxStaticIds),
     +     isOverflow(recNum), lastInBin(nInventoryBins),
     +     lastRecord(maxStaticIds), nStaticIds, numPST,
     +     numericWMOid(recNum), precipAccumICA(recNum),
     +     precipAccumICR(recNum), precipAccumQCA(recNum),
     +     precipAccumQCR(recNum), precipIntensity( maxSensor,
     +     recNum), precipRateQCA(recNum), precipRateQCR(recNum),
     +     precipType( maxSensor, recNum),
     +     pressChange3HourICA(recNum), pressChange3HourICR(recNum),
     +     pressChange3HourQCA(recNum), pressChange3HourQCR(recNum),
     +     pressChangeChar(recNum), prevRecord(recNum),
     +     relHumidityQCA(recNum), relHumidityQCR(recNum),
     +     roadState1(recNum), roadState2(recNum),
     +     roadState3(recNum), roadState4(recNum),
     +     seaLevelPressureICA(recNum), seaLevelPressureICR(recNum),
     +     seaLevelPressureQCA(recNum), seaLevelPressureQCR(recNum),
     +     stationPressureICA(recNum), stationPressureICR(recNum),
     +     stationPressureQCA(recNum), stationPressureQCR(recNum),
     +     temperatureICA(recNum), temperatureICR(recNum),
     +     temperatureQCA(recNum), temperatureQCR(recNum),
     +     visibilityICA(recNum), visibilityICR(recNum),
     +     visibilityQCA(recNum), visibilityQCR(recNum),
     +     windDirICA(recNum), windDirICR(recNum),
     +     windDirQCA(recNum), windDirQCR(recNum),
     +     windSpeedICA(recNum), windSpeedICR(recNum),
     +     windSpeedQCA(recNum), windSpeedQCR(recNum)
      real altimeter(recNum), altimeterQCD( QCcheckNum, recNum),
     +     capPi(recNum), dewpoint(recNum), dewpointQCD( QCcheckNum,
     +     recNum), drySignalDelay(recNum), elevation(recNum),
     +     formalError(recNum), fuelMoisture(recNum),
     +     fuelTemperature(recNum), latitude(recNum),
     +     longitude(recNum), meanWeightedTemperature(recNum),
     +     precipAccum(recNum), precipAccumQCD( QCcheckNum, recNum),
     +     precipRate(recNum), precipRateQCD( QCcheckNum, recNum),
     +     pressChange3Hour(recNum), pressChange3HourQCD( QCcheckNum,
     +     recNum), rawPrecip(recNum), relHumidity(recNum),
     +     relHumidityQCD( QCcheckNum, recNum),
     +     roadLiquidChemFactor1(recNum),
     +     roadLiquidChemFactor2(recNum),
     +     roadLiquidChemFactor3(recNum),
     +     roadLiquidChemFactor4(recNum),
     +     roadLiquidChemPercent1(recNum),
     +     roadLiquidChemPercent2(recNum),
     +     roadLiquidChemPercent3(recNum),
     +     roadLiquidChemPercent4(recNum), roadLiquidDepth1(recNum),
     +     roadLiquidDepth2(recNum), roadLiquidDepth3(recNum),
     +     roadLiquidDepth4(recNum), roadLiquidFreezeTemp1(recNum),
     +     roadLiquidFreezeTemp2(recNum),
     +     roadLiquidFreezeTemp3(recNum),
     +     roadLiquidFreezeTemp4(recNum),
     +     roadLiquidIcePercent1(recNum),
     +     roadLiquidIcePercent2(recNum),
     +     roadLiquidIcePercent3(recNum),
     +     roadLiquidIcePercent4(recNum),
     +     roadSubsurfaceTemp1(recNum), roadSubsurfaceTemp2(recNum),
     +     roadSubsurfaceTemp3(recNum), roadSubsurfaceTemp4(recNum),
     +     roadTemperature1(recNum), roadTemperature2(recNum),
     +     roadTemperature3(recNum), roadTemperature4(recNum),
     +     seaLevelPressure(recNum), seaLevelPressureQCD( QCcheckNum,
     +     recNum), soilMoisture(recNum), soilTemperature(recNum),
     +     solarRadiation(recNum), stationPressure(recNum),
     +     stationPressureQCD( QCcheckNum, recNum),
     +     temperature(recNum), temperatureQCD( QCcheckNum, recNum),
     +     totalColumnPWV(recNum), totalSignalDelay(recNum),
     +     visibility(recNum), visibilityQCD( QCcheckNum, recNum),
     +     wetSignalDelay(recNum), windDir(recNum),
     +     windDirMax(recNum), windDirQCD( QCcheckNum, recNum),
     +     windGust(recNum), windSpeed(recNum), windSpeedQCD(
     +     QCcheckNum, recNum)
      double precision fuelMoistChangeTime(recNum),
     +     fuelTempChangeTime(recNum), observationTime(recNum),
     +     receivedTime(recNum), reportTime(recNum),
     +     rhChangeTime(recNum), solarRadChangeTime(recNum),
     +     stationPressChangeTime(recNum), tempChangeTime(recNum),
     +     timeSinceLastPcp(recNum), windDirChangeTime(recNum),
     +     windGustChangeTime(recNum), windSpeedChangeTime(recNum)
      character*51 stationName(recNum)
      character visibilityDD(recNum)
      character*24 staticIds(maxStaticIds)
      character*11 namePST(maxPSTEntries)
      character*6 handbook5Id(recNum)
      character*11 dataProvider(recNum)
      character*6 stationId(recNum)
      character dewpointDD(recNum)
      character*4 homeWFO(recNum)
      character*60 QCT(QCcheckNum)
      character*12 providerId(recNum)
      character*512 rawMessage(recNum)
      character relHumidityDD(recNum)
      character seaLevelPressureDD(recNum)
      character temperatureDD(recNum)
      character precipAccumDD(recNum)
      character typePST(maxPSTEntries)
      character altimeterDD(recNum)
      character stationPressureDD(recNum)
      character precipRateDD(recNum)
      character*11 stationType(recNum)
      character pressChange3HourDD(recNum)
      character*72 ICT(ICcheckNum)
      character windSpeedDD(recNum)
      character*51 test1(recNum)
      character windDirDD(recNum)


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
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,altimeter)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeter'
      endif
C
C     Variable        NETCDF Long Name
C      altimeterQCD "altimeter setting QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,altimeterQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterQCD'
      endif
C
C     Variable        NETCDF Long Name
C      capPi        "Wet delay mapping function"
C
      nf_status=NF_INQ_VARID(nf_fid,'capPi',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var capPi'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,capPi)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var capPi'
      endif
C
C     Variable        NETCDF Long Name
C      dewpoint     "dew point temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpoint',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dewpoint)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpoint'
      endif
C
C     Variable        NETCDF Long Name
C      dewpointQCD  "dew point QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dewpointQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointQCD'
      endif
C
C     Variable        NETCDF Long Name
C      drySignalDelay"Dry component GPS signal delay"
C
      nf_status=NF_INQ_VARID(nf_fid,'drySignalDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drySignalDelay'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,drySignalDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drySignalDelay'
      endif
C
C     Variable        NETCDF Long Name
C      elevation    "elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
C
C     Variable        NETCDF Long Name
C      formalError  "Formal Error"
C
      nf_status=NF_INQ_VARID(nf_fid,'formalError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var formalError'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,formalError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var formalError'
      endif
C
C     Variable        NETCDF Long Name
C      fuelMoisture "Fuel moisture"
C
      nf_status=NF_INQ_VARID(nf_fid,'fuelMoisture',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelMoisture'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,fuelMoisture)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelMoisture'
      endif
C
C     Variable        NETCDF Long Name
C      fuelTemperature"Fuel temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'fuelTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelTemperature'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,fuelTemperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelTemperature'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
C
C     Variable        NETCDF Long Name
C      meanWeightedTemperature"Mean weighted temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'meanWeightedTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var meanWeightedTemperature'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,meanWeightedTemperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var meanWeightedTemperature'
      endif
C
C     Variable        NETCDF Long Name
C      precipAccum  "precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccum'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccum'
      endif
C
C     Variable        NETCDF Long Name
C      precipAccumQCD"precip amount QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipAccumQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumQCD'
      endif
C
C     Variable        NETCDF Long Name
C      precipRate   "precipitation rate"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRate'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipRate)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRate'
      endif
C
C     Variable        NETCDF Long Name
C      precipRateQCD"precip rate QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipRateQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateQCD'
      endif
C
C     Variable        NETCDF Long Name
C      pressChange3Hour"3 hour pressure change"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3Hour',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3Hour'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressChange3Hour)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3Hour'
      endif
C
C     Variable        NETCDF Long Name
C      pressChange3HourQCD"3h pressure change QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressChange3HourQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourQCD'
      endif
C
C     Variable        NETCDF Long Name
C      rawPrecip    "raw precip accumulation"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawPrecip',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawPrecip'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rawPrecip)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawPrecip'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidity  "relative humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidity'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityQCD"relative humidity QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidityQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityQCD'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemFactor1"Road liquid chem factor  - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemFactor1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemFactor1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor1'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemFactor2"Road liquid chem factor  - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemFactor2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemFactor2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor2'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemFactor3"Road liquid chem factor  - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemFactor3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemFactor3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor3'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemFactor4"Road liquid chem factor  - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemFactor4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemFactor4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemFactor4'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemPercent1"Road liquid chem percent  - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemPercent1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemPercent1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent1'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemPercent2"Road liquid chem percent  - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemPercent2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemPercent2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent2'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemPercent3"Road liquid chem percent  - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemPercent3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemPercent3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent3'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidChemPercent4"Road liquid chem percent  - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidChemPercent4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidChemPercent4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidChemPercent4'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidDepth1"Road liquid depth  - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidDepth1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidDepth1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth1'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidDepth2"Road liquid depth  - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidDepth2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidDepth2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth2'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidDepth3"Road liquid depth  - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidDepth3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidDepth3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth3'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidDepth4"Road liquid depth  - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidDepth4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidDepth4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidDepth4'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidFreezeTemp1"Road liquid freezing temp - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidFreezeTemp1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidFreezeTemp1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp1'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidFreezeTemp2"Road liquid freezing temp - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidFreezeTemp2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidFreezeTemp2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp2'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidFreezeTemp3"Road liquid freezing temp - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidFreezeTemp3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidFreezeTemp3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp3'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidFreezeTemp4"Road liquid freezing temp - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidFreezeTemp4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidFreezeTemp4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidFreezeTemp4'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidIcePercent1"Road liquid ice percent  - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidIcePercent1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidIcePercent1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent1'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidIcePercent2"Road liquid ice percent  - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidIcePercent2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidIcePercent2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent2'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidIcePercent3"Road liquid ice percent  - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidIcePercent3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidIcePercent3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent3'
      endif
C
C     Variable        NETCDF Long Name
C      roadLiquidIcePercent4"Road liquid ice percent  - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadLiquidIcePercent4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadLiquidIcePercent4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadLiquidIcePercent4'
      endif
C
C     Variable        NETCDF Long Name
C      roadSubsurfaceTemp1"Road subsurface temp - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadSubsurfaceTemp1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadSubsurfaceTemp1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp1'
      endif
C
C     Variable        NETCDF Long Name
C      roadSubsurfaceTemp2"Road subsurface temp - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadSubsurfaceTemp2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadSubsurfaceTemp2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp2'
      endif
C
C     Variable        NETCDF Long Name
C      roadSubsurfaceTemp3"Road subsurface temp - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadSubsurfaceTemp3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadSubsurfaceTemp3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp3'
      endif
C
C     Variable        NETCDF Long Name
C      roadSubsurfaceTemp4"Road subsurface temp - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadSubsurfaceTemp4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadSubsurfaceTemp4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadSubsurfaceTemp4'
      endif
C
C     Variable        NETCDF Long Name
C      roadTemperature1"Road temperature - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadTemperature1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature1'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadTemperature1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature1'
      endif
C
C     Variable        NETCDF Long Name
C      roadTemperature2"Road temperature - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadTemperature2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature2'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadTemperature2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature2'
      endif
C
C     Variable        NETCDF Long Name
C      roadTemperature3"Road temperature - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadTemperature3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature3'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadTemperature3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature3'
      endif
C
C     Variable        NETCDF Long Name
C      roadTemperature4"Road temperature - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadTemperature4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature4'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,roadTemperature4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadTemperature4'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressure"sea level pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressure'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressureQCD"sea level pressure QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,seaLevelPressureQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureQCD'
      endif
C
C     Variable        NETCDF Long Name
C      soilMoisture "Soil moisture"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilMoisture',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilMoisture'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilMoisture)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilMoisture'
      endif
C
C     Variable        NETCDF Long Name
C      soilTemperature"Soil temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'soilTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilTemperature'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,soilTemperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var soilTemperature'
      endif
C
C     Variable        NETCDF Long Name
C      solarRadiation"solar radiation"
C
      nf_status=NF_INQ_VARID(nf_fid,'solarRadiation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadiation'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,solarRadiation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadiation'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressure"station pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressure'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressureQCD"station pressure QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,stationPressureQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureQCD'
      endif
C
C     Variable        NETCDF Long Name
C      temperature  "temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureQCD"temperature QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperatureQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureQCD'
      endif
C
C     Variable        NETCDF Long Name
C      totalColumnPWV"Total column precipitable water vapor"
C
      nf_status=NF_INQ_VARID(nf_fid,'totalColumnPWV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var totalColumnPWV'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,totalColumnPWV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var totalColumnPWV'
      endif
C
C     Variable        NETCDF Long Name
C      totalSignalDelay"Total GPS signal delay"
C
      nf_status=NF_INQ_VARID(nf_fid,'totalSignalDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var totalSignalDelay'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,totalSignalDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var totalSignalDelay'
      endif
C
C     Variable        NETCDF Long Name
C      visibility   "visibility"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibility',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,visibility)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibility'
      endif
C
C     Variable        NETCDF Long Name
C      visibilityQCD"visibility QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,visibilityQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityQCD'
      endif
C
C     Variable        NETCDF Long Name
C      wetSignalDelay"Wet component GPS signal delay"
C
      nf_status=NF_INQ_VARID(nf_fid,'wetSignalDelay',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wetSignalDelay'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wetSignalDelay)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wetSignalDelay'
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "wind direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
C
C     Variable        NETCDF Long Name
C      windDirMax   "wind direction at gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirMax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDirMax)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirMax'
      endif
C
C     Variable        NETCDF Long Name
C      windDirQCD   "wind direction QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDirQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirQCD'
      endif
C
C     Variable        NETCDF Long Name
C      windGust     "wind gust"
C
      nf_status=NF_INQ_VARID(nf_fid,'windGust',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windGust)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windGust'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "wind speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedQCD "wind speed QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedQCD'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeedQCD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedQCD'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      altimeterQCA "altimeter setting QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,altimeterQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterQCA'
      endif
C
C     Variable        NETCDF Long Name
C      altimeterQCR "altimeter setting QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'altimeterQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,altimeterQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var altimeterQCR'
      endif
C
C     Variable        NETCDF Long Name
C      code1PST     "precipAccum variable definition"
C
      nf_status=NF_INQ_VARID(nf_fid,'code1PST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var code1PST'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,code1PST)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var code1PST'
      endif
C
C     Variable        NETCDF Long Name
C      code2PST     "solarRadiation variable definition"
C
      nf_status=NF_INQ_VARID(nf_fid,'code2PST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var code2PST'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,code2PST)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var code2PST'
      endif
C
C     Variable        NETCDF Long Name
C      code3PST     "PST code1/code2 usage definition"
C
      nf_status=NF_INQ_VARID(nf_fid,'code3PST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var code3PST'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,code3PST)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var code3PST'
      endif
C
C     Variable        NETCDF Long Name
C      dewpointICA  "dew point IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,dewpointICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointICA'
      endif
C
C     Variable        NETCDF Long Name
C      dewpointICR  "dew point IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,dewpointICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointICR'
      endif
C
C     Variable        NETCDF Long Name
C      dewpointQCA  "dew point QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,dewpointQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointQCA'
      endif
C
C     Variable        NETCDF Long Name
C      dewpointQCR  "dew point QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewpointQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,dewpointQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewpointQCR'
      endif
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
C      firstInBin   
C
      nf_status=NF_INQ_VARID(nf_fid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
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
C      invTime      
C
      nf_status=NF_INQ_VARID(nf_fid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,invTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
C
C     Variable        NETCDF Long Name
C      inventory    
C
      nf_status=NF_INQ_VARID(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,inventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
C
C     Variable        NETCDF Long Name
C      isOverflow   
C
      nf_status=NF_INQ_VARID(nf_fid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,isOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      lastInBin    
C
      nf_status=NF_INQ_VARID(nf_fid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
C
C     Variable        NETCDF Long Name
C      lastRecord   
C
      nf_status=NF_INQ_VARID(nf_fid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
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
C      precipAccumICA"precip amount IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumICA'
      endif
C
C     Variable        NETCDF Long Name
C      precipAccumICR"precip amount IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumICR'
      endif
C
C     Variable        NETCDF Long Name
C      precipAccumQCA"precip amount QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumQCA'
      endif
C
C     Variable        NETCDF Long Name
C      precipAccumQCR"precip amount QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipAccumQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipAccumQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipAccumQCR'
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
C      precipRateQCA"precip rate QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipRateQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateQCA'
      endif
C
C     Variable        NETCDF Long Name
C      precipRateQCR"precip rate QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipRateQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,precipRateQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipRateQCR'
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
C      pressChange3HourICA"3h pressure change IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChange3HourICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourICA'
      endif
C
C     Variable        NETCDF Long Name
C      pressChange3HourICR"3h pressure change IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChange3HourICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourICR'
      endif
C
C     Variable        NETCDF Long Name
C      pressChange3HourQCA"3h pressure change QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChange3HourQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourQCA'
      endif
C
C     Variable        NETCDF Long Name
C      pressChange3HourQCR"3h pressure change QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressChange3HourQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressChange3HourQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressChange3HourQCR'
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
C
C     Variable        NETCDF Long Name
C      prevRecord   
C
      nf_status=NF_INQ_VARID(nf_fid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,prevRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityQCA"relative humidity QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityQCA'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityQCR"relative humidity QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityQCR'
      endif
C
C     Variable        NETCDF Long Name
C      roadState1   "Road State - sensor 1"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadState1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState1'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,roadState1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState1'
      endif
C
C     Variable        NETCDF Long Name
C      roadState2   "Road State - sensor 2"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadState2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState2'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,roadState2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState2'
      endif
C
C     Variable        NETCDF Long Name
C      roadState3   "Road State - sensor 3"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadState3',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState3'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,roadState3)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState3'
      endif
C
C     Variable        NETCDF Long Name
C      roadState4   "Road State - sensor 4"
C
      nf_status=NF_INQ_VARID(nf_fid,'roadState4',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState4'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,roadState4)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var roadState4'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressureICA"sea level pressure IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,seaLevelPressureICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureICA'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressureICR"sea level pressure IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,seaLevelPressureICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureICR'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressureQCA"sea level pressure QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,seaLevelPressureQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureQCA'
      endif
C
C     Variable        NETCDF Long Name
C      seaLevelPressureQCR"sea level pressure QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'seaLevelPressureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,seaLevelPressureQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var seaLevelPressureQCR'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressureICA"station pressure IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,stationPressureICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureICA'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressureICR"station pressure IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,stationPressureICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureICR'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressureQCA"station pressure QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,stationPressureQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureQCA'
      endif
C
C     Variable        NETCDF Long Name
C      stationPressureQCR"station pressure QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationPressureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,stationPressureQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationPressureQCR'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureICA"temperature IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureICA'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureICR"temperature IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureICR'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureQCA"temperature QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureQCA'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureQCR"Temperature QC Results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureQCR'
      endif
C
C     Variable        NETCDF Long Name
C      visibilityICA"visibility QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,visibilityICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityICA'
      endif
C
C     Variable        NETCDF Long Name
C      visibilityICR"visibility IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,visibilityICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityICR'
      endif
C
C     Variable        NETCDF Long Name
C      visibilityQCA"visibility QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,visibilityQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityQCA'
      endif
C
C     Variable        NETCDF Long Name
C      visibilityQCR"visibility QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'visibilityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,visibilityQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var visibilityQCR'
      endif
C
C     Variable        NETCDF Long Name
C      windDirICA   "wind direction IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirICA'
      endif
C
C     Variable        NETCDF Long Name
C      windDirICR   "wind direction IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirICR'
      endif
C
C     Variable        NETCDF Long Name
C      windDirQCA   "wind direction QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirQCA'
      endif
C
C     Variable        NETCDF Long Name
C      windDirQCR   "wind direction QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDirQCR'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedICA "wind speed IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedICA'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedICR "wind speed IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedICR'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedQCA "wind speed QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedQCA'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeedQCR "wind speed QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeedQCR'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      fuelMoistChangeTime"fuel moisture time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'fuelMoistChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelMoistChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,fuelMoistChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelMoistChangeTime'
      endif
C
C     Variable        NETCDF Long Name
C      fuelTempChangeTime"fuel temperature time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'fuelTempChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelTempChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,fuelTempChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fuelTempChangeTime'
      endif
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
C      solarRadChangeTime"solar radiation time of last change"
C
      nf_status=NF_INQ_VARID(nf_fid,'solarRadChangeTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadChangeTime'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,solarRadChangeTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var solarRadChangeTime'
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
C      timeSinceLastPcp"time since last precip"
C
      nf_status=NF_INQ_VARID(nf_fid,'timeSinceLastPcp',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeSinceLastPcp'
      endif
      nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeSinceLastPcp)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeSinceLastPcp'
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
C      ICT          "list of possible IC checks"
C
      nf_status=NF_INQ_VARID(nf_fid,'ICT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ICT'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,ICT)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ICT'
      endif
C
C     Variable        NETCDF Long Name
C      QCT          "list of possible QC checks"
C
      nf_status=NF_INQ_VARID(nf_fid,'QCT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var QCT'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,QCT)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var QCT'
      endif
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
C      handbook5Id  "Handbook5 Id (AFOS or SHEF id)"
C
      nf_status=NF_INQ_VARID(nf_fid,'handbook5Id',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var handbook5Id'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,handbook5Id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var handbook5Id'
      endif
C
C     Variable        NETCDF Long Name
C      homeWFO      "home WFO Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'homeWFO',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var homeWFO'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,homeWFO)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var homeWFO'
      endif
C
C     Variable        NETCDF Long Name
C      namePST      "PST Provider or Subprovider name"
C
      nf_status=NF_INQ_VARID(nf_fid,'namePST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var namePST'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,namePST)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var namePST'
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
C      rawMessage   "raw text LDAD mesonet message"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawMessage',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawMessage'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,rawMessage)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rawMessage'
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
C      typePST      "PST type"
C
      nf_status=NF_INQ_VARID(nf_fid,'typePST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var typePST'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,typePST)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var typePST'
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
