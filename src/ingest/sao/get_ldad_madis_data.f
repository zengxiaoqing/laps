      subroutine get_ldadmadis_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer ICcheckNum, QCcheckNum, maxPSTEntries, maxSensor,
     +     maxStaticIds, nInventoryBins, recNum,nf_fid, nf_vid,
     +     nf_status
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of ICcheckNum
C
      nf_status=NF_INQ_DIMID(nf_fid,'ICcheckNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim ICcheckNum'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,ICcheckNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim ICcheckNum'
      endif
C
C Get size of QCcheckNum
C
      nf_status=NF_INQ_DIMID(nf_fid,'QCcheckNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim QCcheckNum'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,QCcheckNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim QCcheckNum'
      endif
C
C Get size of maxPSTEntries
C
      nf_status=NF_INQ_DIMID(nf_fid,'maxPSTEntries',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxPSTEntries'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,maxPSTEntries)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxPSTEntries'
      endif
C
C Get size of maxSensor
C
      nf_status=NF_INQ_DIMID(nf_fid,'maxSensor',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxSensor'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,maxSensor)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxSensor'
      endif
C
C Get size of maxStaticIds
C
      nf_status=NF_INQ_DIMID(nf_fid,'maxStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxStaticIds'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,maxStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxStaticIds'
      endif
C
C Get size of nInventoryBins
C
      nf_status=NF_INQ_DIMID(nf_fid,'nInventoryBins',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nInventoryBins'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nInventoryBins)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nInventoryBins'
      endif
C
C Get size of recNum
C
      nf_status=NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      call read_ldadmadis_data(nf_fid, ICcheckNum, QCcheckNum,
     +     maxPSTEntries, maxSensor, maxStaticIds, nInventoryBins,
     +     recNum, i4time_sys, ilaps_cycle_time, NX_L, NY_L,
     +     i4time_earliest, i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_ldadmadis_data(nf_fid, ICcheckNum, QCcheckNum,
     +     maxPSTEntries, maxSensor, maxStaticIds, nInventoryBins,
     +     recNum, i4time_sys, ilaps_cycle_time, NX_L, NY_L,
     +     i4time_earliest, i4time_latest, lun_out, istatus)


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
      character dewpointDD(recNum)
      character*6 stationId(recNum)
      character*11 dataProvider(recNum)
      character*4 homeWFO(recNum)
      character*60 QCT(QCcheckNum)
      character*12 providerId(recNum)
      character*512 rawMessage(recNum)
      character relHumidityDD(recNum)
      character temperatureDD(recNum)
      character seaLevelPressureDD(recNum)
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

      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

      call read_ldadmadis_netcdf(nf_fid, ICcheckNum, QCcheckNum, 
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
C The netcdf variables are filled - your lso write call may go here
C
      do iob = 1,recNum
      enddo ! iob
      return
      end
