      subroutine get_ldadmadis_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer maxSensor, recNum,nf_fid, nf_vid, nf_status
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
      call read_ldadmadis_data(nf_fid, maxSensor, recNum, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_ldadmadis_data(nf_fid, maxSensor, recNum,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


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
      character*6 stationId(recNum)
      character dewpointDD(recNum)
      character seaLevelPressureDD(recNum)
      character visibilityDD(recNum)
      character precipAccumDD(recNum)
      character*51 stationName(recNum)
      character*12 providerId(recNum)
      character temperatureDD(recNum)

!     Declarations for 'write_lso' call
      integer iwmostanum(recNum)
      logical l_closest_time, l_closest_time_i
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

      call read_ldadmadis_netcdf(nf_fid, maxSensor, recNum, 
     +     filterSetNum, firstOverflow, globalInventory, nStaticIds, 
     +     numPST, numericWMOid, precipIntensity, precipType, 
     +     pressChangeChar, altimeter(ix), dewpoint(ix), 
     +     elevation(ix), latitude(ix), longitude(ix), 
     +     meanWeightedTemperature(ix), precipAccum(ix), 
     +     precipRate(ix), pressChange3Hour(ix), rawPrecip(ix), 
     +     relHumidity(ix), seaLevelPressure(ix), soilMoisture(ix), 
     +     soilTemperature(ix), solarRadiation(ix), 
     +     stationPressure(ix), temperature(ix), visibility(ix), 
     +     windDir(ix), windDirMax(ix), windGust(ix), windSpeed(ix), 
     +     observationTime, receivedTime, reportTime, rhChangeTime, 
     +     stationPressChangeTime, tempChangeTime, windDirChangeTime, 
     +     windGustChangeTime, windSpeedChangeTime, altimeterDD, 
     +     dataProvider, dewpointDD, precipAccumDD, precipRateDD, 
     +     pressChange3HourDD, providerId, relHumidityDD, 
     +     seaLevelPressureDD, staticIds, stationId, stationName, 
     +     stationPressureDD, stationType, temperatureDD, test1, 
     +     visibilityDD, windDirDD, windSpeedDD)
C
C The netcdf variables are filled - your lso write call may go here
C
!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          iwmostanum(iob) = 0
          if(abs(observationTime(iob)) .le. 1e10)then
              i4time_ob = idint(observationTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

!     c8_obstype = 

      do iob = 1,recNum
          l_closest_time = .true.

      enddo ! iob
      return
      end
