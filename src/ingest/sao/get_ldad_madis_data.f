      subroutine get_ldad_madis_data
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
      call read_ldad_madis_data(nf_fid, maxSensor, recNum, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_ldad_madis_data(nf_fid, maxSensor, recNum,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer maxSensor, recNum,nf_fid, nf_vid, nf_status
      integer altimeterQCR(recNum), dewpointQCR(recNum),
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
     +     seaSurfaceTemp(recNum), soilMoisture(recNum),
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
      character altimeterDD(recNum)
      character*51 stationName(recNum)
      character precipRateDD(recNum)
      character*11 stationType(recNum)
      character stationPressureDD(recNum)
      character seaLevelPressureDD(recNum)
      character windSpeedDD(recNum)

!     Declarations for 'write_lso' call
      integer iwmostanum(recNum)
      logical l_closest_time, l_closest_time_i
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

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

      call read_ldad_madis_netcdf(nf_fid, maxSensor, recNum, 
     +     altimeterQCR(ix), dewpointQCR(ix), firstOverflow, 
     +     globalInventory, nStaticIds, numericWMOid, 
     +     precipAccumQCR(ix), precipIntensity, 
     +     precipRateQCR(ix), precipType, pressChange3HourQCR(ix), 
     +     pressChangeChar, relHumidityQCR(ix), seaLevelPressureQCR(ix),       
     +     stationPressureQCR(ix), temperatureQCR(ix), 
     +     visibilityQCR(ix),       
     +     windDirQCR(ix), windSpeedQCR(ix), altimeter(ix), 
     +     dewpoint(ix), 
     +     elevation(ix), latitude(ix), longitude(ix), 
     +     meanWeightedTemperature(ix), precipAccum(ix), 
     +     precipRate(ix), pressChange3Hour(ix), relHumidity(ix), 
     +     seaLevelPressure(ix), seaSurfaceTemp(ix), 
     +     soilMoisture(ix), soilTemperature(ix), solarRadiation(ix), 
     +     stationPressure(ix), temperature(ix), visibility(ix), 
     +     windDir(ix), windDirMax(ix), windGust(ix), windSpeed(ix), 
     +     altimeterDD(ix), dataProvider(ix), dewpointDD(ix), 
     +     precipAccumDD(ix), precipRateDD(ix), presWeather(ix), 
     +     pressChange3HourDD(ix), providerId(ix), relHumidityDD(ix), 
     +     seaLevelPressureDD(ix), stationId(ix), stationName(ix), 
     +     stationPressureDD(ix), stationType(ix), temperatureDD(ix), 
     +     visibilityDD(ix), windDirDD(ix), windSpeedDD(ix), 
     +     observationTime(ix), receivedTime(ix), reportTime(ix), 
     +     rhChangeTime(ix), stationPressChangeTime(ix), 
     +     tempChangeTime(ix), windDirChangeTime(ix), 
     +     windGustChangeTime(ix), windSpeedChangeTime(ix),badflag)
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
