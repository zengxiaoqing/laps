      subroutine get_wisdom_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     1                   ,i4time_earliest          
     1                   ,i4time_latest           
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer QCcheckNum, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
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
      call read_wisdom_data(nf_fid, QCcheckNum, maxStaticIds,
     +     nInventoryBins, recNum, i4time_sys, ilaps_cycle_time,
     +     NX_L, NY_L, i4time_earliest, i4time_latest, lun_out,
     +     istatus)

      return
      end
C
C
      subroutine read_wisdom_data(nf_fid, QCcheckNum, maxStaticIds,
     +     nInventoryBins, recNum, i4time_sys, ilaps_cycle_time,
     +     NX_L, NY_L, i4time_earliest, i4time_latest, lun_out,
     +     istatus)


      include 'netcdf.inc'
      integer QCcheckNum, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
      integer cutdownFlag(recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     latitudeQCA(recNum), latitudeQCR(recNum),
     +     longitudeQCA(recNum), longitudeQCR(recNum), nStaticIds,
     +     pressureQCA(recNum), pressureQCR(recNum),
     +     prevRecord(recNum), relHumidityQCA(recNum),
     +     relHumidityQCR(recNum), secondsStage1_2(recNum),
     +     temperatureQCA(recNum), temperatureQCR(recNum),
     +     windDirQCA(recNum), windDirQCR(recNum),
     +     windSpeedQCA(recNum), windSpeedQCR(recNum)
      real circularError(recNum), elevation(recNum), latitude(recNum),
     +     latitudeQCD( QCcheckNum, recNum), longitude(recNum),
     +     longitudeQCD( QCcheckNum, recNum), mobileElev(recNum),
     +     mobileLat(recNum), mobileLon(recNum), pressure(recNum),
     +     pressureQCD( QCcheckNum, recNum), relHumidity(recNum),
     +     relHumidityQCD( QCcheckNum, recNum), temperature(recNum),
     +     temperatureQCD( QCcheckNum, recNum), windDir(recNum),
     +     windDirQCD( QCcheckNum, recNum), windSpeed(recNum),
     +     windSpeedQCD( QCcheckNum, recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum)
      character pressureDD(recNum)
      character*11 dataProvider(recNum)
      character longitudeDD(recNum)
      character windDirDD(recNum)
      character*60 QCT(QCcheckNum)
      character*6 handbook5Id(recNum)
      character*51 stationName(recNum)
      character*11 stationType(recNum)
      character latitudeDD(recNum)
      character relHumidityDD(recNum)
      character*512 rawMessage(recNum)
      character temperatureDD(recNum)
      character*24 staticIds(maxStaticIds)
      character windSpeedDD(recNum)
      character*12 providerId(recNum)
      character*6 stationId(recNum)
      character*4 homeWFO(recNum)

!     Declarations for 'write_pin' call
      integer iwmostanum(recNum)
      character a9time_ob_r(recNum)*9
      logical l_closest_time, l_closest_time_i, l_in_domain, l_geoalt
      logical l_debug
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

      call read_wisdom_netcdf(nf_fid, QCcheckNum, maxStaticIds, 
     +     nInventoryBins, recNum, cutdownFlag, firstInBin, 
     +     firstOverflow, globalInventory, invTime, inventory, 
     +     isOverflow, lastInBin, lastRecord, latitudeQCA, 
     +     latitudeQCR, longitudeQCA, longitudeQCR, nStaticIds, 
     +     pressureQCA, pressureQCR, prevRecord, relHumidityQCA, 
     +     relHumidityQCR, secondsStage1_2, temperatureQCA, 
     +     temperatureQCR, windDirQCA, windDirQCR, windSpeedQCA, 
     +     windSpeedQCR, circularError, elevation, latitude, 
     +     latitudeQCD, longitude, longitudeQCD, mobileElev, 
     +     mobileLat, mobileLon, pressure, pressureQCD, relHumidity, 
     +     relHumidityQCD, temperature, temperatureQCD, windDir, 
     +     windDirQCD, windSpeed, windSpeedQCD, QCT, dataProvider, 
     +     handbook5Id, homeWFO, latitudeDD, longitudeDD, pressureDD, 
     +     providerId, rawMessage, relHumidityDD, staticIds, 
     +     stationId, stationName, stationType, temperatureDD, 
     +     windDirDD, windSpeedDD, observationTime, receivedTime, 
     +     reportTime)
C
C The netcdf variables are filled - your pin write call may go here
C
      write(6,*)' # of WISDOM reports read in = ',recNum

!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          iwmostanum(iob) = 0
          if(abs(observationTime(iob)) .le. 1e10)then
              i4time_ob = idint(observationTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      icount_ob_written = 0

      do iob = 1,recNum
          l_geoalt = .true.

!         MADIS QC flag checks can be added here if desired

          if(iob .le. 10)then
              l_debug = .true.
          else
              l_debug = .false.
          endif

          call write_aircraft_sub(lun_out,'pin'
     1                           ,a9time_ob_r(iob),a9time_ob_r(iob)
     1                           ,i4time_sys
     1                           ,i4time_earliest          
     1                           ,i4time_latest            
     1                           ,latitude(iob),longitude(iob)
     1                           ,elevation(iob)
     1                           ,windDir(iob),windSpeed(iob)
     1                           ,temperature(iob),relHumidity(iob)
     1                           ,l_geoalt                          ! I
     1                           ,l_debug                           ! I
     1                           ,istat_ob)                         ! O

          icount_ob_written = icount_ob_written + istat_ob

      enddo ! iob

      write(6,*)' # of WISDOM reports written to PIN file = '
     1          ,icount_ob_written        

      return
      end
C
C  Subroutine to read the file "MADIS - WISDOM" 
C
      subroutine read_wisdom_netcdf(nf_fid, QCcheckNum, maxStaticIds, 
     +     nInventoryBins, recNum, cutdownFlag, firstInBin, 
     +     firstOverflow, globalInventory, invTime, inventory, 
     +     isOverflow, lastInBin, lastRecord, latitudeQCA, 
     +     latitudeQCR, longitudeQCA, longitudeQCR, nStaticIds, 
     +     pressureQCA, pressureQCR, prevRecord, relHumidityQCA, 
     +     relHumidityQCR, secondsStage1_2, temperatureQCA, 
     +     temperatureQCR, windDirQCA, windDirQCR, windSpeedQCA, 
     +     windSpeedQCR, circularError, elevation, latitude, 
     +     latitudeQCD, longitude, longitudeQCD, mobileElev, 
     +     mobileLat, mobileLon, pressure, pressureQCD, relHumidity, 
     +     relHumidityQCD, temperature, temperatureQCD, windDir, 
     +     windDirQCD, windSpeed, windSpeedQCD, QCT, dataProvider, 
     +     handbook5Id, homeWFO, latitudeDD, longitudeDD, pressureDD, 
     +     providerId, rawMessage, relHumidityDD, staticIds, 
     +     stationId, stationName, stationType, temperatureDD, 
     +     windDirDD, windSpeedDD, observationTime, receivedTime, 
     +     reportTime)
C
      include 'netcdf.inc'
      integer QCcheckNum, maxStaticIds, nInventoryBins, recNum,nf_fid, 
     +     nf_vid, nf_status
      integer cutdownFlag(recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     latitudeQCA(recNum), latitudeQCR(recNum),
     +     longitudeQCA(recNum), longitudeQCR(recNum), nStaticIds,
     +     pressureQCA(recNum), pressureQCR(recNum),
     +     prevRecord(recNum), relHumidityQCA(recNum),
     +     relHumidityQCR(recNum), secondsStage1_2(recNum),
     +     temperatureQCA(recNum), temperatureQCR(recNum),
     +     windDirQCA(recNum), windDirQCR(recNum),
     +     windSpeedQCA(recNum), windSpeedQCR(recNum)
      real circularError(recNum), elevation(recNum), latitude(recNum),
     +     latitudeQCD( QCcheckNum, recNum), longitude(recNum),
     +     longitudeQCD( QCcheckNum, recNum), mobileElev(recNum),
     +     mobileLat(recNum), mobileLon(recNum), pressure(recNum),
     +     pressureQCD( QCcheckNum, recNum), relHumidity(recNum),
     +     relHumidityQCD( QCcheckNum, recNum), temperature(recNum),
     +     temperatureQCD( QCcheckNum, recNum), windDir(recNum),
     +     windDirQCD( QCcheckNum, recNum), windSpeed(recNum),
     +     windSpeedQCD( QCcheckNum, recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum)
      character pressureDD(recNum)
      character*11 dataProvider(recNum)
      character longitudeDD(recNum)
      character windDirDD(recNum)
      character*60 QCT(QCcheckNum)
      character*6 handbook5Id(recNum)
      character*51 stationName(recNum)
      character*11 stationType(recNum)
      character latitudeDD(recNum)
      character relHumidityDD(recNum)
      character*512 rawMessage(recNum)
      character*24 staticIds(maxStaticIds)
      character temperatureDD(recNum)
      character*12 providerId(recNum)
      character windSpeedDD(recNum)
      character*4 homeWFO(recNum)
      character*6 stationId(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     circularError "position fix error"
C
      nf_status=NF_INQ_VARID(nf_fid,'circularError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for circularError'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,circularError)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for circularError'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     elevation     "elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for elevation'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for elevation'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitude      "latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitude'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitude'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitudeQCD   "latitude QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitudeQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitudeQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,latitudeQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitudeQCD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitude     "longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitude'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitude'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitudeQCD  "longitude QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitudeQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitudeQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,longitudeQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitudeQCD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     mobileElev    "mobile elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'mobileElev',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mobileElev'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,mobileElev)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mobileElev'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     mobileLat     "mobile latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'mobileLat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mobileLat'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,mobileLat)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mobileLat'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     mobileLon     "mobile longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'mobileLon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mobileLon'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,mobileLon)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mobileLon'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressure      "pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressure'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressure'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressureQCD   "pressure QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressureQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressureQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressureQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressureQCD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidity   "relHumidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidity'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidity)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidity'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityQCD"relHumidity QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidityQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumidityQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidityQCD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperature   "temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperature'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperature'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperatureQCD"temperature QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperatureQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperatureQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperatureQCD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDir       "wind direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDir'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDir)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDir'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirQCD    "wind direction QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windDirQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirQCD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeed     "wind speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeed'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeed'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedQCD  "wind speed QC departures"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeedQCD'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeedQCD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeedQCD'
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     cutdownFlag   "cutdown flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'cutdownFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for cutdownFlag'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,cutdownFlag)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for cutdownFlag'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     firstInBin    
C
      nf_status=NF_INQ_VARID(nf_fid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for firstInBin'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,firstInBin)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for firstInBin'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     firstOverflow 
C
      nf_status=NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for firstOverflow'
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
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for globalInventory'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     invTime       
C
      nf_status=NF_INQ_VARID(nf_fid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for invTime'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,invTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for invTime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     inventory     
C
      nf_status=NF_INQ_VARID(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for inventory'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,inventory)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for inventory'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     isOverflow    
C
      nf_status=NF_INQ_VARID(nf_fid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for isOverflow'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,isOverflow)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for isOverflow'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     lastInBin     
C
      nf_status=NF_INQ_VARID(nf_fid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for lastInBin'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastInBin)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for lastInBin'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     lastRecord    
C
      nf_status=NF_INQ_VARID(nf_fid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for lastRecord'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,lastRecord)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for lastRecord'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitudeQCA   "latitude QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitudeQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitudeQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,latitudeQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitudeQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitudeQCR   "latitude QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitudeQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitudeQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,latitudeQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitudeQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitudeQCA  "longitude QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitudeQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitudeQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,longitudeQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitudeQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitudeQCR  "longitude QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitudeQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitudeQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,longitudeQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitudeQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     nStaticIds    
C
      nf_status=NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for nStaticIds'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for nStaticIds'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressureQCA   "pressure QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressureQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressureQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressureQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressureQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressureQCR   "pressure QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressureQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,pressureQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressureQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     prevRecord    
C
      nf_status=NF_INQ_VARID(nf_fid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for prevRecord'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,prevRecord)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for prevRecord'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityQCA"relHumidity QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidityQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidityQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityQCR"relHumidity QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidityQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidityQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     secondsStage1_2
C
      nf_status=NF_INQ_VARID(nf_fid,'secondsStage1_2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for secondsStage1_2'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,secondsStage1_2)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for secondsStage1_2'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperatureQCA"temperature QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperatureQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperatureQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperatureQCR"temperature QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for temperatureQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,temperatureQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperatureQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirQCA    "wind direction QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirQCR    "wind direction QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windDirQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windDirQCR'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedQCA  "wind speed QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeedQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeedQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windSpeedQCR  "wind speed QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'windSpeedQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windSpeedQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,windSpeedQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for windSpeedQCR'
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
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for observationTime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     receivedTime  "time data was received"
C
      nf_status=NF_INQ_VARID(nf_fid,'receivedTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for receivedTime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receivedTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for receivedTime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     reportTime    "time data was processed by the provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for reportTime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for reportTime'
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     QCT           "list of possible QC checks"
C
      nf_status=NF_INQ_VARID(nf_fid,'QCT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for QCT'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,QCT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for QCT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dataProvider  "Local data provider"
C
      nf_status=NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dataProvider'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dataProvider'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     handbook5Id   "Handbook5 Id (AFOS or SHEF id)"
C
      nf_status=NF_INQ_VARID(nf_fid,'handbook5Id',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for handbook5Id'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,handbook5Id)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for handbook5Id'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     homeWFO       "home WFO Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'homeWFO',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for homeWFO'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,homeWFO)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for homeWFO'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     latitudeDD    "latitude QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'latitudeDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for latitudeDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,latitudeDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for latitudeDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     longitudeDD   "longitude QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'longitudeDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for longitudeDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,longitudeDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for longitudeDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressureDD    "pressure QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressureDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for pressureDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,pressureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for pressureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     providerId    "Data Provider station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for providerId'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for providerId'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rawMessage    "raw text LDAD mesonet message"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawMessage',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rawMessage'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,rawMessage)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rawMessage'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     relHumidityDD "relHumidity QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for relHumidityDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,relHumidityDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for relHumidityDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staticIds     
C
      nf_status=NF_INQ_VARID(nf_fid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staticIds'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,staticIds)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staticIds'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationId     "alphanumeric station Id"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationId'
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
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for stationName'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     stationType   "LDAD station type"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for stationType'
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
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,temperatureDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for temperatureDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     windDirDD     "wind direction QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'windDirDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for windDirDD'
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
