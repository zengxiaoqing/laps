      subroutine get_radiometer_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer ICcheckNum, QCcheckNum, level, maxStaticIds,
     +     nInventoryBins, recNum,nf_fid, nf_vid, nf_status
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
C Get size of level
C
      nf_status=NF_INQ_DIMID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,level)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
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
      call read_radiometer_data(nf_fid, ICcheckNum, QCcheckNum, level,
     +     maxStaticIds, nInventoryBins, recNum, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_radiometer_data(nf_fid, ICcheckNum, QCcheckNum,
     +     level, maxStaticIds, nInventoryBins, recNum, i4time_sys,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer ICcheckNum, QCcheckNum, level, maxStaticIds,
     +     nInventoryBins, recNum,nf_fid, nf_vid, nf_status
      integer cloudBaseTempICA(recNum), cloudBaseTempICR(recNum),
     +     cloudBaseTempQCA(recNum), cloudBaseTempQCR(recNum),
     +     firstInBin(nInventoryBins), firstOverflow,
     +     globalInventory, integratedLiquidICA(recNum),
     +     integratedLiquidICR(recNum), integratedLiquidQCA(recNum),
     +     integratedLiquidQCR(recNum), integratedVaporICA(recNum),
     +     integratedVaporICR(recNum), integratedVaporQCA(recNum),
     +     integratedVaporQCR(recNum), invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     liquidDensityICA( level, recNum), liquidDensityICR( level,
     +     recNum), liquidDensityQCA( level, recNum),
     +     liquidDensityQCR( level, recNum), nStaticIds,
     +     prevRecord(recNum), rainFlag(recNum), relHumidityICA(
     +     level, recNum), relHumidityICR( level, recNum),
     +     relHumidityQCA( level, recNum), relHumidityQCR( level,
     +     recNum), stationType(recNum), temperatureICA( level,
     +     recNum), temperatureICR( level, recNum), temperatureQCA(
     +     level, recNum), temperatureQCR( level, recNum),
     +     vaporDensityICA( level, recNum), vaporDensityICR( level,
     +     recNum), vaporDensityQCA( level, recNum), vaporDensityQCR(
     +     level, recNum), wmoStaNum(recNum)
      real cloudBaseTemp(recNum), elevation(recNum),
     +     integratedLiquid(recNum), integratedVapor(recNum),
     +     latitude(recNum), levels( level, recNum), liquidDensity(
     +     level, recNum), longitude(recNum), pressure(recNum),
     +     relHumidity( level, recNum), relHumiditySfc(recNum),
     +     temperature( level, recNum), temperatureSfc(recNum),
     +     vaporDensity( level, recNum)
      double precision observationTime(recNum)
      character*51 stationName(recNum)
      character*6 providerId(recNum)
      character vaporDensityDD( level, recNum)
      character*30 staticIds(maxStaticIds)
      character temperatureDD( level, recNum)
      character relHumidityDD( level, recNum)
      character*11 dataProvider(recNum)
      character*72 ICT(ICcheckNum)
      character cloudBaseTempDD(recNum)
      character*60 QCT(QCcheckNum)
      character liquidDensityDD( level, recNum)
      character integratedLiquidDD(recNum)
      character integratedVaporDD(recNum)

!     Declarations for 'write_snd' call
      integer iwmostanum(recNum)
      real stalat(level),stalon(level)
      character a9time_ob_r(recNum)*9,a9time_ob_l(level)*9
      character c8_obstype*8
      real height_m(level)
      real pressure_mb(level)
      real temp_c(level)
      real dewpoint_c(level)
      real dir_deg(level)
      real spd_mps(level)

      logical l_closest_time, l_closest_time_i
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

      integer max_lvls
      parameter (max_lvls=200)
      real liquid_a(max_lvls)

      common /write_snd_data/ cloud_base_temp,cloud_integrated_liquid
     1                       ,liquid_a

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

      call read_radiometer_netcdf(nf_fid, ICcheckNum, QCcheckNum, 
     +     level, maxStaticIds, nInventoryBins, recNum, 
     +     cloudBaseTempICA, cloudBaseTempICR, cloudBaseTempQCA, 
     +     cloudBaseTempQCR, firstInBin, firstOverflow, 
     +     globalInventory, integratedLiquidICA, integratedLiquidICR, 
     +     integratedLiquidQCA, integratedLiquidQCR, 
     +     integratedVaporICA, integratedVaporICR, 
     +     integratedVaporQCA, integratedVaporQCR, invTime, 
     +     inventory, isOverflow, lastInBin, lastRecord, 
     +     liquidDensityICA, liquidDensityICR, liquidDensityQCA, 
     +     liquidDensityQCR, nStaticIds, prevRecord, rainFlag, 
     +     relHumidityICA, relHumidityICR, relHumidityQCA, 
     +     relHumidityQCR, stationType, temperatureICA, 
     +     temperatureICR, temperatureQCA, temperatureQCR, 
     +     vaporDensityICA, vaporDensityICR, vaporDensityQCA, 
     +     vaporDensityQCR, wmoStaNum, cloudBaseTemp, elevation, 
     +     integratedLiquid, integratedVapor, latitude, levels, 
     +     liquidDensity, longitude, pressure, relHumidity, 
     +     relHumiditySfc, temperature, temperatureSfc, vaporDensity, 
     +     observationTime, ICT, QCT, cloudBaseTempDD, dataProvider, 
     +     integratedLiquidDD, integratedVaporDD, liquidDensityDD, 
     +     providerId, relHumidityDD, staticIds, stationName, 
     +     temperatureDD, vaporDensityDD)
C
C The netcdf variables are filled - your snd write call may go here
C
!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          read(providerId(iob),'(3x,i2)',err=101)iwmostanum(iob)
          goto 102
101       iwmostanum(iob) = iob
          write(6,*)' Warning: unreadable providerId '
     1             ,trim(providerId(iob))
102       continue

          if(abs(observationTime(iob)) .le. 1e10)then
              i4time_ob = idint(observationTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      c8_obstype = 'RADIOMTR'

      do iob = 1,recNum
          call convert_array(levels(:,iob),height_m,level
     1                      ,'none',r_missing_data,istatus)

          call addcon_miss(height_m,elevation(iob),height_m,level,1)

          stalat = latitude(iob)
          stalon = longitude(iob)

!         Convert arrays for a single sounding
          a9time_ob_l = a9time_ob_r(iob)

          pressure_mb = r_missing_data

          call convert_array(pressure(iob),pressure_mb(1),1
     1                      ,'pa_to_mb',r_missing_data,istatus)

          call convert_array(temperature(:,iob),temp_c,level
     1                      ,'k_to_c',r_missing_data,istatus)

          do ilvl = 1,level
              dewpoint_c(ilvl) =
     1            DWPT_laps(temp_c(ilvl),relHumidity(ilvl,iob))
          enddo

          dir_deg = r_missing_data

          spd_mps = r_missing_data


          call get_nlevels_snd(pressure_mb,height_m,r_missing_data
     +                        ,level,nlevels_snd)

          if(integratedVapor(iob) .gt. .06)nlevels_snd=0

          cloud_base_temp = cloudBaseTemp(iob)
          cloud_integrated_liquid = integratedLiquid(iob)

          l_closest_time = l_closest_time_i(iwmostanum,a9time_ob_r
     1                                     ,recNum,iob,i4time_sys
     1                                     ,istatus)

          if(nlevels_snd .gt. 0 .and. l_closest_time)then
              call filter_string(providerId(iob))
              write(6,*)' valid radiometer near analysis time'
     1                 ,providerId(iob)

              write(6,*)' cloud_base_temp/liq = ',cloud_base_temp
     1                  ,cloud_integrated_liquid

              write(6,*)' liquid density:'
              do ilvl = 1,nlevels_snd
                  write(6,*)ilvl,height_m(ilvl),liquidDensity(ilvl,iob) 
                  liquid_a(ilvl) = liquidDensity(ilvl,iob)       
              enddo ! ilvl

!             call 'write_snd' for a single profile
              call open_ext(lun_out,i4time_sys,'snd',istatus)

              call write_snd(lun_out
     +                      ,1,nlevels_snd,1
     +                      ,iwmostanum
     +                      ,stalat,stalon,elevation(iob)
     +                      ,providerId(iob)
     +                      ,a9time_ob_l,c8_obstype
     +                      ,nlevels_snd
     +                      ,height_m
     +                      ,pressure_mb
     +                      ,temp_c
     +                      ,dewpoint_c
     +                      ,dir_deg
     +                      ,spd_mps
     +                      ,istatus)
          endif ! valid profile

      enddo ! iob
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_radiometer_netcdf(nf_fid, ICcheckNum, 
     +     QCcheckNum, level, maxStaticIds, nInventoryBins, recNum, 
     +     cloudBaseTempICA, cloudBaseTempICR, cloudBaseTempQCA, 
     +     cloudBaseTempQCR, firstInBin, firstOverflow, 
     +     globalInventory, integratedLiquidICA, integratedLiquidICR, 
     +     integratedLiquidQCA, integratedLiquidQCR, 
     +     integratedVaporICA, integratedVaporICR, 
     +     integratedVaporQCA, integratedVaporQCR, invTime, 
     +     inventory, isOverflow, lastInBin, lastRecord, 
     +     liquidDensityICA, liquidDensityICR, liquidDensityQCA, 
     +     liquidDensityQCR, nStaticIds, prevRecord, rainFlag, 
     +     relHumidityICA, relHumidityICR, relHumidityQCA, 
     +     relHumidityQCR, stationType, temperatureICA, 
     +     temperatureICR, temperatureQCA, temperatureQCR, 
     +     vaporDensityICA, vaporDensityICR, vaporDensityQCA, 
     +     vaporDensityQCR, wmoStaNum, cloudBaseTemp, elevation, 
     +     integratedLiquid, integratedVapor, latitude, levels, 
     +     liquidDensity, longitude, pressure, relHumidity, 
     +     relHumiditySfc, temperature, temperatureSfc, vaporDensity, 
     +     observationTime, ICT, QCT, cloudBaseTempDD, dataProvider, 
     +     integratedLiquidDD, integratedVaporDD, liquidDensityDD, 
     +     providerId, relHumidityDD, staticIds, stationName, 
     +     temperatureDD, vaporDensityDD)
C
      include 'netcdf.inc'
      integer ICcheckNum, QCcheckNum, level, maxStaticIds, 
     +     nInventoryBins, recNum,nf_fid, nf_vid, nf_status
      integer cloudBaseTempICA(recNum), cloudBaseTempICR(recNum),
     +     cloudBaseTempQCA(recNum), cloudBaseTempQCR(recNum),
     +     firstInBin(nInventoryBins), firstOverflow,
     +     globalInventory, integratedLiquidICA(recNum),
     +     integratedLiquidICR(recNum), integratedLiquidQCA(recNum),
     +     integratedLiquidQCR(recNum), integratedVaporICA(recNum),
     +     integratedVaporICR(recNum), integratedVaporQCA(recNum),
     +     integratedVaporQCR(recNum), invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     liquidDensityICA( level, recNum), liquidDensityICR( level,
     +     recNum), liquidDensityQCA( level, recNum),
     +     liquidDensityQCR( level, recNum), nStaticIds,
     +     prevRecord(recNum), rainFlag(recNum), relHumidityICA(
     +     level, recNum), relHumidityICR( level, recNum),
     +     relHumidityQCA( level, recNum), relHumidityQCR( level,
     +     recNum), stationType(recNum), temperatureICA( level,
     +     recNum), temperatureICR( level, recNum), temperatureQCA(
     +     level, recNum), temperatureQCR( level, recNum),
     +     vaporDensityICA( level, recNum), vaporDensityICR( level,
     +     recNum), vaporDensityQCA( level, recNum), vaporDensityQCR(
     +     level, recNum), wmoStaNum(recNum)
      real cloudBaseTemp(recNum), elevation(recNum),
     +     integratedLiquid(recNum), integratedVapor(recNum),
     +     latitude(recNum), levels( level, recNum), liquidDensity(
     +     level, recNum), longitude(recNum), pressure(recNum),
     +     relHumidity( level, recNum), relHumiditySfc(recNum),
     +     temperature( level, recNum), temperatureSfc(recNum),
     +     vaporDensity( level, recNum)
      double precision observationTime(recNum)
      character*51 stationName(recNum)
      character*6 providerId(recNum)
      character vaporDensityDD( level, recNum)
      character*30 staticIds(maxStaticIds)
      character temperatureDD( level, recNum)
      character relHumidityDD( level, recNum)
      character*11 dataProvider(recNum)
      character*72 ICT(ICcheckNum)
      character cloudBaseTempDD(recNum)
      character*60 QCT(QCcheckNum)
      character liquidDensityDD( level, recNum)
      character integratedLiquidDD(recNum)
      character integratedVaporDD(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      cloudBaseTemp"Infrared cloud base temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudBaseTemp',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTemp'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,cloudBaseTemp)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTemp'
      endif
C
C     Variable        NETCDF Long Name
C      elevation    "Elevation above MSL"
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
C      integratedLiquid"Integrated liquid"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedLiquid',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquid'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,integratedLiquid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquid'
      endif
C
C     Variable        NETCDF Long Name
C      integratedVapor"Integrated vapor"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedVapor',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVapor'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,integratedVapor)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVapor'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "Station latitude"
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
C      levels       "Instrument Level, height above station"
C
      nf_status=NF_INQ_VARID(nf_fid,'levels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,levels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      endif
C
C     Variable        NETCDF Long Name
C      liquidDensity"Liquid density"
C
      nf_status=NF_INQ_VARID(nf_fid,'liquidDensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensity'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,liquidDensity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensity'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "Station longitude"
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
C      pressure     "Station pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressure'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var pressure'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidity  "Relative humidity"
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
C      relHumiditySfc"Surface relative humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumiditySfc',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumiditySfc'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,relHumiditySfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumiditySfc'
      endif
C
C     Variable        NETCDF Long Name
C      temperature  "Temperature"
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
C      temperatureSfc"Surface temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'temperatureSfc',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureSfc'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,temperatureSfc)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperatureSfc'
      endif
C
C     Variable        NETCDF Long Name
C      vaporDensity "Vapor density"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporDensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensity'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vaporDensity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensity'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      cloudBaseTempICA"Cloud base temp IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudBaseTempICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,cloudBaseTempICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempICA'
      endif
C
C     Variable        NETCDF Long Name
C      cloudBaseTempICR"Cloud base temp IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudBaseTempICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,cloudBaseTempICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempICR'
      endif
C
C     Variable        NETCDF Long Name
C      cloudBaseTempQCA"Cloud base temp QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudBaseTempQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,cloudBaseTempQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempQCA'
      endif
C
C     Variable        NETCDF Long Name
C      cloudBaseTempQCR"Cloud base temp QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudBaseTempQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,cloudBaseTempQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempQCR'
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
C      integratedLiquidICA"Integrated liquid IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedLiquidICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedLiquidICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidICA'
      endif
C
C     Variable        NETCDF Long Name
C      integratedLiquidICR"Integrated liquid IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedLiquidICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedLiquidICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidICR'
      endif
C
C     Variable        NETCDF Long Name
C      integratedLiquidQCA"Integrated liquid QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedLiquidQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedLiquidQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidQCA'
      endif
C
C     Variable        NETCDF Long Name
C      integratedLiquidQCR"Integrated liquid QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedLiquidQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedLiquidQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidQCR'
      endif
C
C     Variable        NETCDF Long Name
C      integratedVaporICA"Integrated vapor IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedVaporICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedVaporICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporICA'
      endif
C
C     Variable        NETCDF Long Name
C      integratedVaporICR"Integrated vapor IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedVaporICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedVaporICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporICR'
      endif
C
C     Variable        NETCDF Long Name
C      integratedVaporQCA"Integrated vapor QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedVaporQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedVaporQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporQCA'
      endif
C
C     Variable        NETCDF Long Name
C      integratedVaporQCR"Integrated vapor QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedVaporQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,integratedVaporQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporQCR'
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
C      liquidDensityICA"Liquid density IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'liquidDensityICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,liquidDensityICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityICA'
      endif
C
C     Variable        NETCDF Long Name
C      liquidDensityICR"Liquid density IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'liquidDensityICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,liquidDensityICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityICR'
      endif
C
C     Variable        NETCDF Long Name
C      liquidDensityQCA"Liquid density QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'liquidDensityQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,liquidDensityQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityQCA'
      endif
C
C     Variable        NETCDF Long Name
C      liquidDensityQCR"Liquid density QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'liquidDensityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,liquidDensityQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityQCR'
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
C      rainFlag     "Rain flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'rainFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rainFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,rainFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rainFlag'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityICA"Relative humidity IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityICA'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityICR"Relative humidity IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'relHumidityICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,relHumidityICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var relHumidityICR'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidityQCA"Relative humidity QC applied word"
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
C      relHumidityQCR"Relative humidity QC results word"
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
C      stationType  "Station type"
C
      nf_status=NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,stationType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
C
C     Variable        NETCDF Long Name
C      temperatureICA"Temperature IC applied word"
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
C      temperatureICR"Temperature IC results word"
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
C      temperatureQCA"Temperature QC applied word"
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
C      temperatureQCR"Temperature QC results word"
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
C      vaporDensityICA"Vapor density IC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporDensityICA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityICA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vaporDensityICA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityICA'
      endif
C
C     Variable        NETCDF Long Name
C      vaporDensityICR"Vapor density IC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporDensityICR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityICR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vaporDensityICR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityICR'
      endif
C
C     Variable        NETCDF Long Name
C      vaporDensityQCA"Vapor density QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporDensityQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityQCA'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vaporDensityQCA)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityQCA'
      endif
C
C     Variable        NETCDF Long Name
C      vaporDensityQCR"Vapor density QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporDensityQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityQCR'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,vaporDensityQCR)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityQCR'
      endif
C
C     Variable        NETCDF Long Name
C      wmoStaNum    "WMO numeric station ID"
C
      nf_status=NF_INQ_VARID(nf_fid,'wmoStaNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wmoStaNum'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wmoStaNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wmoStaNum'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      observationTime"Time of observation"
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
C      cloudBaseTempDD"Cloud base temp QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudBaseTempDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,cloudBaseTempDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var cloudBaseTempDD'
      endif
C
C     Variable        NETCDF Long Name
C      dataProvider "Name of organization responsible for delivering the data"
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
C      integratedLiquidDD"Integrated liquid QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedLiquidDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,integratedLiquidDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedLiquidDD'
      endif
C
C     Variable        NETCDF Long Name
C      integratedVaporDD"Integrated vapor QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'integratedVaporDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,integratedVaporDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var integratedVaporDD'
      endif
C
C     Variable        NETCDF Long Name
C      liquidDensityDD"Liquid density QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'liquidDensityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,liquidDensityDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var liquidDensityDD'
      endif
C
C     Variable        NETCDF Long Name
C      providerId   "Alphanumeric station name"
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
C      relHumidityDD"Relative humidity QC summary value"
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
C      temperatureDD"Temperature QC summary value"
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
C      vaporDensityDD"Vapor density QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporDensityDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityDD'
      endif
      nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,vaporDensityDD)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporDensityDD'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
