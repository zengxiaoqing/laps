      subroutine get_rtamps_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*170 filename

      integer manLevel, maxStaticIds, nInventoryBins, rawLevel,
     +     recNum, stdLevel, termLevel, tropLevel,nf_fid, nf_vid,
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
C Get size of manLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'manLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim manLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,manLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim manLevel'
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
C Get size of rawLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'rawLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim rawLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,rawLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim rawLevel'
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
C
C Get size of stdLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'stdLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim stdLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,stdLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim stdLevel'
      endif
C
C Get size of termLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'termLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim termLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,termLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim termLevel'
      endif
C
C Get size of tropLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'tropLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim tropLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,tropLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim tropLevel'
      endif
      call read_rtamps_data(nf_fid, manLevel, maxStaticIds,
     +     nInventoryBins, rawLevel, recNum, stdLevel, termLevel,
     +     tropLevel, i4time_sys, ilaps_cycle_time, NX_L, NY_L,
     +     i4time_earliest, i4time_latest, lun_out, istatus)

      return
      end
C
C
      subroutine read_rtamps_data(nf_fid, manLevel, maxStaticIds,
     +     nInventoryBins, rawLevel, recNum, stdLevel, termLevel,
     +     tropLevel, i4time_sys, ilaps_cycle_time, NX_L, NY_L,
     +     i4time_earliest, i4time_latest, lun_out, istatus)


      include 'netcdf.inc'
      integer manLevel, maxStaticIds, nInventoryBins, rawLevel,
     +     recNum, stdLevel, termLevel, tropLevel,nf_fid, nf_vid,
     +     nf_status
      integer editFlag( rawLevel, recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, indxRefr( rawLevel,
     +     recNum), invTime(recNum), inventory(maxStaticIds), irMan(
     +     manLevel, recNum), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, oiMan( manLevel, recNum), optIndxRefr(
     +     rawLevel, recNum), prevRecord(recNum), storedObs(recNum)
      real absHumidity( rawLevel, recNum), airDensity( rawLevel,
     +     recNum), baromPressure( rawLevel, recNum), bpMan(
     +     manLevel, recNum), bpStd( stdLevel, recNum), bpTerm(
     +     termLevel, recNum), bpTrop( tropLevel, recNum), dewPt(
     +     rawLevel, recNum), direction( rawLevel, recNum), dpMan(
     +     manLevel, recNum), dpStd( stdLevel, recNum), dpTrop(
     +     tropLevel, recNum), drMan( manLevel, recNum), drStd(
     +     stdLevel, recNum), elevation(recNum), geomHeight(
     +     rawLevel, recNum), geopHeight( rawLevel, recNum), ghMan(
     +     manLevel, recNum), ghTerm( termLevel, recNum), ghTrop(
     +     tropLevel, recNum), gpStd( stdLevel, recNum), gpTerm(
     +     termLevel, recNum), latitude(recNum), longitude(recNum),
     +     precipWater( rawLevel, recNum), relHumidity( rawLevel,
     +     recNum), rhMan( manLevel, recNum), rhStd( stdLevel,
     +     recNum), riseRate( rawLevel, recNum), shear( rawLevel,
     +     recNum), shearDir( rawLevel, recNum), shearMagX( rawLevel,
     +     recNum), shearMagY( rawLevel, recNum), spMan( manLevel,
     +     recNum), spStd( stdLevel, recNum), speed( rawLevel,
     +     recNum), temperature( rawLevel, recNum), tpMan( manLevel,
     +     recNum), tpStd( stdLevel, recNum), tpTrop( tropLevel,
     +     recNum), vaporPressure( rawLevel, recNum), velError(
     +     rawLevel, recNum), velSound( rawLevel, recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum)
      character*12 providerId(recNum)
      character*11 dataProvider(recNum)
      character*51 stationName(recNum)
      character*24 staticIds(maxStaticIds)

!     Declarations for 'write_snd' call
      integer iwmostanum(recNum)
      real stalat(rawLevel),stalon(rawLevel)
      character a9time_ob_r(recNum)*9,a9time_ob_l(rawLevel)*9
      character c8_obstype*8
      real height_m(rawLevel)
      real pressure_mb(rawLevel)
      real temp_c(rawLevel)
      real dewpoint_c(rawLevel)
      real dir_deg(rawLevel)
      real spd_mps(rawLevel)

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

      call read_rtamps_netcdf(nf_fid, manLevel, maxStaticIds, 
     +     nInventoryBins, rawLevel, recNum, stdLevel, termLevel, 
     +     tropLevel, editFlag, firstInBin, firstOverflow, 
     +     globalInventory, indxRefr, invTime, inventory, irMan, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, oiMan, 
     +     optIndxRefr, prevRecord, storedObs, absHumidity, 
     +     airDensity, baromPressure, bpMan, bpStd, bpTerm, bpTrop, 
     +     dewPt, direction, dpMan, dpStd, dpTrop, drMan, drStd, 
     +     elevation, geomHeight, geopHeight, ghMan, ghTerm, ghTrop, 
     +     gpStd, gpTerm, latitude, longitude, precipWater, 
     +     relHumidity, rhMan, rhStd, riseRate, shear, shearDir, 
     +     shearMagX, shearMagY, spMan, spStd, speed, temperature, 
     +     tpMan, tpStd, tpTrop, vaporPressure, velError, velSound, 
     +     observationTime, receivedTime, reportTime, dataProvider, 
     +     providerId, staticIds, stationName)
C
C The netcdf variables are filled - your snd write call may go here
C
!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          read(providerId(iob),*)iwmostanum(iob)
          if(abs(observationTime(iob)) .le. 1e10)then
              i4time_ob = idint(observationTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

      enddo ! iob

      c8_obstype = 'RAOB    '

      do iob = 1,recNum
          call convert_array(geopHeight(:,iob),height_m,rawLevel
     1                      ,'none',r_missing_data,istatus)

          call addcon_miss(height_m,elevation(iob),height_m,rawLevel,1)

          stalat = latitude(iob)
          stalon = longitude(iob)

!         Convert arrays for a single sounding
          a9time_ob_l = a9time_ob_r(iob)

          call convert_array(baromPressure(:,iob),pressure_mb,rawLevel
     1                      ,'none',r_missing_data,istatus)

          call convert_array(temperature(:,iob),temp_c,rawLevel
     1                      ,'k_to_c',r_missing_data,istatus)

          call convert_array(dewPt(:,iob),dewpoint_c,rawLevel
     1                      ,'k_to_c',r_missing_data,istatus)

          call convert_array(direction(:,iob),dir_deg,rawLevel
     1                      ,'none',r_missing_data,istatus)

          call convert_array(speed(:,iob),spd_mps,rawLevel
     1                      ,'none',r_missing_data,istatus)


          call get_nlevels_snd(pressure_mb,height_m,r_missing_data
     +                        ,rawLevel,nlevels_snd)

!         Apply QC editflag
          do i = 1,nlevels_snd
              if(editflag(i,iob) .eq. 3)then ! set wind to missing
                  direction(iob,3) = r_missing_data
                  speed(iob,3) = r_missing_data
              endif
          enddo ! i

          l_closest_time = .true.

          if(nlevels_snd .gt. 0 .and. l_closest_time)then
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
      subroutine read_rtamps_netcdf(nf_fid, manLevel, maxStaticIds, 
     +     nInventoryBins, rawLevel, recNum, stdLevel, termLevel, 
     +     tropLevel, editFlag, firstInBin, firstOverflow, 
     +     globalInventory, indxRefr, invTime, inventory, irMan, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, oiMan, 
     +     optIndxRefr, prevRecord, storedObs, absHumidity, 
     +     airDensity, baromPressure, bpMan, bpStd, bpTerm, bpTrop, 
     +     dewPt, direction, dpMan, dpStd, dpTrop, drMan, drStd, 
     +     elevation, geomHeight, geopHeight, ghMan, ghTerm, ghTrop, 
     +     gpStd, gpTerm, latitude, longitude, precipWater, 
     +     relHumidity, rhMan, rhStd, riseRate, shear, shearDir, 
     +     shearMagX, shearMagY, spMan, spStd, speed, temperature, 
     +     tpMan, tpStd, tpTrop, vaporPressure, velError, velSound, 
     +     observationTime, receivedTime, reportTime, dataProvider, 
     +     providerId, staticIds, stationName)
C
      include 'netcdf.inc'
      integer manLevel, maxStaticIds, nInventoryBins, rawLevel, 
     +     recNum, stdLevel, termLevel, tropLevel,nf_fid, nf_vid, 
     +     nf_status
      integer editFlag( rawLevel, recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, indxRefr( rawLevel,
     +     recNum), invTime(recNum), inventory(maxStaticIds), irMan(
     +     manLevel, recNum), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, oiMan( manLevel, recNum), optIndxRefr(
     +     rawLevel, recNum), prevRecord(recNum), storedObs(recNum)
      real absHumidity( rawLevel, recNum), airDensity( rawLevel,
     +     recNum), baromPressure( rawLevel, recNum), bpMan(
     +     manLevel, recNum), bpStd( stdLevel, recNum), bpTerm(
     +     termLevel, recNum), bpTrop( tropLevel, recNum), dewPt(
     +     rawLevel, recNum), direction( rawLevel, recNum), dpMan(
     +     manLevel, recNum), dpStd( stdLevel, recNum), dpTrop(
     +     tropLevel, recNum), drMan( manLevel, recNum), drStd(
     +     stdLevel, recNum), elevation(recNum), geomHeight(
     +     rawLevel, recNum), geopHeight( rawLevel, recNum), ghMan(
     +     manLevel, recNum), ghTerm( termLevel, recNum), ghTrop(
     +     tropLevel, recNum), gpStd( stdLevel, recNum), gpTerm(
     +     termLevel, recNum), latitude(recNum), longitude(recNum),
     +     precipWater( rawLevel, recNum), relHumidity( rawLevel,
     +     recNum), rhMan( manLevel, recNum), rhStd( stdLevel,
     +     recNum), riseRate( rawLevel, recNum), shear( rawLevel,
     +     recNum), shearDir( rawLevel, recNum), shearMagX( rawLevel,
     +     recNum), shearMagY( rawLevel, recNum), spMan( manLevel,
     +     recNum), spStd( stdLevel, recNum), speed( rawLevel,
     +     recNum), temperature( rawLevel, recNum), tpMan( manLevel,
     +     recNum), tpStd( stdLevel, recNum), tpTrop( tropLevel,
     +     recNum), vaporPressure( rawLevel, recNum), velError(
     +     rawLevel, recNum), velSound( rawLevel, recNum)
      double precision observationTime(recNum), receivedTime(recNum),
     +     reportTime(recNum)
      character*12 providerId(recNum)
      character*11 dataProvider(recNum)
      character*51 stationName(recNum)
      character*24 staticIds(maxStaticIds)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      absHumidity  "Absolute Humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'absHumidity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var absHumidity'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,absHumidity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var absHumidity'
      endif
C
C     Variable        NETCDF Long Name
C      airDensity   "Density of Air"
C
      nf_status=NF_INQ_VARID(nf_fid,'airDensity',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var airDensity'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,airDensity)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var airDensity'
      endif
C
C     Variable        NETCDF Long Name
C      baromPressure"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'baromPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var baromPressure'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,baromPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var baromPressure'
      endif
C
C     Variable        NETCDF Long Name
C      bpMan        "Pressure - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'bpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,bpMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpMan'
      endif
C
C     Variable        NETCDF Long Name
C      bpStd        "Pressure - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'bpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,bpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpStd'
      endif
C
C     Variable        NETCDF Long Name
C      bpTerm       "Pressure - Termination"
C
      nf_status=NF_INQ_VARID(nf_fid,'bpTerm',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTerm'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,bpTerm)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTerm'
      endif
C
C     Variable        NETCDF Long Name
C      bpTrop       "Pressure - Tropopause Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'bpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTrop'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,bpTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var bpTrop'
      endif
C
C     Variable        NETCDF Long Name
C      dewPt        "Dew Point Temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'dewPt',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewPt'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dewPt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dewPt'
      endif
C
C     Variable        NETCDF Long Name
C      direction    "Wind Direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'direction',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var direction'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,direction)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var direction'
      endif
C
C     Variable        NETCDF Long Name
C      dpMan        "Dew Point Temperature - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'dpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dpMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpMan'
      endif
C
C     Variable        NETCDF Long Name
C      dpStd        "Dew Point Temperature - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'dpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpStd'
      endif
C
C     Variable        NETCDF Long Name
C      dpTrop       "Dew Point Temperature - Tropopause Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'dpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpTrop'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,dpTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dpTrop'
      endif
C
C     Variable        NETCDF Long Name
C      drMan        "Wind Direction - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'drMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,drMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drMan'
      endif
C
C     Variable        NETCDF Long Name
C      drStd        "Wind Direction - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'drStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,drStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var drStd'
      endif
C
C     Variable        NETCDF Long Name
C      elevation    "Station Elevation"
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
C      geomHeight   "Geometric Height"
C
      nf_status=NF_INQ_VARID(nf_fid,'geomHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geomHeight'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,geomHeight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geomHeight'
      endif
C
C     Variable        NETCDF Long Name
C      geopHeight   "Geopotential Height"
C
      nf_status=NF_INQ_VARID(nf_fid,'geopHeight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geopHeight'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,geopHeight)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var geopHeight'
      endif
C
C     Variable        NETCDF Long Name
C      ghMan        "Geometric - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'ghMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ghMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghMan'
      endif
C
C     Variable        NETCDF Long Name
C      ghTerm       "Geometric - Termination"
C
      nf_status=NF_INQ_VARID(nf_fid,'ghTerm',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTerm'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ghTerm)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTerm'
      endif
C
C     Variable        NETCDF Long Name
C      ghTrop       "Geometric - Tropopause Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'ghTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTrop'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,ghTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ghTrop'
      endif
C
C     Variable        NETCDF Long Name
C      gpStd        "Geopotential - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'gpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,gpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpStd'
      endif
C
C     Variable        NETCDF Long Name
C      gpTerm       "Geopotential - Termination"
C
      nf_status=NF_INQ_VARID(nf_fid,'gpTerm',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpTerm'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,gpTerm)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var gpTerm'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "Station Latitude"
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
C      longitude    "Station Longitude"
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
C      precipWater  "Precipitable Water"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipWater',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipWater'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipWater)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var precipWater'
      endif
C
C     Variable        NETCDF Long Name
C      relHumidity  "Relative Humidity"
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
C      rhMan        "Relative Humidity - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'rhMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rhMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhMan'
      endif
C
C     Variable        NETCDF Long Name
C      rhStd        "Relative Humidity - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'rhStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,rhStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var rhStd'
      endif
C
C     Variable        NETCDF Long Name
C      riseRate     "Rise Rate"
C
      nf_status=NF_INQ_VARID(nf_fid,'riseRate',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var riseRate'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,riseRate)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var riseRate'
      endif
C
C     Variable        NETCDF Long Name
C      shear        "Shear"
C
      nf_status=NF_INQ_VARID(nf_fid,'shear',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shear'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,shear)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shear'
      endif
C
C     Variable        NETCDF Long Name
C      shearDir     "Shear Direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'shearDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearDir'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,shearDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearDir'
      endif
C
C     Variable        NETCDF Long Name
C      shearMagX    "Shear Magnitude X-direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'shearMagX',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagX'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,shearMagX)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagX'
      endif
C
C     Variable        NETCDF Long Name
C      shearMagY    "Shear Magnitude Y-direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'shearMagY',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagY'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,shearMagY)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var shearMagY'
      endif
C
C     Variable        NETCDF Long Name
C      spMan        "Wind Speed - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'spMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,spMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spMan'
      endif
C
C     Variable        NETCDF Long Name
C      spStd        "Wind Speed - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'spStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,spStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var spStd'
      endif
C
C     Variable        NETCDF Long Name
C      speed        "Wind Speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'speed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var speed'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,speed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var speed'
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
C      tpMan        "Temperature - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'tpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpMan'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tpMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpMan'
      endif
C
C     Variable        NETCDF Long Name
C      tpStd        "Temperature - Standard Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'tpStd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpStd'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tpStd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpStd'
      endif
C
C     Variable        NETCDF Long Name
C      tpTrop       "Temperature - Tropopause Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'tpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpTrop'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tpTrop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tpTrop'
      endif
C
C     Variable        NETCDF Long Name
C      vaporPressure"Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'vaporPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporPressure'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vaporPressure)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vaporPressure'
      endif
C
C     Variable        NETCDF Long Name
C      velError     "Velocity Error"
C
      nf_status=NF_INQ_VARID(nf_fid,'velError',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velError'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,velError)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velError'
      endif
C
C     Variable        NETCDF Long Name
C      velSound     "Velocity of Sound"
C
      nf_status=NF_INQ_VARID(nf_fid,'velSound',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velSound'
      endif
      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,velSound)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var velSound'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      editFlag     "Edit Flag"
C
      nf_status=NF_INQ_VARID(nf_fid,'editFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var editFlag'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,editFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var editFlag'
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
C      indxRefr     "Microwave Index of Refraction"
C
      nf_status=NF_INQ_VARID(nf_fid,'indxRefr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var indxRefr'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,indxRefr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var indxRefr'
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
C      irMan        "Microwave Index of Refraction - Mandatory level"
C
      nf_status=NF_INQ_VARID(nf_fid,'irMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var irMan'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,irMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var irMan'
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
C      oiMan        "Optical Index of Refraction - Mandatory Level"
C
      nf_status=NF_INQ_VARID(nf_fid,'oiMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var oiMan'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,oiMan)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var oiMan'
      endif
C
C     Variable        NETCDF Long Name
C      optIndxRefr  "Optical Index of Refraction"
C
      nf_status=NF_INQ_VARID(nf_fid,'optIndxRefr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var optIndxRefr'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,optIndxRefr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var optIndxRefr'
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
C      storedObs    "Stored \'Raw\' Profile Observations"
C
      nf_status=NF_INQ_VARID(nf_fid,'storedObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var storedObs'
      endif
      nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,storedObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var storedObs'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      observationTime"Observation Time"
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
C      receivedTime "Received Time"
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
C      reportTime   "Report Time"
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


C   Variables of type CHAR
C
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

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
