      subroutine get_poes_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer maxLevels, maxStaticIds, nInventoryBins, recNum,nf_fid,
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
C Get size of maxLevels
C
      nf_status=NF_INQ_DIMID(nf_fid,'maxLevels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxLevels'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,maxLevels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxLevels'
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

      write(6,*)' get_poes_data: number of records is ',recNum
      call read_poes_data(nf_fid, maxLevels, maxStaticIds,
     +     nInventoryBins, recNum, i4time_sys, ilaps_cycle_time,
     +     NX_L, NY_L, i4time_earliest, i4time_latest, lun_out,
     +     istatus)

      return
      end
C
C
      subroutine read_poes_data(nf_fid, maxLevels, maxStaticIds,
     +     nInventoryBins, recNum, i4time_sys, ilaps_cycle_time,
     +     NX_L, NY_L, i4time_earliest, i4time_latest, lun_out,
     +     istatus)


      include 'netcdf.inc'
      integer maxLevels, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
      integer cloudAmount(recNum), dayNight(recNum),
     +     fieldOfViewNum(recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     landSea(recNum), lastInBin(nInventoryBins),
     +     lastRecord(maxStaticIds), mixingRatioQCA( maxLevels,
     +     recNum), mixingRatioQCR( maxLevels, recNum), nStaticIds,
     +     numLevels(recNum), prevRecord(recNum), satProc(recNum),
     +     satelliteID(recNum), superAdiabatic(recNum),
     +     temperatureQCA( maxLevels, recNum), temperatureQCR(
     +     maxLevels, recNum), terrain(recNum)
      real cloudTopPressure(recNum), cloudTopTemperature(recNum),
     +     mixingRatio( maxLevels, recNum), precipWater(recNum),
     +     pressure( maxLevels, recNum), skinTemp(recNum),
     +     solarElev(recNum), staElev(recNum), staLat(recNum),
     +     staLon(recNum), staPress(recNum), temperature( maxLevels,
     +     recNum), zenithAngle(recNum)
      double precision validTime(recNum)
      character temperatureDD( maxLevels, recNum)
      character*8 staticIds(maxStaticIds)
      character*7 staName(recNum)
      character mixingRatioDD( maxLevels, recNum)

!     Declarations for 'write_snd' call
      integer iwmostanum(recNum)
      real stalatl(maxLevels),stalonl(maxLevels)
      character a9time_ob_r(recNum)*9,a9time_ob_l(maxLevels)*9
      character staname_o(recNum)*5
      character c8_obstype*8
      real height_m(maxLevels)
      real pressure_mb(maxLevels)
      real temp_c(maxLevels)
      real dewpoint_c(maxLevels)
      real dir_deg(maxLevels)
      real spd_mps(maxLevels)

      integer iob_tot
      save iob_tot
      data iob_tot /0/

      logical l_closest_time, l_closest_time_i, l_in_domain
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

      call read_poes_netcdf(nf_fid, maxLevels, maxStaticIds, 
     +     nInventoryBins, recNum, cloudAmount, dayNight, 
     +     fieldOfViewNum, firstInBin, firstOverflow, 
     +     globalInventory, invTime, inventory, isOverflow, landSea, 
     +     lastInBin, lastRecord, mixingRatioQCA, mixingRatioQCR, 
     +     nStaticIds, numLevels, prevRecord, satProc, satelliteID, 
     +     superAdiabatic, temperatureQCA, temperatureQCR, terrain, 
     +     cloudTopPressure, cloudTopTemperature, mixingRatio, 
     +     precipWater, pressure, skinTemp, solarElev, staElev, 
     +     staLat, staLon, staPress, temperature, zenithAngle, 
     +     mixingRatioDD, staName, staticIds, temperatureDD, 
     +     validTime)
C
C The netcdf variables are filled - your snd write call may go here
C
!     Initial loop through obs to get times and stanums
      do iob = 1,recNum
          iwmostanum(iob) = 0
          if(abs(validTime(iob)) .le. 1e10)then
              i4time_ob = idint(validTime(iob))+315619200
              call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
          endif

!         Create station name from ob count (try hex later if more needed)
          iob_tot = iob_tot + 1
          write(staname_o(iob),44) iob_tot
 44       format(i5.5)

      enddo ! iob

      c8_obstype = 'POESSND '

      do iob = 1,recNum
          height_m = r_missing_data
          stalatl = staLat(iob)
          stalonl = staLon(iob)

!         Convert arrays for a single sounding
          a9time_ob_l = a9time_ob_r(iob)

          call convert_array(pressure(:,iob),pressure_mb,maxLevels
     1                      ,'pa_to_mb',r_missing_data,istatus)

          call convert_array(temperature(:,iob),temp_c,maxLevels
     1                      ,'k_to_c',r_missing_data,istatus)

          dewpoint_c = r_missing_data
          do ilvl = 1,maxLevels
              if(        mixingRatio(ilvl,iob) .gt. 0.
     1             .and. mixingRatio(ilvl,iob) .le. 1.)then
                  dewpoint_c(ilvl) =
     1              tmr(mixingRatio(ilvl,iob)*1000.,pressure_mb(ilvl))
              endif
          enddo

          dir_deg = r_missing_data

          spd_mps = r_missing_data


          call get_nlevels_snd(pressure_mb,height_m,r_missing_data
     +                        ,maxLevels,nlevels_snd)

          l_closest_time = .true.

          if(   stalatl(1) .le. rnorth .and. stalatl(1) .ge. south
     1    .and. stalonl(1) .le. east   .and. stalonl(1) .ge. west
     1                                                       )then
              l_in_domain = .true.
          else
              l_in_domain = .false.
          endif

          if(nlevels_snd .gt. 0 .and. l_closest_time
     1                           .and. l_in_domain)then
!             call 'write_snd' for a single profile
              call open_ext(lun_out,i4time_sys,'snd',istatus)

              call write_snd(lun_out
     +                      ,1,nlevels_snd,1
     +                      ,iwmostanum
     +                      ,stalatl,stalonl,staElev(iob)
     +                      ,staName_o(iob)
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
      subroutine read_poes_netcdf(nf_fid, maxLevels, maxStaticIds, 
     +     nInventoryBins, recNum, cloudAmount, dayNight, 
     +     fieldOfViewNum, firstInBin, firstOverflow, 
     +     globalInventory, invTime, inventory, isOverflow, landSea, 
     +     lastInBin, lastRecord, mixingRatioQCA, mixingRatioQCR, 
     +     nStaticIds, numLevels, prevRecord, satProc, satelliteID, 
     +     superAdiabatic, temperatureQCA, temperatureQCR, terrain, 
     +     cloudTopPressure, cloudTopTemperature, mixingRatio, 
     +     precipWater, pressure, skinTemp, solarElev, staElev, 
     +     staLat, staLon, staPress, temperature, zenithAngle, 
     +     mixingRatioDD, staName, staticIds, temperatureDD, 
     +     validTime)
C
      include 'netcdf.inc'
      integer maxLevels, maxStaticIds, nInventoryBins, recNum,nf_fid, 
     +     nf_vid, nf_status
      integer cloudAmount(recNum), dayNight(recNum),
     +     fieldOfViewNum(recNum), firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     landSea(recNum), lastInBin(nInventoryBins),
     +     lastRecord(maxStaticIds), mixingRatioQCA( maxLevels,
     +     recNum), mixingRatioQCR( maxLevels, recNum), nStaticIds,
     +     numLevels(recNum), prevRecord(recNum), satProc(recNum),
     +     satelliteID(recNum), superAdiabatic(recNum),
     +     temperatureQCA( maxLevels, recNum), temperatureQCR(
     +     maxLevels, recNum), terrain(recNum)
      real cloudTopPressure(recNum), cloudTopTemperature(recNum),
     +     mixingRatio( maxLevels, recNum), precipWater(recNum),
     +     pressure( maxLevels, recNum), skinTemp(recNum),
     +     solarElev(recNum), staElev(recNum), staLat(recNum),
     +     staLon(recNum), staPress(recNum), temperature( maxLevels,
     +     recNum), zenithAngle(recNum)
      double precision validTime(recNum)
      character temperatureDD( maxLevels, recNum)
      character*8 staticIds(maxStaticIds)
      character*7 staName(recNum)
      character mixingRatioDD( maxLevels, recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     cloudTopPressure"Cloud Top Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudTopPressure',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for cloudTopPressure'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,cloudTopPressure)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for cloudTopPressure'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     cloudTopTemperature"Cloud Top Temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudTopTemperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for cloudTopTemperature'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,cloudTopTemperature)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for cloudTopTemperature'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     mixingRatio   "Mixing Ratio"
C
      nf_status=NF_INQ_VARID(nf_fid,'mixingRatio',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mixingRatio'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,mixingRatio)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mixingRatio'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     precipWater   "Total Precipitable Water"
C
      nf_status=NF_INQ_VARID(nf_fid,'precipWater',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for precipWater'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,precipWater)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for precipWater'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     pressure      "Pressure"
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
C     skinTemp      "Skin Temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'skinTemp',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for skinTemp'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,skinTemp)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for skinTemp'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     solarElev     "Grid Point Solar Elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'solarElev',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for solarElev'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,solarElev)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for solarElev'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staElev       "Grid Point Elevation"
C
      nf_status=NF_INQ_VARID(nf_fid,'staElev',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staElev'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,staElev)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staElev'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staLat        "Grid Point Latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'staLat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staLat'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,staLat)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staLat'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staLon        "Grid Point Longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'staLon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staLon'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,staLon)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staLon'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staPress      "Grid Point Pressure"
C
      nf_status=NF_INQ_VARID(nf_fid,'staPress',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staPress'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,staPress)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staPress'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     temperature   "Retrieved Temperature"
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
C     zenithAngle   "Satellite Zenith Angle"
C
      nf_status=NF_INQ_VARID(nf_fid,'zenithAngle',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for zenithAngle'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,zenithAngle)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for zenithAngle'
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     cloudAmount   "Cloud Amount"
C
      nf_status=NF_INQ_VARID(nf_fid,'cloudAmount',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for cloudAmount'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,cloudAmount)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for cloudAmount'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     dayNight      "Day/Night Qualifier"
C
      nf_status=NF_INQ_VARID(nf_fid,'dayNight',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for dayNight'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,dayNight)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for dayNight'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     fieldOfViewNum"Field of View Number"
C
      nf_status=NF_INQ_VARID(nf_fid,'fieldOfViewNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for fieldOfViewNum'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,fieldOfViewNum)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for fieldOfViewNum'
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
C     landSea       "land/sea mask"
C
      nf_status=NF_INQ_VARID(nf_fid,'landSea',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for landSea'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,landSea)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for landSea'
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
C     mixingRatioQCA"mixingRatio QC applied word"
C
      nf_status=NF_INQ_VARID(nf_fid,'mixingRatioQCA',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mixingRatioQCA'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,mixingRatioQCA)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mixingRatioQCA'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     mixingRatioQCR"mixingRatio QC results word"
C
      nf_status=NF_INQ_VARID(nf_fid,'mixingRatioQCR',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mixingRatioQCR'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,mixingRatioQCR)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mixingRatioQCR'
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
C     numLevels     "Number of Sounding Levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'numLevels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numLevels'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numLevels)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numLevels'
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
C     satProc       "Satellite Processing Technique Used"
C
      nf_status=NF_INQ_VARID(nf_fid,'satProc',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for satProc'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,satProc)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for satProc'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     satelliteID   "Satellite Identifier"
C
      nf_status=NF_INQ_VARID(nf_fid,'satelliteID',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for satelliteID'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,satelliteID)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for satelliteID'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     superAdiabatic"Superadiabatic Indicator"
C
      nf_status=NF_INQ_VARID(nf_fid,'superAdiabatic',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for superAdiabatic'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,superAdiabatic)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for superAdiabatic'
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
C     terrain       "Terrain Type"
C
      nf_status=NF_INQ_VARID(nf_fid,'terrain',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for terrain'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,terrain)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for terrain'
       endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     validTime     "Sounding Valid Time"
C
      nf_status=NF_INQ_VARID(nf_fid,'validTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for validTime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,validTime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for validTime'
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     mixingRatioDD "mixingRatio QC summary value"
C
      nf_status=NF_INQ_VARID(nf_fid,'mixingRatioDD',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for mixingRatioDD'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,mixingRatioDD)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for mixingRatioDD'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     staName       "Grid Point Station Identifier"
C
      nf_status=NF_INQ_VARID(nf_fid,'staName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for staName'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,staName)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for staName'
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

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
