
      subroutine get_pirep_data(i4time_sys,ilaps_cycle_time,filename
     1                                                     ,ext
     1                                                     ,NX_L,NY_L
     1                                                     ,istatus)

      character*170 filename
      character*(*) ext

!.............................................................................

      include 'netcdf.inc'
      integer maxIcingLvls, maxLocs, maxSkyLvls, maxTurbElements,
     +     maxTurbLvls, recNum,nf_fid, nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename
      endif
C
C  Fill all dimension values
C
C
C Get size of maxIcingLvls
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxIcingLvls',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxIcingLvls'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxIcingLvls)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxIcingLvls'
      endif
C
C Get size of maxLocs
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxLocs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxLocs'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxLocs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxLocs'
      endif
C
C Get size of maxSkyLvls
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxSkyLvls',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxSkyLvls'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxSkyLvls)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxSkyLvls'
      endif
C
C Get size of maxTurbElements
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxTurbElements',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxTurbElements'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxTurbElements)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxTurbElements'
      endif
C
C Get size of maxTurbLvls
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxTurbLvls',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxTurbLvls'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxTurbLvls)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxTurbLvls'
      endif
C
C Get size of recNum
C
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
      endif
      call pireps_sub(nf_fid, maxIcingLvls, maxLocs, maxSkyLvls,
     +     maxTurbElements, maxTurbLvls, recNum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,NX_L,NY_L,istatus)
      return
!.............................................................................
      end
C
C
      subroutine pireps_sub(nf_fid, maxIcingLvls, maxLocs, maxSkyLvls,
     +     maxTurbElements, maxTurbLvls, recNum,
!.............................................................................
     1     ext,i4time_sys,ilaps_cycle_time,NX_L,NY_L,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer maxIcingLvls, maxLocs, maxSkyLvls, maxTurbElements,
     +     maxTurbLvls, recNum,nf_fid, nf_vid, nf_status
      integer icingIntens( maxIcingLvls, recNum),
     +     lowLvlWndShr(recNum), temperature(recNum), turbIntens(
     +     maxTurbElements,  maxTurbLvls, recNum), vis(recNum),
     +     windDir(recNum), windSpd(recNum)
      real fltLvlBottom(recNum), fltLvlTop(recNum), icingBottom(
     +     maxIcingLvls, recNum), icingTop( maxIcingLvls, recNum),
     +     lat( maxLocs, recNum), lon( maxLocs, recNum),
     +     skyCvrBottom( maxSkyLvls, recNum), skyCvrTop( maxSkyLvls,
     +     recNum), turbBottom( maxTurbLvls, recNum), turbTop(
     +     maxTurbLvls, recNum)
      double precision recptTime(recNum), timeObs(recNum)
      character*4 icingHtInd( maxIcingLvls, recNum)
      character*8 skyCvrAmt( maxSkyLvls, recNum)
      character*4 turbHtInd( maxTurbLvls, recNum)
      character*16 locStr( maxLocs, recNum)
      character*4 skyCvrHtInd( maxSkyLvls, recNum)
      character*4 turbFreq( maxTurbElements,  maxTurbLvls, recNum)
      character*33 fltWthr(recNum)
      character*5 icingType( maxIcingLvls, recNum)
      character*201 origReport(recNum)
      character*5 aircraftType(recNum)
      character*4 reportType(recNum)
      character*4 fltLvlInd(recNum)
      character*4 collSite(recNum)
      character*5 turbType( maxTurbLvls, recNum)

!.............................................................................

      character*170 filename
      character*(*) ext
      character*9 a9_timeObs,a9_recptTime 
      character*7 c7_skycover
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)

!.............................................................................

      call read_pireps_netcdf(nf_fid, maxIcingLvls, maxLocs, maxSkyLvls,     
     +     maxTurbElements, maxTurbLvls, recNum, icingIntens, 
     +     lowLvlWndShr, temperature, turbIntens, vis, windDir, 
     +     windSpd, fltLvlBottom, fltLvlTop, icingBottom, icingTop, 
     +     lat, lon, skyCvrBottom, skyCvrTop, turbBottom, turbTop, 
     +     recptTime, timeObs, aircraftType, collSite, fltLvlInd, 
     +     fltWthr, icingHtInd, icingType, locStr, origReport, 
     +     reportType, skyCvrAmt, skyCvrHtInd, turbFreq, turbHtInd, 
     +     turbType)
C
C The netcdf variables are filled - your code goes here
C
!.............................................................................

      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

!     Write All Pireps to LAPS PIN file
      altitude = 0.

      r_nc_missing_data = 1e20

      write(6,*)' recNum = ',recNum

      num_pireps = recNum

      do i = 1,num_pireps

          write(6,*)
          write(6,*)' Pirep #',i

          if(lat(1,i) .ge. r_nc_missing_data)then
              write(6,*)' Missing first latitude',i
              goto 999
          endif
          if(lon(1,i) .ge. r_nc_missing_data)then
              write(6,*)' Missing first longitude',i
              goto 999
          endif
          if(timeObs(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing ob time',i
              goto 999
          endif
          if(recptTime(i) .ge. r_nc_missing_data)then
              write(6,*)' Missing received time',i
              goto 999
          endif

!         Test to see how many lat/lons are present
          n_latitude_present = 0
          do j = 1,4
            if(lat(j,i) .lt. r_nc_missing_data)then
              write(6,*)' lat/lon = ',lat(j,i),lon(j,i)
              n_latitude_present = n_latitude_present + 1
            endif
          enddo ! j

          write(6,*)' Num locations = ',n_latitude_present       

          if(n_latitude_present .gt. 1)then
              write(6,*)' Multiple locations, reject ob',i
              goto 999
          endif

          if(lat(1,i) .le. rnorth .and. lat(1,i) .ge. south .and.
     1       lon(1,i) .ge. west   .and. lon(1,i) .le. east      )then        
              write(6,*)' Pirep is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' Outside domain lat/lon perimeter - reject'
              goto 999
          endif

!         Write Individual Pirep to LAPS PIN file
          call c_time2fname(nint(timeObs(i)),a9_timeObs)
          call c_time2fname(nint(recptTime(i)),a9_recptTime)

          call cv_asc_i4time(a9_timeObs,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then
              write(6,*)' Outside time window - reject '
     1                              ,a9_timeObs,i4_resid
              goto 999        
          endif

          call open_ext(31,i4time_sys,ext(1:3),istatus)

          write(6,1)a9_timeObs,a9_recptTime 
          write(31,1)a9_timeObs,a9_recptTime 
 1        format(' Time - prp/rcvd:'/1x,a9,2x,a9) 

          write(6,2)lat(1,i),lon(1,i),altitude
          write(31,2)lat(1,i),lon(1,i),altitude
 2        format(' Lat, lon, altitude'/f8.3,f10.3,f8.0)  

          write(6,33)
          write(31,33)
 33       format(' Cloud layer')

!         Write out cloud base/top in feet and cloud amount in eighths
          do ilyr = 1,3

!             Use -1000. for missing value of cloud base/top
              if(skyCvrBottom(ilyr,i) .ge. 1e10)then
                  rbase = -1000.
              else
                  rbase = skyCvrBottom(ilyr,i) * 3.281
              endif
              if(skyCvrTop(ilyr,i) .ge. 1e10)then
                  rtop = -1000.
              else
                  rtop = skyCvrTop(ilyr,i) * 3.281
              endif

              c7_skycover(1:7) = skyCvrAmt(ilyr,i)(1:7)

              call skycover_to_frac(c7_skycover,fraction,istatus)

!             Test for missing or unusable cloud cover
              if(istatus .ne. 1)then 
                  ieighths = -999
              else 
                  ieighths = nint(fraction * 8.0)
              endif                  

              write(6,3)rbase,rtop,ieighths,skyCvrAmt(ilyr,i)
     1                                     ,fraction
 3            format(' Cloud layer',2f8.0,i5,1x,a8,f5.2)

!             if(istatus .eq. 1)then
!                 write(6,*)' Above layer written to PIN file'
                  write(31,3)rbase,rtop,ieighths
!             endif

          enddo ! ilyr

 999      continue

      enddo ! i

!.............................................................................

      return
      end
C
C  Subroutine to read the file "PIREP data : selected by receipt time : time range from 885926400 to 885926700" 
C
      subroutine read_pireps_netcdf(nf_fid, maxIcingLvls, maxLocs, 
     +     maxSkyLvls, maxTurbElements, maxTurbLvls, recNum, 
     +     icingIntens, lowLvlWndShr, temperature, turbIntens, vis, 
     +     windDir, windSpd, fltLvlBottom, fltLvlTop, icingBottom, 
     +     icingTop, lat, lon, skyCvrBottom, skyCvrTop, turbBottom, 
     +     turbTop, recptTime, timeObs, aircraftType, collSite, 
     +     fltLvlInd, fltWthr, icingHtInd, icingType, locStr, 
     +     origReport, reportType, skyCvrAmt, skyCvrHtInd, turbFreq, 
     +     turbHtInd, turbType)
C
      include 'netcdf.inc'
      integer maxIcingLvls, maxLocs, maxSkyLvls, maxTurbElements, 
     +     maxTurbLvls, recNum,nf_fid, nf_vid, nf_status
      integer icingIntens( maxIcingLvls, recNum),
     +     lowLvlWndShr(recNum), temperature(recNum), turbIntens(
     +     maxTurbElements,  maxTurbLvls, recNum), vis(recNum),
     +     windDir(recNum), windSpd(recNum)
      real fltLvlBottom(recNum), fltLvlTop(recNum), icingBottom(
     +     maxIcingLvls, recNum), icingTop( maxIcingLvls, recNum),
     +     lat( maxLocs, recNum), lon( maxLocs, recNum),
     +     skyCvrBottom( maxSkyLvls, recNum), skyCvrTop( maxSkyLvls,
     +     recNum), turbBottom( maxTurbLvls, recNum), turbTop(
     +     maxTurbLvls, recNum)
      double precision recptTime(recNum), timeObs(recNum)
      character*4 icingHtInd( maxIcingLvls, recNum)
      character*8 skyCvrAmt( maxSkyLvls, recNum)
      character*4 turbHtInd( maxTurbLvls, recNum)
      character*16 locStr( maxLocs, recNum)
      character*4 skyCvrHtInd( maxSkyLvls, recNum)
      character*33 fltWthr(recNum)
      character*4 turbFreq( maxTurbElements,  maxTurbLvls, recNum)
      character*5 icingType( maxIcingLvls, recNum)
      character*201 origReport(recNum)
      character*5 aircraftType(recNum)
      character*4 reportType(recNum)
      character*4 fltLvlInd(recNum)
      character*4 collSite(recNum)
      character*5 turbType( maxTurbLvls, recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      fltLvlBottom "Bottom of PIREPs flight level"
C
        nf_status = NF_INQ_VARID(nf_fid,'fltLvlBottom',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltLvlBottom'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,fltLvlBottom)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltLvlBottom'
      endif
C
C     Variable        NETCDF Long Name
C      fltLvlTop    "Top of PIREPs flight level"
C
        nf_status = NF_INQ_VARID(nf_fid,'fltLvlTop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltLvlTop'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,fltLvlTop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltLvlTop'
      endif
C
C     Variable        NETCDF Long Name
C      icingBottom  "Bottom of Icing Layer"
C
        nf_status = NF_INQ_VARID(nf_fid,'icingBottom',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingBottom'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,icingBottom)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingBottom'
      endif
C
C     Variable        NETCDF Long Name
C      icingTop     "Top of Icing Layer"
C
        nf_status = NF_INQ_VARID(nf_fid,'icingTop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingTop'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,icingTop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingTop'
      endif
C
C     Variable        NETCDF Long Name
C      lat          "Latitude of report"
C
        nf_status = NF_INQ_VARID(nf_fid,'lat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lat'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,lat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lat'
      endif
C
C     Variable        NETCDF Long Name
C      lon          "Longitude of report"
C
        nf_status = NF_INQ_VARID(nf_fid,'lon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lon'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,lon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lon'
      endif
C
C     Variable        NETCDF Long Name
C      skyCvrBottom "Bottom of Sky Cover Layer"
C
        nf_status = NF_INQ_VARID(nf_fid,'skyCvrBottom',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrBottom'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,skyCvrBottom)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrBottom'
      endif
C
C     Variable        NETCDF Long Name
C      skyCvrTop    "Bottom of Sky Cover Layer"
C
        nf_status = NF_INQ_VARID(nf_fid,'skyCvrTop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrTop'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,skyCvrTop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrTop'
      endif
C
C     Variable        NETCDF Long Name
C      turbBottom   "Bottom of Turbulence Layer"
C
        nf_status = NF_INQ_VARID(nf_fid,'turbBottom',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbBottom'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,turbBottom)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbBottom'
      endif
C
C     Variable        NETCDF Long Name
C      turbTop      "Top of Turbulence Layer"
C
        nf_status = NF_INQ_VARID(nf_fid,'turbTop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbTop'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,turbTop)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbTop'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      icingIntens  "Icing intensity"
C
        nf_status = NF_INQ_VARID(nf_fid,'icingIntens',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingIntens'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,icingIntens)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingIntens'
      endif
C
C     Variable        NETCDF Long Name
C      lowLvlWndShr "Low Level Wind Shear Indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'lowLvlWndShr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lowLvlWndShr'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,lowLvlWndShr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lowLvlWndShr'
      endif
C
C     Variable        NETCDF Long Name
C      temperature  "Temperature"
C
        nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      turbIntens   "Turbulence Intensity"
C
        nf_status = NF_INQ_VARID(nf_fid,'turbIntens',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbIntens'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,turbIntens)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbIntens'
      endif
C
C     Variable        NETCDF Long Name
C      vis          "Visibility"
C
        nf_status = NF_INQ_VARID(nf_fid,'vis',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vis'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,vis)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vis'
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "Wind Direction"
C
        nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,windDir)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      endif
C
C     Variable        NETCDF Long Name
C      windSpd      "Wind Speed"
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpd'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,windSpd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpd'
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      recptTime    "Time report was received"
C
        nf_status = NF_INQ_VARID(nf_fid,'recptTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var recptTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,recptTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var recptTime'
      endif
C
C     Variable        NETCDF Long Name
C      timeObs      "Time of report"
C
        nf_status = NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var timeObs'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      aircraftType "Type of Aircraft issuing report"
C
        nf_status = NF_INQ_VARID(nf_fid,'aircraftType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var aircraftType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,aircraftType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var aircraftType'
      endif
C
C     Variable        NETCDF Long Name
C      collSite     "Data Collection Site"
C
        nf_status = NF_INQ_VARID(nf_fid,'collSite',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var collSite'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,collSite)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var collSite'
      endif
C
C     Variable        NETCDF Long Name
C      fltLvlInd    "Flight Level Indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'fltLvlInd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltLvlInd'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,fltLvlInd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltLvlInd'
      endif
C
C     Variable        NETCDF Long Name
C      fltWthr      "Flight Weather"
C
        nf_status = NF_INQ_VARID(nf_fid,'fltWthr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltWthr'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,fltWthr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var fltWthr'
      endif
C
C     Variable        NETCDF Long Name
C      icingHtInd   "Icing Height Indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'icingHtInd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingHtInd'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,icingHtInd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingHtInd'
      endif
C
C     Variable        NETCDF Long Name
C      icingType    "Type of Icing observed"
C
        nf_status = NF_INQ_VARID(nf_fid,'icingType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,icingType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var icingType'
      endif
C
C     Variable        NETCDF Long Name
C      locStr       "Location String for each derived lat/lon"
C
        nf_status = NF_INQ_VARID(nf_fid,'locStr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var locStr'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,locStr)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var locStr'
      endif
C
C     Variable        NETCDF Long Name
C      origReport   "Original PIREPs report"
C
        nf_status = NF_INQ_VARID(nf_fid,'origReport',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origReport'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,origReport)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origReport'
      endif
C
C     Variable        NETCDF Long Name
C      reportType   "Report Type"
C
        nf_status = NF_INQ_VARID(nf_fid,'reportType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,reportType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportType'
      endif
C
C     Variable        NETCDF Long Name
C      skyCvrAmt    "Amount of Sky Cover"
C
        nf_status = NF_INQ_VARID(nf_fid,'skyCvrAmt',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrAmt'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,skyCvrAmt)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrAmt'
      endif
C
C     Variable        NETCDF Long Name
C      skyCvrHtInd  "Sky Cover Height Indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'skyCvrHtInd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrHtInd'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,skyCvrHtInd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var skyCvrHtInd'
      endif
C
C     Variable        NETCDF Long Name
C      turbFreq     "Weather qualifier"
C
        nf_status = NF_INQ_VARID(nf_fid,'turbFreq',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbFreq'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,turbFreq)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbFreq'
      endif
C
C     Variable        NETCDF Long Name
C      turbHtInd    "Turbulence Height Indicator"
C
        nf_status = NF_INQ_VARID(nf_fid,'turbHtInd',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbHtInd'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,turbHtInd)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbHtInd'
      endif
C
C     Variable        NETCDF Long Name
C      turbType     "Type of turbulence observed"
C
        nf_status = NF_INQ_VARID(nf_fid,'turbType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,turbType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var turbType'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
