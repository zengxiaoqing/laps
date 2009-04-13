      subroutine get_drpsnd_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_drpsnd_earliest,i4time_drpsnd_latest       
     +                   ,filename,isnd_staname
     +                   ,lun_out,l_fill_ht
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      logical l_fill_ht

      integer manLevel, recNum, sigTLevel, sigWLevel,
     +     tropLevel,nf_fid, nf_vid, nf_status, nlvl_out
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
C Get size of sigTLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'sigTLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigTLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,sigTLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigTLevel'
      endif
C
C Get size of sigWLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'sigWLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigWLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,sigWLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigWLevel'
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
      nlvl_out = manLevel + sigTLevel + sigWLevel
      call read_drpsnd_data(nf_fid, manLevel, recNum, sigTLevel,
     +     sigWLevel, tropLevel, i4time_sys, ilaps_cycle_time, NX_L,
     +     NY_L, i4time_drpsnd_earliest, i4time_drpsnd_latest, 
     +     nlvl_out, lun_out, l_fill_ht, isnd_staname, istatus)

      return
      end
C
C
      subroutine read_drpsnd_data(nf_fid, manLevel, recNum, sigTLevel,
     +     sigWLevel, tropLevel, i4time_sys, ilaps_cycle_time, NX_L,
     +     NY_L, i4time_drpsnd_earliest, i4time_drpsnd_latest, nlvl_out,
     +     lun_out, l_fill_ht, isnd_staname, istatus)


      include 'netcdf.inc'
      integer manLevel, recNum, sigTLevel, sigWLevel,
     +     tropLevel,nf_fid, nf_vid, nf_status
      integer marsdenSquare(recNum), numMand(recNum), numSigT(recNum),
     +     numSigW(recNum), numTrop(recNum), prMan( manLevel,
     +     recNum), prTrop( tropLevel, recNum), wdMan( manLevel,
     +     recNum), wdSigW( sigWLevel, recNum), wdTrop( tropLevel,
     +     recNum)
      real htMan( manLevel, recNum), latitude(recNum),
     +     longitude(recNum), prSigT( sigTLevel, recNum), prSigW(
     +     sigWLevel, recNum), tdMan( manLevel, recNum), tdSigT(
     +     sigTLevel, recNum), tdTrop( tropLevel, recNum), tpMan(
     +     manLevel, recNum), tpSigT( sigTLevel, recNum), tpTrop(
     +     tropLevel, recNum), wsMan( manLevel, recNum), wsSigW(
     +     sigWLevel, recNum), wsTrop( tropLevel, recNum)
      double precision timeNominal(recNum), timeObs(recNum)
      character*12 DropSondeLocation(recNum)
      character*2048 rawDropSonde(recNum)

      integer wmostanum(recNum)
      real staelev(recNum)
      CHARACTER*1 staName(6,recNum)
      character*6 c6_staname
      real wdMan_r(manLevel,recNum)
      real prMan_r(manLevel,recNum)
      real wdSigW_r(sigWLevel,recNum)
      real htSigW(sigWLevel, recNum)

      logical l_closest_time, l_closest_time_i, l_in_domain, l_fill_ht       
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)
      character*8 c8_obstype
      character*9 a9time_sys,a9time_release,a9time_syn,a9time_drpsnd       

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

      call read_drpsnd_netcdf(nf_fid, manLevel, recNum, sigTLevel, 
     +     sigWLevel, tropLevel, marsdenSquare, numMand, numSigT, 
     +     numSigW, numTrop, prMan, prTrop, wdMan, wdSigW, wdTrop, 
     +     htMan, latitude, longitude, prSigT, prSigW, tdMan, tdSigT, 
     +     tdTrop, tpMan, tpSigT, tpTrop, wsMan, wsSigW, wsTrop, 
     +     DropSondeLocation, rawDropSonde, timeNominal, timeObs)
C
C The netcdf variables are filled - your snd write call may go here
C
!     Write All Dropsondes to LAPS SND file

      r_nc_missing_data = 1e20

      htSigW = r_missing_data

      n_snd = recNum

      do isnd = 1,n_snd

!         QC and write out the sounding
          i4time_drpsnd = 0

          if(abs(timeObs(isnd)) .lt. 1e10 .and.
     1       abs(timeObs(isnd)) .gt.    0.      )then
              i4time_release = idint(timeObs(isnd))+315619200

              i4time_diff = i4time_release - i4time_sys

              if(abs(i4time_diff) .gt. 20000)then
                  write(6,*)' Warning: i4time_release is not '
     1                     ,'consistent with i4time_diff'
     1                     ,i4time_release,i4time_sys
              endif

!             Correction for balloon fall time to mid-sounding
              i4time_drpsnd = i4time_release + 100

          else
              i4time_release = 0
              i4time_diff = 0

          endif

          write(6,*)
          write(6,*)' Drpsnd #',isnd,i4time_sys,i4time_release
     1                         ,i4time_diff,i4time_syn

          call make_fnam_lp(i4time_sys    , a9time_sys    , istatus)
          call make_fnam_lp(i4time_release, a9time_release, istatus)
          call make_fnam_lp(i4time_syn    , a9time_syn    , istatus)
          call make_fnam_lp(i4time_drpsnd , a9time_drpsnd , istatus)

          write(6,*)' times - sys/release/syn/drpsnd: '
     1             ,a9time_sys,' ',a9time_release,' '
     1             ,a9time_syn,' ',a9time_drpsnd

          if(latitude(isnd) .ge. r_nc_missing_data)then
              write(6,*)' Missing first latitude',isnd
              goto 999
          endif

          if(longitude(isnd) .ge. r_nc_missing_data)then
              write(6,*)' Missing first longitude',isnd
              goto 999
          endif

          if(latitude(isnd) .le. rnorth .and. latitude(isnd) .ge. south 
     1                                .AND.      
     1       longitude(isnd) .ge. west   .and. longitude(isnd) .le. east      
     1                                                            )then       

!         if(.true.)then      ! for testing

              write(6,*)' Drpsnd is inside domain lat/lon perimeter'
          else
              write(6,*)
     1            ' Outside domain lat/lon perimeter - reject'
              goto 999
          endif

          if(i4time_drpsnd .ne. 0)then ! test window
              if(i4time_drpsnd .ge. i4time_drpsnd_earliest .and.
     1           i4time_drpsnd .le. i4time_drpsnd_latest)then
                  write(6,*)' Inside time window'
              else
                  write(6,*)' Outside time window - reject'
                  goto 999
              endif
          endif

          staelev(isnd) = -999.
          c8_obstype = 'DROPSND'

          isnd_staname = isnd_staname + 1
          write(c6_staname,1)isnd_staname
 1        format('D',i4.4,1x)
          do i = 1,6
              staname(i,isnd) = c6_staname(i:i)
          enddo 

!         Convert arrays from integer to real
          do ilev = 1,nummand(isnd) 
              wdMan_r(ilev,isnd) = wdMan(ilev,isnd)
              prMan_r(ilev,isnd) = prMan(ilev,isnd)
          enddo ! l

          do ilev = 1,numsigw(isnd) 
              wdSigW_r(ilev,isnd) = wdSigW(ilev,isnd)
          enddo ! l

          call sort_and_write(i4time_sys,lun_out,l_fill_ht
     1                       ,recNum,isnd,r_missing_data,a9time_drpsnd
     1                       ,c8_obstype
     1                       ,wmostanum,staname,latitude,longitude
     1                       ,staelev
     1                       ,nummand,htman,prman_r,tpman,tdman      
     1                       ,wdMan_r,wsMan
     1                       ,numsigt,prsigt,tpsigt,tdsigt
     1                       ,numsigw,prsigw,htsigw,wdSigW_r,wssigw
     1                       ,nlvl_out 
     1                       ,manLevel,sigTLevel,sigWLevel,istatus)

          go to 999

 998      write(6,*)' Error writing out DRPSND'

 999      continue

      enddo ! i

      return
      end
C
C  Subroutine to read the file "Drop Sonde data : selected by ob time : time range from 1201190400 to 1201201200" 
C
      subroutine read_drpsnd_netcdf(nf_fid, manLevel, recNum, 
     +     sigTLevel, sigWLevel, tropLevel, marsdenSquare, numMand, 
     +     numSigT, numSigW, numTrop, prMan, prTrop, wdMan, wdSigW, 
     +     wdTrop, htMan, latitude, longitude, prSigT, prSigW, tdMan, 
     +     tdSigT, tdTrop, tpMan, tpSigT, tpTrop, wsMan, wsSigW, 
     +     wsTrop, DropSondeLocation, rawDropSonde, timeNominal, 
     +     timeObs)
C
      include 'netcdf.inc'
      integer manLevel, recNum, sigTLevel, sigWLevel, 
     +     tropLevel,nf_fid, nf_vid, nf_status
      integer marsdenSquare(recNum), numMand(recNum), numSigT(recNum),
     +     numSigW(recNum), numTrop(recNum), prMan( manLevel,
     +     recNum), prTrop( tropLevel, recNum), wdMan( manLevel,
     +     recNum), wdSigW( sigWLevel, recNum), wdTrop( tropLevel,
     +     recNum)
      real htMan( manLevel, recNum), latitude(recNum),
     +     longitude(recNum), prSigT( sigTLevel, recNum), prSigW(
     +     sigWLevel, recNum), tdMan( manLevel, recNum), tdSigT(
     +     sigTLevel, recNum), tdTrop( tropLevel, recNum), tpMan(
     +     manLevel, recNum), tpSigT( sigTLevel, recNum), tpTrop(
     +     tropLevel, recNum), wsMan( manLevel, recNum), wsSigW(
     +     sigWLevel, recNum), wsTrop( tropLevel, recNum)
      double precision timeNominal(recNum), timeObs(recNum)
      character*12 DropSondeLocation(recNum)
      character*2048 rawDropSonde(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     htMan         "Geopotential - Mandatory level"
C
      nf_status=NF_INQ_VARID(nf_fid,'htMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for htMan'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,htMan)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for htMan'
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
C     prSigT        "Pressure - Significant level wrt T"
C
      nf_status=NF_INQ_VARID(nf_fid,'prSigT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for prSigT'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,prSigT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for prSigT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     prSigW        "Pressure - Significant level wrt W"
C
      nf_status=NF_INQ_VARID(nf_fid,'prSigW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for prSigW'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,prSigW)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for prSigW'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tdMan         "Dew Point Depression - Mandatory level"
C
      nf_status=NF_INQ_VARID(nf_fid,'tdMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tdMan'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tdMan)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tdMan'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tdSigT        "Dew Point Depression - Significant level wrt T"
C
      nf_status=NF_INQ_VARID(nf_fid,'tdSigT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tdSigT'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tdSigT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tdSigT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tdTrop        "Dew Point Depression - Tropopause levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'tdTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tdTrop'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tdTrop)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tdTrop'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tpMan         "Temperature - Mandatory level"
C
      nf_status=NF_INQ_VARID(nf_fid,'tpMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tpMan'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tpMan)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tpMan'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tpSigT        "Temperature - Significant level wrt T"
C
      nf_status=NF_INQ_VARID(nf_fid,'tpSigT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tpSigT'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tpSigT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tpSigT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     tpTrop        "Temperature - Tropopause level"
C
      nf_status=NF_INQ_VARID(nf_fid,'tpTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for tpTrop'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,tpTrop)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for tpTrop'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wsMan         "wind speed"
C
      nf_status=NF_INQ_VARID(nf_fid,'wsMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wsMan'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wsMan)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wsMan'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wsSigW        "wind speed - Significant Level wrt W"
C
      nf_status=NF_INQ_VARID(nf_fid,'wsSigW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wsSigW'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wsSigW)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wsSigW'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wsTrop        "wind speed - Tropopause levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'wsTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wsTrop'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,wsTrop)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wsTrop'
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     marsdenSquare "Marsden Square"
C
      nf_status=NF_INQ_VARID(nf_fid,'marsdenSquare',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for marsdenSquare'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,marsdenSquare)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for marsdenSquare'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numMand       "Number of Mandatory Levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'numMand',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numMand'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numMand)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numMand'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numSigT       "Number of Significant Levels wrt T"
C
      nf_status=NF_INQ_VARID(nf_fid,'numSigT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numSigT'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numSigT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numSigT'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numSigW       "Number of Significant Levels wrt W"
C
      nf_status=NF_INQ_VARID(nf_fid,'numSigW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numSigW'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numSigW)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numSigW'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     numTrop       "Number of Tropopause Levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'numTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for numTrop'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,numTrop)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for numTrop'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     prMan         "Pressure - Mandatory level"
C
      nf_status=NF_INQ_VARID(nf_fid,'prMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for prMan'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,prMan)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for prMan'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     prTrop        "Pressure - Tropopause level"
C
      nf_status=NF_INQ_VARID(nf_fid,'prTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for prTrop'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,prTrop)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for prTrop'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wdMan         "wind direction"
C
      nf_status=NF_INQ_VARID(nf_fid,'wdMan',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wdMan'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wdMan)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wdMan'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wdSigW        "wind direction - Significant Level wrt W"
C
      nf_status=NF_INQ_VARID(nf_fid,'wdSigW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wdSigW'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wdSigW)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wdSigW'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     wdTrop        "wind direction - Tropopause levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'wdTrop',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for wdTrop'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,wdTrop)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for wdTrop'
       endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     timeNominal   "drop sonde data hour"
C
      nf_status=NF_INQ_VARID(nf_fid,'timeNominal',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for timeNominal'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeNominal)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for timeNominal'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     timeObs       "Observation Time"
C
      nf_status=NF_INQ_VARID(nf_fid,'timeObs',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for timeObs'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,timeObs)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for timeObs'
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     DropSondeLocation"DropSonde Location or Observation Type (from 62626)"
C
      nf_status=NF_INQ_VARID(nf_fid,'DropSondeLocation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for DropSondeLocation'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,DropSondeLocation)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for DropSondeLocation'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     rawDropSonde  "raw Drop Sonde ASCII message"
C
      nf_status=NF_INQ_VARID(nf_fid,'rawDropSonde',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for rawDropSonde'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,rawDropSonde)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for rawDropSonde'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
