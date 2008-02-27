      subroutine get_raob_data_a(nx,ny,i4time_sys,i4time_raob_earliest,
     1                           i4time_raob_latest,filename,maxRaob,
     1                           maxM,maxT,maxW,lat,lon,timeSyn,timeRel,
     1                           numSigT, numSigW, wmoStaNum,staname,
     1                           typeW,typeT,prSigT, tSigT,tdSigT,
     1                           htSigT, htSigW, wdSigW,wsSigW, 
     1                           staLat, staLon, staElev,
     1                           max_ht_m_proc,min_pres_mb_proc,numRaob, 
     1                           n_raobs_avail,verif_missing_data, 
     1                           raob_missing_data,istatus)

      implicit none
      include 'netcdf.inc'

!     Input Variables
      integer     nx,ny,i4time_sys,i4time_raob_earliest,
     1              i4time_raob_latest,maxRaob,maxM,maxT,maxW
      character*(*) filename
      real        lat(nx,ny),lon(nx,ny)

!     Variables filled for output
      integer     timeSyn(maxRaob),timeRel(maxRaob),
     1              numSigT(maxRaob),numSigW(maxRaob),
     1              istatus
      integer       numRaob, n_raobs_avail
      real	    verif_missing_data,raob_missing_data
      character*6   staName(maxRaob)
      integer	    wmoStaNum(maxRaob),wmoStaNum_f
      character*1   typeW(maxW,maxRaob),typeT(maxT,maxRaob)
      real        prSigT(maxT,maxRaob),tSigT(maxT,maxRaob),
     1              tdSigT(maxT,maxRaob),
     1              htSigT(maxT,maxRaob),
     1              htSigW(maxW,maxRaob),wdSigW(maxW,maxRaob),
     1              wsSigW(maxW,maxRaob),staElev(maxRaob), 
     1              max_ht_m_proc, min_pres_mb_proc,
     1              staLat(maxRaob), staLon(maxRaob) 

!     Other Variables
      real        htMan(maxM,maxRaob), prMan(maxM,maxRaob),
     1              tMan(maxM,maxRaob), tdMan(maxM,maxRaob),
     1              wsMan(maxM,maxRaob), wdMan(maxM,maxRaob),
     1              prSigTI(maxT,maxRaob),tSigTI(maxT,maxRaob),
     1              tdSigTI(maxT,maxRaob),
     1              htSigWI(maxW,maxRaob),wdSigWI(maxW,maxRaob),
     1              wsSigWI(maxW,maxRaob)

      double precision d_timeRel, d_timeSyn
      character*9   a9time_syn, a9time_release, a9time_raob, 
     1              a9time_sys
      character*8   c8_project
      character*6   staName_f
      integer     i, j, recNum, sigTLevel, sigWLevel, staNameLen, 
     1              staNameLen_f, nf_fid, nf_vid, nf_status,
     1              nObs, i4time_syn, i4time_release, i4time_diff,
     1              i4time_raob, manLevel, numMan(maxRaob),
     1              itime_delay,status
      real        staLat_f, staLon_f, staElev_f, ri, rj,
     1              r_nc_missing_data
      integer       index_1(1), start(2), count(2)

!.............................................................................
! BEGIN

      n_raobs_avail = 0

      istatus = 1    ! assume a good return
      staNameLen = len(staName(1))

C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        write(6,*) NF_STRERROR(nf_status)
        write(6,*)'NF_OPEN ',filename
        istatus = 0
        return
      else
        write(6,*) 'Opening file ',filename
      endif
C
C  Read dimension values recNum, manLevel, sigTLevel, sigWLevel, staNameLen
C
C
C Get size of recNum
C
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim recNum'
        istatus = 0
        return
      endif
C
C Get size of manLevel
C
      nf_status = NF_INQ_DIMID(nf_fid,'manLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim manLevel'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,manLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim manLevel'
        istatus = 0
        return
      endif
C
C Get size of sigTLevel
C
      nf_status = NF_INQ_DIMID(nf_fid,'sigTLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigTLevel'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,sigTLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigTLevel'
        istatus = 0
        return
      endif
C
C Get size of sigWLevel
C
      nf_status = NF_INQ_DIMID(nf_fid,'sigWLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigWLevel'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,sigWLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim sigWLevel'
        istatus = 0
        return
      endif
C
C Get size of staNameLen
C
      nf_status = NF_INQ_DIMID(nf_fid,'staNameLen',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim staNameLen'
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,staNameLen_f)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim staNameLen'
        istatus = 0
        return
      endif

      if (staNameLen_f .gt. staNameLen) then
        print *, 'staName truncated to ',staNameLen,' characters.'
        staNameLen_f = staNameLen
      endif
C
C     Read missing data value for staLat and staLon
C
        nf_status = NF_INQ_VARID(nf_fid,'staLat',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staLat'
          istatus = 0
          return
        endif
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',
     1                              r_nc_missing_data)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in att for staLat'
          istatus = 0
          return
        endif

! Read data from netCDF file
      nObs = 0   !number of RAOBs stored to returning arrays
      start(1) = 1
      count(2) = 1

      do i = 1, recNum
        index_1(1) = i
        start(2) = i
C
C       Read staName      "Station Identifier"
C
        nf_status = NF_INQ_VARID(nf_fid,'staName',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staName'
        endif

        count(1) = staNameLen_f
        nf_status = NF_GET_VARA_TEXT(nf_fid,nf_vid,start,
     1                               count,staName_f)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staName'
        endif
C
C       Read wmoStaNum    "WMO Station Number"
C
        nf_status = NF_INQ_VARID(nf_fid,'wmoStaNum',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wmoStaNum'
        endif
        nf_status = NF_GET_VAR1_INT(nf_fid,nf_vid,index_1,wmoStaNum_f)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wmoStaNum'
        endif
C
C       Read staLat       "Station Latitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'staLat',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staLat'
        endif
        nf_status = NF_GET_VAR1_REAL(nf_fid,nf_vid,index_1,staLat_f)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staLat'
        endif
C
C       Read staLon       "Station Longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'staLon',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staLon'
        endif
        nf_status = NF_GET_VAR1_REAL(nf_fid,nf_vid,index_1,staLon_f)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staLon'
        endif

        if(staLat_f .ge. r_nc_missing_data)then
          goto 888
        endif

        if(stalon_f .ge. r_nc_missing_data)then
          goto 888
        endif

        call latlon_to_rlapsgrid(staLat_f,staLon_f,lat,lon,nx,ny,
     1                           ri,rj,status)
        if (status .ne. 1) then   !raob is not in laps domain
          goto 888
        else
          write(6,*)
          write(6,*) 'Raob ',wmoStaNum_f,' is in Laps domain',
     1               staLat_f,' ',staLon_f
          n_raobs_avail = n_raobs_avail + 1
        endif

!       Read synTime and relTime and see if raob is in time window
C
C       Read relTime      "Sounding Release Time"

        nf_status = NF_INQ_VARID(nf_fid,'relTime',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var relTime'
        endif
        nf_status = NF_GET_VAR1_DOUBLE(nf_fid,nf_vid,index_1,d_timeRel)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var relTime'
        endif
C
C       Read synTime      "Synoptic Time"
C
        nf_status = NF_INQ_VARID(nf_fid,'synTime',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var synTime'
        endif
        nf_status = NF_GET_VAR1_DOUBLE(nf_fid,nf_vid,index_1,d_timeSyn)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var synTime'
        endif

        i4time_raob = 0

        if(abs(d_timeSyn) .lt. 1e10)then
          i4time_syn  = idint(d_timeSyn)+315619200
          i4time_raob = i4time_syn
        else
          i4time_syn = 0
        endif
c
c Somewhat bogus but allows synch of delays with time of data
c
        call get_c8_project(c8_project,istatus)
        call upcase(c8_project,c8_project)
        if(c8_project.eq.'AIRDROP')itime_delay=0  !2*3600

        if(abs(d_timeRel) .lt. 1e10)then
          i4time_release = idint(d_timeRel)+315619200

          i4time_diff = i4time_release - i4time_sys

          if(abs(i4time_diff) .gt. 20000)then
            write(6,*)' Warning: i4time_release is not '
     1               ,'consistent with i4time_diff'
     1               ,i4time_release,i4time_sys
            call make_fnam_lp(i4time_sys,a9time_sys,istatus)
!            if (a9time_sys(6:9) .eq. '0000') then
!              i4time_release = i4time_release - 210600
!            else
!              i4time_release = i4time_release - 166200
!            endif
          endif

!         Correction for balloon rise time to mid-troposphere
!         and time delay due to cron and sched inputs.
          i4time_raob = i4time_release + 1800 + itime_delay
        else
          i4time_release = 0
          i4time_diff = 0
        endif

        if(i4time_raob .ne. 0)then ! test window
          if(i4time_raob .ge. i4time_raob_earliest .and.
     1       i4time_raob .le. i4time_raob_latest)then
c.or.(abs(i4time_raob-i4time_sys).gt.10800))then
          else   !Outside time window - reject
            write(6,*) 'Raob ',wmoStaNum_f,' is outside of time window '
            write(6,*) 'Raob:',i4time_raob,' window:',
     1                 i4time_raob_earliest, '-',i4time_raob_latest
            n_raobs_avail = n_raobs_avail - 1
            goto 888
          endif
        elseif (i4time_raob .eq. 0) then  !not a valid time
          goto 888
        endif
        i4time_diff = abs(i4time_raob-i4time_sys)
        if(i4time_diff.gt.10800)then
           print*,'Raob not valid verification for this cycle'
           print*,'i4time_raob/i4time_sys: ',i4time_raob,i4time_sys
           n_raobs_avail = n_raobs_avail - 1
           goto 888
        endif
! If you get to here, RAOB is in domain and within time window, save it
! Write what have to return arrays: staName, staLat, staLon, timeSyn, timeRel

        nObs = nObs + 1
        staName(nObs) = staName_f
        wmoStaNum(nObs) = wmoStaNum_f
        staLat(nObs) = staLat_f
        staLon(nObs) = staLon_f
        timeSyn(nObs) = i4time_syn
        timeRel(nObs) = i4time_release

! read staElev, numMan, numSigT, numSigW, prSigT, tpSigT, tdSigT, htSigW, 
!               wdSigW, wsSigW, prMan, htMan, tMan, tdMan, wsMan, wdMan 
!   into return arrays
C
C       Read staElev      "Station Elevation"
C
        nf_status = NF_INQ_VARID(nf_fid,'staElev',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staElev'
        endif
        nf_status = NF_GET_VAR1_REAL(nf_fid,nf_vid,index_1,
     1                               staElev(nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var staElev'
        endif
C
C       Read numMan      "Number of Mandatory levels"
C
        nf_status = NF_INQ_VARID(nf_fid,'numMand',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var numMand'
        endif
        nf_status = NF_GET_VAR1_INT(nf_fid,nf_vid,index_1,numMan(nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var numMan'
        endif

!       Read in prMan and tpMan
        count(1) = numMan(nObs)
C
C       Read prMan       "Pressure - Mandatory level"
C
        nf_status = NF_INQ_VARID(nf_fid,'prMan',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var prMan' 
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               prMan(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var prMan'
        endif
C
C       Read missing for prMan...use as missing for rest of data too
C
        nf_status = NF_GET_ATT_REAL(nf_fid,nf_vid,'_FillValue',
     1                              raob_missing_data)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in _FillValue for prMan'
          istatus = 0
          return
        endif
C
C       Read tpMan       "Temperature - Mandatory level"
C
        nf_status = NF_INQ_VARID(nf_fid,'tpMan',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tpMan'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               tMan(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tpMan'
        endif
C
C       Read tdMan       "Dewpoint - Mandatory level"
C
        nf_status = NF_INQ_VARID(nf_fid,'tdMan',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tcMan'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               tdMan(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tdMan'
        endif
C
C       Read htMan       "Geopotential - Mandatory level"
C
        nf_status = NF_INQ_VARID(nf_fid,'htMan',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var htMan'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               htMan(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var htMan'
        endif
C
C       Read wdMan       "Wind Direction - Mandatory level"
C
        nf_status = NF_INQ_VARID(nf_fid,'wdMan',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wdMan'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               wdMan(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wdMan'
        endif
C
C       Read wsMan       "Wind Speed - Mandatory level"
C
        nf_status = NF_INQ_VARID(nf_fid,'wsMan',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wsMan'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               wsMan(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wsMan'
        endif

C
C       Read numSigT      "Number of Significant Levels wrt T"
C
        nf_status = NF_INQ_VARID(nf_fid,'numSigT',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var numSigT'
        endif
        nf_status = NF_GET_VAR1_INT(nf_fid,nf_vid,index_1,numSigT(nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var numSigT'
        endif

!       Read in prSigT and tpSigT
        count(1) = numSigT(nObs)
C
C       Read prSigT       "Pressure - Significant level wrt T"
C
        nf_status = NF_INQ_VARID(nf_fid,'prSigT',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var prSigT' 
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               prSigTI(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var prSigT'
        endif
C
C       Read tpSigT       "Temperature - Significant level wrt T"
C
        nf_status = NF_INQ_VARID(nf_fid,'tpSigT',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tpSigT'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               tSigTI(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tpSigT'
        endif

C
C       Read tdSigT       "Dewpoint - Significant level wrt T"
C
        nf_status = NF_INQ_VARID(nf_fid,'tdSigT',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tdSigT'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               tdSigTI(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var tdSigT'
        endif

C
C       Read numSigW      "Number of Significant Levels wrt W"
C
        nf_status = NF_INQ_VARID(nf_fid,'numSigW',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var numSigW'
        endif
        nf_status = NF_GET_VAR1_INT(nf_fid,nf_vid,index_1,numSigW(nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var numSigW'
        endif

!       Read in htSigW, wdSigW and wsSigW
        count(1) = numSigW(nObs)
C
C       Read htSigW       "Geopotential - Significant level wrt W"
C
        nf_status = NF_INQ_VARID(nf_fid,'htSigW',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var htSigW'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               htSigWI(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var htSigW'
        endif
C
C       Read wdSigW       "Wind Direction - Significant level wrt W"
C
        nf_status = NF_INQ_VARID(nf_fid,'wdSigW',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wdSigW'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               wdSigWI(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wdSigW'
        endif
C
C       Read wsSigW       "Wind Speed - Significant level wrt W"
C
        nf_status = NF_INQ_VARID(nf_fid,'wsSigW',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wsSigW'
        endif
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,start,count,
     1                               wsSigWI(1,nObs))
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status),'in var wsSigW'
        endif

888     continue
      enddo
      numRaob = nObs
C
C     Close netCDF file
C
      nf_status = NF_CLOSE(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_CLOSE ',filename
      endif

C     convert dewpoint depression to dewpoint
      do i = 1, numRaob
        do j = 1, numMan(i)
          if ((tdMan(j,i) .ne. raob_missing_data) .and.
     1        (tMan(j,i) .ne. raob_missing_data))
     1      tdMan(j,i) = tMan(j,i) - tdMan(j,i)
        enddo
 
        do j = 1, numSigT(i)
          if ((tdSigTI(j,i) .ne. raob_missing_data) .and.
     1        (tSigTI(j,i) .ne. raob_missing_data))
     1      tdSigTI(j,i) = tSigTI(j,i) - tdSigTI(j,i)
        enddo
      enddo

          

C
C     interleave Man obs into SigT and SigW obs
C
      call mingleManSig(maxM,maxT,maxW,maxRaob,numRaob,
     1                  max_ht_m_proc, min_pres_mb_proc,
     1                  numSigT,numSigW,numMan,typeW,typeT,
     1                  htMan,prMan,tMan,tdMan,wsMan,wdMan,
     1                  staElev,raob_missing_data, 
     1                  prSigTI,tSigTI,tdSigTI,htSigWI,
     1                  wsSigWI,wdSigWI,prSigT,tSigT,tdSigT,
     1                  htSigT,htSigW,wsSigW,wdSigW,
     1                  verif_missing_data,istatus)

      if (istatus .ne. 1) then
        write(6,*) 'Error mingling Man and Sig obs '
        istatus = 0
      endif

      return
      end
!1............................................................................
      subroutine mingleManSig(maxM,maxT,maxW,maxRaob,numRaob,
     1                  max_ht_m_proc, min_pres_mb_proc,
     1                  numSigT,numSigW,numMan,typeW,typeT,
     1                  htMan,prMan,tMan,tdMan,wsMan,wdMan,
     1                  staElev,raob_missing_data,
     1                  prSigTI,tSigTI,tdSigTI,htSigWI,
     1                  wsSigWI,wdSigWI,prSigT,tSigT,tdSigT,
     1                  htSigT,htSigW,wsSigW,wdSigW,
     1                  verif_missing_data,istatus)

      implicit none

      integer     maxM,maxT,maxW,maxRaob,numRaob,
     1              numSigT(maxRaob),numSigW(maxRaob),
     1              numMan(maxRaob)
      real	    max_ht_m_proc, min_pres_mb_proc
      character*1   typeW(maxW,maxRaob),typeT(maxT,maxRaob)
      real        htMan(maxM,maxRaob), prMan(maxM,maxRaob),
     1              tMan(maxM,maxRaob), tdMan(maxM,maxRaob),
     1              wsMan(maxM,maxRaob), wdMan(maxM,maxRaob),
     1              staElev(maxRaob), raob_missing_data,
     1              prSigTI(maxT,maxRaob),tSigTI(maxT,maxRaob),
     1              tdSigTI(maxT,maxRaob),
     1              htSigWI(maxW,maxRaob),wdSigWI(maxW,maxRaob),
     1              wsSigWI(maxW,maxRaob)
      real        prSigT(maxT,maxRaob),tSigT(maxT,maxRaob),
     1              tdSigT(maxT,maxRaob),
     1              htSigT(maxT,maxRaob),
     1              htSigW(maxW,maxRaob),wdSigW(maxW,maxRaob),
     1              wsSigW(maxW,maxRaob), htSfc, prSfc,
     1              verif_missing_data
      integer	    istatus

      integer     mPtr, sPtr, jPtr, i, j

C     BEGIN

C     loop through raobs to interleave Man with SigT and Man with SigW
      do i = 1, numRaob

C       DEBUG  print out data
        if(.true.)then

        write(6,*) 'Raob ',i
        write(6,*) 'htMan, wsMan, wdMan'
        do j = 1, numMan(i)
          write(6,*) j,htMan(j,i),wsMan(j,i),wdMan(j,i)
        enddo

        write(6,*) 'htSigWI, wsSigWI, wdSigWI'
        do j = 1, numSigW(i)
          write(6,*) j,htSigWI(j,i),wsSigWI(j,i),wdSigWI(j,i)
        enddo

        write(6,*) 'prMan, tMan, tdMan'
        do j = 1, numMan(i)
          write(6,*) j,prMan(j,i),tMan(j,i),tdMan(j,i)
        enddo

        write(6,*) 'prSigTI, tSigTI, tdSigTI'
        do j = 1, numSigT(i)
          write(6,*) j,prSigTI(j,i),tSigTI(j,i),tdSigTI(j,i)
        enddo
        write(6,*)

        endif

C       get surface ht and pr setup
        jPtr = 1
        mPtr = 1
        sPtr = 1

C       see if first level of Man is missing
        if (htMan(1,i) .eq. raob_missing_data) then  !see if it's in htSigW
          if ((htSigWI(1,i) .eq. raob_missing_data) .or.
     1        ((htSigWI(1,i) .eq. 0.0) .and.
     1         (wsSigWI(1,i) .eq. raob_missing_data) .and.
     1         (wdSigWI(1,i) .eq. raob_missing_data))) then  !skip surface level
            htSfc = raob_missing_data 
            mPtr = 2
            sPtr = 2
          else
            htSfc = htSigWI(1,i)
            htSigW(jPtr,i) = htSigWI(sPtr,i)
            wsSigW(jPtr,i) = wsSigWI(sPtr,i)
            wdSigW(jPtr,i) = wdSigWI(sPtr,i)
            typeW(jPtr,i) = 'W'
            mPtr = mPtr + 1  !skip sfc in Man
            sPtr = sPtr + 1
            jPtr = jPtr + 1
          endif
        else
          htSfc = htMan(1,i)
          htSigW(jPtr,i) = htMan(mPtr,i)
          wsSigW(jPtr,i) = wsMan(mPtr,i)
          wdSigW(jPtr,i) = wdMan(mPtr,i)
          typeW(jPtr,i) = 'M'
          mPtr = mPtr + 1
          sPtr = sPtr + 1  !skip sfc in sigW
          jPtr = jPtr + 1
        endif

C       set mPtr and sPtr so ht .gt. htSfc 
        do while ((mPtr .le. numMan(i)) .and.
     1            (htMan(mPtr,i) .le. htSfc))
          mPtr = mPtr + 1
        enddo
        do while ((sPtr .le. numSigW(i)) .and.
     1            (htSigWI(sPtr,i) .le. htSfc))
          sPtr = sPtr + 1
        enddo

C       mingle Man with SigW - Go until SigW levels exhausted
        do while (sPtr .le. numSigW(i))  

C         make sure mPtr is past any missing htMan
          do while ((mPtr .le. numMan(i)) .and.
     1          ((htMan(mPtr,i) .gt. max_ht_m_proc)
     1      .or. (htMan(mPtr,i) .eq. raob_missing_data)
     1      .or. (wsMan(mPtr,i) .eq. raob_missing_data) 
     1      .or. (wdMan(mPtr,i) .eq. raob_missing_data)))
            mPtr = mPtr + 1
          enddo
       
C         make sure sPtr is past any missing htSigW
          do while ((sPtr .le. numSigW(i)) .and.
     1          ((htSigWI(sPtr,i) .gt. max_ht_m_proc)
     1      .or. (htSigWI(sPtr,i) .eq. raob_missing_data)
     1      .or. (wsSigWI(mPtr,i) .eq. raob_missing_data) 
     1      .or. (wdSigWI(mPtr,i) .eq. raob_missing_data)))
            sPtr = sPtr + 1
          enddo

          if (sPtr .gt. numSigW(i)) goto 777  !no sigW obs to mingle
       
          if (mPtr .le. numMan(i)) then
            if (htSigWI(sPtr,i) .lt. htMan(mPtr,i)) then
              htSigW(jPtr,i) = htSigWI(sPtr,i)
              wsSigW(jPtr,i) = wsSigWI(sPtr,i)
              wdSigW(jPtr,i) = wdSigWI(sPtr,i)
              typeW(jPtr,i) = 'W'
              sPtr = sPtr + 1
              jPtr = jPtr + 1
            else 
              if (htMan(mPtr,i) .lt. htSigWI(sPtr,i)) then
                htSigW(jPtr,i) = htMan(mPtr,i)
                wsSigW(jPtr,i) = wsMan(mPtr,i)
                wdSigW(jPtr,i) = wdMan(mPtr,i)
                typeW(jPtr,i) = 'M'
                mPtr = mPtr + 1
                jPtr = jPtr + 1
              else
                if (htMan(mPtr,i) .eq. htSigWI(sPtr,i)) then
                  htSigW(jPtr,i) = htMan(mPtr,i)
                  wsSigW(jPtr,i) = wsMan(mPtr,i)
                  wdSigW(jPtr,i) = wdMan(mPtr,i)
                  typeW(jPtr,i) = 'M'
                  mPtr = mPtr + 1
                  sPtr = sPtr + 1
                  jPtr = jPtr + 1
                endif
              endif
            endif
          else
            if (htSigW(jPtr,i) .le. max_ht_m_proc) then
              htSigW(jPtr,i) = htSigWI(sPtr,i)
              wsSigW(jPtr,i) = wsSigWI(sPtr,i)
              wdSigW(jPtr,i) = wdSigWI(sPtr,i)
              typeW(jPtr,i) = 'W'
              sPtr = sPtr + 1
              jPtr = jPtr + 1
            endif
          endif

777       continue 

        enddo

C       Finish rest of man if any are left
        do while (mPtr .le. numMan(i)) 
          if ((htMan(mPtr,i) .ne. raob_missing_data) 
     1      .and. (htMan(mPtr,i) .le. max_ht_m_proc) 
     1      .and. (wsMan(mPtr,i) .ne. raob_missing_data) 
     1      .and. (wdMan(mPtr,i) .ne. raob_missing_data))
     1      then
            htSigW(jPtr,i) = htMan(mPtr,i)
            wsSigW(jPtr,i) = wsMan(mPtr,i)
            wdSigW(jPtr,i) = wdMan(mPtr,i)
            typeW(jPtr,i) = 'M'
            jPtr = jPtr + 1
          endif
          mPtr = mPtr + 1
        enddo

C       set numSigW to mingled value
        numSigW(i) = jPtr - 1

        write(6,*) numSigW(i)
        write(6,*) 'htSigW, typeW, wsSigW, wdSigW '        
        do j = 1, numSigW(i)
          write(6,*) j, htSigW(j,i), typeW(j,i), wsSigW(j,i), 
     1               wdSigW(j,i)         
        enddo

C       mingle Man with SigT
C       get surface pr setup
        jPtr = 1
        mPtr = 1
        sPtr = 1

C       see if first level of Man is missing
        if (prMan(1,i) .eq. raob_missing_data) then  !see if it's in prSigT
          if ((prSigTI(1,i) .eq. raob_missing_data) .or.
     1        (tSigTI(1,i) .eq. raob_missing_data))then  !skip surface level
            prSfc = raob_missing_data 
            mPtr = 2
            sPtr = 2
          else
            prSfc = prSigTI(1,i)
            prSigT(jPtr,i) = prSigTI(sPtr,i)
            tSigT(jPtr,i) = tSigTI(sPtr,i)
            tdSigT(jPtr,i) = tdSigTI(sPtr,i)
            htSigT(jPtr,i) = verif_missing_data
            typeT(jPtr,i) = 'T'
            mPtr = mPtr + 1  !skip sfc in Man
            sPtr = sPtr + 1
            jPtr = jPtr + 1
          endif
        else
          prSfc = prMan(1,i)
          prSigT(jPtr,i) = prMan(mPtr,i)
          tSigT(jPtr,i) = tMan(mPtr,i)
          tdSigT(jPtr,i) = tdMan(mPtr,i)
          htSigT(jPtr,i) = htMan(mPtr,i)
          typeT(jPtr,i) = 'M'
          mPtr = mPtr + 1
          sPtr = sPtr + 1  !skip sfc in sigT
          jPtr = jPtr + 1
        endif

C       set mPtr and sPtr so pr .lt. prSfc
        do while ((mPtr .le. numMan(i)) .and.
     1            (prMan(mPtr,i) .ge. prSfc))
          mPtr = mPtr + 1
        enddo
        do while ((sPtr .le. numSigT(i)) .and.
     1            (prSigTI(sPtr,i) .ge. prSfc))
          sPtr = sPtr + 1
        enddo

C       mingle Man with SigT - Go until SigT levels exhausted
        do while (sPtr .le. numSigT(i))  

C         make sure mPtr is past any missing prMan
          do while ((mPtr .le. numMan(i)) .and.
     1          ((prMan(mPtr,i) .lt. min_pres_mb_proc)
     1      .or. (prMan(mPtr,i) .eq. raob_missing_data)
     1      .or. (tMan(mPtr,i) .eq. raob_missing_data)))
            mPtr = mPtr + 1
          enddo
       
C         make sure sPtr is past any missing prSigT
          do while ((sPtr .le. numSigT(i)) .and.
     1         ((prSigTI(sPtr,i) .lt. min_pres_mb_proc)
     1      .or.(prSigTI(sPtr,i) .eq. raob_missing_data)
     1      .or.(tSigTI(mPtr,i) .eq. raob_missing_data)))
            sPtr = sPtr + 1
          enddo

          if (sPtr .gt. numSigT(i)) goto 778

          if (mPtr .le. numMan(i)) then
            if (prSigTI(sPtr,i) .gt. prMan(mPtr,i)) then
              prSigT(jPtr,i) = prSigTI(sPtr,i)
              tSigT(jPtr,i) = tSigTI(sPtr,i)
              tdSigT(jPtr,i) = tdSigTI(sPtr,i)
              htSigT(jPtr,i) = verif_missing_data
              typeT(jPtr,i) = 'T'
              sPtr = sPtr + 1
              jPtr = jPtr + 1
            else 
              if (prMan(mPtr,i) .gt. prSigTI(sPtr,i)) then
                prSigT(jPtr,i) = prMan(mPtr,i)
                tSigT(jPtr,i) = tMan(mPtr,i)
                tdSigT(jPtr,i) = tdMan(mPtr,i)
                htSigT(jPtr,i) = htMan(mPtr,i)
                typeT(jPtr,i) = 'M'
                mPtr = mPtr + 1
                jPtr = jPtr + 1
              else
                if (prMan(mPtr,i) .eq. prSigTI(sPtr,i)) then
                  prSigT(jPtr,i) = prMan(mPtr,i)
                  tSigT(jPtr,i) = tMan(mPtr,i)
                  tdSigT(jPtr,i) = tdMan(mPtr,i)
                  htSigT(jPtr,i) = htMan(mPtr,i)
                  typeT(jPtr,i) = 'M'
                  mPtr = mPtr + 1
                  sPtr = sPtr + 1
                  jPtr = jPtr + 1
                endif
              endif
            endif
          else
            prSigT(jPtr,i) = prSigTI(sPtr,i)
            tSigT(jPtr,i) = tSigTI(sPtr,i)
            tdSigT(jPtr,i) = tdSigTI(sPtr,i)
            htSigT(jPtr,i) = verif_missing_data
            typeT(jPtr,i) = 'T'
            sPtr = sPtr + 1
            jPtr = jPtr + 1
          endif

778       continue

        enddo

C       Finish rest of man if any are left
        if (mPtr .le. numMan(i)) then
          do while (mPtr .le. numMan(i)) 
            if ((prMan(mPtr,i) .ne. raob_missing_data).and. 
     1          (tMan(mPtr,i) .ne. raob_missing_data)) then
              prSigT(jPtr,i) = prMan(mPtr,i)
              tSigT(jPtr,i) = tMan(mPtr,i)
              tdSigT(jPtr,i) = tdMan(mPtr,i)
              htSigT(jPtr,i) = htMan(mPtr,i)
              typeT(jPtr,i) = 'M'
              jPtr = jPtr + 1
            endif
            mPtr = mPtr + 1
          enddo
        endif

C       set numSigT to mingled value
        numSigT(i) = jPtr - 1

        write(6,*) 'Raob ',i,' after mingle'
        write(6,*) numSigT(i)
        write(6,*) 'prSigT, typeT tSigT, tdSigT htSigT'        
        do j = 1, numSigT(i)
          write(6,*) j,prSigT(j,i),typeT(j,i),tSigT(j,i), 
     1               tdSigT(j,i), htSigT(j,i) 
        enddo

      enddo

      return
      end
!2............................................................................
      subroutine calc_domain_perim(nx, ny, lat, lon, north, south,
     1                             east,west,r_buffer)

      integer     nx, ny, i, j
      real        north, south, east, west, r_buffer,
     1              lat(nx,ny), lon(nx,ny)

!     calculate domain perimeter
        north = -90.
        south  = +90.
        west = +1000.
        east = -1000.

        do i = 1,nx
          north = max(north,lat(i,1),lat(i,ny))
          south = min(south,lat(i,1),lat(i,ny))
          east  = max(east ,lon(i,1),lon(i,ny))
          west  = min(west ,lon(i,1),lon(i,ny))
        enddo ! i

        do j = 1,ny
          north = max(rnorth,lat(1,j),lat(nx,j))
          south = min(south ,lat(1,j),lat(nx,j))
          east  = max(east  ,lon(1,j),lon(nx,j))
          west  = min(west  ,lon(1,j),lon(nx,j))
        enddo ! j

        north = north + r_buffer
        south = south - r_buffer
        east  = east  + r_buffer
        west  = west  - r_buffer

       write(6,101)north,south,east,west
101     format(1x,' Box around LAPS grid - NSEW ',4f9.2)

       return
       end
!.............................................................................
