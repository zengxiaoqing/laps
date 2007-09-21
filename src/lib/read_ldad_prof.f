
      subroutine read_ldad_prof(i4time_sys,i4_prof_window
     1                                    ,NX_L,NY_L
     1                                    ,ext,lun
     1                                    ,filename,n_good_obs,istatus)

      character*(*) filename,ext

!.............................................................................

      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status

      write(6,*)
      write(6,*)' Subroutine read_ldad_prof...'

C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN filename'
      endif
C
C  Fill all dimension values
C
C
C Get size of level
C
      nf_status = NF_INQ_DIMID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,level)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim level'
      endif
C
C Get size of maxStaticIds
C
      nf_status = NF_INQ_DIMID(nf_fid,'maxStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxStaticIds'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim maxStaticIds'
      endif
C
C Get size of nInventoryBins
C
      nf_status = NF_INQ_DIMID(nf_fid,'nInventoryBins',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nInventoryBins'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nInventoryBins)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nInventoryBins'
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
      call read_prof(nf_fid, level, maxStaticIds, nInventoryBins,
     +     recNum,ext,lun,
!.............................................................................
     1     i4time_sys,i4_prof_window,NX_L,NY_L,n_good_obs,istatus)

      return
!.............................................................................

      end
C
C
      subroutine read_prof(nf_fid, level, maxStaticIds, nInventoryBins,       
     +     recNum,ext,lun,
!.............................................................................
     1     i4time_sys,i4_prof_window,NX_L,NY_L,n_good_obs,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
      integer assetId(recNum), averageMinutes(recNum), 
     +     firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, prevRecord(recNum), tempQcFlag( level,
     +     recNum), wdQcFlag( level, recNum), windDir( level,
     +     recNum), wsQcFlag( level, recNum)
      real elevation(recNum), latitude(recNum), levels( level,
     +     recNum), longitude(recNum), temperature( level, recNum),
     +     windSpeed( level, recNum)
      double precision observationTime(recNum), receiptTime(recNum),
     +     reportTime(recNum)
      character*6 stationId(recNum)
      character*4 homeWFO(recNum)
      character*24 dataProvider(recNum)
      character*6 providerId(recNum)
      integer stationType(recNum)
      character*30 staticIds(maxStaticIds)
      character*51 stationName(recNum)
!.............................................................................

      character*9 a9_timeObs,a9_recptTime,a9_closest,a9time_ob
      character*8 c8_project,c8_format
      character*6 provider_ref
      character*(*)ext
      real lat_a(NX_L,NY_L)
      real lon_a(NX_L,NY_L)
      real topo_a(NX_L,NY_L)
      real ht_out(200),di_out(200),sp_out(200),temp_out(200)
      integer iqc1_out(200),iqc2_out(200)
      integer assetId_ref
      real mspkt
      data mspkt/.518/

!............................................................................

!     Initialize arrays
      averageMinutes = -999
      stationId = 'UNK   '

      call get_c8_project(c8_project,istatus)

      if(c8_project .eq. 'RSA')then
          c8_format = 'LDAD'
      else
          c8_format = 'MADIS'
      endif

      write(6,*)' Reading profiler/RASS data, format = ',c8_format

      call read_ldad_prof_netcdf(nf_fid, level, maxStaticIds, 
     +     nInventoryBins,       
     +     recNum, assetId, averageMinutes, firstInBin, firstOverflow, 
     +     globalInventory, invTime, inventory, isOverflow, 
     +     lastInBin, lastRecord, nStaticIds, prevRecord, tempQcFlag, 
     +     wdQcFlag, windDir, wsQcFlag, elevation, latitude, levels, 
     +     longitude, temperature, windSpeed, observationTime, 
     +     receiptTime, reportTime, dataProvider, homeWFO, 
     +     providerId, staticIds, stationId, stationName, stationType,
     +     c8_format)
C
C The netcdf variables are filled - your code goes here
C
!............................................................................

      call get_latlon_perimeter(NX_L,NY_L,1.0
     1                         ,lat_a,lon_a,topo_a
     1                         ,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_latlon_perimeter'
          return
      endif

      write(6,*)' # of profilers = ',nStaticIds
      write(6,*)' # of records = ',recnum

      do i_sta = 1,nStaticIds
        ilast_rec = lastRecord(i_sta) + 1 ! Offset going from C to FORTRAN
        provider_ref = providerId(ilast_rec)

        call s_len(staticIds(i_sta),len_id)

        write(6,*)
        write(6,*)' Looping for profiler ',i_sta
     1                                    ,staticIds(i_sta)(1:len_id)

        write(6,*)' Station / last record / provider '
     1           ,i_sta,ilast_rec,provider_ref         

        i4_resid_closest = 999999
        i4_avg_window = -999
        a9_closest = '---------'

        do irec = 1,recnum

          call s_len(providerId(irec),len_prov)

          if(providerId(irec) .eq. provider_ref)then

              write(6,*)
              write(6,*)' ProviderId Match = ',irec
     1                 ,providerId(irec)(1:len_prov)

              rlat = latitude(irec)
              rlon = longitude(irec)

              if(rlat .le. rnorth .and. rlat .ge. south .and.
     1           rlon .ge. west   .and. rlon .le. east            )then        

                  write(6,*)irec,providerId(irec)(1:len_prov)
     1                     ,' is in box'       

                  elev = elevation(irec)

!                 Convert u_std, v_std to rms

!                 Test observation time
                  if(abs(observationTime(irec))      .lt. 3d9)then
                      ictime_ob = nint(observationTime(irec))

                      if(averageMinutes(irec) .ne. -999)then
                          i4_avg_window = averageMinutes(irec) * 60
                          i4_avg_window_half = i4_avg_window / 2
                          ictime_ob = ictime_ob - i4_avg_window_half
                      endif

                      call c_time2fname(ictime_ob,a9_timeObs)

                  else
                      write(6,*)' Bad observation time - reject record'         
     1                           ,observationTime(irec)
                      goto 300

                  endif

!                 Test number of good levels
                  n_good_levels = 0

                  if(ext(1:3) .eq. 'pro')then ! test wind profile
                      do i = 1,level
                          if(windDir(i,irec) .ge. 0      .and.
     1                       windDir(i,irec) .le. 360    .and.
     1                 iqc_rsa(wdQcFlag(i,irec)) .ne. -1 .and.       
     1                 iqc_rsa(wsQcFlag(i,irec)) .ne. -1     
     1                                                          )then ! Good QC
                              n_good_levels = n_good_levels + 1
                          endif
                      enddo ! i

                  elseif(ext(1:3) .eq. 'lrs')then ! test temperature profile
                      do i = 1,level
                          if(iqc_rsa(tempQcFlag(i,irec)) .ne. -1
     1                            .and.
     1                        temperature(i,irec) .gt. 200.
     1                            .and.
     1                        temperature(i,irec) .lt. 400.
     1                                                   )then ! Good QC
                              n_good_levels = n_good_levels + 1
                          endif
                      enddo ! i
                  endif

                  if(n_good_levels .le. 0)then
                      write(6,*)' No good levels - reject record'         
     1                           ,observationTime(irec)
                      goto 300

                  else ! good levels detected
                      write(6,*)' Good levels = ',n_good_levels
     1                           ,observationTime(irec)
                  endif

!                 Determine if this report is closest to the analysis time
                  call cv_asc_i4time(a9_timeObs,i4time_ob)
                  i4_resid = abs(i4time_ob - i4time_sys)
                  if(i4_resid .lt. i4_resid_closest)then
                      i4_resid_closest = i4_resid
                      i_pr_cl = irec
                      a9_closest = a9_timeobs
                  endif
                  write(6,*)'i4_resid/closest/avg_window = '
     1                      ,i4_resid,i4_resid_closest,i4_avg_window       

              else !
                  write(6,*)irec,providerId(irec)(1:len_prov)
     1                     ,' is outside of domain'

                  go to 900 ! loop back to next station

              endif ! in box

!         else
!             write(6,*)' ProviderId = ',irec
!    1                 ,providerId(irec)(1:len_prov)
 
          endif ! correct provider ID

 300      continue

        enddo ! irec 

        write(6,*)
        write(6,*)' Evaluating profile closest in time'

        if(i4_resid_closest .gt. i4_prof_window)then ! outside time window
            write(6,*)' outside time window - reject '
     1               ,a9_closest,i4_resid,i4_prof_window
        
        else
!           lun=1 ! for both 'pro' and 'lrs'
C
C           Open intermediate output file.
C
            call open_ext(lun,i4time_sys,ext(1:3),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error opening product file',ext
                goto980
            endif

            n_good_obs = n_good_obs + 1

            a9time_ob = a9_closest

            n_good_levels = 0

            call filter_string(stationId(i_pr_cl))

            if(ext(1:3) .eq. 'pro')then

                rms = 1.0

                write(6,*)'i/windDir/wdQcFlag/wsQcFlag, i_pr_cl='
     1                                                 ,i_pr_cl     

                do i = 1,level
                    write(6,*)i,windDir(i,i_pr_cl)
     1                         ,wdQcFlag(i,i_pr_cl),wsQcFlag(i,i_pr_cl)       
                    if(windDir(i,i_pr_cl) .ge. 0      .and.
     1                 windDir(i,i_pr_cl) .le. 360    .and.
     1                 iqc_rsa(wdQcFlag(i,i_pr_cl)) .ne. -1 .and.       
     1                 iqc_rsa(wsQcFlag(i,i_pr_cl)) .ne. -1     
     1                                                          )then ! Good QC
                        n_good_levels = n_good_levels + 1
                        ht_out(n_good_levels) = levels(i,i_pr_cl)
                        di_out(n_good_levels) = windDir(i,i_pr_cl)
                        sp_out(n_good_levels) = windSpeed(i,i_pr_cl)
     1                                        * mspkt
                        iqc1_out(n_good_levels) = wdQcFlag(i,i_pr_cl)     
                        iqc2_out(n_good_levels) = wsQcFlag(i,i_pr_cl)     
                    endif
                enddo ! i

                call filter_string(provider_ref)

                if(c8_format .eq. 'LDAD')then
                    write(6,401)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationId(i_pr_cl)(1:6)
     1                         ,a9time_ob,'PROFILER'
                    write(lun,401)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationId(i_pr_cl)(1:6)
     1                         ,a9time_ob,'PROFILER'
                else ! MADIS
                    write(6,402)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:6)
     1                         ,a9time_ob,'PROFILER'
                    write(lun,402)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:6)
     1                         ,a9time_ob,'PROFILER'
                endif

401             format(a12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)
402             format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

                do i = 1,n_good_levels
                    write(lun,411,err=421)ht_out(i)
     1                                   ,di_out(i),sp_out(i)
     1                                   ,rms
                    write(6  ,411,err=421)ht_out(i)
     1                                   ,di_out(i),sp_out(i)
     1                                   ,rms
     1                                   ,iqc1_out(i)
     1                                   ,iqc2_out(i)
411                 format(1x,f6.0,f6.0,2f6.1,2i7)
421                 continue
                enddo ! i

            elseif(ext(1:3) .eq. 'lrs')then

                rms = 1.0
                iqc = 1

                do i = 1,level
                    if(iqc_rsa(tempQcFlag(i,i_pr_cl)) .ne. -1
     1                            .and.
     1                 temperature(i,i_pr_cl) .gt. 200.
     1                            .and.
     1                 temperature(i,i_pr_cl) .lt. 400.
     1                                                   )then ! Good QC
                        n_good_levels = n_good_levels + 1
                        ht_out(n_good_levels) = levels(i,i_pr_cl)
                        temp_out(n_good_levels) = temperature(i,i_pr_cl)
                        iqc1_out(n_good_levels) = tempQcFlag(i,i_pr_cl)
                    endif
                enddo ! i

                call filter_string(provider_ref)

                if(c8_format .eq. 'LDAD')then
                    write(6,501)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationId(i_pr_cl)(1:5)
     1                         ,a9time_ob,'RASS    '
                    write(lun,501)provider_ref
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,stationId(i_pr_cl)(1:5)
     1                         ,a9time_ob,'RASS    '
                else ! MADIS
                    write(6,502)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:5)
     1                         ,a9time_ob,'RASS    '
                    write(lun,502)i_sta
     1                         ,n_good_levels
     1                         ,rlat,rlon,elev,provider_ref(1:5)
     1                         ,a9time_ob,'RASS    '
                endif

501             format(a12,i12,f11.3,f15.3,f15.0,5x,a5,3x,a9,1x,a8)
502             format(i12,i12,f11.3,f15.3,f15.0,5x,a5,3x,a9,1x,a8)

                do i = 1,n_good_levels
                    write(lun,511,err=521)ht_out(i)
     1                                   ,temp_out(i)
     1                                   ,iqc
     1                                   ,rms
                    write(6  ,511,err=521)ht_out(i)
     1                                   ,temp_out(i)
     1                                   ,iqc1_out(i)
     1                                   ,rms
511                 format(1x,f6.0,f6.1,i6,f6.1)
521                 continue
                enddo ! i

            endif ! ext

        endif ! in time window

980     continue

900   enddo ! ista

!............................................................................
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_ldad_prof_netcdf(nf_fid, level, maxStaticIds, 
     +     nInventoryBins, recNum, assetId, averageMinutes, firstInBin,        
     +     firstOverflow, globalInventory, invTime, inventory, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, prevRecord, 
     +     tempQcFlag, wdQcFlag, windDir, wsQcFlag, elevation, 
     +     latitude, levels, longitude, temperature, windSpeed, 
     +     observationTime, receiptTime, reportTime, dataProvider, 
     +     homeWFO, providerId, staticIds, stationId, stationName, 
     +     stationType,c8_format)
C
      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid, 
     +     nf_vid, nf_status
      integer assetId(recNum), averageMinutes(recNum), 
     +     firstInBin(nInventoryBins),
     +     firstOverflow, globalInventory, invTime(recNum),
     +     inventory(maxStaticIds), isOverflow(recNum),
     +     lastInBin(nInventoryBins), lastRecord(maxStaticIds),
     +     nStaticIds, prevRecord(recNum), tempQcFlag( level,
     +     recNum), wdQcFlag( level, recNum), windDir( level,
     +     recNum), wsQcFlag( level, recNum)
      real elevation(recNum), latitude(recNum), levels( level,
     +     recNum), longitude(recNum), temperature( level, recNum),
     +     windSpeed( level, recNum)
      double precision observationTime(recNum), receiptTime(recNum),
     +     reportTime(recNum)
      character*6 stationId(recNum)
      character*4 homeWFO(recNum)
      character*24 dataProvider(recNum)
      character*6 providerId(recNum)
      integer stationType(recNum)
      character*30 staticIds(maxStaticIds)
      character*51 stationName(recNum)

      character*8 c8_format


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      elevation    "Elevation above MSL"
C
      nf_status = NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var elevation'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "Station latitude"
C
      nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var latitude'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      levels       "Instrument Level, height above station"
C
      nf_status = NF_INQ_VARID(nf_fid,'levels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,levels)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var levels'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "Station longitude"
C
      nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var longitude'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      temperature  "Temperature"
C
      nf_status = NF_INQ_VARID(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var temperature'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "Wind Speed (scalar)"
C
      nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var windSpeed'
        endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      assetId      "RSA Asset Identifier"
C
!       nf_status = NF_INQ_VARID(nf_fid,'assetId',nf_vid)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var assetId'
!     endif
!       nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,assetId)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var assetId'
!     endif
C
C     Variable        NETCDF Long Name
C      firstInBin   
C
      nf_status = NF_INQ_VARID(nf_fid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,firstInBin)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var firstInBin'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      averageMinutes
C
      nf_status = NF_INQ_VARID(nf_fid,'averageMinutes',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var averageMinutes'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,averageMinutes)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var averageMinutes'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      firstOverflow
C
      nf_status = NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var firstOverflow'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      globalInventory
C
      nf_status = NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var globalInventory'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      invTime      
C
      nf_status = NF_INQ_VARID(nf_fid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,invTime)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var invTime'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      inventory    
C
      nf_status = NF_INQ_VARID(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,inventory)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var inventory'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      isOverflow   
C
      nf_status = NF_INQ_VARID(nf_fid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,isOverflow)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var isOverflow'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      lastInBin    
C
      nf_status = NF_INQ_VARID(nf_fid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,lastInBin)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var lastInBin'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      lastRecord   
C
      nf_status = NF_INQ_VARID(nf_fid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,lastRecord)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var lastRecord'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      nStaticIds   
C
      nf_status = NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var nStaticIds'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      prevRecord   
C
      nf_status = NF_INQ_VARID(nf_fid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,prevRecord)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var prevRecord'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      tempQcFlag   "RSA Temperature quality control flag"
C
      nf_status = NF_INQ_VARID(nf_fid,'tempQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempQcFlag'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,tempQcFlag)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var tempQcFlag'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      wdQcFlag     "RSA Wind direction quality control flag"
C
      nf_status = NF_INQ_VARID(nf_fid,'wdQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wdQcFlag'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,wdQcFlag)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var wdQcFlag'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "Wind Direction (scalar)"
C
      nf_status = NF_INQ_VARID(nf_fid,'windDir',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windDir'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,windDir)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var windDir'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      wsQcFlag     "RSA Wind speed quality control flag"
C
      nf_status = NF_INQ_VARID(nf_fid,'wsQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wsQcFlag'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,wsQcFlag)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var wsQcFlag'
        endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      observationTime"Time of observation"
C
      nf_status = NF_INQ_VARID(nf_fid,'observationTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      else
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var observationTime'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      receiptTime  "File time stamp (time file was received)"
C
!       nf_status = NF_INQ_VARID(nf_fid,'receiptTime',nf_vid)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var receiptTime'
!     endif
!       nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receiptTime)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var receiptTime'
!     endif
C
C     Variable        NETCDF Long Name
C      reportTime   "Time of observation"
C
!       nf_status = NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var reportTime'
!     endif
!       nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var reportTime'
!      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      dataProvider "Name of organization responsible for delivering the data"
C
      nf_status = NF_INQ_VARID(nf_fid,'dataProvider',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      else
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var dataProvider'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      homeWFO      "home WFO Id"
C
!       nf_status = NF_INQ_VARID(nf_fid,'homeWFO',nf_vid)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var homeWFO'
!     endif
!       nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,homeWFO)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var homeWFO'
!     endif
C
C     Variable        NETCDF Long Name
C      providerId   "Alphanumeric station name"
C
      nf_status = NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      else
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var providerId'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      staticIds    
C
      nf_status = NF_INQ_VARID(nf_fid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      else
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,staticIds)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var staticIds'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      stationId    "alphanumeric station Id"
C
      nf_status = NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      else
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var stationId'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      stationName  "alphanumeric station name"
C
!       nf_status = NF_INQ_VARID(nf_fid,'stationName',nf_vid)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var stationName'
!     endif
!       nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
!     if(nf_status.ne.NF_NOERR) then
!       print *, NF_STRERROR(nf_status)
!       print *,'in var stationName'
!     endif


!     For RSA, this variable doesn't contain any useful information.
!     This could be useful if 'c8_format' is MADIS, if we read it in 
!     given its "short" declaration as an INT variable.

C
C     Variable        NETCDF Long Name
C      stationType  "LDAD station type"
C
      if(c8_format .eq. 'MADIS')then
        write(6,*)' read stationType as short/integer'
        nf_status = NF_INQ_VARID(nf_fid,'stationType',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var stationType'
        else
          nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,stationType)
          if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
            print *,'in var stationType'
          endif
        endif

      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end


      function iqc_rsa(iflag)

!     Possible Outputs
!     -1 Bad Data
!      0 Unknown quality or unknown input flag
!     +1 Good Data

      iqc_rsa = 0

      if(iflag .eq.     0)iqc_rsa = +1 ! OK
      if(iflag .eq.     1)iqc_rsa = -1 ! Out of range
      if(iflag .eq.     2)iqc_rsa = -1 ! Questionable
      if(iflag .eq.     3)iqc_rsa =  0 ! Not tested
      if(iflag .eq.     4)iqc_rsa = -1 ! Missing Data
      if(iflag .eq.     5)iqc_rsa = +1 ! MFFG Auto Mode Algorithm good data
      if(iflag .eq.     6)iqc_rsa = +1 ! MFFG Manual Mode Algorithm good data 
                                       ! Operator no action 
      if(iflag .eq.     7)iqc_rsa = -1 ! MFFG Auto Mode Algorithm bad data
      if(iflag .eq.     8)iqc_rsa = -1 ! MFFG Auto Mode Algorithm bad data
                                       ! Operator no action 
      if(iflag .eq.     9)iqc_rsa =  0 ! reserved for MFFG but not used
      if(iflag .eq.    10)iqc_rsa = -1 ! MFFG Manual Mode Algorithm good data
                                       ! Operator flagged as bad data
      if(iflag .eq.    11)iqc_rsa =  0 ! reserved for MFFG but not used
      if(iflag .eq.    12)iqc_rsa = -1 ! MFFG Manual Mode Algorithm bad data 
                                       ! Operator flagged as bad data
      if(iflag .eq.    13)iqc_rsa =  0 ! reserved for MFFG but not used
      if(iflag .eq.    14)iqc_rsa = +1 ! MFFG Manual Mode Algorithm good data 
                                       ! Operator no action
      if(iflag .eq.    15)iqc_rsa =  0 ! reserved for MFFG but not used
      if(iflag .eq.    16)iqc_rsa = -1 ! MFFG Manual Mode Algorithm bad data 
                                       ! Operator no action
      if(iflag .eq.    17)iqc_rsa =  0 ! reserved for MFFG but not used
      if(iflag .eq.    18)iqc_rsa = -1 ! MFFG Manual Mode Algorithm good data 
                                       ! Operator flagged as bad data
      if(iflag .eq.    19)iqc_rsa =  0 ! reserved for MFFG but not used
      if(iflag .eq.    20)iqc_rsa = -1 ! MFFG Manual Mode Algorithm bad data 
                                       ! Operator flagged as bad data
      if(iflag .eq. -9999)iqc_rsa = -1 ! Missing Data

      return
      end
