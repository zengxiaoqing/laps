
      subroutine read_ldad_prof(i4time_sys,i4_prof_window
     1                                    ,NX_L,NY_L
     1                                    ,ext
     1                                    ,filename,istatus)

      character*(*) filename,ext

!.............................................................................

      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
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
     +     recNum,ext,
!.............................................................................
     1              i4time_sys,i4_prof_window,NX_L,NY_L,istatus)

      return
!.............................................................................

      end
C
C
      subroutine read_prof(nf_fid, level, maxStaticIds, nInventoryBins,       
     +     recNum,
!.............................................................................
     1              ext,i4time_sys,i4_prof_window,NX_L,NY_L,istatus)
!.............................................................................

      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid,
     +     nf_vid, nf_status
      integer assetId(recNum), firstInBin(nInventoryBins),
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
      character*11 stationType(recNum)
      character*30 staticIds(maxStaticIds)
      character*51 stationName(recNum)
!.............................................................................

      character*9 a9_timeObs,a9_recptTime,a9_closest,a9time_ob
      character*(*)ext
      real*4 lat_a(NX_L,NY_L)
      real*4 lon_a(NX_L,NY_L)
      real*4 topo_a(NX_L,NY_L)
      real*4 ht_out(200),di_out(200),sp_out(200)
      integer assetId_ref
!............................................................................

      call read_ldad_prof_netcdf(nf_fid, level, maxStaticIds, 
     +     nInventoryBins,       
     +     recNum, assetId, firstInBin, firstOverflow, 
     +     globalInventory, invTime, inventory, isOverflow, 
     +     lastInBin, lastRecord, nStaticIds, prevRecord, tempQcFlag, 
     +     wdQcFlag, windDir, wsQcFlag, elevation, latitude, levels, 
     +     longitude, temperature, windSpeed, observationTime, 
     +     receiptTime, reportTime, dataProvider, homeWFO, 
     +     providerId, staticIds, stationId, stationName, stationType)
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
        i_pr_ref = lastRecord(i_sta) + 1
        assetId_ref = assetId(i_pr_ref)

        write(6,*)' Looping for profiler ',i_sta,staticIds(i_sta)
        write(6,*)' Station / last record ',i_sta,i_pr_ref,assetId_ref         

        i4_resid_closest = 999999

        do irec = 1,recnum

          if(assetId(irec) .eq. assetId_ref)then

              rlat = latitude(irec)
              rlon = longitude(irec)

              if(rlat .le. rnorth .and. rlat .ge. south .and.
     1           rlon .ge. west   .and. rlon .le. east            )then        

                  write(6,*)staticIds(i_sta),' is in box'

                  elev = elevation(irec)

                  write(6,*)staticIds(i_sta)

!                 Convert u_std, v_std to rms

                  write(6,*)

                  if(abs(observationTime(irec))      .lt. 3d9)then
                      call c_time2fname(nint(observationTime(irec))
     1                                        ,a9_timeObs)
                  else
                      write(6,*)' Bad observation time - reject'       
     1                           ,observationTime(irec)
                      goto 900
                  endif

                  call cv_asc_i4time(a9_timeObs,i4time_ob)
                  i4_resid = abs(i4time_ob - i4time_sys)
                  if(i4_resid .lt. i4_resid_closest)then
                      i4_resid_closest = i4_resid
                      i_pr_cl = irec
                      a9_closest = a9_timeobs
                  endif
                  write(6,*)'i4_resid/closest = '
     1                      ,i4_resid,i4_resid_closest       

              else !
                  write(6,*)staticIds(i_sta),' is outside of domain'

              endif ! in box

          endif ! correct asset ID

        enddo ! irec 

        if(i4_resid_closest .gt. i4_prof_window)then ! outside time window
            write(6,*)' time - reject '
     1               ,a9_closest,i4_resid,i4_prof_window
        
        else
            if(ext(1:3) .eq. 'pro')lun=1
            if(ext(1:3) .eq. 'lrs')lun=1
C
C           Open intermediate output file.
C
            call open_ext(lun,i4time_sys,ext(1:3),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error opening product file',ext
                goto980
            endif

            a9time_ob = a9_closest

            n_good_levels = 0

            call filter_string(stationId(i_pr_cl))

            if(ext(1:3) .eq. 'pro')then

                rms = 1.0

                do i = 1,level
                    if(windDir(i,i_pr_cl) .ge. 0.   .and.
     1                 windDir(i,i_pr_cl) .le. 360.       )then ! Good QC
                        n_good_levels = n_good_levels + 1
                        ht_out(n_good_levels) = levels(i,i_pr_cl)
                        di_out(n_good_levels) = windDir(i,i_pr_cl)
                        sp_out(n_good_levels) = windSpeed(i,i_pr_cl)
                    endif
                enddo ! i

                write(6,401)assetId_ref
     1                     ,n_good_levels
     1                     ,rlat,rlon,elev,stationId(i_pr_cl)(1:6)
     1                     ,a9time_ob,'PROFILER'
                write(lun,401)assetId_ref
     1                     ,n_good_levels
     1                     ,rlat,rlon,elev,stationId(i_pr_cl)(1:6)
     1                     ,a9time_ob,'PROFILER'
401             format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

                do i = 1,n_good_levels
                    write(lun,301,err=303)ht_out(i)
     1                                   ,di_out(i),sp_out(i)
     1                                   ,rms
                    write(6  ,301,err=303)ht_out(i)
     1                                   ,di_out(i),sp_out(i)
     1                                   ,rms
301                 format(1x,f6.0,f6.0,2f6.1)
303                 continue
                enddo ! i

            elseif(ext(1:3) .eq. 'lrs')then

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
     +     nInventoryBins, recNum, assetId, firstInBin, 
     +     firstOverflow, globalInventory, invTime, inventory, 
     +     isOverflow, lastInBin, lastRecord, nStaticIds, prevRecord, 
     +     tempQcFlag, wdQcFlag, windDir, wsQcFlag, elevation, 
     +     latitude, levels, longitude, temperature, windSpeed, 
     +     observationTime, receiptTime, reportTime, dataProvider, 
     +     homeWFO, providerId, staticIds, stationId, stationName, 
     +     stationType)
C
      include 'netcdf.inc'
      integer level, maxStaticIds, nInventoryBins, recNum,nf_fid, 
     +     nf_vid, nf_status
      integer assetId(recNum), firstInBin(nInventoryBins),
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
      character*11 stationType(recNum)
      character*30 staticIds(maxStaticIds)
      character*51 stationName(recNum)


C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      elevation    "Elevation above MSL"
C
        nf_status = NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,elevation)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var elevation'
      endif
C
C     Variable        NETCDF Long Name
C      latitude     "Station latitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,latitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var latitude'
      endif
C
C     Variable        NETCDF Long Name
C      levels       "Instrument Level, height above station"
C
        nf_status = NF_INQ_VARID(nf_fid,'levels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,levels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var levels'
      endif
C
C     Variable        NETCDF Long Name
C      longitude    "Station longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,longitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var longitude'
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
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,temperature)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var temperature'
      endif
C
C     Variable        NETCDF Long Name
C      windSpeed    "Wind Speed (scalar)"
C
        nf_status = NF_INQ_VARID(nf_fid,'windSpeed',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,windSpeed)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var windSpeed'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      assetId      "RSA Asset Identifier"
C
        nf_status = NF_INQ_VARID(nf_fid,'assetId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var assetId'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,assetId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var assetId'
      endif
C
C     Variable        NETCDF Long Name
C      firstInBin   
C
        nf_status = NF_INQ_VARID(nf_fid,'firstInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,firstInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstInBin'
      endif
C
C     Variable        NETCDF Long Name
C      firstOverflow
C
        nf_status = NF_INQ_VARID(nf_fid,'firstOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,firstOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var firstOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      globalInventory
C
        nf_status = NF_INQ_VARID(nf_fid,'globalInventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,globalInventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var globalInventory'
      endif
C
C     Variable        NETCDF Long Name
C      invTime      
C
        nf_status = NF_INQ_VARID(nf_fid,'invTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,invTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var invTime'
      endif
C
C     Variable        NETCDF Long Name
C      inventory    
C
        nf_status = NF_INQ_VARID(nf_fid,'inventory',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,inventory)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var inventory'
      endif
C
C     Variable        NETCDF Long Name
C      isOverflow   
C
        nf_status = NF_INQ_VARID(nf_fid,'isOverflow',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,isOverflow)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isOverflow'
      endif
C
C     Variable        NETCDF Long Name
C      lastInBin    
C
        nf_status = NF_INQ_VARID(nf_fid,'lastInBin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,lastInBin)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastInBin'
      endif
C
C     Variable        NETCDF Long Name
C      lastRecord   
C
        nf_status = NF_INQ_VARID(nf_fid,'lastRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,lastRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lastRecord'
      endif
C
C     Variable        NETCDF Long Name
C      nStaticIds   
C
        nf_status = NF_INQ_VARID(nf_fid,'nStaticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,nStaticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var nStaticIds'
      endif
C
C     Variable        NETCDF Long Name
C      prevRecord   
C
        nf_status = NF_INQ_VARID(nf_fid,'prevRecord',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,prevRecord)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var prevRecord'
      endif
C
C     Variable        NETCDF Long Name
C      tempQcFlag   "RSA Temperature quality control flag"
C
        nf_status = NF_INQ_VARID(nf_fid,'tempQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempQcFlag'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,tempQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var tempQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      wdQcFlag     "RSA Wind direction quality control flag"
C
        nf_status = NF_INQ_VARID(nf_fid,'wdQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wdQcFlag'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,wdQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wdQcFlag'
      endif
C
C     Variable        NETCDF Long Name
C      windDir      "Wind Direction (scalar)"
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
C      wsQcFlag     "RSA Wind speed quality control flag"
C
        nf_status = NF_INQ_VARID(nf_fid,'wsQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wsQcFlag'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,wsQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wsQcFlag'
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
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,observationTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var observationTime'
      endif
C
C     Variable        NETCDF Long Name
C      receiptTime  "File time stamp (time file was received)"
C
        nf_status = NF_INQ_VARID(nf_fid,'receiptTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receiptTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,receiptTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var receiptTime'
      endif
C
C     Variable        NETCDF Long Name
C      reportTime   "Time of observation"
C
        nf_status = NF_INQ_VARID(nf_fid,'reportTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif
        nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reportTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reportTime'
      endif


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
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,dataProvider)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var dataProvider'
      endif
C
C     Variable        NETCDF Long Name
C      homeWFO      "home WFO Id"
C
        nf_status = NF_INQ_VARID(nf_fid,'homeWFO',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var homeWFO'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,homeWFO)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var homeWFO'
      endif
C
C     Variable        NETCDF Long Name
C      providerId   "Alphanumeric station name"
C
        nf_status = NF_INQ_VARID(nf_fid,'providerId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,providerId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var providerId'
      endif
C
C     Variable        NETCDF Long Name
C      staticIds    
C
        nf_status = NF_INQ_VARID(nf_fid,'staticIds',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,staticIds)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var staticIds'
      endif
C
C     Variable        NETCDF Long Name
C      stationId    "alphanumeric station Id"
C
        nf_status = NF_INQ_VARID(nf_fid,'stationId',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationId)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationId'
      endif
C
C     Variable        NETCDF Long Name
C      stationName  "alphanumeric station name"
C
        nf_status = NF_INQ_VARID(nf_fid,'stationName',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationName)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationName'
      endif
C
C     Variable        NETCDF Long Name
C      stationType  "LDAD station type"
C
        nf_status = NF_INQ_VARID(nf_fid,'stationType',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,stationType)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var stationType'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
