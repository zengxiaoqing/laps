c
        subroutine get_local_towerobs(maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 lat,lon,ni,nj, nobs, istatus)

c
c.....  Input variables/arrays
c
        integer maxsta ! processed stations for SND file
        integer maxlvls ! raw/processed stations for SND file

        parameter (maxlvls=100)

        character*(*) path_to_local_data, local_format

        real    lat(ni,nj), lon(ni,nj)
c
c.....  Local variables/arrays
c
        double precision d_timeobs

!       Obs arrays (raw files)
        integer     nobs,nlvls
	real*4      lats(maxsta), lons(maxsta)
        real*4      lvls_m(maxlvls,maxsta)
        real*4      fill_stationP, fill_lvls, stationP
        integer*4   i4time
        character*51  stationName
        character*6   c_staId
	real*4      dd(maxlvls,maxsta), ff(maxlvls,maxsta)
        character*9 a9time_a(maxlvls,maxsta)

        integer*4  i4time_ob_a(maxsta), before, after
c
c.....  Variables used by write_snd
c
        integer    nsnd_all ! combined # of obs over multiple files
	integer*4  wmoid(maxsta)
        real*4     stalat(maxsta,maxlvls),stalon(maxsta,maxlvls)
        real*4     staelev(maxsta)
        character  c5_staid(maxsta)*5, a9time_ob(maxsta,maxlvls)*9
        character  c8_obstype(maxsta)*8
        real*4     height_m(maxsta,maxlvls), pressure_pa(maxsta,maxlvls)
        real*4     temp_c(maxsta,maxlvls), dewpoint_c(maxsta,maxlvls)
        real*4     dir_deg(maxsta,maxlvls),spd_mps(maxsta,maxlvls)
c
c.....  Output arrays.
c
	character  stations(maxsta)*20
	character  reptype(maxsta)*6, atype(maxsta)*6

c.....  Unknown vars.
	character  stname(maxsta)*6, save_stn(maxsta)*6

	integer    rtime
	integer    recNum, nf_fid, nf_vid, nf_status
	character  timech*9, time*4
        character*13 filename13, cvt_i4time_wfo_fname13
        character*150 data_file 
c
c.....  Start.
c
c       Get r_missing_data
        call get_r_missing_data(r_missing_data,istatus)
c       
c.....	Set istatus flag for the local data to bad until we find otherwise.
c
	istatus = -1

        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

!       call get_box_size(box_size,istatus)
!       if(istatus .ne. 1)return

        box_size = 1.0

c
c.....  Figure out the size of the "box" in gridpoints.  User defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c

        call get_grid_spacing_cen(grid_spacing_m,istatus)
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing_m / 1000.) !in km
c
c.....	Zero out the counters.
c
        nobs = 0	        ! # of local obs in the laps grid
c
c.....  Get the data from the NetCDF file(s).  First, open the file(s).
c.....  If none are there, return to tower_driver.
c
        ix = 1
c
c.....  Set up the time window.
c
	i4time_before = i4time_sys - itime_before
	i4time_after  = i4time_sys + itime_after

!       Ob times contained in each file
        i4_contains_early = 0 
        i4_contains_late = 3599

        call get_filetime_range(i4time_before,i4time_after                
     1                         ,i4_contains_early,i4_contains_late       
     1                         ,3600                                     
     1                         ,i4time_file_before,i4time_file_after)

        write(6,*)' i4time_file_before,i4time_file_after:'
     1             ,i4time_file_before,i4time_file_after

        do i4time_file = i4time_file_before, i4time_file_after, +3600       

            call s_len(path_to_local_data,len_path)
            filename13= cvt_i4time_wfo_fname13(i4time_file)

            if(len_path .lt. 1)goto590
 	    data_file = path_to_local_data(1:len_path)//filename13

            write(6,*)' LDAD tower file = ',data_file(1:len_path+13)
c
c.....  Call the read routine.
c
	    call read_local_tower(data_file, maxsta, maxlvls,     ! I
     &         r_missing_data,                                    ! I
     &         nsnd_file, nlvls, lvls_m(1,ix),                    ! O
     &         staelev(ix), stalat(1,ix), stalon(1,ix),           ! O
     &         temp_c(1,ix), dewpoint_c(1,ix),                    ! O
     &         dd(1,ix), ff(1,ix),                                ! O
     &         a9time_ob(1,ix), stname(ix), wmoid(ix),            ! O
     &         istatus)                                           ! O
 
	    if(istatus .ne. 1)then
                write(6,*)
     1          '     Warning: bad status return from read_local_tower'       
                nsnd_file = 0

            else
                write(6,*)'     nsnd_file = ',nsnd_file

            endif

            ix = ix + nsnd_file

590     enddo ! i4time_file

        nsnd_all = ix - 1
        write(6,*)' nsnd_all = ',nsnd_all
c
        max_write = 100
      
c.....  Call the routine to write the SND file.
c

        print *
	print *,'  Appending SND file, # of obs (in grid) = ',nobs

        call write_snd(    lun_out                         ! I
     1                    ,maxsta,maxlvl,nsnd_all          ! I
     1                    ,iwmostanum                      ! I
     1                    ,stalat,stalon,staelev           ! I
     1                    ,stname,a9time_ob,c8_obstype     ! I
     1                    ,nlvl                            ! I
     1                    ,height_m                        ! I
     1                    ,pressure_pa                     ! I
     1                    ,temp_c                          ! I
     1                    ,dewpoint_c                      ! I
     1                    ,dir_deg                         ! I
     1                    ,spd_mps                         ! I
     1                    ,istatus)                        ! O


         istatus = 1            ! everything's ok...
         return
c
 990     continue               ! no data available
         istatus = 0
         print *,' No data available from GET_LOCAL_TOWEROBS'
         return
c
         end

         subroutine read_local_tower(data_file, maxsta, maxlvls,  ! I
     &         r_missing_data,                                    ! I
     &         nobs, nlvls, lvls_m,                               ! O
     &         staelev, stalat, stalon,                           ! O
     &         temp_c, dewpoint_c, dir_deg, spd_mps               ! O
     &         a9time_ob, stname, wmoid,                          ! O
     &         istatus)                                           ! O

      character*(*) data_file 
      integer       maxsta ! processed stations for SND file
      integer       maxlvls ! raw/processed stations for SND file
      real*4        r_missing_data 
      integer       nobs,nlvls
      real*4        lats(maxsta), lons(maxsta)
      real*4        lvls_m(maxlvls,maxsta)
      real*4        staelev(maxsta)
      real*4        stalat(maxsta,maxlvls),stalon(maxsta,maxlvls)
      real*4        dd(maxlvls,maxsta), ff(maxlvls,maxsta)
      real*4        temp_k, rh_pct,stationP,ws,wd
      real*4        height_m(maxsta,maxlvls)
      real*4        pressure_pa(maxsta,maxlvls)      
      real*4        temp_c(maxsta,maxlvls), dewpoint_c(maxsta,maxlvls)
      real*4        dir_deg(maxsta,maxlvls),spd_mps(maxsta,maxlvls)
      real*4        sp_fill,levels_fill
      integer*4     wmoid(maxsta),nlvl(maxsta),sp_id
      integer       istatus
      integer       sta_id,sn_id,wd_id,ws_id,rh_id,temp_id,lev_id,ot_id
      integer       index_1(1), index_2(2),start(2),count(2)
      integer       pi_len, sn_len
      character     stname(maxsta)*6, stationName*51, c_staid*6
      character     a9time_ob(maxsta,maxlvls)*9

c     open data_file
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error on NF_OPEN - open file: ',filename
      endif

c     read dim recNum -> nobs
      nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error finding dim recNum' 
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nobs)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error reading dim recNum'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

c     verify nobs .lt. maxsta
      if (nobs .gt. maxsta) then
        print *,'nobs is greater than maxsta: ',nobs,' ',maxsta
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

c     read dim level -> nlvls
      nf_status = NF_INQ_DIMID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error finding dim level' 
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nlvls)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error reading dim level'
      endif

c     verify nlvls .lt. maxlvls
      if (nlvls .gt. maxlvls) then
        print *,'nlvls is greater than maxlvls: ',nlvls,' ',maxlvls
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

c     read dim level staNameLen
      nf_status = NF_INQ_DIMID(nf_fid,'staNameLen',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error finding dim staNameLen' 
        print *, 'set sn_len = 51'
        sn_len = 51
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,sn_len)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error reading dim staNameLen'
        print *, 'set sn_len = 51'
        sn_len = 51
      endif

c     read dim level providerIDLen
      nf_status = NF_INQ_DIMID(nf_fid,'providerIDLen',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error finding dim providerIDLen' 
        print *, 'set pi_len = 6'
        pi_len = 6
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,pi_len)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error reading dim providerIDLen'
        print *, 'set pi_len = 6'
        pi_len = 6
      endif

c     read var elevation(recNum) -> staelev(maxsta)
      nf_status = NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var elevation'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,staelev)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var elevation'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read var latitude(recNum -> stalat(maxsta)
      nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var latitude'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,stalat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var latitude'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read var longitude(recNum) -> stalon(maxsta)
      nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var longitude'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,stalon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var longitude'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read var levels(recNum,level) -> lvls_m(maxlvls,maxsta)
      nf_status = NF_INQ_VARID(nf_fid,'levels',lev_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var levels'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,lvls_m(1,recnum))
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var levels'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read _fillValue for levels
      nf_status = NF_INQ_ATTLEN(nf_fid, nf_vid,'_FillValue',
     1                          sp_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading stationPressure _FillValue attribute'
        levels_fill = 1.0e+38
      endif 

c     get varids for stationId, stationName, temperature,stationPressure
c       relHumidity, windSpeed, windDir, observationTime

      nf_status = NF_INQ_VARID(nf_fid,'stationId',sta_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var stationId'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'stationName',sta_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var stationName'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'temperature',sta_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var temperature'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'stationPressure',sp_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var stationPressure'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

c     read _fillValue for stationPressure
      nf_status = NF_INQ_ATTLEN(nf_fid, sp_id,'_FillValue',
     1                          sp_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading stationPressure _FillValue attribute'
        sp_fill = 3.4028e+38
      endif 

      nf_status = NF_INQ_VARID(nf_fid,'relHumidity',rh_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var relHumidity'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'windSpeed',ws_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var windSpeed'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'windDir',wd_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var windDir'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'observationTime',ot_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var observationTime'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      do obno = 1, nobs

        index_1(1) = obno
        index_2(1) = obno
        start(1) = obno
        start(2) = 1
        count(1) = 1

c       read var observationTime(recNum) -> d_timeobs
        nf_status = NF_GET_VAR1_double(nf_fid,ot_id,index_1,d_timeobs)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading var levels'
          print *, 'Aborting read'
          nf_status = NF_CLOSE(nf_fid)
          istatus = 0
          return
        endif 

        i4_tim = int(d_timeobs - 315619200)
        call make_fnam_lp(i4_tim,a9time_ob(obno,1),istatus)

c       read var stationName(obno,staNameLen) -> stationName
        count(2) = sn_len 
        nf_status = NF_GET_VARA_TEXT(nf_fid,sn_id,start,count,
     1                               stationName)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading var levels'
          print *, 'Aborting read'
          nf_status = NF_CLOSE(nf_fid)
          istatus = 0
          return
        endif

c       truncate stname(obno)  = stationName(1:5)
        stname(obno)  = stationName(1:5)

c       read var stationId(recNum,providerIDLen) -> c_staid 
        count(2) = pi_len 
        nf_status = NF_GET_VARA_TEXT(nf_fid,sta_id,start,count,
     1                               c_staid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading var stationId'
          print *, 'Aborting read'
          nf_status = NF_CLOSE(nf_fid)
          istatus = 0
          return
        endif

c       NEEDS DOING
c       convert string to iwmostanum(maxsta) (cvt S to 0 and N to 1)
       
!       Replace blank/UNK station names with wmoid if possible, else set blank
        iblank = 0
        call s_len(stname(obno),lensta)
!       if(lensta .eq. 0 .or. stations(i)(1:3) .eq. 'UNK')then
        if(.false.)then
          if(wmoid(i) .ne. ibadflag .and. wmoid(i) .ne. 0)then
!           write(stations(i),511,err=512)wmoid(i)
511         format(i8)
512         continue
          else
!           stations(i) = 'UNK                 '
            iblank = iblank + 1
          endif
        endif

        do lvl = 1, nlvls  
          if(lvls_m(lvl,recnum) .ne. _fillValue)then
            index_2(2) = lvl
c           read stationPressure
            nf_status = NF_GET_VAR1_REAL(nf_fid,sp_id,index_1,stationP)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var stationPressure'
            endif 

c           check stationPressure for _FillValue
            if (stationP .eq. sp_fill) stationP = r_missing_data
            if (lvl.eq.1) then
              pressure_pa(obno,lvl) = stationP
            else
              pressure_pa(obno,lvl) = r_missing_data
            endif

c           read var temperature(recNum,lvl) -> temp
            nf_status = NF_GET_VAR1_REAL(nf_fid,temp_id,index_2,temp_k)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var temperature'
            endif 

c           NEEDS DOING?
            temp_c(lvl,maxsta) = k_to_c(temp_k) ! to convert kelvin to celsius

c           read var relHumidity(recNum,level) -> rh
            nf_status = NF_GET_VAR1_REAL(nf_fid,rh_id,index_2,rh_pct)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var relHumidity'
            endif 

c           NEEDS DOING?
            dewpoint_c(lvl,maxsta) = dwpt(temp_c(lvl,maxsta), rh_pct) ! celsius

c           read var windSpeed(recNum,level) -> ws
            nf_status = NF_GET_VAR1_REAL(nf_fid,ws_id,index_2,ws)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var windSpeed'
            endif 
            spd_mps(obno,lvl) = ws

c           read var windDir(recNum,level) -> dd(lvl,obno)
            nf_status = NF_GET_VAR1_REAL(nf_fid,wd_id,index_2,wd)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var windDir'
            endif 
            dir_deg(obno,lvl) = wd

c           NEEDS DOING?
c           duplicate a9_time_obs(obno,1) for each level
            a9time_ob(obno,:) = a9time_ob(obno,1)
          else
            nlvl(obno) = lvl - 1
          endif
        enddo
      enddo

!     Final QC check
      call get_ibadflag(ibadflag,istatus)
      if(istatus .ne. 1)return

      if(iblank .gt. 0)then
        write(6,*)' Warning: number of UNK stanames = ',iblank
      endif


      return
      end
