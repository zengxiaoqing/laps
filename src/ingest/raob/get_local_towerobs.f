c
        subroutine get_local_towerobs(maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 lat,lon,ni,nj, nsta, istatus)

c
c.....  Input variables/arrays
c
        integer maxlvls ! raw/processed stations for SND file
        integer maxobs ! raw stations in NetCDF files
        integer maxsta ! processed stations for SND file 

        parameter (maxlvls=100)
        parameter (maxobs=1000) ! Raw stations in NetCDf files

        character*(*) path_to_local_data, local_format

        real    lat(ni,nj), lon(ni,nj)
c
c.....  Local variables/arrays
c
        double precision d_timeobs

!       Obs arrays (raw files)
        integer     nobs,nlvls,nlvl(maxobs)
	real*4      lats(maxobs), lons(maxobs)
        real*4      lvls_m(maxlvls,maxobs)
        real*4      fill_stationP, fill_lvls, stationP
        integer*4   i4time
        character*51  stationName
        character*6   c_staId
	real*4      dd(maxlvls,maxobs), ff(maxlvls,maxobs)

        integer*4  i4time_ob_a(maxobs), before, after
c
c.....  Variables used by write_snd - all obs not yet whittled down
c
        integer    nsnd_all ! combined # of obs over multiple files
	integer*4  wmoid(maxobs)
        real*4     stalat(maxobs,maxlvls),stalon(maxobs,maxlvls)
        real*4     staelev(maxobs)
        character  c5_staid(maxobs)*5, a9time_ob(maxobs)*9
!       character  c8_obstype(maxobs)*8
        real*4     height_m(maxobs,maxlvls), pressure_mb(maxobs,maxlvls)
        real*4     temp_c(maxobs,maxlvls), dewpoint_c(maxobs,maxlvls)
        real*4     dir_deg(maxobs,maxlvls),spd_mps(maxobs,maxlvls)
	character  stname(maxobs)*6
c
c.....  Output arrays.
c
        integer    nsnd_all_s ! combined # of obs over multiple files
        integer    nlvl_s(maxsta)
	integer*4  wmoid_s(maxsta)
        real*4     stalat_s(maxsta,maxlvls),stalon_s(maxsta,maxlvls)
        real*4     staelev_s(maxsta)
        character  c5_staid_s(maxsta)*5, a9time_ob_s(maxsta,maxlvls)*9
        character  c8_obstype_s(maxsta)*8
        real*4     height_m_s(maxsta,maxlvls)
        real*4     pressure_mb_s(maxsta,maxlvls)
        real*4     temp_c_s(maxsta,maxlvls)
        real*4     dewpoint_c_s(maxsta,maxlvls)      
        real*4     dir_deg_s(maxsta,maxlvls),spd_mps_s(maxsta,maxlvls)       
	character  stname_s(maxsta)*6

c.....  Unknown vars.
	character  save_stn(maxsta)*6

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
	    call read_local_tower(data_file,len_path+13,          ! I 
     &         maxobs, maxlvls,                                   ! I
     &         r_missing_data,                                    ! I
     &         nsnd_file, nlvl, lvls_m(1,ix),                     ! O
     &         staelev(ix), stalat(ix,1), stalon(ix,1),           ! O
     &         temp_c(ix,1), dewpoint_c(ix,1),                    ! O
     &         height_m(ix,1),                                    ! O
c    &         pressure_mb(ix,1),                                 ! O
c    &         dir_deg(ix,1), spd_mps(ix,1),                      ! O
     &         a9time_ob(ix), stname(ix), wmoid(ix),              ! O
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
c.....  Call the routine to write the SND file.
c

        print *
	print *,'  Appending SND file, # of obs (in grid) = ',nsnd_all       

        pressure_mb_s = r_missing_data

        nsta = 0
        do i = 1,nsnd_all
            if(nlvl(i) .gt. 0)then ! Valid sounding - use for output
                nsta = nsta + 1
                write(6,*)
     1           ' Valid sounding - transferring to output arrays ',nsta      
                stalat_s(nsta,:) = stalat(i,:)           
                stalon_s(nsta,:) = stalon(i,:)           
                staelev_s(nsta) = staelev(i)           
                stname_s(nsta) = stname(i)           
                a9time_ob_s(nsta,:) = a9time_ob(i)           
                c8_obstype_s(nsta) = 'TOWER   '
                nlvl_s(nsta) = nlvl(i)           
!               height_m_s(nsta,:) = height_m(i,:)           
                height_m_s(nsta,:) = lvls_m(:,i) + staelev_s(nsta)           
                temp_c_s(nsta,:) = temp_c(i,:)           
                dewpoint_c_s(nsta,:) = dewpoint_c(i,:)           
                dir_deg_s(nsta,:) = dir_deg(i,:)           
                spd_mps_s(nsta,:) = spd_mps(i,:)           
            endif
        enddo ! i

        lun_out = 11
        if(nsta .gt. 0)then
            call open_ext(lun_out,i4time_sys,'snd',istatus)
        endif

        call write_snd(    lun_out                               ! I
     1                    ,maxsta,maxlvl,nsta                    ! I
     1                    ,wmoid_s                               ! I
     1                    ,stalat_s,stalon_s,staelev_s           ! I
     1                    ,stname_s,a9time_ob_s,c8_obstype_s     ! I
     1                    ,nlvl_s                                ! I
     1                    ,height_m_s                            ! I
     1                    ,pressure_mb_s                         ! I
     1                    ,temp_c_s                              ! I
     1                    ,dewpoint_c_s                          ! I
     1                    ,dir_deg_s                             ! I
     1                    ,spd_mps_s                             ! I
     1                    ,istatus)                              ! O

        if(istatus .ne. 1)then
            write(6,*)
     1       ' get_local_towerobs: Bad status returned from write_snd'       
        endif

        return
c
 990    continue               ! no data available
        istatus = 0
        print *,' No data available from GET_LOCAL_TOWEROBS'
        return
c
        end

         subroutine read_local_tower(filename,fn_len,             ! I 
     &         maxobs, maxlvls,                                   ! I
     &         r_missing_data,                                    ! I
     &         nobs, nlvl, lvls_m,                                ! O
     &         staelev, stalat, stalon,                           ! O
     &         temp_c, dewpoint_c, height_m,                      ! O
c    &         pressure_pa, dir_deg, spd_mps                      ! O
     &         a9time_ob, stname, wmoid,                          ! O
     &         istatus)                                           ! O

      character*(*) filename 
      integer       maxobs ! raw stations for SND file
      integer       maxlvls ! raw/processed stations for SND file
      real*4        r_missing_data 
      integer       nobs,nlvls,lev_set
      real*4        lats(maxobs), lons(maxobs)
      real*4        lvls_m(maxlvls,maxobs)
      real*4        staelev(maxobs)
      real*4        stalat(maxobs,maxlvls),stalon(maxobs,maxlvls)
      real*4        dd(maxlvls,maxobs), ff(maxlvls,maxobs)
      real*4        temp_k, rh_pct,stationP,ws,wd
      real*4        height_m(maxobs,maxlvls)
      real*4        pressure_pa(maxobs,maxlvls)      
      real*4        temp_c(maxobs,maxlvls), dewpoint_c(maxobs,maxlvls)
      real*4        dir_deg(maxobs,maxlvls),spd_mps(maxobs,maxlvls)
      real*4        sp_fill,levels_fill
      real*4        fill_t, fill_ws, fill_rh
      integer       fill_wd
      integer       wmoid(maxobs),nlvl(maxobs),sp_id
      integer       istatus
      integer       sta_id,sn_id,wd_id,ws_id,rh_id,temp_id,lev_id,ot_id
      integer       index_1(1), index_2(2),start(2),count(2)
      integer       start1(1),count1(1)
      integer       pi_len, sn_len, fn_len, obno
      character     stname(maxobs)*6, stationName*51, c_staid*6
      character     a9time_ob(maxobs)*9
      double precision d_timeobs

c     open data_file
      nf_status = NF_OPEN(filename(1:fn_len),NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error on NF_OPEN - open file: ',filename(1:fn_len)
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

c     verify nobs .lt. maxobs

      print *, 'Number of records in file: ',nobs
      if (nobs .gt. maxobs) then
        print *,
     1'WARNING: nobs is greater than maxobs: ',
     1'nobs / maxobs',nobs,' / ',maxobs
        nobs = maxobs
        print *, 'WARNING: Reading only the first ',
     1nobs,' records'
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

c     read dim level staNamLen
      nf_status = NF_INQ_DIMID(nf_fid,'staNamLen',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error finding dim staNamLen' 
        print *, 'set sn_len = 51'
        sn_len = 51
      else
        nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,sn_len)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'Error reading dim staNamLen'
          print *, 'set sn_len = 51'
          sn_len = 51
        endif
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

c     setup start and count
      start1(1) = 1
      count1(1) = nobs

c     read var elevation(recNum) -> staelev(maxobs)
      nf_status = NF_INQ_VARID(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var elevation'
      endif
      nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,
     1 start1,count1,staelev)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var elevation'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read var latitude(recNum -> stalat(maxobs)
      nf_status = NF_INQ_VARID(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var latitude'
      endif
      nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,
     1 start1,count1,stalat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var latitude'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read var longitude(recNum) -> stalon(maxobs)
      nf_status = NF_INQ_VARID(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var longitude'
      endif
      nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,
     1 start1,count1,stalon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var longitude'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read var levels(recNum,level) -> lvls_m(maxlvls,maxobs)
      start(1) = 1
      start(2) = 1
      count(1) = nlvls
      count(2) = nobs
      nf_status = NF_INQ_VARID(nf_fid,'levels',lev_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var levels'
      endif
      nf_status = NF_GET_VARA_REAL(nf_fid,lev_id,
     1 start,count,lvls_m)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading var levels'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif 

c     read _fillValue for levels
      nf_status = NF_INQ_ATTLEN(nf_fid, lev_id,'_FillValue',
     1                          levels_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading levels _FillValue attribute'
        levels_fill = 1.0e+38
      endif 

c     LW _FillValue variable returns 1.4012985E-45...hard wire for now
      levels_fill = 9.9999997E+37

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

      nf_status = NF_INQ_VARID(nf_fid,'stationName',sn_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var stationName'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'temperature',temp_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'finding var temperature'
        print *, 'Aborting read'
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

c     read _fillValue for temperature
      nf_status = NF_INQ_ATTLEN(nf_fid, temp_id,'_FillValue',
     1                          t_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading temperature _FillValue attribute'
        t_fill = 1.e+38
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

c     read _fillValue for relHumidity
      nf_status = NF_INQ_ATTLEN(nf_fid, rh_id,'_FillValue',
     1                          rh_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading relHumidity _FillValue attribute'
        rh_fill = 3.4028e+38
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

c     read _fillValue for windSpeed
      nf_status = NF_INQ_ATTLEN(nf_fid, ws_id,'_FillValue',
     1                          ws_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading windSpeed _FillValue attribute'
        ws_fill = 1.e+38
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

c     read _fillValue for windDir
      nf_status = NF_INQ_ATTLEN(nf_fid, wd_id,'_FillValue',
     1                          wd_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'reading windDir _FillValue attribute'
        wd_fill = 2147483647
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

      print *, 'LW nobs = ',nobs
      lev_set = 0
      do obno = 1, nobs

        lev_set = 0

        index_1(1) = obno
        index_2(2) = obno
        start(2) = obno
        start(1) = 1
        count(2) = 1

        do lno = 1, 10
          print *, 'LW level ',lno,'= ',lvls_m(lno,obno)
        enddo

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

        i4_tim = int(d_timeobs + 315619200)
        call make_fnam_lp(i4_tim,a9time_ob(obno),istatus)

c       read var stationName(obno,staNamLen) -> stationName
        count(1) = sn_len 
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
        count(1) = pi_len 
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
c       convert string to iwmostanum(maxobs) (cvt S to 0 and N to 1)
        call s_len(c_staid,lensta)
        do ic = 1, lensta 
          chr = ichar(c_staid(ic:ic))
          if ((chr.ge.48).and.(chr.le.57)) then 
            ichr = chr - 48
          else
            if (chr.eq.83) then
              ichr = 0
            else
              if (chr.eq.78) then
                ichr = 1
              else
                ichr = 0
              endif
            endif 
          endif
          if (ic.eq.1) ichr1 = ichr
          if (ic.eq.2) ichr2 = ichr
          if (ic.eq.3) ichr3 = ichr
          if (ic.eq.4) ichr4 = ichr
          if (ic.eq.5) ichr5 = ichr
          if (ic.eq.6) ichr6 = ichr
        enddo
        if (lensta .eq.1) wmoid(obno) = ichr1
        if (lensta .eq.2) 
     1    wmoid(obno) = ichr1*10 + ichr2
        if (lensta .eq.3) 
     1    wmoid(obno) = ichr1*100 + ichr2*10 + ichr3
        if (lensta .eq.4) 
     1    wmoid(obno) = ichr1*1000 + ichr2*100 + ichr3*10 + 
     1                  ichr4
        if (lensta .eq.5) 
     1    wmoid(obno) = ichr1*10000 + ichr2*1000 + ichr3*100 + 
     1                  ichr4*10 +ichr5
        if (lensta .eq.6) 
     1    wmoid(obno) = ichr1*100000 + ichr2*10000 + ichr3*1000 + 
     1                  ichr4*100 +ichr5*10 + ichr6

        write(6,*) 'LW c_staid wmoid >',c_staid(1:lensta),
     1'<  >',wmoid(obno),'<'
        
        print *,'LW obno nobs nlvls ',obno,nobs,nlvls

        lvl = 1
        do while ((lvl.le.nlvls).and.(lev_set.eq.0))

          print *, 'LW levels_fill lvls_m ',
     1levels_fill, lvls_m(lvl,obno)

          if(lvls_m(lvl,obno) .ne. levels_fill)then

            height_m(obno,lvl) = lvls_m(lvl,obno)
            index_2(1) = lvl
c           read stationPressure
            nf_status = NF_GET_VAR1_REAL(nf_fid,sp_id,index_1,stationP)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var stationPressure'
            endif 

            write(6,*) 'LW o l stationP ',obno,lvl,stationP

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

c           Convert temp_k to temp_c
            if ((temp_k .eq. -9999.).or.(temp_k .eq. t_fill)) then
              temp_c(obno,lvl) = r_missing_data
            else
              if ((temp_k .ge. 227.50).and.(temp_k .le. 328.15)) then
                temp_c(obno,lvl) = temp_k - 273.15
              else
                temp_c(obno,lvl) = r_missing_data
              endif
            endif

      write(6,*) 'LW temp_k temp_c',temp_k,'   ',temp_c(obno,lvl)

c           read var relHumidity(recNum,level) -> rh
            nf_status = NF_GET_VAR1_REAL(nf_fid,rh_id,index_2,rh_pct)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var relHumidity'
            endif 

c           Convert rh to dewpoint
            if ((rh_pct .eq. -9999.).or.(rh_pct .eq. rh_fill)) then
              dewpoint_c(obno,lvl) = r_missing_data
            else
              if ((rh_pct .ge. 0).and.(rh_pct .le. 100)
     1            .and.(temp_c(obno,lvl).ne.r_missing_data)) then
                dewpoint_c(obno,lvl)=dwpt(temp_c(obno,lvl), rh_pct) ! celsius
              else
                dewpoint_c(obno,lvl) = r_missing_data
              endif
            endif

      write(6,*) 'LW rh_pct dpt_c ',rh_pct,'   ',dewpoint_c(obno,lvl)

c           read var windSpeed(recNum,level) -> ws
            nf_status = NF_GET_VAR1_REAL(nf_fid,ws_id,index_2,ws)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var windSpeed'
            endif 

            if ((ws .eq. -9999.).or.(ws .eq. ws_fill)) then
              spd_mps(obno,lvl) = r_missing_data
            else
              if ((ws .ge. 0.).and.(ws .le. 150.)) then
                spd_mps(obno,lvl) = ws
              else
                spd_mps(obno,lvl) = r_missing_data
              endif
            endif

            write(6,*) 'LW spd_mps      ',spd_mps(obno,lvl)

c           read var windDir(recNum,level) -> dd(lvl,obno)
            nf_status = NF_GET_VAR1_REAL(nf_fid,wd_id,index_2,wd)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var windDir'
            endif 

            if ((wd .eq. -9999.).or.(wd .eq. wd_fill)) then
              dir_deg(obno,lvl) = r_missing_data
            else
              if ((ws .ge. 0).and.(ws .le. 360)) then
                dir_deg(obno,lvl) = wd
              else
                dir_deg(obno,lvl) = r_missing_data
              endif
            endif

            write(6,*) 'LW dir_deg      ',dir_deg(obno,lvl)

          else
            nlvl(obno) = lvl - 1
            lev_set = 1
            print *, 'LW lvl nlvl(obno) ',lvl, nlvl(obno)
          endif
          lvl = lvl + 1
        enddo
      enddo

      write(6,*) 'End of read_local_tower'

!     Final QC check
      call get_ibadflag(ibadflag,istatus)
      if(istatus .ne. 1)return

      if(iblank .gt. 0)then
        write(6,*)' Warning: number of UNK stanames = ',iblank
      endif


      return
      end
