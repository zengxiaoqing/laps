c
        subroutine get_local_towerobs(maxsta,maxlvls,                    ! I
     &                 i4time_sys,lun_out,
     &                 path_to_local_data,tower_format,ext,
     &                 itime_before,itime_after,
     &                 lat,lon,ni,nj,                                    ! I
     &                 nsta,                                             ! O
     &                 stalat_s,stalon_s,staelev_s,                      ! O
     &                 stname_s,                                         ! O
     &                 soilmoist_p,                                      ! O
     &                 istatus)

c       Read tower data from either RSA or NIMBUS netCDF CDLs
c
c.....  Input variables/arrays
c
        integer maxlvls ! raw/processed stations for SND file
        integer maxobs ! raw stations in NetCDF files
        integer maxsta ! processed stations for SND file 

        parameter (maxobs=10000) ! Raw stations in NetCDf files

        character*(*) path_to_local_data, tower_format, ext

        real    lat(ni,nj), lon(ni,nj)
c
c.....  Local variables/arrays
c
        double precision d_timeobs

!       Obs arrays (raw files)
        integer     nobs,nlvl(maxobs)
        real      lvls_m(maxlvls,maxobs)
        real      fill_stationP, fill_lvls, stationP
        integer   i4time
        character*51  stationName
        character*6   c_staId
	real      dd(maxlvls,maxobs), ff(maxlvls,maxobs)

        integer  i4time_ob_a(maxobs), before, after
c
c.....  Variables returned from 'read_local_tower'
c
        integer    nsnd_all ! combined # of obs over multiple files
	integer    wmoid(maxobs)
        real     stalat(maxobs),stalon(maxobs)
        real     staelev(maxobs)
        real     soilMoisture(maxobs)
        character  c5_staid(maxobs)*5, a9time_ob(maxobs)*9
        character  a9time*9
!       character  c8_obstype(maxobs)*8
        real     height_m(maxobs,maxlvls), pressure_mb(maxobs,maxlvls)
        real     temp_c(maxobs,maxlvls), dewpoint_c(maxobs,maxlvls)
        real     dir_deg(maxobs,maxlvls),spd_mps(maxobs,maxlvls)
	character  stname(maxobs)*6
        integer    tempQcFlag(maxobs,maxlvls)
        integer    prsQcFlag(maxobs,maxlvls)
        integer    rhQcFlag(maxobs,maxlvls)
        integer    wsQcFlag(maxobs,maxlvls)
        integer    wdQcFlag(maxobs,maxlvls)
        integer    smQcFlag(maxobs)

        logical    l_closest_time(maxobs), l_closest_time_i
c
c.....  Output arrays used by 'write_snd' 
c
        integer    nsnd_all_s ! combined # of obs over multiple files
        integer    nlvl_s(maxsta)
	integer  wmoid_s(maxsta)
        real     stalat_s(maxsta,maxlvls),stalon_s(maxsta,maxlvls)
        real     staelev_s(maxsta)
        character  c5_staid_s(maxsta)*5, a9time_ob_s(maxsta,maxlvls)*9
        character  c8_obstype_s(maxsta)*8
        real     height_m_s(maxsta,maxlvls)
        real     pressure_mb_s(maxsta,maxlvls)
        real     temp_c_s(maxsta,maxlvls)
        real     dewpoint_c_s(maxsta,maxlvls)      
        real     dir_deg_s(maxsta,maxlvls),spd_mps_s(maxsta,maxlvls)       
        real     soilmoist_p(maxsta)       
	character  stname_s(maxsta)*5

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

        write(6,*)' Tower format = ',tower_format

        do i4time_file = i4time_file_before, i4time_file_after, +3600       

            call s_len(path_to_local_data,len_path)

            if(tower_format .eq. 'WFO' .or. tower_format .eq. 'RSA')then
                filename13= cvt_i4time_wfo_fname13(i4time_file)
                if(len_path .lt. 1)goto590
 	        data_file = path_to_local_data(1:len_path)//'/'
     1                                                    //filename13
            else
                call make_fnam_lp(i4time_file,a9time,istatus)
                if(len_path .lt. 1)goto590
 	        data_file = path_to_local_data(1:len_path)//'/'
     1                                       //a9time//'0100o'       
            endif

            write(6,*)' LDAD tower file = ',trim(data_file)
            call s_len(data_file,lenf)
c
c.....  Call the read routine.
c
	    call read_local_tower(data_file,lenf,                 ! I 
     &         tower_format,                                      ! I
     &         maxobs, maxlvls,                                   ! I
     &         r_missing_data,                                    ! I
     &         nsnd_file, nlvl(ix), lvls_m(1,ix),                 ! O
     &         staelev(ix), stalat(ix), stalon(ix),               ! O
     &         soilMoisture(ix),                                  ! O
     &         temp_c(ix,1), dewpoint_c(ix,1),                    ! O
     &         height_m(ix,1),                                    ! O
     &         dir_deg(ix,1), spd_mps(ix,1),                      ! O
c    &         pressure_mb(ix,1), prsQcFlag(ix,1),                ! O
     &         tempQcFlag(ix,1), rhQcFlag(ix,1),                  ! O
     &         wsQcFlag(ix,1), wdQcFlag(ix,1),                    ! O
     &         smQcFlag(ix),                                      ! O
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

        if(nsnd_all .gt. maxobs)then
            write(6,*)' ERROR: nsnd_all > maxobs ',nsnd_all,maxobs
            nsta = 0
            istatus = 0
            return
        endif
c
c.....  Post process the soundings.....
c
        print *
	print *,'  Appending SND file, # of obs (in grid) = ',nsnd_all       

        pressure_mb_s = r_missing_data

!       Flag those reports that are the closest times for the station
        do i = 1,nsnd_all
            l_closest_time(i) = l_closest_time_i(wmoid,a9time_ob
     1                                          ,nsnd_all,i,i4time_sys       
     1                                          ,istatus)       
        enddo ! i

!       Transfer arrays (with various QC steps)

        nsta = 0
        do i = 1,nsnd_all
            if(nlvl(i) .gt. 0 .and. nlvl(i) .le. maxlvls 
     1                        .and. l_closest_time(i)     )then 

!               Valid sounding - use for output
                nsta = nsta + 1
                write(6,601)nsta,i,stname(i),wmoid(i),a9time_ob(i)
 601            format(' Valid sounding - transferring to output arrays'      
     1                ,2i5,2x,a7,2x,i10,2x,a10)
                stalat_s(nsta,:) = stalat(i)           
                stalon_s(nsta,:) = stalon(i)           
                staelev_s(nsta) = staelev(i)           
                stname_s(nsta) = stname(i)(1:5)           

                wmoid_s(nsta) = wmoid(i)           
                a9time_ob_s(nsta,:) = a9time_ob(i)           
                c8_obstype_s(nsta) = 'TOWER   '

!               Single level data
                soilmoist_p(nsta) = soilMoisture(i)

                nlvl_s(nsta) = nlvl(i)

!               QC the levels 
                do il = nlvl(i),1,-1
                    if(     lvls_m(il,i) .ge. 1e10 
!    1                 .or. lvls_m(il,i) .le. 0.          
     1                                             )then
                        write(6,*)' ERROR: invalid lvls_m',i,wmoid(i)
     1                           ,il,lvls_m(il,i)
                        nlvl_s(nsta) = il-1
                    else
                        if(il .ge. 2)then
                            if(lvls_m(il,i) .le. lvls_m(il-1,i))then
                                write(6,*)' ERROR: levels out of order'
     1                                   ,i,lvls_m(il-1,i),lvls_m(il,i)      
                                go to 1500
                            endif 
                        endif
                    endif
                enddo ! il

                do il = 1,nlvl_s(nsta)
                    height_m_s(nsta,il) = lvls_m(il,i) + staelev_s(nsta)       

                    if(iqc_rsa(tempQcFlag(nsta,il)) .ne. -1)then
                        temp_c_s(nsta,il) = temp_c(i,il)           
                    else
                        temp_c_s(nsta,il) = r_missing_data
                    endif

                    if(iqc_rsa(rhQcFlag(nsta,il)) .ne. -1)then
                        dewpoint_c_s(nsta,il) = dewpoint_c(i,il)           
                    else
                        dewpoint_c_s(nsta,il) = r_missing_data
                    endif

                    if(iqc_rsa(wdQcFlag(nsta,il)) .ne. -1 .and.
     1                 iqc_rsa(wsQcFlag(nsta,il)) .ne. -1      )then
                        dir_deg_s(nsta,il) = dir_deg(i,il)           
                        spd_mps_s(nsta,il) = spd_mps(i,il)           
                    else
                        dir_deg_s(nsta,il) = r_missing_data
                        spd_mps_s(nsta,il) = r_missing_data
                    endif
                enddo

                go to 1600 

 1500           write(6,*)' Sounding rejected: ' ,nsta

                nsta = nsta - 1

 1600           continue ! Normal status

            endif
        enddo ! i
c
c.....  Call the routine to write the SND file.
c
        if(ext(1:3) .eq. 'snd')then

            if(nsta .gt. 0)then
                call open_ext(lun_out,i4time_sys,'snd',istatus)
            endif

            call write_snd(lun_out                               ! I
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

        elseif(ext(1:3) .eq. 'lso')then
            write(6,*)
     1       ' LSO option, pass data back instead of writing data out'       

        else
            write(6,*)' Error - unknown ext in get_local_towerobs',ext
            nsta = 0
            istatus = 0
            return

        endif
c
 990    continue               ! no data available
        istatus = 0
        print *,' No data available from GET_LOCAL_TOWEROBS'
        return
c
        end

         subroutine read_local_tower(filename,fn_len,             ! I 
     &         tower_format,                                      ! I
     &         maxobs, maxlvls,                                   ! I
     &         r_missing_data,                                    ! I
     &         nobs, nlvl, lvls_m,                                ! O
     &         staelev, stalat, stalon,                           ! O
     &         soilMoisture,                                      ! O
     &         temp_c, dewpoint_c, height_m,                      ! O
     &         dir_deg, spd_mps,                                  ! O
c    &         pressure_pa, prsQcFlag,                            ! O 
     &         tempQcFlag, rhQcFlag, wsQcFlag, wdQcFlag,          ! O
     &         smQcFlag,                                          ! O
     &         a9time_ob, stname, wmoid,                          ! O
     &         istatus)                                           ! O

      include 'netcdf.inc'

      character*(*) filename, tower_format 
      integer       maxobs ! raw stations for SND file
      integer       maxlvls ! raw/processed stations for SND file
      real        r_missing_data 
      integer       nobs,nlvls,lev_set
      real        lvls_m(maxlvls,maxobs)
      real        staelev(maxobs)
      real        stalat(maxobs),stalon(maxobs)
      real        soilMoisture(maxobs)
      real        dd(maxlvls,maxobs), ff(maxlvls,maxobs)
      real        temp_k, rh_pct,stationP,ws,wd
      integer       tempQcFlag(maxobs,maxlvls)
      integer       prsQcFlag(maxobs,maxlvls)
      integer       rhQcFlag(maxobs,maxlvls)
      integer       wsQcFlag(maxobs,maxlvls)
      integer       wdQcFlag(maxobs,maxlvls)
      integer       smQcFlag(maxobs)
      integer       tempQF, prsQF, rhQF, wsQF, wdQF
      real        height_m(maxobs,maxlvls)
      real        pressure_pa(maxobs,maxlvls)      
      real        temp_c(maxobs,maxlvls), dewpoint_c(maxobs,maxlvls)
      real        dir_deg(maxobs,maxlvls),spd_mps(maxobs,maxlvls)
      real        sp_fill,levels_fill
      real        fill_t, fill_ws, fill_rh
      integer       fill_wd
      integer       wmoid(maxobs),nlvl(maxobs),sp_id
      integer       istatus
      integer       sta_id,sn_id,wd_id,ws_id,rh_id,temp_id,lev_id,ot_id
      integer       tempQc_id,prsQc_id,rhQc_id,wsQc_id,wdQc_id
      integer       index_1(1), index_2(2),start(2),count(2)
      integer       start1(1),count1(1)
      integer       pi_len, sn_len, fn_len, obno
      character     stname(maxobs)*6, stationName*51, c_staid*6
      character     a9time_ob(maxobs)*9, rh_var*30
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

c LW
      print *, 'nlvls / nobs ',nlvls,' / ',nobs

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

c     read var soilMoisture(recNum) -> soilMoisture(maxobs)
      nf_status = NF_INQ_VARID(nf_fid,'soilMoisture',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Warning: could not find var soilMoisture'
      else
        nf_status = NF_GET_VARA_REAL(nf_fid,nf_vid,
     1                               start1,count1,soilMoisture)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading var soilMoisture'
          print *, 'Aborting read'
          nf_status = NF_CLOSE(nf_fid)
          istatus = 0
          return
        endif 
      endif
      
c     read dim smQcFlag -> smQcFlag
      nf_status = NF_INQ_DIMID(nf_fid,'smQcFlag',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'could not find dim smQcFlag' 
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,smQcFlag)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'could not read dim smQcFlag'
      endif

c LW
      print *, 'b4 read levels: nlvls / nobs ',nlvls,' / ',nobs

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

      nf_status = NF_INQ_VARID(nf_fid,'tempQcFlag',tempQc_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Warning: could not find var tempQcFlag'
!       nf_status = NF_CLOSE(nf_fid)
!       istatus = 0
!       return
        istat_tempQcFlag = 0
      else
        istat_tempQcFlag = 1
      endif

c     read _fillValue for temperature
      nf_status = NF_INQ_ATTLEN(nf_fid, temp_id,'_FillValue',
     1                          t_fill)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'use guess for filling temperature _FillValue attribute'
        t_fill = 1.e+38
      endif 

      nf_status = NF_INQ_VARID(nf_fid,'stationPressure',sp_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Warning: could not find var stationPressure'
!       nf_status = NF_CLOSE(nf_fid)
!       istatus = 0
!       return
        istat_stationPressure = 0
      else
        istat_stationPressure = 1
      endif

      if(istat_stationPressure .eq. 1)then
        nf_status = NF_INQ_VARID(nf_fid,'prsQcFlag',spQc_id)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'finding var spQcFlag'
          print *, 'Aborting read'
          nf_status = NF_CLOSE(nf_fid)
          istatus = 0
          return
        endif

c       read _fillValue for stationPressure
        nf_status = NF_INQ_ATTLEN(nf_fid, sp_id,'_FillValue',
     1                            sp_fill)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading stationPressure _FillValue attribute'
          sp_fill = 3.4028e+38
        endif 
      endif

      if(.true.)then ! NIMBUS
          rh_var = 'relativeHumidity'
      else ! RSA
          rh_var = 'relHumidity'
      endif

      nf_status = NF_INQ_VARID(nf_fid,trim(rh_var),rh_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Warning: could not find var ',rh_var     
!       nf_status = NF_CLOSE(nf_fid)
!       istatus = 0
!       return
        istat_relHumidity = 0
      else
        istat_relHumidity = 1
      endif

      if(istat_relHumidity .eq. 1)then
        nf_status = NF_INQ_VARID(nf_fid,'rhQcFlag',rhQc_id)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'finding var rhQcFlag'
!         print *, 'Aborting read'
!         nf_status = NF_CLOSE(nf_fid)
!         istatus = 0
!         return
        endif

c       read _fillValue for relHumidity
        nf_status = NF_INQ_ATTLEN(nf_fid, rh_id,'_FillValue',
     1                            rh_fill)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'reading relHumidity _FillValue attribute'
          rh_fill = 3.4028e+38
        endif 
      endif

      nf_status = NF_INQ_VARID(nf_fid,'windSpeed',ws_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Error, could not find var windSpeed',NF_NOERR,nf_status
        print *, 'Aborting read ',nf_fid
        nf_status = NF_CLOSE(nf_fid)
        istatus = 0
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'wsQcFlag',wsQc_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'could not find var wsQcFlag'
!       print *, 'Aborting read'
!       nf_status = NF_CLOSE(nf_fid)
!       istatus = 0
!       return
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

      nf_status = NF_INQ_VARID(nf_fid,'wdQcFlag',wdQc_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'could not find var wdQcFlag'
!       print *, 'Aborting read'
!       nf_status = NF_CLOSE(nf_fid)
!       istatus = 0
!       return
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
        print *,'could not find observationTime'
        print *, 'try for timeObs '

        nf_status = NF_INQ_VARID(nf_fid,'timeObs',ot_id)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'could not find timeObs'
          print *, 'Aborting routine '
          nf_status = NF_CLOSE(nf_fid)
          istatus = 0
          return
        endif

      endif

      print *, 'LW nobs = ',nobs
      lev_set = 0
      do obno = 1, nobs

        if(obno .le. 50
     1     .or. (obno .eq. (obno/100 * 100) )     
     1     .or. (obno .ge. 1700 .and. obno .le. 1800)      
     1                                                 )then
            id = 1
        else
            id = 0
        endif

        if(id.eq.1)write(6,*)
        if(id.eq.1)write(6,*)' SA: obno,stalat,stalon,staelev'
     1                 ,obno,stalat(obno),stalon(obno),staelev(obno)       

        lev_set = 0

        index_1(1) = obno
        index_2(2) = obno
        start(2) = obno
        start(1) = 1
        count(2) = 1

        do lno = 1, 10
          if(id.eq.1)print *, 'LW level ',lno,'= ',lvls_m(lno,obno)
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

c       truncate stname(obno)  
        call left_justify(stationName)
        call remove_blanks(stationName)
        stname(obno) = stationName(1:6)

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
        call s_len(c_staid,len_sta)
        stname(obno) = '      '
        stname(obno)(1:len_sta) = c_staid(1:len_sta)          

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

        if(id.eq.1)write(6,*) 'LW c_staid wmoid >',c_staid(1:lensta),
     1'<  >',wmoid(obno),'<'
        
        if(id.eq.1)print *,'LW obno nobs nlvls a9time staname'
     1            ,obno,nobs,nlvls,a9time_ob(obno),' ',stationName(1:20)       
     1            ,' ',stname(obno)

        lvl = 1
        do while ((lvl.le.nlvls).and.(lev_set.eq.0))

          if(id.eq.1)print *, 'LW levels_fill lvls_m ',
     1               levels_fill, lvls_m(lvl,obno)

          if(lvls_m(lvl,obno) .ne. levels_fill)then

            height_m(obno,lvl) = lvls_m(lvl,obno)
            index_2(1) = lvl
c           read stationPressure
            if(istat_stationPressure .eq. 1)then
              nf_status = NF_GET_VAR1_REAL(nf_fid,sp_id,index_1
     1                                                 ,stationP)
              if(nf_status.ne.NF_NOERR) then
                print *, NF_STRERROR(nf_status)
                print *,'error reading var stationPressure'
              endif 

c             read prsQcFlag
              nf_status=NF_GET_VAR1_REAL(nf_fid,spQc_id,index_1,prsQF)    
              if(nf_status.ne.NF_NOERR) then
                print *, NF_STRERROR(nf_status)
                print *,'reading var prsQcFlag'
              endif 

c             write prsQcFlag
              prsQcFlag(obno,lvl) = prsQF

              if(id.eq.1)write(6,*) 'LW o l stationP ',obno,lvl,stationP       

c             check stationPressure for _FillValue
              if (stationP .eq. sp_fill) stationP = r_missing_data
              if (lvl.eq.1) then
                pressure_pa(obno,lvl) = stationP
              else
                pressure_pa(obno,lvl) = r_missing_data
              endif
            else
                pressure_pa(obno,lvl) = r_missing_data
            endif

c           read var temperature(recNum,lvl) -> temp
            nf_status = NF_GET_VAR1_REAL(nf_fid,temp_id,index_2,temp_k)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var temperature'
            endif 

c           read var tempQcFlag(recNum,lvl) -> tempQcFlag
            if(istat_tempQcFlag .eq. 1)then
                nf_status = NF_GET_VAR1_REAL(nf_fid,tempQc_id,index_2,
     1                                       tempQF)
                if(nf_status.ne.NF_NOERR) then
                  print *, NF_STRERROR(nf_status)
                  print *,'reading var tempQcflag'
                endif 

c               write tempQcFlag
                tempQcFlag(obno,lvl) = tempQF
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

            if(id.eq.1)
     1      write(6,*) 'LW temp_k temp_c',temp_k,'   ',temp_c(obno,lvl)

c           read var relHumidity(recNum,level) -> rh
            nf_status = NF_GET_VAR1_REAL(nf_fid,rh_id,index_2,rh_pct)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var relHumidity'
            endif 

c           read var rhQcFlag(recNum,level) -> rhQcFlag
            nf_status=NF_GET_VAR1_REAL(nf_fid,rhQc_id,index_2,rhQF)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var rhQcFlag'
            endif 

c           write rhQcFlag
            rhQcFlag(obno,lvl) = rhQF

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

            if(id.eq.1)write(6,*)'LW rh_pct dpt_c ',rh_pct,'   '
     1                ,dewpoint_c(obno,lvl)       

c           read var windSpeed(recNum,level) -> ws
            nf_status = NF_GET_VAR1_REAL(nf_fid,ws_id,index_2,ws)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var windSpeed'
            endif 

c           read var wsQcFlag(recNum,level) -> wsQcFlag
            nf_status=NF_GET_VAR1_REAL(nf_fid,wsQc_id,index_2,wsQF)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var wsQcFlag'
            endif 

c           write wsQcFlag
            wsQcFlag(obno,lvl) = wsQF

            if ((ws .eq. -9999.).or.(ws .eq. ws_fill)) then
              spd_mps(obno,lvl) = r_missing_data
            else
              if ((ws .ge. 0.).and.(ws .le. 150.)) then
                spd_mps(obno,lvl) = ws
              else
                spd_mps(obno,lvl) = r_missing_data
              endif
            endif

            if(id.eq.1)write(6,*) 'LW spd_mps      ',spd_mps(obno,lvl)

c           read var windDir(recNum,level) -> dd(lvl,obno)
            nf_status = NF_GET_VAR1_REAL(nf_fid,wd_id,index_2,wd)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var windDir'
            endif 

c           read var wdQcFlag(recNum,level) -> wdQcFlag(lvl,obno)
            nf_status=NF_GET_VAR1_REAL(nf_fid,wdQc_id,index_2,wdQF)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'reading var wdQcFlag'
            endif 

c           write wdQcFlag
            wdQcFlag(obno,lvl) = wdQF

            if ((wd .eq. -9999.).or.(wd .eq. wd_fill)) then
              dir_deg(obno,lvl) = r_missing_data
            else
              if ((ws .ge. 0).and.(ws .le. 360)) then
                dir_deg(obno,lvl) = wd
              else
                dir_deg(obno,lvl) = r_missing_data
              endif
            endif

            if(id.eq.1)write(6,*) 'LW dir_deg      ',dir_deg(obno,lvl)

          else
            nlvl(obno) = lvl - 1
            lev_set = 1
            if(id.eq.1)print *, 'LW lvl nlvl(obno) ',lvl, nlvl(obno)
          endif
          lvl = lvl + 1
        enddo
      enddo

      write(6,*) 'End of read_local_tower'

      return
      end

      subroutine remove_blanks(string)

      character*(*)string

      len1 = len(string)
      call s_len2(string,len2)

      len_loop = min(len2,len1-1)

      do i = 1,len_loop
          if(string(i:i) .eq. ' ')then ! blank found
!             ip1 = i+1
!             len1m1 = len1-1
              string(i:len1-1) = string(i+1:len1)
              string(len1:len1) = ' '
          endif
      enddo ! i

      return
      end


