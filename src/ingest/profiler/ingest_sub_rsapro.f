
        subroutine ingest_rsapro(i4time_sys,NX_L,NY_L,istatus)

        integer cdfid,status,MAX_PROFILES,MAX_LEVELS,file_n_prof       

        parameter (MAX_PROFILES = 1000)
	parameter (MAX_LEVELS = 300)

        real ht_out(max_levels)
        real di_out(max_levels)
        real sp_out(max_levels)

        character*6 prof_name(MAX_PROFILES)

        integer i4_mid_window_pr(MAX_PROFILES)
        integer wmo_id(MAX_PROFILES)

        real lat_pr(MAX_PROFILES)
        real lon_pr(MAX_PROFILES)
        real elev_m_pr(MAX_PROFILES)
        real n_lvls_pr(MAX_PROFILES)
        real ht_m_pr(MAX_PROFILES,MAX_LEVELS)
        real dir_dg_pr(MAX_PROFILES,MAX_LEVELS)
        real spd_ms_pr(MAX_PROFILES,MAX_LEVELS)
        real u_std_ms_pr(MAX_PROFILES,MAX_LEVELS)
        real v_std_ms_pr(MAX_PROFILES,MAX_LEVELS)

        integer bad,missing, start(2), count(2), staNamLen
        integer start_time(1), count_time(1)
        parameter (bad = 12)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'G','B','M'/
        integer*4 byte_to_i4

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec
        integer error_code
        data error_code/1/
        logical l_in_box
        data l_in_box/.true./

        integer varid
        include 'netcdf.inc'
        character*(MAXNCNAM) dimname 

        character*13 filename13,outfile
        character*9 asc9_tim,a9time_ob

        character*31    ext
        integer*4       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*9 a9_timeObs
        integer*4 timeObs

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           return
        endif
 
        r_mspkt = .518

        call get_latlon_perimeter(NX_L,NY_L,1.0
     1                           ,lat,lon,topo
     1                           ,rnorth,south,east,west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS perimeter'
            return
        endif

        write(6,*)' return from ingest_rsapro'
        return

        outfile = filename13(i4time_sys,'pro')
        asc9_tim = outfile(1:9)

!       dir_in = path_to_raw_blpprofiler

        c_vars_req = 'path_to_raw_blpprofiler'
        call get_static_info(c_vars_req,c_values_req,1,istatus)
        if(istatus .eq. 1)then
            write(6,*)c_vars_req(1:30),' = ',c_values_req
            dir_in = c_values_req
        else
            write(6,*)' Error getting ',c_vars_req
            return
        endif

        call s_len(dir_in,len_dir_in)

C       READ IN THE RAW PROFILER DATA

!       call cvt_wfo_fname()
!       fnam_in = dir_in(1:len_dir_in)//asc9_tim//'0100o'
        call s_len(fnam_in,len_fnam_in)
        write(6,*)fnam_in(1:len_fnam_in)

!       call read_prof_rsa(fnam_in(1:len_fnam_in)                      ! I
!     1                   ,MAX_PROFILES,MAX_LEVELS                     ! I
!     1                   ,n_profiles                                  ! O
!     1                   ,n_lvls_pr                                   ! O
!     1                   ,prof_name,wmo_id                            ! O
!     1                   ,lat_pr,lon_pr,elev_m_pr                     ! O
!     1                   ,ht_m_pr,di_dg_pr,sp_ms_pr                   ! O
!     1                   ,u_std_ms_pr,v_std_ms_pr                     ! O
!     1                   ,i4_mid_window_pr,istatus)                   ! O
        istatus = 0
        if(istatus.ne.1)then
            write(6,*)' Warning: bad status on read_prof_rsa ',status       
            return
        endif
C
C       Open intermediate output file.
C
        ext = 'pro'
        call open_lapsprd_file_append(1,i4time_sys,ext,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error opening product file',ext
            return
        endif

        do i_pr = 1,n_profiles

          write(6,*)
          write(6,*)' Looping for profiler ',prof_name(i_pr)
C
C
C       Get the surface pressure.  This time we'll use the
C       5-character site name (plus a terminating blank) to select the station.
C
          CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0,'pressure',2,prs
     1                      ,status)
          if(status.ne.0)then
            if(status .eq. -3)then
                write(6,*)prof_name(i_pr),' not found'
                goto 900
            else
                write(6,*)' Warning: bad pressure read ',status
                return
            endif
          endif
          write(6,*)
          write(6,*)prof_name(i_pr),' pressure ',prs

          rlat = lat_pr(i_pr)
          rlon = lon_pr(i_pr)

          if(rlat .le. rnorth .and. rlat .ge. south .and.
     1       rlon .ge. west   .and. rlon .le. east            )then

            write(6,*)prof_name(i_pr),' is in box'

            elev = elev_m_pr(i_pr)

            CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0
     1                         ,'windSpeedSfc',2,sp_sfc,istatus1)
            CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0
     1                         ,'windDirSfc',2,di_sfc,istatus2)

            if(abs(sp_sfc) .gt. 500.)istatus1 = 1
            if(abs(di_sfc) .gt. 500.)istatus2 = 1

            if(istatus1 .eq. 0 .and. istatus2 .eq. 0)then
                n_good_sfc = 1
            else
                n_good_sfc = 0
            endif

C
C           Get the array of profiler winds for the profiler station
C           at Platteville.  For this call, we'll use the WMO
C           identifier to tell the subroutine what station we want.
C
            CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0,'uComponent',2,
     1                                             u,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad uComponent read ',status
                return
            endif

            CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0,'vComponent',2,
     1                              v,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad vComponent read ',status
                return
            endif
C
C           Get the associated quality control flags.
C
            CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0,'uvQualityCode'
     $                     ,0,qc_flag,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad qualityCode read ',status
                return
            endif
C
C           Get the associated levels for this profiler
C
            CALL PROF_CDF_READ(cdfid,prof_name(i_pr),0,'levels',2,
     $                     level,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad level read ',status
                return
            endif
C
            write(6,*)prof_name(i_pr)

!           Convert u_std, v_std to rms

            write(6,*)

            i4time_ob = i4_mid_window_pr(i_pr)

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)

            write(*,401)wmo_id(i_pr),n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(i_pr),a9time_ob,'PROFILER'
            write(1,401)wmo_id(i_pr),n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(i_pr),a9time_ob,'PROFILER'
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

            rms = 1.0

            if(n_good_sfc .eq. 1)then
!               write surface winds as first level
                write(1,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
                write(6,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
            endif

            n_good_levels = n_lvls_pr(i_pr)

            do i = 1,n_good_levels
!           do i = n_good_levels, 1, -1
                write(1,301,err=303)ht_out(i),di_out(i),sp_out(i),rms
                write(6,301,err=303)ht_out(i),di_out(i),sp_out(i),rms
301             format(1x,f6.0,f6.0,2f6.1)
303             continue
            enddo ! i

          else !
             write(6,*)prof_name(i_pr),' is outside of domain'

          endif ! l_in_box

          write(6,*)prof_name(i_pr)

900     enddo ! i_pr

        close(1)

        return
        end


