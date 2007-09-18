
        subroutine ingest_rsapro(i4time_sys,NX_L,NY_L,lun_out,istatus)       

        integer cdfid,status,MAX_PROFILES,MAX_LEVELS,file_n_prof       

!       parameter (MAX_PROFILES = 1000)
!       parameter (MAX_LEVELS = 300)
        parameter (MAX_SUBDIRS = 3)

!       real ht_out(max_levels)
!       real di_out(max_levels)
!       real sp_out(max_levels)

!       character*6 prof_name(MAX_PROFILES)

        character*255 prof_subdirs(MAX_SUBDIRS)

!       integer i4_mid_window_pr(MAX_PROFILES)
!       integer wmo_id(MAX_PROFILES)

!       real lat_pr(MAX_PROFILES)
!       real lon_pr(MAX_PROFILES)
!       real elev_m_pr(MAX_PROFILES)
!       real n_lvls_pr(MAX_PROFILES)
!       real ht_m_pr(MAX_PROFILES,MAX_LEVELS)
!       real dir_dg_pr(MAX_PROFILES,MAX_LEVELS)
!       real spd_ms_pr(MAX_PROFILES,MAX_LEVELS)
!       real u_std_ms_pr(MAX_PROFILES,MAX_LEVELS)
!       real v_std_ms_pr(MAX_PROFILES,MAX_LEVELS)

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec
        logical l_exist

        character*13 a13_time,filename13,cvt_i4time_wfo_fname13,outfile       
        character*9 asc9_tim,a9time_ob

        character*31    ext
        integer       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*9 a9_timeObs
        integer timeObs

        real lat(NX_L,NY_L),lon(NX_L,NY_L)
        real topo(NX_L,NY_L)

        lun_out = 1

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

        prof_subdirs(1) = '50mhz'
        prof_subdirs(2) = '915mhz'
        prof_subdirs(3) = 'minisodar'

        ext = 'pro'

        i4_prof_window = 1800 ! could be reset to laps_cycle_time

        do idir = 1,MAX_SUBDIRS
            call s_len(prof_subdirs(idir),len_subdir)
 
C           READ IN THE RAW PROFILER DATA
            a13_time = cvt_i4time_wfo_fname13(i4time_sys)
            fnam_in = dir_in(1:len_dir_in)
     1                //prof_subdirs(idir)(1:len_subdir)
     1                //'/netCDF/'//a13_time
            call s_len(fnam_in,len_fnam_in)
            write(6,*)' file = ',fnam_in(1:len_fnam_in)

            inquire(file=fnam_in(1:len_fnam_in),exist=l_exist)

            if(l_exist)then

                istatus = 0

                n_good_obs = 0

                if(idir .eq. 1 .or. idir .eq. 2)then
                    call read_ldad_prof(i4time_sys,i4_prof_window      ! I
     1                                    ,NX_L,NY_L                   ! I
     1                                    ,ext                         ! I
     1                                    ,lun_out                     ! I
     1                                    ,fnam_in(1:len_fnam_in)      ! I
     1                                    ,n_good_obs                  ! I/O
     1                                    ,istatus)                    ! O

                    write(6,*)' n_good_obs = ',n_good_obs

                else ! idir .eq. 3
                    i4time_earliest = i4time_sys - i4_prof_window
                    i4time_latest   = i4time_sys + i4_prof_window

                    call get_sodar_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L        ! I
     +                   ,i4time_earliest,i4time_latest                ! I
     +                   ,fnam_in(1:len_fnam_in)                       ! I
     +                   ,lun_out                                      ! I
     +                   ,istatus)                                     ! O

                endif ! idir

                if(istatus.ne.1)then
                    write(6,*)' Warning: bad status on '
     1                       ,prof_subdirs(idir),istatus           
                endif

            else
                write(6,*)' Warning: cannot find file '
     1                   ,fnam_in(1:len_fnam_in)

            endif ! file exists

            write(6,*)

        enddo ! idir

        return
        end


