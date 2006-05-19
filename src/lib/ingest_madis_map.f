
        subroutine ingest_madis_map(i4time_sys,NX_L,NY_L,ext,istatus)

        integer status,MAX_PROFILES,MAX_LEVELS

        parameter (MAX_SUBDIRS = 1)

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec
        logical l_exist

        character*13 a13_time,filename13,cvt_i4time_wfo_fname13,outfile       
        character*9 asc9_tim,a9time_ob

        character*(*)   ext
        integer*4       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*9 a9_timeObs
        integer*4 timeObs

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

        write(6,*)' Start ingest_madis_map, ext = ',ext

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

        outfile = filename13(i4time_sys,ext)
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

        i4_prof_window = 1800 ! could be reset to laps_cycle_time

C       READ IN THE RAW PROFILER/RASS DATA
        a13_time = cvt_i4time_wfo_fname13(i4time_sys)
        fnam_in = dir_in(1:len_dir_in)//a13_time
        call s_len(fnam_in,len_fnam_in)
        write(6,*)' file = ',fnam_in(1:len_fnam_in)

        inquire(file=fnam_in(1:len_fnam_in),exist=l_exist)

        if(l_exist)then

            istatus = 0

            call read_ldad_prof(i4time_sys,i4_prof_window              ! I
     1                                    ,NX_L,NY_L                   ! I
     1                                    ,ext                         ! I
     1                                    ,fnam_in(1:len_fnam_in)      ! I
     1                                    ,istatus)                    ! O

            if(istatus.ne.1)then
                write(6,*)' Warning: bad status in ingest_madis_map'
                goto980
            endif

 980        continue

        else
                write(6,*)' Warning: cannot find file '
     1                   ,fnam_in(1:len_fnam_in)

        endif ! file exists

        write(6,*)

        return
        end


