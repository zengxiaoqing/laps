cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 



        subroutine ingest_pro(i4time_sys,NX_L,NY_L,istatus)

C       Michael Barth           12-Aug-1993
C       Steve Albers               Nov-1993         Reworked for LAPS ingest
C                                  Oct-1994         Improve QC
C       Steve Albers               Sep-1996         WFO database compatability
!       Ken Dritz                3-Jul-1997  Added NX_L, NY_L as dummy
!                                            arguments.
!       Ken Dritz                3-Jul-1997  Changed include of lapsparms.for
!                                            to laps_static_parameters.inc.
!       Ken Dritz                3-Jul-1997  Added call to get_r_missing_data.
!       Ken Dritz                8-Jul-1997  Replaced LAPS_DOMAIN_FILE with
!                                            'nest7grid' and removed include
!                                            of laps_static_parameters.inc.
C
C       This file shows examples of how the use PROF_CDF subroutines to read
C       WPDN 60-minute RASS data in netCDF files.
C
C       NOTE: Profiler winds are written out in KNOTS, and are sorted by HEIGHT

        integer cdfid,status,i,j,max_levels,max_levels_out,n_profilers
     1         ,file_n_prof
        parameter (max_levels = 72)
        parameter (max_levels_out = 72)
! number of profilers upped to 100 from 33 to accomodate RSA project
!       parameter (n_profilers = 33)
        parameter (n_profilers = 100)

        real u(max_levels),v(max_levels),prs
        real ht_out(max_levels_out),di_out(max_levels_out)
     1                             ,sp_out(max_levels_out)

        character*1 c1_qc_flag(max_levels)        ! for /public
        integer*4 i4_qc_flag(max_levels)          ! for WFO

        real*4 level(max_levels)
        integer good,bad,missing, start(2), count(2), staNamLen
        parameter (good = 0)
        parameter (bad = 12)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'G','B','M'/
        integer*4 byte_to_i4

        character*100 fnam_in
        character*80 dir_in
        character*255 c_filespec
        integer max_files

        parameter(max_files = 3000)
        character*255 c_filenames(max_files)

        integer wsmr_wmo_id
        data wsmr_wmo_id/74533/
        integer error_code
        data error_code/1/
        logical l_in_box
        data l_in_box/.true./
        character*8 c8_project
        character*5 c5_data_interval
        character*1 c1_char

        integer varid
        integer n_levels 
        include 'netcdf.inc'
        character*(MAXNCNAM) dimname 
C
C       Set error handling mode.  Note that you don't have to do this, if this
C       call isn't made, default error processing will occur:
C
C       ERROR_CODE              Meaning
C
C       0                       Return status codes -- THIS IS THE DEFAULT.
C       1                       Return status codes and write an error message
C                               to standard output (SYS$OUTPUT on VMS).
C       2                       Write an error message to standard output and
C                               exit the program.
C

        character*13 filename13,outfile,asc13_tim,fname9_to_wfo_fname13       
        character*9 asc9_tim,a9time_ob

        character*31    ext
        integer*4       len_dir_in

        character*40 c_vars_req
        character*100 c_values_req

        character*6 prof_name(n_profilers)
        character*9 a9_timeObs
        double precision timeObs

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           return
        endif

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write(6,*)'Error getting laps_cycle_time'
           return
        else
           write(6,*)'laps_cycle_time = ',laps_cycle_time
        endif

        r_mspkt = .518

        call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat,lon,    
     1                  topo,1.0,rnorth,south,east,west,istatus)

        CALL PROF_CDF_SET_ERROR(error_code,status)
        if(status.ne.0)then
                write(6,*)'bad set_error ',status
                return
        endif

        outfile = filename13(i4time_sys,'pro')
        asc9_tim = outfile(1:9)
C
C       Open a 60-minute RASS netCDF file for 20:00:00.00 on Julian date 217,
C       1993.  Both 6-minute and 60-minute files have the same filename
C       convention (you're supposed to use different directories to hold the
C       different resolutions).
C
C       YYJJJHHMMhhmmO
C
C       YY   = last 2 digits of year
C       JJJ  = Julian date
C       HHMM = hour, minute (UTC) of the data
C       hhmm = hour, minute (UTC) of the data
C              (One of these is supposed to be observation time, one receipt
C               time, but both are identical for RASS files created from
C               Demonstration Division tapes.)
C       O    = ASCII "Oh" character, not "zero" -- stands for observation.
C
C       For VMS systems, you must explicitly put in a period after the oh if
C       you don't want a file extension, thus the string used in the open call.
C

        c_vars_req = 'path_to_raw_profiler'
        call get_static_info(c_vars_req,c_values_req,1,istatus)
        if(istatus .eq. 1)then
            write(6,*)c_vars_req(1:30),' = ',c_values_req
            dir_in = c_values_req
        else
            write(6,*)' Error getting ',c_vars_req
            return
        endif

        call s_len(dir_in,len_dir_in) ! Slash should be on the end of dir_in

!       Determine whether we are using /public or WFO Advanced filenames...
        if(.false.)then
            call get_c8_project(c8_project,istatus)
            if (istatus .ne. 1) then
               write(6,*)'Error getting c8_project'
               return
            else
               write(6,*)'c8_project = ',c8_project
            endif

        else
C           Determine file format by looking at the file name convention
            call get_file_names(dir_in(1:len_dir_in),numoffiles
     1                         ,c_filenames,max_files,istatus)
            if(istatus .ne. 1 .or. numoffiles .eq. 0)then
                write(6,*)' Error calling get_file_names'
                istatus = 0
                return
            endif
            ipos = len_dir_in + 9
            c1_char = c_filenames(1)(ipos:ipos)
            if(c1_char .eq. '_')then
                c8_project = 'WFO' 
            else
                c8_project = 'NIMBUS'
            endif
            write(6,*)' 9th character of filename: ',c1_char
            write(6,*)' Setting c8_project parameter: ',c8_project

        endif

        if(c8_project(1:6) .eq. 'NIMBUS')then
            write(6,*)' Assumming /public filename format'
        else
            write(6,*)' Assumming WFO filename format'
        endif

!       Do we want hourly or 6 minute profiler data?
        if(c8_project(1:6) .eq. 'NIMBUS')then
            if(laps_cycle_time .le. 1800)then
!           if(.true.)then
                c5_data_interval = '0006o'
                write(6,*)' Using 6 minute data'
            else
                c5_data_interval = '0100o'
                write(6,*)' Using hourly data'
            endif
            c_filespec = dir_in(1:len_dir_in)//'*'//c5_data_interval

        else ! WFO
            c_filespec = dir_in(1:len_dir_in)//'*'

        endif

        write(6,*)c_filespec(1:80)

C       Wait for the data
        i4time_desired = i4time_sys
        i4_check_interval = 10
        i4_thresh_age = 3600
        i4time_now = i4time_now_gg()
        i4_total_wait = min(300,i4time_desired+25*60 - i4time_now)

        open(31,file='zzzz', status = 'old', err=10)
        read(31,*,err=10)i4_check_interval
        read(31,*,err=10)i4_total_wait
        read(31,*,err=10)i4_thresh_age
 10     continue
        close(31)

        if(i4_total_wait .gt. 0)then
            call wait_for_data(c_filespec,i4time_desired
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age       ! Only loop through the waiting
                                          ! if data is younger than this
                                          ! threshold
     1               ,istatus)

            if(istatus .ne. 1)then
                write(6,*)' No recent data'
                return        ! Normal action
!               continue      ! Do this for testing on the WFO
            endif
        endif

C       READ IN THE RAW PROFILER DATA
        if(c8_project(1:6) .eq. 'NIMBUS')then
            fnam_in = dir_in(1:len_dir_in)//asc9_tim//c5_data_interval
        else ! WFO
!           Convert from asc9_tim to asc13_tim
!           asc13_tim = '19960903_2200'                 ! Hardwired for testing.
            asc13_tim = fname9_to_wfo_fname13(asc9_tim) ! John Smart's routine
            fnam_in = dir_in(1:len_dir_in)//asc13_tim
        endif

        call s_len(fnam_in,len_fnam_in)
        write(6,*)fnam_in(1:len_fnam_in)
        CALL PROF_CDF_OPEN(fnam_in(1:len_fnam_in),cdfid,status)
C
C       Status of 0 means success.  Positive status codes and -1 are returned
C       from netCDF routines.  Any netCDF errors encountered will cause an
C       error message to be written to standard output (SYS$OUTPUT in VMS).
C       Errors -2 through -6 are from PROF_CDF routines and are explained in
C       the documentation for each PROF_CDF routine.
C
        if(status.ne.0)then
            write(6,*)' Warning: bad open ',status
            return
        endif

! added to read from file and set lag_time LW 8-27-98
C 	read global attribute avgTimePeriod from input file and set lag_time
        call prof_i4_avg_wdw(i4_avg_wdw_sec,cdfid,istatus)
        if(istatus .eq. 1)then
            lag_time = i4_avg_wdw_sec/2
        else
            write(6,*)' ingest_sub_pro: '
     1               ,'Error obtaining i4_avg_wdw_sec'
            return
        endif
C
C       Open an output file.
        ext = 'pro'
        call open_lapsprd_file(1,i4time_sys,ext,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error opening product file',ext
            return
        endif

!       Get the number of levels from the NetCDF file
        varid = NCDID(cdfid,'level',status)
        CALL NCDINQ(cdfid,varid,dimname,n_levels,status)

        write(6,*)' # of levels = ',n_levels

        if(n_levels .gt. max_levels)then
            write(6,*)' Error: too many levels in the data'
            istatus = 0
            return
        endif

!       Get the number of profilers from the NetCDF file
        varid = NCDID(cdfid,'recNum',status)
        CALL NCDINQ(cdfid,varid,dimname,file_n_prof,status)
        write(6,*)' # of profilers = ',file_n_prof
        if (file_n_prof .gt. n_profilers) then
          write(6,*)' Too many profilers to process'
          istatus = 0
          return
        endif

!	Get the number of characters in the station name
        varid = NCDID(cdfid,'staNamLen',status)
        CALL NCDINQ(cdfid,varid,dimname,staNamLen,status)

        start(1) = 1
        count(2) = 1
        count(1) = staNamLen

        do ista = 1,file_n_prof

!       fill profiler name into prof_name(ista)
          varid = NCVID(cdfid, 'staName', status)
          start(2) = ista
          CALL NCVGTC(cdfid, varid, start, count, 
     1                prof_name(ista), staNamLen, status)
          prof_name(ista)(6:6) = ' '
C
C         Get the surface pressure for Platteville.  This time we'll use the
C         5-character site name (plus a terminating blank) to select the station.
C
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'pressure',0,prs
     1                    ,status)
          if(status.ne.0)then
            if(status .eq. -3)then
                write(6,*)prof_name(ista),' not found'
                goto 900
            else
                write(6,*)' Warning: bad pressure read ',status
                return
            endif
          endif
          write(6,*)
          write(6,*)prof_name(ista),' pressure ',prs

!         Get the latitude
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'staLat',0,rlat
     1                      ,status)
          if(status.ne.0)then
              write(6,*)' Warning: bad lat read ',status
              return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' Lat ',rlat

!         Get the longitude
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'staLon',0,rlon
     1                      ,status)
          if(status.ne.0)then
              write(6,*)' Warning: bad lon read ',status
              return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' Lon ',rlon

!         Get the observation time
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'timeObs',0,timeObs
     1                      ,status)

          if(status.ne.0)then
              write(6,*)' Warning: bad timeObs read ',status,timeObs

          elseif(abs(timeObs) .gt. 3d9)then
              write(6,*)' Warning: Bad observation time',timeObs

          else
              call c_time2fname(nint(timeObs),a9_timeObs)

              write(6,*)
              write(6,*)' timeObs ',a9_timeObs

              call cv_asc_i4time(a9_timeObs,i4_timeObs)
              i4_resid = abs(i4_timeObs - i4time_sys)
!             if(i4_resid .gt. (ilaps_cycle_time / 2) )then ! outside time window
              if(i4_resid .gt. 0)then
                  write(6,*)' Warning, time is suspect '
     1                     ,a9_timeObs,i4_resid       
              endif

          endif

          if(rlat .le. rnorth .and. rlat .ge. south .and.
     1       rlon .ge. west   .and. rlon .le. east            )then       

            write(6,*)prof_name(ista),' is in box'

            CALL PROF_CDF_READ
     1          (cdfid,prof_name(ista),0,'staElev',0,elev,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad elev read ',status
                return
            endif
            write(6,*)
            write(6,*)prof_name(ista),' elev ',elev

            CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     1                         ,'windSpeedSfc',0,sp_sfc,istatus1)

            if(c8_project(1:6) .eq. 'NIMBUS')then
                CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     1                         ,'windDirSfc',0,di_sfc,istatus2)
            else
                CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     1                         ,'windDirSfc',0,i4_di_sfc,istatus2)

                di_sfc = i4_di_sfc
            endif

!           Test whether di and sp lie within valid range
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
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'uComponent',0,
     1                                             u,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad uComponent read ',status
                return
            endif

            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'vComponent',0,
     1                              v,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad vComponent read ',status
                return
            endif
C
C           Get the associated quality control flags.
C
            if(c8_project(1:6) .eq. 'NIMBUS')then
                CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     $                     ,'uvQualityCode',0,c1_qc_flag,status)
                if(status.ne.0)then
                    write(6,*)' Warning: bad qualityCode read ',status       
                    return
                endif

            else
                CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     $                     ,'uvQualityCode',0,i4_qc_flag,status)
                if(status.ne.0)then
                    write(6,*)' Warning: bad qualityCode read ',status
                    return
                endif

            endif

C
C           Get the associated levels for this profiler (also works with
C           global data statement on NIMBUS)
C
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'levels',0,
     $                     level,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad level read ',status
                return
            endif
C
            write(6,*)prof_name(ista)

            n_good_levels = 0

            do i = n_levels, 1, -1

                if(c8_project(1:6) .eq. 'NIMBUS')then
                    iqc_flag = byte_to_i4(c1_qc_flag(i))
                else
                    iqc_flag = i4_qc_flag(i)
                endif

                if(iqc_flag.eq.good)then
                        j = 1
                        iqc = 1
                else if(iqc_flag.eq.bad)then
                        j = 2
                        iqc = 0
                        write(6,*)prof_name(ista),' ',i,iqc_flag,' Bad'
                else if(iqc_flag.eq.missing)then
                        j = 3
                        iqc = 0
                        write(6,*)prof_name(ista),' ',i,iqc_flag,' Msg'
                else
                        iqc = 0
                        write(6,*)prof_name(ista),' ',i,iqc_flag,' QCd'
                endif
                height_msl = level(i) + elev
                write(6,*)i,height_msl,u(i),v(i),' ',iqc,iqc_flag ! ,iqc_flag2

                mode_flag = 1

                if(n_good_levels .ge. 1)then
                    if(height_msl .ge. ht_out(n_good_levels))then
                        mode_flag = 0 ! This is a low mode level repeating the
                                      ! wind from the high mode
                    endif
                endif


                if(  (.not.
     1    (u(i) .gt. r_missing_data .or. v(i) .gt. r_missing_data)  )
     1                          .and.
     1                      iqc .eq. 1
     1                          .and.
     1                           mode_flag .eq. 1
     1                                                          )then

                    n_good_levels = n_good_levels + 1
                    ht_out(n_good_levels) = height_msl
                    call uv_to_disp(u(i),v(i)
     1         ,di_out(n_good_levels),sp_out(n_good_levels))
                endif
            enddo ! i

            i4time_ob = i4_timeObs - lag_time ! i4time_sys - lag_time 

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)

            write(6,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'PROFILER'       
            write(1,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'PROFILER'       
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

            rms = 1.0

            if(n_good_sfc .eq. 1)then
!               write surface winds as first level
                write(1,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
                write(6,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
            endif

            do i = n_good_levels, 1, -1
                write(1,301,err=303)ht_out(i),di_out(i),sp_out(i),rms 
                write(6,301,err=303)ht_out(i),di_out(i),sp_out(i),rms 
301             format(1x,f6.0,f6.0,2f6.1)
303             continue
            enddo ! i

        else !
             write(6,*)prof_name(ista),' is outside of domain'

        endif ! l_in_box

        write(6,*)prof_name(ista)

900     enddo ! stations

        close(1)


C       Close the netCDF file.  This isn't necessary, but here's a sample of
C       the call.  A program could have up to 4 profiler netCDF files open at
C       a time.  The CDFID's are what indicate which file is which.
C
        CALL PROF_CDF_CLOSE(cdfid,status)
        if(status.ne.0)then
                write(6,*)' Warning: bad close ',status
        endif
C
        return
        end

