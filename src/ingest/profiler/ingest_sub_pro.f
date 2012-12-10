cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   



        subroutine ingest_pro(i4time_sys,NX_L,NY_L,lun_out,istatus)

C       Michael Barth           12-Aug-1993
C       Steve Albers               Nov-1993         Reworked for LAPS ingest
C                                  Oct-1994         Improve QC
C       Steve Albers               Sep-1996         WFO database compatability
!       Ken Dritz                3-Jul-1997  Added NX_L, NY_L as dummy
!                                            arguments.
!       Ken Dritz                3-Jul-1997  Changed include of lapsparms.for
!                                            to laps_static_parameters.inc.
!       Ken Dritz                3-Jul-1997  Added call to get_r_missing_data.
C
C       This file shows examples of how the use PROF_CDF subroutines to read
C       WPDN 60-minute RASS data in netCDF files.
C
C       NOTE: Profiler winds are written out in KNOTS, and are sorted by HEIGHT

        integer cdfid,status,i,j,max_levels,max_levels_out,n_profilers
     1         ,file_n_prof
        parameter (max_levels = 72)
        parameter (max_levels_out = 72)
        parameter (n_profilers = 200) ! Accomodates RSA

        real u(max_levels),v(max_levels),prs
        real ht_out(max_levels_out),di_out(max_levels_out)
     1                             ,sp_out(max_levels_out)

        character*1 c1_qc_flag(max_levels)        ! for /public
        integer i4_qc_flag(max_levels)          ! for WFO

        real level(max_levels)
        integer good,bad,missing, start(2), count(2), staNamLen
        integer start_time(1), count_time(1)
        parameter (good = 0)
        parameter (bad = 12)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'G','B','M'/
        integer byte_to_i4

        character*200 fnam_in
        character*180 dir_in
        character*255 c_filespec

        include 'lapsparms.for'

        integer max_files

        parameter(max_files = MAX_INGEST_FILES)
        character*255 c_filenames(max_files)

        integer wsmr_wmo_id
        data wsmr_wmo_id/0/
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
        character*9 asc9_tim,a9time_ob,a9time_infile

        character*31    ext
        integer       len_dir_in

        character*40 c_vars_req
        character*180 c_values_req

        character*6 prof_name(n_profilers)
        character*9 a9_timeObs
        integer timeObs

        real lat(NX_L,NY_L),lon(NX_L,NY_L)
        real topo(NX_L,NY_L)

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

        call get_latlon_perimeter(NX_L,NY_L,1.0
     1                           ,lat,lon,topo
     1                           ,rnorth,south,east,west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS perimeter'
            return
        endif

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
        if(.false.)then !! WNI-BLS ... changed to false to force file name determination
            call get_c8_project(c8_project,istatus)
            if (istatus .ne. 1) then
               write(6,*)'Error getting c8_project'
               return
            else
               write(6,*)'c8_project = ',c8_project
            endif

        else
C           Determine file format by looking at the file name convention
            call get_file_names(dir_in(1:len_dir_in),numoffiles_ret
     1                         ,c_filenames,max_files,istatus)

!           Note that the GFN call may not work unless we also call
!           'filter_nonnumeric_fnames'.
            call Filter_non_numeric_fnames(c_filenames,
     1                   numoffiles_ret,
     1                   numoffiles,
     1                   max_files,
     1                   istatus)

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
                i4time_desired = (i4time_sys / 360) * 360
            else
                c5_data_interval = '0100o'
                write(6,*)' Using hourly data'
                i4time_desired = i4time_sys
            endif
            c_filespec = dir_in(1:len_dir_in)//'*'//c5_data_interval

        else ! WFO
            c_filespec = dir_in(1:len_dir_in)//'*'
            i4time_desired = i4time_sys

        endif

        call make_fnam_lp(i4time_desired,a9time_infile,istatus)

        write(6,*)c_filespec(1:80)

C       Wait for the data
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
            fnam_in = dir_in(1:len_dir_in)//a9time_infile
     1                                    //c5_data_interval
        else ! WFO
!           Convert from asc9_tim to asc13_tim
!           asc13_tim = '19960903_2200'                      ! Hardwired for testing.
            asc13_tim = fname9_to_wfo_fname13(a9time_infile) ! John Smart's routine
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

!       Get the number of levels from the NetCDF file
        status = NF_INQ_DIMID(cdfid,'level',varid)
        status = NF_INQ_DIM(cdfid,varid,dimname,n_levels)

        write(6,*)' # of levels = ',n_levels

        if(n_levels .gt. max_levels)then
            write(6,*)' Error: too many levels in the data'
            istatus = 0
            return
        endif

!       Get the number of profilers from the NetCDF file
        status = NF_INQ_DIMID(cdfid,'recNum',varid)
        status = NF_INQ_DIM(cdfid,varid,dimname,file_n_prof)
        write(6,*)' # of profilers = ',file_n_prof
        if (file_n_prof .gt. n_profilers) then
          write(6,*)' Too many profilers to process'
          istatus = 0
          return
        endif

!	Get the number of characters in the station name
        status = NF_INQ_DIMID(cdfid,'staNamLen',varid)
        status = NF_INQ_DIM(cdfid,varid,dimname,staNamLen)

        start(1) = 1
        count(2) = 1
        count(1) = staNamLen

        do ista = 1,file_n_prof

!       fill profiler name into prof_name(ista)
          status = NF_INQ_VARID(cdfid,'staName',varid)
          start(2) = ista
         
          status = NF_GET_VARA_TEXT(cdfid, varid, start, count,
     1                prof_name(ista))
          prof_name(ista)(6:6) = ' '
C
C         Get the surface pressure for Platteville.  This time we'll use the
C         5-character site name (plus a terminating blank) to select the station.
C
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'pressure',2,prs
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
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'staLat',2,rlat
     1                      ,status)
          if(status.ne.0)then
              write(6,*)' Warning: bad lat read ',status
              return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' Lat ',rlat

!         Get the longitude
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'staLon',2,rlon
     1                      ,status)
          if(status.ne.0)then
              write(6,*)' Warning: bad lon read ',status
              return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' Lon ',rlon

!         Get the observation time
          status = NF_INQ_VARID(cdfid,'timeObs',varid)
          start_time(1) = ista
          count_time(1) = 1
          status = NF_GET_VARA_INT(cdfid, varid, start_time, 
     1                             count_time,timeObs)
          if(status.ne.NF_NOERR)then
              write(6,*)' Warning: bad timeObs read ',status,timeObs

          elseif(abs(timeObs) .gt. 3d9)then
              write(6,*)' Warning: Bad observation time',timeObs

          else
              call c_time2fname(timeObs,a9_timeObs)

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
     1          (cdfid,prof_name(ista),0,'staElev',2,elev,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad elev read ',status
                return
            endif
            write(6,*)
            write(6,*)prof_name(ista),' elev ',elev

            call prof_sfcob_read(cdfid,prof_name(ista)          ! I
     1                          ,r_missing_data                 ! I
     1                          ,di_sfc                         ! O
     1                          ,sp_sfc                         ! O
     1                          ,p_sfc_hpa                      ! O
     1                          ,t_sfc_k                        ! O
     1                          ,rh_sfc_pct                     ! O
     1                          ,status)                        ! O

!           Test whether di and sp lie within valid range
            if(abs(sp_sfc) .gt. 500.)status = 1
            if(abs(di_sfc) .gt. 500.)status = 1

            if(status .eq. 0)then
                n_good_sfc = 1
            else
                n_good_sfc = 0
            endif
C
C           Get the array of profiler winds for the profiler station
C           at Platteville.  For this call, we'll use the WMO
C           identifier to tell the subroutine what station we want.
C
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'uComponent',2,
     1                                             u,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad uComponent read ',status
                return
            endif

            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'vComponent',2,
     1                              v,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad vComponent read ',status
                return
            endif
C
C           Get the associated quality control flags.
C
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     $                 ,'uvQualityCode',0,i4_qc_flag,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad qualityCode read ',status       
                return
            endif
C
C           Get the associated levels for this profiler (also works with
C           global data statement on NIMBUS)
C
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'levels',2,
     $                     level,status)
            if(status.ne.0)then
                write(6,*)' Warning: bad level read ',status
                return
            endif
C
            write(6,*)prof_name(ista)

            n_good_levels = 0

            do i = n_levels, 1, -1

                iqc_flag = i4_qc_flag(i)

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
C
C           Open an output file if needed
            ext = 'pro'
            call open_ext(lun_out,i4time_sys,ext,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error opening product file',ext
                return
            endif

            write(6,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'PROFILER'       
            write(1,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob,'PROFILER'       
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9,1x,a8)

            rms = 1.0

            if(n_good_sfc .eq. 1)then
!               write surface data as first level
                write(1,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
     1                          ,p_sfc_hpa 
     1                          ,t_sfc_k   
     1                          ,rh_sfc_pct
                write(6,301)elev,di_sfc,sp_sfc,rms ! /r_mspkt
     1                          ,p_sfc_hpa 
     1                          ,t_sfc_k   
     1                          ,rh_sfc_pct
            endif

            do i = n_good_levels, 1, -1
                write(1,301,err=303)ht_out(i),di_out(i),sp_out(i),rms 
                write(6,301,err=303)ht_out(i),di_out(i),sp_out(i),rms 
301             format(1x,f6.0,f6.0,2f6.1,3f7.1)
303             continue
            enddo ! i

        else !
             write(6,*)prof_name(ista),' is outside of domain'

        endif ! l_in_box

        write(6,*)prof_name(ista)

900     enddo ! stations


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

        subroutine prof_sfcob_read(cdfid,prof_name              ! I
     1                          ,r_missing_data                 ! I
     1                          ,di_sfc                         ! O
     1                          ,sp_sfc                         ! O
     1                          ,p_sfc_hpa                      ! O
     1                          ,t_sfc_k                        ! O
     1                          ,rh_sfc_pct                     ! O
     1                          ,istatus)                       ! O

        integer istatus   ! A value of zero indicates good dir & speed

        integer cdfid
        character*8 c8_project
        character*(*) prof_name

        di_sfc = r_missing_data
        sp_sfc = r_missing_data

        call get_sfc_badflag(sfc_badflag,istatus)
        if(istatus .ne. 1)return

        t_sfc_k = sfc_badflag
        rh_sfc_pct = sfc_badflag
        p_sfc_hpa = sfc_badflag

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)return

        CALL PROF_CDF_READ(cdfid,prof_name,0
     1                         ,'windSpeedSfc',2,sp_sfc,istat_sp)

        if(c8_project(1:6) .eq. 'NIMBUS')then
            CALL PROF_CDF_READ(cdfid,prof_name,0
     1                         ,'windDirSfc',2,di_sfc,istat_di)
        else
            CALL PROF_CDF_READ(cdfid,prof_name,0
     1                         ,'windDirSfc',0,i4_di_sfc,istat_di)

            di_sfc = i4_di_sfc
        endif

        CALL PROF_CDF_READ(cdfid,prof_name,0
     1                         ,'pressure',2,p_sfc_hpa,istatus2)

        CALL PROF_CDF_READ(cdfid,prof_name,0
     1                         ,'temperature',2,t_sfc_k,istatus2)

        CALL PROF_CDF_READ(cdfid,prof_name,0
     1                         ,'relHumidity',2,rh_sfc_pct,istatus2)

        istatus = istat_di * istat_sp

        if(abs(p_sfc_hpa) .gt. 5000.)p_sfc_hpa  = sfc_badflag
        if(abs(t_sfc_k) .gt. 500.)   t_sfc_k    = sfc_badflag
        if(abs(rh_sfc_pct) .gt. 500.)rh_sfc_pct = sfc_badflag

        return
        end
