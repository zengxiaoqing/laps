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
cdis



        subroutine ingest_lrs(i4time_sys,NX_L,NY_L,istatus)

C       Michael Barth           12-Aug-1993
C       Steve Albers               Nov-1993         Reworked for LAPS ingest
!       Ken Dritz                1-Jul-1997  Added NX_L and NY_L as dummy
!                                            arguments.
!       Ken Dritz                1-Jul-1997  Changed include of lapsparms.for
!                                            to laps_static_parameters.inc.
!       Ken Dritz                1-Jul-1997  Added call to get_r_missing_data.
!       Ken Dritz                8-Jul-1997  Replaced LAPS_DOMAIN_FILE by
!                                            'nest7grid' and removed include
!                                            of laps_static_parameters.inc.
C
C       This file shows examples of how the use PROF_CDF subroutines to read
C       WPDN 60-minute RASS data in netCDF files.
C
        integer cdfid,status,i,j,max_levels,max_stations
        parameter (max_levels = 100)
        parameter (max_stations = 200)

        real temp(max_levels),prs
c       character*1 qc_flag(max_levels)
        integer i_qc_flag(max_levels)

        integer level(max_levels), start(1), count(1)
        integer good,bad,missing
        parameter (good = 0)
        parameter (bad = 8)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'G','B','M'/
        character*6 staname
        character*200 fnam_in
        character*150 line
        character*180 dir_in
        character*255 c_filespec
        character*5 c5_data_interval

        integer wsmr_wmo_id
        integer wsmr_wmo_id_a(max_stations)
        integer error_code
        data error_code/1/
        integer byte_to_i4
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

        integer varid
        include 'netcdf.inc'
        character*(MAXNCNAM) dimname

        character*13 filename13,c13_dum
        character*9 asc9_tim,a9time_ob

        character*31    ext

        character*40 c_vars_req
        character*180 c_values_req

        character*9 a9_timeObs
        integer timeObs

        real lat(NX_L,NY_L),lon(NX_L,NY_L)
        real topo(NX_L,NY_L)

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           return
        endif
           
        CALL PROF_CDF_SET_ERROR(error_code,status)
        if(status.ne.0)then
                write(*,*)'bad set_error ',status
                return
        endif

        c13_dum = filename13(i4time_sys,'lrs')
        asc9_tim = c13_dum(1:9)
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

        c_vars_req = 'path_to_raw_rass'
        call get_static_info(c_vars_req,c_values_req,1,istatus)
        if(istatus .eq. 1)then
            write(6,*)c_vars_req(1:30),' = ',c_values_req
            dir_in = c_values_req
        else
            write(6,*)' Error getting ',c_vars_req
            return
        endif

        call s_len(dir_in,len_dir_in)

        call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat,lon,
     1            topo,1.0,rnorth,south,east,west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading domain perimeter'
            return
        endif

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write(6,*)'Error getting laps_cycle_time'
           return
        else
           write(6,*)'laps_cycle_time = ',laps_cycle_time
        endif

!       Do we want hourly or 6 minute profiler data?
        if(laps_cycle_time .le. 1800)then
!       if(.true.)then
            c5_data_interval = '0006o'
            write(6,*)' Using 6 minute data'
            i4time_desired = (i4time_sys / 360) * 360
            i4_avg_wdw_sec = 360 ! default value
        else
            c5_data_interval = '0100o'
            write(6,*)' Using hourly data'
            i4time_desired = i4time_sys
            i4_avg_wdw_sec = 3600 ! default value
        endif

        if(len_dir_in .gt. 0)then
            c_filespec = dir_in(1:len_dir_in)//'*'//c5_data_interval
        else
            write(6,*)' path_to_raw_rass has zero length'
            istatus = 0
            return
        endif

C       Wait for the data
        write(6,*)c_filespec(1:80)

        i4_check_interval = 10
        i4_total_wait = 300
        i4_thresh_age = 3600

        open(31,file='zzzz', status = 'old', err=10)
        read(31,*,err=10)i4_check_interval
        read(31,*,err=10)i4_total_wait
        read(31,*,err=10)i4_thresh_age
 10     continue
        close(31)

        call wait_for_data(c_filespec,i4time_desired
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age       ! Only loop through the waiting
                                          ! if data is younger than this
                                          ! threshold
     1               ,istatus)

        if(istatus .ne. 1)then
            write(6,*)' No recent data'
            return
        endif


C       READ IN THE RAW RASS DATA

        fnam_in = dir_in(1:len_dir_in)//asc9_tim//c5_data_interval
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
                write(*,*)'bad open ',status
                return
        endif

! added to read from file and set lag_time LW 8-27-98
C 	read global attribute avgTimePeriod from input file and set lag_time
        call prof_i4_avg_wdw(i4_avg_wdw_sec,cdfid,istatus)
        if(istatus .eq. 1)then
            write(6,*)' i4_avg_wdw_sec from file = ',i4_avg_wdw_sec
        else
            write(6,*)' ingest_sub_lrs: '
     1               ,'Warning: could not obtain i4_avg_wdw_sec'
            write(6,*)' Assuming i4_avg_wdw_sec = ',i4_avg_wdw_sec
        endif

        lag_time = i4_avg_wdw_sec/2
C
C       Open an output file.
C
        ext = 'lrs'
        call open_lapsprd_file(1,i4time_sys,ext(1:3),istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error opening output file'
            istatus = 0
            return
        endif

!       Get the number of levels from the NetCDF file
        status = NF_INQ_DIMID(cdfid,'level',varid)
        status = NF_INQ_DIM(cdfid,varid,dimname,n_levels)

        write(6,*)' # of levels = ',n_levels

        if(n_levels .gt. max_levels)then
            write(6,*)' Too many levels in the data'
            istatus = 0
            return
        endif

!       Get the number of profilers from the NetCDF file
        status = NF_INQ_DIMID(cdfid,'recNum',varid)
        status = NF_INQ_DIM(cdfid,varid,dimname,n_profilers)
        write(6,*)' # of RASSs = ',n_profilers
        if (n_profilers .gt. max_stations) then
          write(6,*)' Too many profilers to process'
     1             ,n_profilers,max_stations
          istatus = 0
          return
        endif

!       fill array wsmr_wmo_id_a with wmo station numbers
        status = NF_INQ_VARID(cdfid,'wmoStaNum',varid)
        start(1) = 1
        count(1) = n_profilers

        status = NF_GET_VARA_INT(cdfid, varid, start, count,
     !                           wsmr_wmo_id_a)
        if(status.ne.NF_NOERR)then
              write(*,*)'bad wmo_id array read ',status
              return
        endif

        do ista = 1, n_profilers

C         Here is where we'd like to try to read in the wsmr_wmo_id_a(ista)
C         (or the station name) for the particular station. We really want
C         both. Once we can get one via a direct NetCDF read, we should be
C         able to get the other by calling PROF_CDF_READ.
C
C         Get the station name for White Sands.  This shows how to get a
C         character string:  use the number of characters as the number of
C         array elements.
C
          wsmr_wmo_id = wsmr_wmo_id_a(ista) 

!         Using the wmo_id as input, attempt to retrieve the station name.
!         It is not essential we have the station name for the program to work.
          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'staName',1,
     $                     staname,status)
          call filter_string(staname)
          if(status.ne.0)then
                write(*,*)'bad staName read ',status
          else
                write(6,*)
                write(6,*)'staName ',staname
          endif
C
C         Get the surface pressure.  This time we'll use the
C         5-character site name (plus a terminating blank) to select the station.
C
          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'pressure',2,prs       
     1                      ,status)
          if(status.ne.0)then
                write(*,*)'bad pressure read ',status
                prs = r_missing_data
!               return
          endif
          write(6,*)
          write(6,*)'pressure ',prs

          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'temperature',2
     1                      ,t_sfc,status)
          i_qc_sfc = 1
          if(status.ne.0)then
                write(*,*)'bad t_sfc read ',status
                i_qc_sfc = 0
          endif
          write(6,*)
          write(6,*)'t_sfc ',t_sfc,i_qc_sfc

          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'relHumidity',2
     1                      ,rh_sfc,status)
          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'windSpeedSfc',2
     1                      ,sp_sfc,status)
          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'windDirSfc',2
     1                      ,di_sfc,status)

          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'staLat',2,rlat
     1                      ,status)
          if(status.ne.0)then
                write(*,*)'bad lat read ',status
                return
          endif
          write(6,*)
          write(6,*)'Lat ',rlat

          CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'staLon',2,rlon
     1                      ,status)
          if(status.ne.0)then
                write(*,*)'bad lon read ',status
                return
          endif
          write(6,*)'Lon ',rlon

!         Get the observation time
          status = NF_INQ_VARID(cdfid,'timeObs',varid)
          start(1) = ista
          count(1) = 1
          status = NF_GET_VARA_INT(cdfid, varid, start, count,
     1                             timeObs)
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

            write(6,*)
            write(6,*)'wmo id ', wsmr_wmo_id,' is in box'


            CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'staElev',2
     1          ,elev,status)
            if(status.ne.0)then
                write(*,*)'bad elev read ',status
                return
            endif
            write(6,*)
            write(6,*)'elev ',elev

!           i4time_ob = i4time_sys - lag_time
            i4time_ob = i4_timeObs - lag_time 

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)

            write(6,401)wsmr_wmo_id/100,n_levels+1,rlat,rlon,elev
     1                 ,staname(1:5),a9time_ob,'RASS    '
            write(1,401)wsmr_wmo_id/100,n_levels+1,rlat,rlon,elev
     1                 ,staname(1:5),a9time_ob,'RASS    '
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a5,3x,a9,1x,a8)
C
C           Get the array of RASS virtual temperatures for the profiler station.
C           For this call, we'll use the WMO
C           identifier to tell the subroutine what station we want.
C
            CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'virtualTemp',
     $                     2,temp,status)
            if(status.ne.0)then
                write(*,*)'bad virtualTemp read ',status
                return
            endif
C
C           Get the associated quality control flags.
C
            CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'qualityCode',       
     $                     0,i_qc_flag,status)
            if(status.ne.0)then
                write(*,*)'bad qualityCode read ',status
                return
            endif
C
C           Get the associated levels
C
            CALL PROF_CDF_READ(cdfid,'      ',wsmr_wmo_id,'level',0,
     $                     level,status)
            if(status.ne.0)then
                write(*,*)'bad level read ',status
                return
            endif
C
            write(6,*)'virtualTemp'
c           write(1,*)'virtualTemp'

            rms = 1.0

!           write surface temperature (and other data) as first level
            write(line,*)elev,t_sfc,' ',i_qc_sfc,rms
     1                                 ,rh_sfc,di_sfc,sp_sfc,prs      

            write(1,11)line
            write(6,11)line
 11         format(1x,a150)

            do i = 1, n_levels
                iqc_flag  = i_qc_flag(i)
                iqc_flag2 = i_qc_flag(i)

                if(iqc_flag.eq.good)then
                        j = 1
                        iqc = 1
                else if(iqc_flag.eq.bad)then
                        j = 2
                        iqc = 0
                else if(iqc_flag.eq.missing)then
                        j = 3
                        iqc = 0
                else ! Cover other cases too
                        iqc = 0
                endif

                height_msl = float(level(i)) + elev

                if(temp(i) .gt. r_missing_data)then
                    temp_out = r_missing_data
                else
                    temp_out = temp(i)
                endif

                write(1,*)height_msl,temp_out,' ',iqc,rms
                write(6,*)height_msl,temp_out,' ',iqc,rms
     1                                           ,iqc_flag,iqc_flag2       
            enddo ! i (level)

          else !
            write(6,*)'wmo id ', wsmr_wmo_id,' is outside of domain'

          endif ! l_in_box

        enddo ! ista

        close(1)
C
C       Close the netCDF file.  This isn't necessary, but here's a sample of
C       the call.  A program could have up to 4 profiler netCDF files open at
C       a time.  The CDFID's are what indicate which file is which.
C
        CALL PROF_CDF_CLOSE(cdfid,status)
        if(status.ne.0)then
            write(*,*)'bad close ',status
        endif
C
        return
        end
