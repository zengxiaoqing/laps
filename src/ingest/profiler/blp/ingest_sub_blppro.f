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



        subroutine ingest_blppro(i4time,NX_L,NY_L,istatus)

C       Michael Barth           12-Aug-1993
C       Steve Albers               Nov-1993         Reworked for LAPS ingest
C                                  Oct-1994         Improve QC
C       Steve Albers                   1996         BL Profilers
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
C       NOTE: Profiler winds are written out in KNOTS

        integer cdfid,status,i,j,n_levels,n_profilers,file_n_prof
        parameter (n_levels = 150)
        parameter (n_profilers = 1000)

	parameter (max_modes = 3)
	parameter (max_gates = 50)
        integer*2 nmodes_short             
	integer nmodes

        real u(max_modes,max_gates)
        real v(max_modes,max_gates), prs

        character*1 qc_flag(max_modes,max_gates)
        character*4 c4_qc
        real*4 level(max_modes,max_gates)

        integer*2 ngates_short(max_modes)
        integer   ngates      (max_modes)

        real ht_out(n_levels)
        real di_out(n_levels)
        real sp_out(n_levels)

        integer good,bad,missing, start(2), count(2), staNamLen
        parameter (good = 0)
        parameter (bad = 12)
        parameter (missing = -1)
        character*1 qc_char(3)
        data qc_char/'G','B','M'/
        integer*4 byte_to_i4

        character*6 staname
        character*1 submode
        character*6 pltc2_name
        data pltc2_name/'PLTC2 '/

        character*100 fnam_in
        character*80 dir_in
        character*255 c_filespec
        integer wsmr_wmo_id
        data wsmr_wmo_id/74533/
        integer error_code
        data error_code/1/
        logical l_in_box
        data l_in_box/.true./

        integer varid
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

        character*13 filename13,outfile
        character*9 asc9_tim,a9time_ob

        character*31    ext
        integer*4       len_dir,len_dir_in

        character*40 c_vars_req
        character*100 c_values_req

        character*6 prof_name(n_profilers)

!       data prof_name /'ASTOR ','BRICO ','BSNOR ','FRACO ','GBYCO '
!    1                 ,'MEDOR ','MCGQU ','FTMCO ','NPTOR ','PLTC2 '
!    1                 ,'OLYAT'/

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           return
        endif
 
        r_mspkt = .518

        call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat,lon,
     1                  topo,1.0,rnorth,south,east,west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS perimeter'
            return
        endif

        CALL PROF_CDF_SET_ERROR(error_code,status)
        if(status.ne.0)then
            write(6,*)'bad set_error ',status
            return
        endif
C
C       Open an output file.
C
        ext = 'pro'
        call open_lapsprd_file(1,i4time,ext,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error opening product file',ext
            return
        endif

!       Read to end of file so we can append to it
        do i = 1,1000000
            read(1,*,end=2)
        enddo ! i
 2      continue

        outfile = filename13(i4time,'pro')
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

!       dir_in = '/public/data/wpdn/netcdf/'
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

C       Wait for the data

!       len_dir_in = 25
!       c_filespec = dir_in(1:len_dir_in)//'*0100o'
        c_filespec = dir_in(1:len_dir_in)//'*0??0o'

        write(6,*)c_filespec(1:80)

        i4time_desired = i4time

        i4time_stop_waiting = i4time + 25 * 60
        i4time_now = i4time_now_gg()
        i4_wait_period = i4time_stop_waiting - i4time_now

        i4_check_interval = 10
        i4_total_wait = min(300,i4_wait_period)
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

C       READ IN THE RAW PROFILER DATA

        fnam_in = dir_in(1:len_dir_in)//asc9_tim//'0100o'
!       fnam_in = dir_in(1:len_dir_in)//asc9_tim//'0??0o'
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
            write(6,*)'bad open ',status
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

!       do ista = 1,n_profilers

        do ista = 1,file_n_prof

!       fill profiler name into prof_name(ista)
          varid = NCVID(cdfid, 'staName', status)
          start(2) = ista
          CALL NCVGTC(cdfid, varid, start, count, 
     1                prof_name(ista), staNamLen, status)
          prof_name(ista)(6:6) = ' '

          write(6,*)
          write(6,*)' Looping for profiler ',prof_name(ista)
C
	  CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'numModesUsed',
     1                     0,nmodes_short,
     $			   status)
	  nmodes = nmodes_short		! INT*2 in file for compactness.

          if(status .ne. 0)then
            write(6,*)prof_name(ista)
     1               ,' bad status reading numModesUsed',status
            go to 900
          else
            write(6,*)'numModesused =',nmodes
          endif

	  CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'numLevelsUsed',
     1                     0,ngates_short,
     $			   status)

	  do i = 1, max_modes
		ngates(i) = ngates_short(i)
	  enddo

          if(status .ne. 0)then
            write(6,*)'bad status reading numLevelsUsed'
            go to 900
          endif

C
C       Get the surface pressure.  This time we'll use the
C       5-character site name (plus a terminating blank) to select the station.
C
          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'pressure',0,prs,st
     1atus)
          if(status.ne.0)then
            if(status .eq. -3)then
                write(6,*)prof_name(ista),' not found'
                goto 900
            else
                write(6,*)'bad pressure read ',status
                return
            endif
          endif
          write(6,*)
          write(6,*)prof_name(ista),' pressure ',prs

          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'staLat',0,rlat,sta       
     1tus)
          if(status.ne.0)then
                write(6,*)'bad lat read ',status
                return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' Lat ',rlat

          CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'staLon',0,rlon,sta
     1tus)
          if(status.ne.0)then
                write(6,*)'bad lon read ',status
                return
          endif
          write(6,*)
          write(6,*)prof_name(ista),' Lon ',rlon

          if(rlat .le. rnorth .and. rlat .ge. south .and.
     1       rlon .ge. west   .and. rlon .le. east            )then

            write(6,*)prof_name(ista),' is in box'

            CALL PROF_CDF_READ
     1          (cdfid,prof_name(ista),0,'staElev',0,elev,status)
            if(status.ne.0)then
                write(6,*)'bad elev read ',status
                return
            endif
            write(6,*)
            write(6,*)prof_name(ista),' elev ',elev

            CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     1                         ,'windSpeedSfc',0,sp_sfc,istatus1)
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0
     1                         ,'windDirSfc',0,di_sfc,istatus2)

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
                write(6,*)'bad uComponent read ',status
                return
            endif

            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'vComponent',0,
     1                              v,status)
            if(status.ne.0)then
                write(6,*)'bad vComponent read ',status
                return
            endif
C
C           Get the associated quality control flags.
C
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'uvQualityCode'
     $                     ,0,qc_flag,status)
            if(status.ne.0)then
                write(6,*)'bad qualityCode read ',status
                return
            endif
C
C           Get the associated levels
C
            CALL PROF_CDF_READ(cdfid,prof_name(ista),0,'levels',0,
     $                     level,status)
            if(status.ne.0)then
                write(6,*)'bad level read ',status
                return
            endif
C
            write(6,*)prof_name(ista)

            n_good_levels = 0
!           do i = 1, max_gates

            do im = 2,1,-1

              write(6,*)
              write(6,*)' Checking mode ',im,' for ',prof_name(ista)

              nlevels = ngates(im)

              do i = 1, nlevels
                iqc_flag = byte_to_i4(qc_flag(im,i))
!               iqc_flag2 = qc_flag(i)

                if(iqc_flag.eq.good)then
                        j = 1
                        iqc = 1
                        c4_qc = 'Good'
                else if(iqc_flag.eq.bad)then
                        j = 2
                        iqc = 0
                        c4_qc = 'Bad'
                else if(iqc_flag.eq.missing)then
                        j = 3
                        iqc = 0
                        c4_qc = 'Msg'
                else
                        iqc = 0
                        c4_qc = 'QCd'
                endif

                height_msl = level(im,i) + elev

                mode_flag = 1

!               Use mode 1 (100 km) only if it's above the highest good level 
!               from mode 2
                if(n_good_levels .ge. 1 .and. im .eq. 1)then
                    if(height_msl .le. ht_out(n_good_levels))then
                        mode_flag = 0 
                    endif
                endif

                write(6,*)i,height_msl,u(im,i),v(im,i),' ',iqc_flag,iqc
     1                   ,mode_flag,' ',c4_qc

                if(  (.not.
     1                    (u(im,i) .gt. r_missing_data      .or. 
     1                     v(im,i) .gt. r_missing_data)  )
     1                          .and.
     1                      iqc .eq. 1
     1                          .and.
     1                           mode_flag .eq. 1
     1                                                          )then

                    n_good_levels = n_good_levels + 1
                    ht_out(n_good_levels) = height_msl
                    call uv_to_disp(u(im,i),v(im,i)
     1                ,di_out(n_good_levels),sp_out(n_good_levels))
                endif
              enddo ! i
            enddo ! im

            lag_time = 1800
            i4time_ob = i4time - lag_time

            call make_fnam_lp(i4time_ob,a9time_ob,istatus)

            write(6,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob
            write(1,401)wsmr_wmo_id,n_good_levels+n_good_sfc,rlat
     1                 ,rlon,elev,prof_name(ista),a9time_ob
401         format(i12,i12,f11.3,f15.3,f15.0,5x,a6,3x,a9)

            if(n_good_sfc .eq. 1)then
!               write surface winds as first level
                write(1,301)elev,di_sfc,sp_sfc ! /r_mspkt
                write(6,301)elev,di_sfc,sp_sfc ! /r_mspkt
            endif

            do i = 1,n_good_levels
!           do i = n_good_levels, 1, -1
                write(1,301,err=303)ht_out(i),di_out(i),sp_out(i) ! /r_mspkt
                write(6,301,err=303)ht_out(i),di_out(i),sp_out(i) ! /r_mspkt
301             format(1x,f6.0,f6.0,f6.1)
303             continue
            enddo ! i

          else !
             write(6,*)prof_name(ista),' is outside of domain'

          endif ! l_in_box

          write(6,*)prof_name(ista)

900     enddo ! ista

        close(1)

C       Close the netCDF file.  This isn't necessary, but here's a sample of
C       the call.  A program could have up to 4 profiler netCDF files open at
C       a time.  The CDFID's are what indicate which file is which.
C
        CALL PROF_CDF_CLOSE(cdfid,status)
        if(status.ne.0)then
                write(6,*)'bad close ',status
        endif
C
        return
        end


