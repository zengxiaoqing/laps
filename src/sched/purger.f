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
	program purger

        character*31 ext
        character*9 a9_time
        character*200 fname
        include 'lapsparms.for'

        call get_directory('etc',fname,len_fname)
        open(11,file=fname(1:len_fname)//'systime.dat',status='old')
c        open(11,file=fname(1:len_fname),status='old')
        read(11,*,err=1)i4time_now
        read(11,2,err=1)a9_time
2       format(1x,a9)
        close(11)

        write(6,*)' Systime = ',a9_time

!       i4time_now = i4time_now_gg()/3600 * 3600

1       nfiles = 24
 
        ntime_min = 720           ! How far back to keep files
                                   ! 1440 on ||, 2160 on operational

        ntime_min_lga = 720       ! How far back to keep lga files (minutes)
                                   ! 1440 on ||, 2160 on operational

        ifcst_time_max_lga = 1800  ! Maximum possible forward forecast (minutes)

!       LAPS DATA INGEST PRODUCTS

        ext = 'lga'
        call purge_lga(ext,ntime_min_lga,ifcst_time_max_lga,i4time_now)   

        ext = 'pmt'
        call purge_rad(ext,nfiles,ntime_min,i4time_now)

        ext = 'dmt'
        call purge_rad(ext,nfiles,ntime_min,i4time_now)

        ext = 'nhd'
        call purge_nowrad(ext,12,60,i4time_now)

        ext = 'lm1'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lm2'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'ram'
        call purge(ext,12,720,i4time_now)

        ext = 'rsf'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lve'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lh1'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lh2'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'ln3'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lso'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lvd'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'pin'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'pro'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lqo'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lrs'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'snd'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'vrc'
        call purge(ext,nfiles,ntime_min,i4time_now)

!       Purge the Doppler radar directories (v01, v02, etc.)
        do i = 1,max_radars
            if(i .lt. 10)then
                write(ext,91)i
 91             format('v0',i1)
            else
                write(ext,92)i
 92             format('v',i2)
            endif
            call purge(ext,12,120,i4time_now)

            if(i .lt. 10)then
                write(ext,93)i
 93             format('d0',i1)
            else
                write(ext,94)i
 94             format('d',i2)
            endif
            call purge(ext,12,120,i4time_now)

        enddo


!       LAPS ANALYSIS PRODUCTS


!       CLOUD PRODUCTS
        ext = 'lc3'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lwc'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lil'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lcb'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lct'
        call purge(ext,nfiles,2880,     i4time_now)

        ext = 'lcv'
        call purge(ext,nfiles,2880,     i4time_now)

        ext = 'lmd'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lps'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lco'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lrp'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lty'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lcp'
        call purge(ext,nfiles,ntime_min,i4time_now)

!       SPECIFIC HUMIDITY PRODUCTS
        ext = 'lh3'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lh4'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lq3'
        call purge(ext,nfiles,ntime_min,i4time_now)

!       SURFACE PACKAGE PRODUCTS
        ext = 'lsx'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lt1'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lgs'
        call purge(ext,2,120,i4time_now)


!       WIND ANALYSIS GRIDDED PRODUCTS
        ext = 'rdc'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lwm'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lhe'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lmt'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lf1'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lw3'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'liw'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'lmr'
        call purge(ext,nfiles,ntime_min,i4time_now)

!       WIND ANALYSIS GRAPHICAL PRODUCTS
        ext = 'sag'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'prg'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'msg'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'pig'
        call purge(ext,nfiles,ntime_min,i4time_now)

        ext = 'vdr'
        call purge(ext,nfiles,ntime_min,i4time_now) 

!       BALANCE PACKAGE PRODUCT
        ext = 'lba'
        call purge(ext,nfiles,ntime_min,i4time_now)

!       PRECIP PACKAGE PRODUCT
        ext = 'l1s'
        call purge(ext,48,2880,i4time_now)

        end

        subroutine purge(ext,nfiles,ntime_min,i4time_now)

!       Keeps number of files according to nfiles or time span according to
!       ntime_min, whichever is greater

        integer MAX_FILES
        parameter (MAX_FILES = 1000)

        character*9 asc_tim_9
        character*31 ext
        character*255 c_filespec
        character c_fnames(MAX_FILES)*80

        include 'lapsparms.for'

        call get_directory(ext,c_filespec,len)
         
        c_filespec = c_filespec(1:len)//'/*.'//ext(1:3)

        write(6,*)c_filespec

        call    Get_file_names(  c_filespec,
     1			 i_nbr_files_ret,
     1			 c_fnames,MAX_FILES,
     1			 i_status )

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
            write(6,*)i_nbr_files_ret,' file(s) in directory'
        else ! Error Condition
            write(6,*)' No files in directory'
            istatus = 0
            return
        endif

        ntime_sec = ntime_min * 60

10      do i=1,i_nbr_files_ret-nfiles ! Loop through excess versions
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
            if(i4time_now - i4time_file .gt. ntime_sec)then ! File is too old

!               Delete the file
                call rm_file(c_fnames(i)(1:lenf+13),istatus)


            endif
        enddo

        return
        end


        subroutine purge_rad(ext,nfiles,ntime_min,i4time_now)

!       Keeps number of files according to nfiles or time span according to
!       ntime_min, whichever is greater

        integer MAX_FILES
        parameter (MAX_FILES = 1000)

        character*9 asc_tim_9
        character*31 ext
        character*255 c_filespec
        character c_fnames(MAX_FILES)*80

        include 'lapsparms.for'

        call get_directory('radioprd',c_filespec,len)
         
        c_filespec = c_filespec(1:len)//'*.'//ext(1:3)
        call s_len(c_filespec,len)

        write(6,*)c_filespec

        call    Get_file_names(  c_filespec,
     1			 i_nbr_files_ret,
     1			 c_fnames,MAX_FILES,
     1			 i_status )

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
            write(6,*)i_nbr_files_ret,' file(s) in directory'
        else ! Error Condition
            write(6,*)' No files in directory'
            istatus = 0
            return
        endif

        ntime_sec = ntime_min * 60

10      do i=1,i_nbr_files_ret-nfiles ! Loop through excess versions
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
            if(i4time_now - i4time_file .gt. ntime_sec)then ! File is too old

!               Delete the file
                call rm_file(c_fnames(i)(1:lenf+13),istatus)


            endif
        enddo

        return
        end


        subroutine purge_nowrad(ext,nfiles,ntime_min,i4time_now)

!       Keeps number of files according to nfiles or time span according to
!       ntime_min, whichever is greater

        integer MAX_FILES
        parameter (MAX_FILES = 1000)

        character*9 asc_tim_9
        character*31 ext
        character*255 c_filespec
        character c_fnames(MAX_FILES)*80

        include 'lapsparms.for'

        call get_directory('nowrad6',c_filespec,len)
        c_filespec = c_filespec(1:len)//'*.'//ext(1:3)

        write(6,*)c_filespec

        call    Get_file_names(  c_filespec,
     1			 i_nbr_files_ret,
     1			 c_fnames,MAX_FILES,
     1			 i_status )

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
            write(6,*)i_nbr_files_ret,' file(s) in directory'
        else ! Error Condition
            write(6,*)' No files in directory'
            istatus = 0
            return
        endif

        ntime_sec = ntime_min * 60

10      do i=1,i_nbr_files_ret-nfiles ! Loop through excess versions
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
            if(i4time_now - i4time_file .gt. ntime_sec)then ! File is too old

!               Delete the file
                call rm_file(c_fnames(i)(1:lenf+13),istatus)


            endif
        enddo

        return
        end






        subroutine purge_lga(ext_a,ntime_min,ifcst_time_max_lga
     1                                                  ,i4time_now)

        include 'lapsparms.for' ! laps_cycle_time

!       Keeps number of files according to nfiles or time span according to
!       ntime_min, whichever is greater

        integer MAX_FILES
        parameter (MAX_FILES = 1000)

        character*9 asc_tim_9,a9_filename
        character*13 a_filename
        character*31 ext_a
        character*255 c_filespec
        character c_fnames(MAX_FILES)*80
        integer ikeep(MAX_FILES)

        character*50 directory

        ntime_sec = ntime_min * 60
        ifcst_sec_max = ifcst_time_max_lga * 60
!       laps_cycle_time = 3600

!***************** START NEW SECTION ******************************************

        do i = 1,MAX_FILES
            ikeep(i) = 0
        enddo ! i

        call get_directory(ext_a,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext_a(1:3)

!       Obtain list of analysis/forecast filenames
        call    Get_file_names(c_filespec,
     1                         i_nbr_files_ret,
     1                         c_fnames,
     1                         max_files,
     1                         istatus )

        call get_laps_config('nest7grid',istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling get_laps_config'
            return
        endif

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        do i4_vt = i4time_now + ifcst_sec_max, 
     1             i4time_now - ntime_sec, 
     1            -ilaps_cycle_time

!           Determine which file having the proper valid time has the
!           most recent initialization time.

            i_best_file = 0
            i4_fcst_time_min = 9999999

            call make_fnam_lp(i4_vt,asc_tim_9,istatus)

            do i=1,i_nbr_files_ret
                call get_directory_length(c_fnames(i),lend)
                call get_time_length(c_fnames(i),lenf)
                a_filename = c_fnames(i)(lend+1:lenf)
                call get_fcst_times(a_filename,i4_initial,i4_valid
     1                                        ,i4_fn)      
                if(i4_valid .eq. i4_vt)then
                    i4_fcst_time = i4_valid - i4_initial

                    if(i4_fcst_time .lt. i4_fcst_time_min)then
                        i4_fcst_time_min = i4_fcst_time
                        i_best_file = i
                    endif ! Smallest forecast time?
                endif ! Correct valid time
            enddo ! i

      

            if(i_best_file .gt. 0)then ! Keep this file
                i = i_best_file
                ikeep(i) = 1

                write(6,*)' Keeping file: ',c_fnames(i)(lend+1:lenf)
     1                   ,' Valid at: ',asc_tim_9

            else ! i_best_file = 0
                write(6,*)' No file with valid time: ',asc_tim_9         

            endif

        enddo ! Loop through (needed) valid times

!**************************END NEW SECTION************************************


        do i=1,i_nbr_files_ret ! Purge appropriate files
            a_filename = c_fnames(i)(lend+1:lenf)

            call get_fcst_times(a_filename,i4_initial,i4_valid,i4_fn)
            i4_age = i4time_now - i4_valid

            write(6,*)i,ikeep(i),i4_age,' ',c_fnames(i)(lend+1:lenf)

          ! Test whether file valid time is too old
            if(i4_age .gt. ntime_sec 
     1                        .OR.
     1                   ikeep(i) .eq. 0
     1                                                         )then 

!               Delete the file
!               call rm_file(c_fnames(i)(1:lenf+13),istatus)
                call rm_file(c_fnames(i),istatus)

            endif
        enddo

        return
        end


        subroutine rm_file(c_filename,istatus)

        character*(*) c_filename

        integer istatus

        lun = 151

        write(6,*)' rm_file ',c_filename

        open(lun,file=c_filename,status='unknown')

        close(lun,status='delete')

        istatus = 1
 
        return
        end

