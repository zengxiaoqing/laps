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
        subroutine get_modelfg_3d(i4time,var_2d,imax,jmax,kmax
     1                          ,field_3d_laps,istatus)

!       This routine reads the model background from RAMS or LGA/F and
!       interpolates in time to the requested i4time. The time interpolation
!       occurs for LGA/F, but the i4time must be a whole hour for a successful
!       read of the .RAM files to occur.

!          ~1990  Steve Albers - Original Version
!       Sep 1994      "        - Read in .RAM file first, then default over
!                                to LGA/F.
!       Jun 1995      "        - Slight cleanup and added comments
!       Jul 1995      "        - Now handles lga/f files that come as frequently
!                                as the LAPS cycle. This is backwards compatable
!                                as it still can do the time interpolation with
!                                three hourly forecast files.
!                     "        - New subroutine 'get_fcst_filename', modularizes
!                                the handling of forecast file naming convention.
!       Jan 1996               - Now will handle 9 or 13 character versions
!                                of lga filenames
!       Feb 1997               - Added entry get_modelfg_3d.
!       Oct 1998 Linda Wharton - removed variables never used: a9_filename,
!                                ext_f, field_3d_maps, l_fill, field_3d_maps_1,
!                                field_3d_maps_2
!
!       Feb 1999 John Smart    - Add subroutine get_modelfg_3d_sub which allows
!                                specific model extension.
!
!       Jun 2000     "         - Add subroutine get_best_fcst so that this code
!                                can be use elsewhere (eg., lga, dprep, laps_sfc, osse)
!          "         "         - Incorporate namelist parameter fdda_model_source to
!                                determine the subdirectory in which fdda model bckgnd
!                                exist.  Remove 'ram' as fdda extension.

        include 'bgdata.inc'

        real*4 field_3d_laps(imax,jmax,kmax)       ! Output array

        character*3  var_2d
        character*9  a9_time
        character*31 ext_a(maxbgmodels)
        character*31 subdir(maxbgmodels)
        character*15 fdda_model_source(maxbgmodels)

!       ****************** RAMS SECTION ***************************************

!       RESTRICTIONS:

!       1) No time interpolation is performed
!       2) FDDA/LGA file must be available valid at i4time/i4time_needed

        i4time_needed = i4time

        call get_fdda_model_source(fdda_model_source,n_fdda_models
     +,istatus)
        if(istatus .ne. 1)then
           print*,'Error getting fdda_model_source'
           return
        endif

        if(n_fdda_models .eq. 0)then
           n_fdda_models = 1
           subdir(n_fdda_models)='lga'
           ext_a(n_fdda_models) ='lga'
        else
           do i=1,n_fdda_models
              subdir(i)=fdda_model_source(i)
              ext_a(i) ='fua'
           enddo
           if(n_fdda_models.lt.maxbgmodels)then
              n_fdda_models = n_fdda_models + 1
              subdir(n_fdda_models)='lga'
              ext_a(n_fdda_models) ='lga'
           else
              print*,'*** WARNING *** '
              print*,'Cannot add lga to model background list in'
     .,'get_modelfg_3d'
           endif
        endif
           



        do isource = 1,n_fdda_models

           call make_fnam_lp(i4time_needed,a9_time,istatus)

           write(6,*)
           write(6,*)' Searching for model background valid at: '
     1                        ,a9_time,' ',ext_a(isource)(1:6),var_2d

           call get_modelfg_3d_sub(i4time_needed,var_2d,subdir(isource),
     1                      ext_a(isource),imax,jmax,kmax,field_3d_laps,
     1                      istatus)

           if(istatus.eq.1)then
              print*,'file obtained in get_modelfg_3d_sub - return'
              return
           endif

        enddo

        return
        end

!***************** START NEW SECTION ******************************************
! Start subroutine

        subroutine get_modelfg_3d_sub(i4time_needed,var_2d,subdir,ext_a
     1                         ,imax,jmax,kmax,field_3d_laps,istatus)
!
!
        real*4 field_3d_laps(imax,jmax,kmax)       ! Output array

        character*(*) var_2d
        character*(*) subdir
        character*9 a9_time
        character*13 a_filename

        character*125 comment_2d
        character*10 units_2d
        character*31 ext,ext_a
        character*150  directory
        character*255 c_filespec

        integer*4 MAX_FILES
        parameter (MAX_FILES = 300)
        character c_fnames(MAX_FILES)*180


        call get_directory(ext_a,directory,lend)
        if(ext_a.ne.'lga')then
           call s_len(subdir,ld)
           directory = directory(1:lend)//subdir(1:ld)//'/'
           lend=lend+ld+1
        endif
        c_filespec=directory(1:lend)

        c_filespec=c_filespec(1:lend)//'*.'//ext_a(1:3)

!       Obtain list of analysis/forecast filenames
        call Get_file_names(c_filespec,
     1                      i_nbr_files_ret,
     1                      c_fnames,
     1                      max_files,
     1                      istatus )

!       Determine which file having the proper valid time has the
!       most recent initialization time.

        call get_best_fcst(max_files,i4time_needed,i_nbr_files_ret
     1,c_fnames,i_best_file)

        if(i_best_file .gt. 0)then ! File for this ext exists with proper
           i = i_best_file
           call get_directory_length(c_fnames(i),lend)
           call get_time_length(c_fnames(i),lenf)                               ! valid time.
           ext = ext_a
           a_filename = c_fnames(i)(lend+1:lenf)

           write(6,*)' Found file for: ',c_fnames(i)(lend+1:lenf)
     1                                       ,' ',ext(1:6),var_2d

           call get_fcst_times(a_filename,i4_initial,i4_valid
     1                        ,i4_fn)

           if(lenf - lend .eq. 9)then

c             call get_laps_3d(i4_fn,imax,jmax,kmax,ext,var_2d
c    1                ,units_2d,comment_2d,field_3d_laps,istatus)
c
c
c new
c
              call get_3d_dir_time(directory,i4_fn
     1                      ,EXT,var_2d,units_2d,comment_2d
     1                      ,imax,jmax,kmax,field_3d_laps,istatus)

           elseif(lenf - lend .eq. 13)then
c
c new: changed variable name "ext" to "directory".
               call get_lapsdata_3d(i4_initial,i4_valid,imax
     1                 ,jmax,kmax,directory,var_2d
     1                 ,units_2d,comment_2d,field_3d_laps,istatus)

           else
               write(6,*)' Error, illegal length of lga filename'
     1                 ,lend,lenf
               istatus = 0

           endif

           if(istatus .ne. 1)then
              write(6,*)'get_modelfg_3d_sub: Warning - could not read'
     1                 ,' 3-D file'
           else ! istatus = 1
               call qc_field_3d(var_2d,field_3d_laps
     1                         ,imax,jmax,kmax,istatus)            
           endif

       else ! i_best_file = 0
           write(6,*)' No file with proper valid time'
           istatus = 0

       endif

       if(istatus .eq. 1)then
          write(6,*)' Successfully obtained: '
     1                        ,c_fnames(i)(lend+1:lenf),' '
     1                        ,ext(1:6),var_2d
          write(6,*)' Exiting get_modelfg_3d_sub'
          write(6,*)
          return
       else
          write(6,*)' No good ',ext_a(1:6),' files available.'          
       endif

       write(6,*)
       write(6,*)' No Good Files: exiting get_modelfg_3d_sub'

       istatus = 0
       return
       end


        subroutine get_fcst_times(a_filename,i4_initial,i4_valid,
     1                            i4_filename)

!       a_filename    (Character)                                         I
!       i4_initial                                                        O
!       i4_valid                                                          O
!       i4_filename   (This is what you might pass into read_laps_data)   O

!       This routine deals with either c9 or c13 forecast file naming
!       convention

        character*(*) a_filename
        character*2 c2_hr, c2_mn

        call s_len(a_filename,i_len)

        if(i_len .eq. 9)then
            call cv_asc_i4time(a_filename,i4_filename)
            i4_initial = (i4_filename/3600) * 3600
            i4_fcst = (i4_filename - i4_initial) * 60
            i4_valid = i4_initial + i4_fcst

        elseif(i_len .eq. 13)then
            call cv_asc_i4time(a_filename(1:9),i4_filename)
            i4_initial = i4_filename
            c2_hr = a_filename(10:11)
            c2_mn = a_filename(12:13)
            read(c2_hr,1)i4_fcst_hr
            read(c2_mn,1)i4_fcst_mn
 1          format(i2)
            i4_fcst =  i4_fcst_hr * 3600 + i4_fcst_mn * 60
            i4_valid = i4_initial + i4_fcst

        else
            write(6,*)' get_fcst_time: illegal value of i_len',i_len

        endif

        return
        end

c --------------------------------------------------------------------------
        subroutine get_best_fcst(maxfiles,i4time_needed
     1,i_nbr_files,c_fnames,i_best_file)

c
c determine the best file that matches the i4time_needed input
c J. Smart 6-20-00:  Pulled this section of software out of
c                    get_modelfg_3d_sub for use elsewhere in LAPS.
c
        implicit  none

        integer i4time_needed
        integer i_best_file,i
        integer i_nbr_files
        integer i4_fcst_time_min
        integer i4_valid,i4_fcst_time,i4_fn
        integer i4_initial
        integer lend,lenf
        integer maxfiles

        character*(*)  c_fnames(maxfiles)
        
        i_best_file = 0
        i4_fcst_time_min = 9999999

        do i=1,i_nbr_files
            call get_directory_length(c_fnames(i),lend)
            call get_time_length(c_fnames(i),lenf)
            call get_fcst_times(c_fnames(i)(lend+1:lenf)
     1                 ,i4_initial,i4_valid,i4_fn)
            if(i4_valid .eq. i4time_needed)then
               i4_fcst_time = i4_valid - i4_initial

               if(i4_fcst_time .lt. i4_fcst_time_min)then

                  i4_fcst_time = i4_fcst_time_min
                  i_best_file = i

               endif ! Smallest forecast time?
            endif ! Correct valid time
        enddo ! i

        return
        end
