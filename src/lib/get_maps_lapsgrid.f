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
!       Jan 2001     "         - Obtain 2d fields from fsf or lgb when kmax = 1.
!                                New jacket routine get_modelfg_2d.

        include 'bgdata.inc'

        real field_3d_laps(imax,jmax,kmax)       ! Output array

        logical      lgab

        character*3  var_2d
        character*3  c_bkgd_ext(2)
        character*9  a9_time
        character*31 ext_a(maxbgmodels)
        character*31 subdir(maxbgmodels)
        character*9 fdda_model_source(maxbgmodels)
        character*5 bgmodelnames(maxbgmodels)
 
        integer     nbgm

!       ****************** Model Selection Section *********************

!       RESTRICTIONS:

!       1) No time interpolation is performed
!       2) FDDA/LGA file must be available valid at i4time/i4time_needed

        call bgmodel_name(maxbgmodels,nbgm,bgmodelnames,istatus)
!       do i=1,nbgm
!        print*,'bg model derived from from cmodel = ',bgmodelnames(i)
!       enddo

        if(kmax.eq.1)then
           c_bkgd_ext(1)='fsf'
           c_bkgd_ext(2)='lgb'
        else
           c_bkgd_ext(1)='fua'
           c_bkgd_ext(2)='lga'
        endif

        i4time_needed = i4time

        call get_fdda_model_source(fdda_model_source,n_fdda_models
     +,istatus)
        if(istatus .ne. 1)then
           print*,'Error getting fdda_model_source'
           return
        endif

        lgab=.false.
        if(n_fdda_models .eq. 0)then
           n_fdda_models = 1
c          subdir(n_fdda_models)=bgmodelnames(1:3) #replaces line below (needs mod) when lga subdirs is activated
           subdir(n_fdda_models)=c_bkgd_ext(2)
           ext_a(n_fdda_models) =c_bkgd_ext(2)
           lgab=.true.

c if lga is the first in the list then we'll force this to be the
c only possible background returned from get_modelfg_2/3d..
c removed functionality 12-06-01 (JRS)
c       elseif(fdda_model_source(1).eq.'lga')then
c          n_fdda_models = 1
c          subdir(n_fdda_models)=c_bkgd_ext(2)
c          ext_a(n_fdda_models) =c_bkgd_ext(2)
c          lgab=.true.

        else
           do i=1,n_fdda_models
              subdir(i)=fdda_model_source(i)
              if(subdir(i).eq.'lga')then        !.or.subdir(i).eq.'lgb')then
                 ext_a(i) =c_bkgd_ext(2)
                 lgab=.true.
              else
                 ext_a(i) =c_bkgd_ext(1)
              endif
           enddo
        endif
c this part adds lga/b to the list since it isn't already  part of it.
        if(.not.lgab)then
           if(n_fdda_models.lt.maxbgmodels)then
              n_fdda_models = n_fdda_models + 1
c             subdir(n_fdda_models)=bgmodelnames(1:3) #replaces line below (needs mod) when lga subdirs is activated
              subdir(n_fdda_models)=c_bkgd_ext(2)
              ext_a(n_fdda_models) =c_bkgd_ext(2)
           else
              print*,'*** WARNING *** '
              print*,'Cannot add lga/b to model background list in'
     .,'get_modelfg_3d, n_fdda_models > maxbgmodels'
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
        real field_3d_laps(imax,jmax,kmax)       ! Output array

        character*(*) var_2d
        character*(*) subdir
        character*9 a9_time
        character*14 a_filename

        character*125 comment_2d
        character*10 units_2d
        character*31 ext,ext_a
        character*150  directory
        character*255 c_filespec

        integer MAX_FILES
        parameter (MAX_FILES = 20000)
        character c_fnames(MAX_FILES)*180


        call get_directory(ext_a,directory,lend)
c Here -> reverse the directory order from "model/fua" to "fua/model"
        if(ext_a.ne.'lga'.and.ext_a.ne.'lgb')then
           call s_len(subdir,ld)
           if(ld .le. 0)then
               write(6,*)' subdir has zero length in get_modelfg_3d_sub'
               istatus = 0
               return
           endif
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

           elseif(lenf - lend .eq. 13 .OR. lenf - lend .eq. 14)then       

               if(kmax.gt.1)then
c
c new: changed variable name "ext" to "directory".
                  call get_lapsdata_3d(i4_initial,i4_valid,imax
     1                 ,jmax,kmax,directory,var_2d
     1                 ,units_2d,comment_2d,field_3d_laps,istatus)

                else

                  call get_lapsdata_2d(i4_initial,i4_valid,directory
     1              ,var_2d,units_2d,comment_2d,imax,jmax
     1              ,field_3d_laps(1,1,1),istatus)

c                 call get_2dgrid_dname(directory
c    1         ,i4_fn,0,i4time_nearest,EXT,var_2d,units_2d
c    1         ,comment_2d,imax,jmax,field_3d_laps(1,1,1),0,istatus)


                endif

           else
               write(6,*)' Error, illegal length of bckgd filename'
     1                 ,lend,lenf
               istatus = 0

           endif

           if(istatus .ne. 1)then
              write(6,*)'get_modelfg_3d_sub: Warning - could not read'
     1                 ,' model file'
           elseif(kmax.gt.1)then ! istatus = 1
               call qc_field_3d(var_2d,field_3d_laps
     1                         ,imax,jmax,kmax,istatus)            
           else
               call qc_field_2d(var_2d,field_3d_laps(1,1,1)
     1                         ,imax,jmax,istatus)
           endif

       else ! i_best_file = 0
           write(6,*)' No file with proper valid time'
           istatus = 0

       endif

       if(istatus .eq. 1)then
          call s_len(comment_2d,lenc)
          lenc = max(lenc,1)
          write(6,*)' Successfully obtained: '
     1                        ,c_fnames(i)(lend+1:lenf),' '
     1                        ,ext(1:6),var_2d,comment_2d(1:lenc)
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
        character*3 c3_hr

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

        elseif(i_len .eq. 14)then
            call cv_asc_i4time(a_filename(1:9),i4_filename)
            i4_initial = i4_filename
            c3_hr = a_filename(10:12)
            c2_mn = a_filename(13:14)
            read(c3_hr,3)i4_fcst_hr
            read(c2_mn,1)i4_fcst_mn
 3          format(i3)
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
            if(lenf.le.0)then
               call get_fname_length(c_fnames(i),lenf)
               lenf=lend+lenf
            endif
            if(lenf.eq.lend)then
               print*,'no filenames in c_fnames array: get_best_fcst'
               return
            endif

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
c
c-------------------------------------------------------------------------
        subroutine get_modelfg_2d(i4time,var_2d,imax,jmax
     1                          ,field_2d_laps,istatus)

c
c routine requires kmax = 1
c
        include 'bgdata.inc'

        integer imax,jmax

        real field_2d_laps(imax,jmax)       ! Output array

        character*3  var_2d
        character*9  a9_time
        character*9  fdda_model_source(maxbgmodels)
        character*31 ext_a(maxbgmodels)
        character*31 subdir(maxbgmodels)

        call get_modelfg_3d(i4time,var_2d,imax,jmax,1
     1                          ,field_2d_laps,istatus)

        if(istatus.ne.1)then
           print*,'warning: no 2d field from get_modelfg_2d'
           return
        endif

        return
        end
