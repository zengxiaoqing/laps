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

        subroutine get_laps_2d(i4time,EXT,var_2d,units_2d,
     1                  comment_2d,imax,jmax,field_2d,istatus)

cdoc    Used to read in a surface grid with inputs of time and ext

        character*9 asc9_tim
        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer LVL_2d
        character*4 LVL_COORD_2d

        real field_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_2d: bad istatus, return'
            return
        endif

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)trim(directory),asc9_tim,ext,var_2d
11      format(' Reading 2d ',a,1x,a,1x,a,1x,a)

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!       Check for missing data
        do j = 1,jmax
        do i = 1,imax
            if(istatus .eq. 1)then
                if(field_2d(i,j) .eq. r_missing_data)then
                    write(6,*)' Missing Data Value Detected in 2D Field'
                    istatus = -1
                endif
            endif
        enddo ! i
        enddo ! j

        return
        end

        subroutine get_lapsdata_2d(i4time,i4_valid,EXT
     1,var_2d,units_2d,comment_2d,imax,jmax,field_2d,istatus)

!       Returns a 2D laps grid
!       i4time              Input      Desired i4time initial
!       i4_valid            Input      i4time for valid data time
!       ext                 Input      3 character file extension
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       imax,jmax           Input      LAPS grid dimensions
!       field_2d            Output     2D grid
!       istatus             Output     status

!       Steve Albers            1996

cdoc    This routine can be used to read in a surface grid of known time
cdoc    by calling the new READ_LAPS routine

        include      'bgdata.inc'

        character*9 asc9_tim
        character*150 DIRECTORY
        character*(*) EXT
        character*31  EXT_INT

        character*125 comment_2d
        character*10 units_2d
        character*15 fdda_model_source(maxbgmodels)
        character*3 var_2d
        integer   LVL_2d
        integer   nfdda
        character*4 LVL_COORD_2d

        real field_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_lapsdata_2d: bad istatus, return'
            return
        endif

c       call get_fdda_model_source(fdda_model_source,nfdda,istatus)

        call s_len(ext,len)

        if(len.gt.3)then
           directory = ext
           call s_len(directory,len_dir)
c
c code section below needs more when we reverse the fdda subdirectory
c order (ie., from fsf/"model" to "model"/fsf)
c
           if(directory(len_dir:len_dir).eq.'/')then
              len_dir=len_dir-1
              directory=directory(1:len_dir)
           endif
           call get_directory_length(directory,len)
           if(directory(len+1:len_dir).eq.'lgb')then
              ext_int='lgb'

           else ! maybe fsf (determine extension from directory name)

c             do l=1,nfdda
c                if(directory(len+1:len_dir).eq.
c    &fdda_model_source(l))then

              mark=0
              search_dir: do l=len_dir-1,1,-1
                 if(directory(l:l).eq.'/')then 
                    mark=l
                    exit search_dir
                 endif
              enddo search_dir
                    
              ext_int=directory(mark-3:mark-1)

c                endif
c             enddo
           endif ! directory is lgb

c          if(ext_int.eq.' ')then
c             print*,'Unknown lapsprd extension'
c    &,directory(1:len_dir)
c             istatus = 0
c             return
c          endif

           len=3
           len_dir=len_dir+1
           directory(len_dir:len_dir)='/'
        else
           call get_directory(ext,directory,len_dir)
           ext_int=ext
        endif

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)trim(directory),ext_int(1:len),var_2d
11      format(' Reading 2d ',a,1x,a5,1x,a3)

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL READ_LAPS(I4TIME,i4_valid,DIRECTORY,EXT_INT
     1,imax,jmax,1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D
     1, COMMENT_2D,field_2d,ISTATUS)

!       Check for missing data
        do j = 1,jmax
        do i = 1,imax
            if(istatus .eq. 1)then
                if(field_2d(i,j) .eq. r_missing_data)then
                    write(6,*)' Missing Data Value Detected in 2D Field'
                    istatus = -1
                endif
            endif
        enddo ! i
        enddo ! j

        return
        end

        subroutine get_laps_2dgrid(i4time_needed,i4tol,i4time_nearest
     1         ,ext_in,var_2d,units_2d
     1         ,comment_2d,imax,jmax,field_2d,ilevel,istatus)

cdoc    Returns a 2-D grid. Inputs include the extension and time window.

        character*9 asc9_tim

        character*150 DIRECTORY
        character*(*) EXT_IN
        character*31 ext

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer LVL_2d
        character*4 LVL_COORD_2d

        real field_2d(imax,jmax)

        character*255 c_filespec

        logical ltest_vertical_grid

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_2dgrid: bad istatus, return'
            return
        endif

        call get_directory(ext_in,directory,len_dir)

        if(ext_in .eq. 'balance')then
            directory = directory(1:len_dir)//'lsx/'
            len_dir = len_dir + 4
            ext = 'lsx'
        else
            ext = ext_in
        endif

        call s_len(ext,lenext)

        c_filespec = directory(1:len_dir)//'*.'//ext(1:lenext)
        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .le. i4tol)then
            if(ilevel .ne. 0)then
                if(ltest_vertical_grid('HEIGHT'))then
                    lvl_2d = zcoord_of_level(k)/10
                    lvl_coord_2d = 'MSL'
                elseif(ltest_vertical_grid('PRESSURE'))then
                    lvl_2d = ilevel
                    lvl_coord_2d = 'MB'
                elseif(var_2d .eq. 'LSM' .and. ilevel .ge. -3 
     1                                   .and. ilevel .le. -1)then
                    lvl_2d = ilevel
                    lvl_coord_2d = 'HPA'
                else
                    write(6,*)' Error, vertical grid not supported,'
     1                      ,' this routine supports PRESSURE or HEIGHT'       
                    istatus = 0
                    return
                endif

            else
                lvl_2d = 0
                lvl_coord_2d = 'MSL'

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            write(6,11)trim(directory),asc9_tim,ext(1:5),var_2d
11          format(' Reading 2d ',a,1x,a9,1x,a5,1x,a3)

            CALL READ_LAPS_DATA(I4TIME_nearest,DIRECTORY,EXT,imax,jmax,
     1                     1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!           Check for missing data
            do j = 1,jmax
            do i = 1,imax
            if(istatus .eq. 1)then
                if(field_2d(i,j) .eq. r_missing_data)then
                    write(6,*)' Missing Data Value Detected in 2D Field'
                    istatus = -1
                endif
            endif
            enddo ! i
            enddo ! j

        else
            write(6,*)' No field found within window ',ext(1:10)
            istatus = 0

        endif

        return
        end
c
        subroutine get_2dgrid_dname(directory
     1         ,i4time_needed,i4tol,i4time_nearest
     1         ,EXT,var_2d,units_2d
     1         ,comment_2d,imax,jmax,field_2d,ilevel,istatus)

cdoc    Returns a 2-D grid. Inputs include the directory, ext, and time window.

!       Steve Albers            1990
!           J Smart             1998

        character*9 asc9_tim

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer LVL_2d
        character*4 LVL_COORD_2d

        real field_2d(imax,jmax)

        character*255 c_filespec

        logical ltest_vertical_grid

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_2dgrid_dname: bad istatus, return'
            return
        endif

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lenext = i-1

        len_dir=index(directory,' ')-1
        c_filespec = directory(1:len_dir)//'*.'//ext(1:lenext)
        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .le. i4tol)then
            if(ilevel .ne. 0)then
                if(ltest_vertical_grid('HEIGHT'))then
                    lvl_2d = zcoord_of_level(k)/10
                    lvl_coord_2d = 'MSL'
                elseif(ltest_vertical_grid('PRESSURE'))then
                    lvl_2d = ilevel
                    lvl_coord_2d = 'MB'
                else
                    write(6,*)' Error, vertical grid not supported,'
     1                      ,' this routine supports PRESSURE or HEIGHT'       
                    istatus = 0
                    return
                endif

            else
                lvl_2d = 0
                lvl_coord_2d = 'MSL'

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            write(6,11)trim(directory),asc9_tim,ext(1:5),var_2d
11          format(' Reading 2d ',a,1x,a9,1x,a5,1x,a3)

            CALL READ_LAPS_DATA(I4TIME_nearest,DIRECTORY,EXT,imax,jmax,
     1                     1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!           Check for missing data
            do j = 1,jmax
            do i = 1,imax
            if(istatus .eq. 1)then
                if(field_2d(i,j) .eq. r_missing_data)then
                    write(6,*)' Missing Data Value Detected in 2D Field'
                    istatus = -1
                endif
            endif
            enddo ! i
            enddo ! j

        else
            write(6,*)' No field found within window ',ext(1:10)
            istatus = 0

        endif

        return
        end




        subroutine get_laps_2dvar(i4time_needed,i4tol,i4time_nearest
     1         ,lat,lon
     1         ,subpoint_lat_clo,subpoint_lon_clo      ! O 
     1         ,EXT,var_2d,units_2d
     1         ,comment_2d,imax,jmax,field_2d,ilevel,istatus)

!       Steve Albers            1996
cdoc    This routine tries to read in the desired variable from all files
cdoc    having the proper extension, picking the closest one within the
cdoc    specified time window.
!
!       J Smart                 1998
cdoc    added lvd subdirectory flexibility. Only one 2d satellite field returned.

!       Steve Albers            2011
!       Pass in lat/lon so that it can potentally be used for satellite mosaicing

        character*9 asc9_tim

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer LVL_2d
        character*4 LVL_COORD_2d

        real field_2d(imax,jmax)
        real lat(imax,jmax)
        real lon(imax,jmax)

        integer max_files
        parameter (max_files = 20000)
        character*255 c_filespec
        character*120 c_fnames(max_files)
        integer i4times(max_files)
        integer i_selected(max_files)

        real    subpoint_lat_clo(imax,jmax)
        real    subpoint_lon_clo(imax,jmax)

        include 'satellite_dims_lvd.inc'
        include 'satellite_common_lvd.inc'

        logical ltest_vertical_grid

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_2dvar: bad istatus, return'
            return
        endif

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lenext = i-1

c
c lvd switch
c
        if(ext(1:lenext).eq.'lvd')then
           call config_satellite_lvd(istatus)
           if(istatus.ne.1)then
              return
           endif

           call get_laps_sat(maxsat,c_sat_id,isats
     1     ,i4time_needed,i4tol,i4time_nearest
     1     ,var_2d,units_2d,comment_2d,imax,jmax
     1     ,subpoint_lat_clo,subpoint_lon_clo      ! O 
     1     ,lat,lon,field_2d,istatus)

           if(istatus.ne.1)then
              write(6,*)'No data returned from get_laps_sat'
              return
           else
              return
           endif
        endif

        call get_directory(ext,directory,len_dir)

        do i = 1,max_files
            i_selected(i) = 0
        enddo ! i

        c_filespec = directory(1:len_dir)//'*.'//ext(1:lenext)

        call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)
        if(istatus .ne. 1)then
            write(6,*)'get_laps_2dvar: Bad status returned '
     1               ,'from get_file_times'
            return
        endif

50      i4_diff_min = 999999999
        do i = 1,i_nbr_files_ret
            i4_diff = abs(i4times(i) - i4time_needed)
            if(i_selected(i) .eq. 0)then
                i4_diff_min = min(i4_diff,i4_diff_min)
            endif
        enddo ! i

        if(i4_diff_min .gt. i4tol)then
            write(6,*)' No remaining files found within ',i4tol
     1               ,' sec time window ',ext(1:5),var_2d
            istatus = 0
            return
        endif

        do i = 1,i_nbr_files_ret
            i4_diff = abs(i4times(i) - i4time_needed)
            if(i4_diff .eq. i4_diff_min .and. i_selected(i) .eq. 0)then
                i_selected(i) = 1

                if(ilevel .ne. 0)then
                    if(ltest_vertical_grid('HEIGHT'))then
                        lvl_2d = zcoord_of_level(k)/10
                        lvl_coord_2d = 'MSL'
                    elseif(ltest_vertical_grid('PRESSURE'))then
                        lvl_2d = ilevel
                        lvl_coord_2d = 'MB'
                    else
                        write(6,*)' Error, vertical grid not supported,'
     1                      ,' this routine supports PRESSURE or HEIGHT'
                        istatus = 0
                        return
                    endif

                else
                    lvl_2d = 0
                    lvl_coord_2d = 'MSL'

                endif

                call make_fnam_lp(i4times(i),asc9_tim,istatus)

                write(6,11)trim(directory),asc9_tim,ext(1:5),var_2d
11              format(' Reading 2d ',a,1x,a9,1x,a5,1x,a3)

                CALL READ_LAPS_DATA(i4times(i),DIRECTORY,EXT,imax
     1            ,jmax,1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D
     1            ,COMMENT_2D,field_2d,ISTATUS)

                if(istatus .ne. 1)then
                    write(6,*)' No field found at ',ext(1:10)
     1                       ,var_2d,' ',asc9_tim
                    go to 50

                else   !  istatus = 1, check for missing data
                    do il = 1,imax
                    do jl = 1,jmax
                        if(field_2d(il,jl) .eq. r_missing_data)then
                            write(6,*)il,jl,
     1                        ' Missing Data Value Detected in 2D Field'
                            istatus = -1
                            return
                        endif
                    enddo ! j
                    enddo ! i

                endif

                return ! istatus = 1 and no missing data

            endif ! File is closest unread file to desired time
        enddo ! ith file

        end

        subroutine get_laps_3d(i4time,imax,jmax,kmax
     1  ,EXT,var_2d,units_2d,comment_2d,field_3d,istatus)

cdoc    Returns a 3-D grid. Inputs include the extension and time.

!       i4time              Input      Desired i4time
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       ext                 Input      3 character file extension
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

!       Steve Albers            1990

        character*150 DIRECTORY
        character*(*) EXT, var_2d

        character*125 comment_3d(kmax),comment_2d
        character*10 units_3d(kmax),units_2d
        character*3 var_3d(kmax)
        integer LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real field_3d(imax,jmax,kmax)

        write(6,*)' Subroutine get_laps_3d...'

        call get_directory(ext,directory,len_dir)

        call get_3d_dir_time(directory,i4time
     1                      ,EXT,var_2d,units_2d
     1                      ,comment_2d
     1                      ,imax,jmax,kmax,field_3d,istatus)

        return
        end

        subroutine get_lapsdata_3d(i4time,i4_valid
     1,imax,jmax,kmax,EXT,var_2d,units_2d,comment_2d
     1,field_3d,istatus)

cdoc    Returns a 3-D fcst grid. Inputs include directory, initial and valid time.

!       i4time              Input      Desired i4time initial
!       i4_valid            Input      i4time for valid data time
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       ext                 Input      3 character file extension or model member
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

!       Steve Albers            1990

        use mem_namelist, ONLY: lvl_coord_cdf

        include      'bgdata.inc'

        character*150 DIRECTORY
        character*(*) EXT
        character*31  EXT_INT

        character*125 comment_3d(kmax),comment_2d
        character*10 units_3d(kmax),units_2d
        character*15 fdda_model_source(maxbgmodels)
        character*3  var_3d(kmax),var_2d

        integer     LVL_3d(kmax)
        integer     nfdda,l
        character*4 LVL_COORD_3d(kmax)

        real field_3d(imax,jmax,kmax)
        real sigma_1d(kmax)
        real ht_1d(kmax)

        logical ltest_vertical_grid

        call get_fdda_model_source(fdda_model_source,nfdda,istatus)

        call s_len(ext,len)
        if(len.gt.3)then    !this tests whether ext is truely an extension or not.
           directory = ext
           call s_len(directory,len_dir)
c
c code section below needs more when we reverse the fdda subdirectory
c order (ie., from fsf/"model" to "model"/fsf)
c
           if(directory(len_dir:len_dir).eq.'/')then
              len_dir=len_dir-1
              directory=directory(1:len_dir)
           endif
           call get_directory_length(directory,len)
           if(directory(len+1:len_dir).eq.'lga')then
              ext_int='lga'
           else

c             do l=1,nfdda
c                if(directory(len+1:len_dir).eq.
c    &fdda_model_source(l))then

                    ext_int='fua'

c                endif
c             enddo

           endif
c          if(ext_int.eq.' ')then
c             print*,'Unknown lapsprd extension'
c    &,directory(1:len_dir)
c             istatus = 0
c             return
c          endif

           len=3
           len_dir=len_dir+1
           directory(len_dir:len_dir)='/'
        else
           call get_directory(ext,directory,len_dir)
           ext_int=ext
        endif

        write(6,11)directory(1:len_dir),ext_int(1:len)
     1,var_2d
11      format(' Reading 3d ',a,1x,a5,1x,a3)

        if(kmax.gt.1)then
          do k = 1,kmax
            units_3d(k)   = units_2d
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                if(lvl_coord_cdf(1:3) .eq. 'HPA')then
                    lvl_3d(k) = nint(zcoord_of_level(k))/100
                    lvl_coord_3d(k) = 'MB' ! informational
                elseif(lvl_coord_cdf(1:3) .eq. 'PA')then
                    lvl_3d(k) = nint(zcoord_of_level(k))
                    lvl_coord_3d(k) = 'PA' ! informational
                else
                    write(6,*)
     1               ' Error, pressure units not supported '
     1               ,lvl_coord_cdf
                    istatus = 0
                    return
                endif
            elseif(ltest_vertical_grid('SIGMA_P'))then
                if(k .eq. 1)then
                    write(6,*)' Reading Sigma P Levels'
                    call get_sigma_1d(kmax,sigma_1d,istatus)
                    if(istatus .ne. 1)then
                        return
                    endif
                endif

                lvl_3d(k) = nint(sigma_1d(k) * 1000.)
                lvl_coord_3d(k) = '  ' ! informational
            elseif(ltest_vertical_grid('SIGMA_HT'))then
                if(k .eq. 1)then
                    write(6,*)' Reading Sigma Ht Levels'
                    call get_ht_1d(kmax,ht_1d,istatus)
                    if(istatus .ne. 1)then
                        return
                    endif
                endif

                lvl_3d(k) = nint(ht_1d(k))
                lvl_coord_3d(k) = '  ' ! informational
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' in subroutine get_lapsdata_3d'                
                istatus = 0
                return
            endif


            var_3d(k) = var_2d

          enddo ! k


          CALL READ_LAPS(I4TIME,I4_VALID,DIRECTORY,EXT_INT,imax,jmax,
     1  kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,field_3d,ISTATUS)

          comment_2d=comment_3d(1)
          units_2d=units_3d(1)

        else


          CALL READ_LAPS(I4TIME,i4_valid,DIRECTORY,EXT_INT
     1,imax,jmax,1,1,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_2D
     1, COMMENT_2D,field_3d,ISTATUS)


        endif
        return
        end

        subroutine get_laps_3dgrid(i4time_needed,i4tol,i4time_nearest,
     1          imax,jmax,kmax,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,istatus)

cdoc    Returns a 3-D analysis grid. Inputs include the extension and time window.

!       i4time_needed       Input      Desired i4time
!       i4tol               Input      Tolerance of accepted file times
!       i4time_nearest      Output     Actual File time of returned data
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       var_2d              Input      Which Variable Do you Want?
!       ext                 Input      3 character file extension
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

!       Steve Albers            1990

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d

        real field_3d(imax,jmax,kmax)

        character*255 c_filespec

        call get_directory(ext,directory,len_dir)

        call s_len(ext,lenext)

        c_filespec = directory(1:len_dir)//'*.'//ext(1:lenext)
        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .le. i4tol)then
            call get_laps_3d(i4time_nearest,imax,jmax,kmax
     1                  ,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,istatus)
        else
            write(6,*)' No field found within window ',ext(1:10)
            istatus = 0
        endif

        return
        end

        subroutine put_laps_2d(i4time,EXT,var_2d,units_2d,
     1                  comment_2d,imax,jmax,field_2d,istatus)

cdoc    Writes a 2-D grid. Inputs include the extension and time.

        character*150 DIRECTORY
cc        character*31 EXT
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
cc        character*3 var_2d
        character*(*) var_2d
        integer LVL_2d
        character*4 LVL_COORD_2d

        real field_2d(imax,jmax)

        call get_directory(ext,directory,len_dir)
 
        call s_len(ext,len_ext)

        write(6,11)directory(1:len_dir),ext(1:len_ext),var_2d
11      format(' Writing 2d ',a,1x,a,1x,a)

        call check_nan2(field_2d,imax,jmax,istatus)
        if(istatus .ne. 1)then
            write(6,*)' ERROR: Nan Detected in this 2D field'
            return
        endif

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

        return
        end

        subroutine put_laps_3d(i4time,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,ni,nj,nk)

cdoc    Writes a 3-D grid. Inputs include the extension and time.

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(nk),comment_2d
        character*10 units_3d(nk),units_2d
        character*3 var_3d(nk),var_2d
        integer LVL_3d(nk)
        character*4 LVL_COORD_3d(nk)
        real ht_1d(nk)

        real field_3d(ni,nj,nk)

        call get_directory(ext,directory,len_dir)

        write(6,11)directory(1:min(len_dir,50)),ext,var_2d
11      format(' Writing 3d ',a50,1x,a5,1x,a3)

        do k = 1,nk
            units_3d(k)   = units_2d
            comment_3d(k) = comment_2d
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'HPA'
            elseif(ltest_vertical_grid('SIGMA_HT'))then
                if(k .eq. 1)then
                    write(6,*)' Reading Sigma Ht Levels'
                    call get_ht_1d(nk,ht_1d,istatus)
                    if(istatus .ne. 1)then
                        return
                    endif
                endif
                lvl_3d(k) = nint(ht_1d(k))
                lvl_coord_3d(k) = '  ' ! informational
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(k) = var_2d

        enddo ! k

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,ni,nj,
     1  nk,nk,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,field_3d,ISTATUS)

        return
        end

        subroutine put_laps_multi_3d(i4time,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,ni,nj,nk,nf,istatus)

cdoc    Writes multiple 3-D grids. Inputs include the extension and time.

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(nk*nf),comment_2d(nf)
        character*10 units_3d(nk*nf),units_2d(nf)
        character*3 var_3d(nk*nf),var_2d(nf)
        integer LVL_3d(nk*nf)
        character*4 LVL_COORD_3d(nk*nf)
        real ht_1d(nk)

        real field_3d(ni,nj,nk,nf)

        istatus = 0

        call get_directory(ext,directory,len_dir)

        do l = 1,nf
            write(6,11)directory(1:min(len_dir,50)),ext,var_2d(l)       
11          format(' Writing 3d ',a50,1x,a5,1x,a3)
        enddo ! l

        do l = 1,nf
          do k = 1,nk

            iscript_3d = (l-1) * nk + k

            units_3d(iscript_3d)   = units_2d(l)
            comment_3d(iscript_3d) = comment_2d(l)
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(iscript_3d) = zcoord_of_level(k)/10
                lvl_coord_3d(iscript_3d) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(iscript_3d) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(iscript_3d) = 'HPA'
            elseif(ltest_vertical_grid('SIGMA_HT'))then
                if(k .eq. 1)then
                    write(6,*)' Reading Sigma Ht Levels'
                    call get_ht_1d(nk,ht_1d,istatus)
                    if(istatus .ne. 1)then
                        return
                    endif
                endif
! for other vertical_grid, lvl_3d(iscript_3d). Hongli Jiang 11/2/2011
                lvl_3d(iscript_3d) = nint(ht_1d(k))
                lvl_coord_3d(iscript_3d) = 'M' ! informational
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(iscript_3d) = var_2d(l)

          enddo ! k
        enddo ! l

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,ni,nj,
     1  nk*nf,nk*nf,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,field_3d,ISTATUS)

        if(istatus .ne. 1)return

        istatus = 1

        return
        end

        subroutine put_compressed_multi_3d(i4time,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,ni,nj,nk,nf,istatus)       

cdoc    Writes multiple 3-D compressed grids. Inputs include the extension and
cdoc    time.

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(nk*nf),comment_2d(nf)
        character*10 units_3d(nk*nf),units_2d(nf)
        character*3 var_3d(nk*nf),var_2d(nf)
        integer LVL_3d(nk*nf)
        character*4 LVL_COORD_3d(nk*nf)
        real ht_1d(nk)

        real field_3d(ni,nj,nk,nf)

        istatus = 0

        call get_directory(ext,directory,len_dir)

        do l = 1,nf
            write(6,11)directory(1:min(len_dir,50)),ext,var_2d(l)       
11          format(' Writing 3d ',a50,1x,a5,1x,a3)
        enddo ! l

        do l = 1,nf
          do k = 1,nk

            iscript_3d = (l-1) * nk + k

            units_3d(iscript_3d)   = units_2d(l)
            comment_3d(iscript_3d) = comment_2d(l)
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(iscript_3d) = zcoord_of_level(k)/10
                lvl_coord_3d(iscript_3d) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(iscript_3d) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(iscript_3d) = 'HPA'
            elseif(ltest_vertical_grid('SIGMA_HT'))then
                if(k .eq. 1)then
                    write(6,*)' Reading Sigma Ht Levels'
                    call get_ht_1d(nk,ht_1d,istatus)
                    if(istatus .ne. 1)then
                        return
                    endif
                endif
                lvl_3d(k) = nint(ht_1d(k))
                lvl_coord_3d(k) = '  ' ! informational
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(iscript_3d) = var_2d(l)

          enddo ! k
        enddo ! l

        CALL WRITE_LAPS_COMPRESSED(I4TIME,DIRECTORY,EXT,ni,nj,
     1  nk*nf,nk*nf,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,field_3d,ISTATUS)

        if(istatus .ne. 1)return

        istatus = 1

        return
        end

        subroutine put_laps_multi_2d(i4time,EXT,var_a,units_a,
     1                          comment_a,field_2d,ni,nj,nf,istatus)       

cdoc    Writes multiple 2-D grids. Inputs include the extension and time.

        integer MAX_FIELDS
        parameter (MAX_FIELDS = 10)

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_a(nf) ! ,comment
        character*10 units_a(nf)
        character*3 var_a(nf)
        integer LVL_2d(MAX_FIELDS)
        character*4 LVL_COORD_2d(MAX_FIELDS)

        real field_2d(ni,nj,nf)

        istatus = 0

        if(nf .gt. MAX_FIELDS)then
            write(6,*)' Too many fields in put_laps_multi_2d'
        endif

        call get_directory(ext,directory,len_dir)
        call s_len(ext,len_ext)

        do l = 1,nf
            write(6,11)directory(1:len_dir),ext(1:len_ext),var_a(l)
11          format(' Writing 2d ',a,1x,a,1x,a)

            call check_nan2(field_2d(1,1,l),ni,nj,istatus)
            if(istatus .ne. 1)then
                write(6,*)' ERROR: Nan Detected in above 2D field'
                return
            endif
        enddo ! l

        do l = 1,nf
            lvl_2d (l)= 0
            lvl_coord_2d(l) = 'MSL'
!           comment_a(l) = comment
        enddo ! l

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,ni,nj,
     1  nf,nf,VAR_A,LVL_2D,LVL_COORD_2D,UNITS_A,
     1                     COMMENT_A,field_2d,ISTATUS)

        if(istatus .ne. 1)return

        istatus = 1

        return
        end


        subroutine sort_integer(i_array,i_dim,istatus)

        integer i_array(i_dim)

 10     i_switch = 0
        do i = 2,i_dim
            if(i_array(i) .lt. i_array(i-1))then ! Bubble sort exchange
                izzz = i_array(i-1)
                i_array(i-1) = i_array(i)
                i_array(i) = izzz
                i_switch = 1
            endif
        enddo

        if(i_switch .eq. 1)goto 10

        return
        end

        subroutine lapsprd_file_exist(i4time,ext,l_exist,istatus)

        character*(*)    ext
        character*150    directory
        character*13 filename13

        logical l_exist

!       Test for proper length of extension
        call s_len(ext,len_ext)
        if(len_ext .ne. 3)then
           write(6,*)' Error in open_lapsprd_file: ext has length'
     1               ,len_ext
           istatus = 0
           return
        endif

        call get_directory(ext,directory,len_dir)

        inquire(file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1         ,exist=l_exist)

 999    istatus = 1
        return
        end

        subroutine open_lapsprd_file(lun,i4time,ext,istatus)

!       1997   Steve Albers

        character*(*)    ext
        character*150    directory
        character*13 filename13
        character*14 filename14

!       Test for proper length of extension
        call s_len(ext,len_ext)
        if(len_ext .ne. 3 .and. len_ext .ne. 4)then
           write(6,*)' Error in open_lapsprd_file: ext has length'
     1               ,len_ext
           istatus = 0
           return
        endif

        call get_directory(ext,directory,len_dir)

        if(len_ext .eq. 3)then
         open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1          ,status='unknown',err=998)
        else ! length of 4
         open(lun,file=directory(1:len_dir)//filename14(i4time,ext(1:4))
     1          ,status='unknown',err=998)
        endif
        go to 999

 998    write(6,*)' Error in open_lapsprd_file, cannot open product: '   
     1            ,ext       
        istatus = 0
        return

 999    istatus = 1
        return
        end

        subroutine open_lapsprd_file_read(lun,i4time,ext,istatus)

!       1997   Steve Albers

        character*(*)    ext
        character*150    directory
        character*13 filename13

!       Test for proper length of extension
        call s_len(ext,len_ext)
        if(len_ext .ne. 3 .and. len_ext .ne. 4)then
           write(6,*)' Error in open_lapsprd_file_read: ext has length'       
     1               ,len_ext
           istatus = 0
           return
        endif

        call get_directory(ext,directory,len_dir)

        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))       
     1          ,status='old',err=998)
        go to 999

 998    write(6,*)' Warning in open_lapsprd_file_read: '
     1           ,'cannot open the product ',ext     
        istatus = 0
        return

 999    istatus = 1
        return
        end


        subroutine open_lapsprd_file_append(lun,i4time,ext,istatus)

!       1997   Steve Albers

        character*(*)    ext
        character*150    directory
        character*13 filename13
        integer istatus
!       Test for proper length of extension
        call s_len(ext,len_ext)
        if(len_ext .ne. 3 .and. len_ext .ne. 4)then
           write(6,*)' Error in open_lapsprd_file: ext has length'
     1               ,len_ext
           istatus = 0
           return
        endif

        call get_directory(ext,directory,len_dir)
        
        call open_append(lun,directory(1:len_dir)//
     +       filename13(i4time,ext(1:3)),'unknown',istatus)

        if(istatus.eq.0) then
           write(6,*)
     1          ' Error in open_lapsprd_file_append: ',
     2          'cannot open the file ',
     3          directory(1:len_dir)//filename13(i4time,ext(1:3))       
        else
           istatus = 1
        endif

        return
        end ! open_lapsprd_file_append

 
      subroutine open_ext(lun,i4time_sys,ext,istatus)

      integer iopen
      save iopen
      data iopen /0/

      character*3 ext

!     Open output lapsprd file
      if(iopen .eq. 0)then
          call open_lapsprd_file(lun,i4time_sys,ext,istatus)
          if(istatus .ne. 1)then
              write(6,*)' ERROR: could not open lapsprd file ',ext
              stop
          else
              write(6,*)' Successfully Opened lapsprd file for ',ext
          endif

          iopen = 1
      endif

      return
      end

        subroutine get_filespec(ext,mode,c_filespec,istatus)

!       1997   Steve Albers

        character*(*)    ext         ! input
        integer          mode        ! input
        character*(*)    c_filespec  ! output

        character*150    directory

        istatus = 1

        call s_len(ext,len_ext)
        call get_directory(ext,directory,len_dir)

        if(mode .eq. 1)then ! short form
            if(len(c_filespec) .lt. len_dir)then
                istatus = 0
                write(6,*)' Error in get_filespec, c_filespec too short'       
            else
                c_filespec = directory(1:len_dir)
            endif
        else                ! long form
            if(len(c_filespec) .lt. len_dir+len_ext+2)then
                istatus = 0
                write(6,*)' Error in get_filespec, c_filespec too short'
            else
                c_filespec = directory(1:len_dir)//'*.'//ext(1:len_ext)       
            endif
        endif

        return
        end
c
c J Smart 3/98.
c
        subroutine get_laps_sat(maxsat,c_sat_id,isats
     1     ,i4time_needed,i4tol,i4time_nearest
     1     ,var_2d,units_2d,comment_2d,imax,jmax
     1     ,subpoint_lat_clo,subpoint_lon_clo         ! O
     1     ,lat,lon,field_2d,istatus)
c
c    J. Smart 2-98.
cdoc This routine acquires satellite data from the lvd subdirectories
cdoc and makes decisions about what the best data is to return as the
cdoc 2d field (for var_2d).

        use mem_namelist, ONLY: l_mosaic_sat, r_missing_data
c
        implicit none

        integer imax,jmax,maxsat,min,max,iskip

        real    field_2d_lvd(imax,jmax,maxsat)
        real    field_2d(imax,jmax)
        real    lat(imax,jmax)
        real    lon(imax,jmax)
        real    subpoint_lat_clo(imax,jmax)
        real    subpoint_lon_clo(imax,jmax)
        real    subpoint_lon_sat(maxsat)
        real    cost(imax,jmax,maxsat)
        real    costmin(imax,jmax), fracp, time_diff_sec
        real    satm_time_wt,satm_subp_wt

        integer   isats(maxsat)
        integer   nsats
        integer   npts(maxsat)

        integer   i4time_needed
        integer   i4time_nearest
        integer   i4tol
        integer   i4timedata(maxsat)
        integer   i4time_min
        integer   i4time_sys
        integer   min_i4time

        integer   i,j,isat
        integer   imn
        integer   istatus 
        integer   jstatus

        character comment_2d*125
        character units_2d*10
        character asc9_tim*9
        character var_2d*3
        character c_sat_id(maxsat)*6   !satellite id's known to system
        character csatid(maxsat)*6   !satellite id's returned from routine

        integer   i4time_first
c       data      i4time_first/0/
c       save      i4time_first

        call get_systime_i4(i4time_sys,istatus)
        istatus=1   !data was found
        nsats=0
        npts=0
        i4time_first=0
        subpoint_lon_sat = -75.            ! fill satellite array (default)
        subpoint_lat_clo = 0.
        satm_time_wt = 1./60.
        satm_subp_wt = 1.0

        do i=1,maxsat
           if(isats(i).eq.1)then
              call get_laps_lvd(c_sat_id(i),
     &                 i4time_needed,i4tol,i4time_nearest,
     &                 i4time_sys,
     &                 var_2d,units_2d,comment_2d,
     &                 imax,jmax,field_2d,jstatus)

              if(jstatus.ne.1)then
                 write(6,*)'No data returned from get_laps_lvd',
     &                     ' for ',c_sat_id(i)
              else
                 nsats=nsats+1
                 csatid(nsats)=c_sat_id(i)
                 i4timedata(nsats)=i4time_nearest
                 call move(field_2d,field_2d_lvd(1,1,nsats),imax,jmax)
                 write(6,*)' Comment is:',trim(comment_2d)
                 read(comment_2d,1)subpoint_lon_sat(nsats)
 1               format(83x,e17.8)
                 write(6,*)' Sublon for satellite: ',nsats
     1                      ,subpoint_lon_sat(nsats)      

              endif
           endif
        enddo
c
c this section can make decisions about which satellite data
c to return in the event there is more than 1 2d field.
c
        write(6,*)' nsats / l_mosaic_sat = ',nsats,l_mosaic_sat

        if(nsats.gt.1 .and. l_mosaic_sat)then
           write(6,*)
           write(6,*)'Mosaicing ',nsats,' satellites'
           subpoint_lon_clo = r_missing_data  ! fill grid 

           cost = 9999. ! initialize array

           do isat=1,nsats

              time_diff_sec = i4timedata(isat) - i4time_sys   

              write(6,*)' isat/time_diff = ',isat,time_diff_sec

!             Determine cost of satellite based on valid data, time, and subpoint
              where(field_2d_lvd(:,:,isat) .ne. r_missing_data)
                 cost(:,:,isat) 
     1           = abs(lon(:,:) - subpoint_lon_sat(isat)) * satm_subp_wt
     1           + abs(time_diff_sec) * satm_time_wt
              end where

           enddo ! isat

           do i = 1,imax
           do j = 1,jmax
               costmin(i,j) = minval(cost(i,j,:))
           enddo ! j
           enddo ! i

           write(6,*)' Summary of cost function'
           j = jmax/2
           iskip = max(((imax/100)/5)*5,5)
           do i = 1,imax,iskip
               write(6,*)i,lon(i,j),cost(i,j,1:nsats)                  
           enddo ! i                      

           do isat=1,nsats
              call make_fnam_lp(i4timedata(isat),asc9_tim,istatus)

!             Add satellite to mosaic if it has the lowest cost
              do i = 1,imax
              do j = 1,jmax
                  if(field_2d_lvd(i,j,isat) .ne. r_missing_data 
     1                            .AND.
     1              cost(i,j,isat) .eq. costmin(i,j)
     1                                                           )then       
                    npts(isat) = npts(isat) + 1
                    subpoint_lon_clo(i,j) = subpoint_lon_sat(isat) 
                    field_2d(i,j) = field_2d_lvd(i,j,isat)
                  endif
              enddo ! j
              enddo ! i

              fracp = float(npts(isat)) / float(imax*jmax)

              write(6,11)var_2d,csatid(isat),asc9_tim 
     &                 ,subpoint_lon_sat(isat),npts(isat),fracp  
11            format('Mosaicing ',a,' for ',a,1x,a,f8.2
     1                           ,' npts/frac=',i8,f10.6)

           enddo ! isat
           return

        elseif(nsats .gt.1 .and. (.not. l_mosaic_sat))then
           write(6,*)
           write(6,*)'Found data for ',nsats,' satellites'
           do i=1,nsats
              if(i4timedata(i).eq.i4time_first)then
                 call make_fnam_lp(i4timedata(i),asc9_tim,istatus)
                 call move(field_2d_lvd(1,1,i),field_2d,imax,jmax)
                 i4time_nearest=i4timedata(i)
                 subpoint_lon_clo = subpoint_lon_sat(i)
                 write(6,*)'Returning ',var_2d,' for ',csatid(i),
     &                     ' ',asc9_tim,' 1'
                 return
              endif
           enddo

           min_i4time=i4time_sys
           do i=1,nsats
              i4time_min=abs(i4timedata(i)-i4time_sys)
              if(i4time_min.lt.min_i4time)then
                 imn=i
                 min_i4time=i4time_min
              endif
           enddo

           if(imn.ne.0)then
              i4time_first=i4timedata(imn)
              call make_fnam_lp(i4timedata(imn),asc9_tim,istatus)
              call move(field_2d_lvd(1,1,imn),field_2d,imax,jmax)
              i4time_nearest=i4timedata(imn)
              subpoint_lon_clo = subpoint_lon_sat(imn)
              write(6,*)'Returning ',var_2d,' for ',csatid(imn),
     &                  ' ',asc9_tim,' 2'
              return

           else
c default
              i4time_first=i4timedata(1)
              call make_fnam_lp(i4time_first,asc9_tim,istatus)
              call move(field_2d_lvd(1,1,1),field_2d,imax,jmax)
              i4time_nearest=i4timedata(1)
              subpoint_lon_clo = subpoint_lon_sat(1)
              write(6,*)'Returning ',var_2d,' for ',csatid(1),
     &                  ' ',asc9_tim,' 3'
              return
           endif

        elseif(nsats.eq.1)then
           i4time_first=i4timedata(1)
           call make_fnam_lp(i4time_first,asc9_tim,istatus)
           call move(field_2d_lvd(1,1,1),field_2d,imax,jmax)
           i4time_nearest=i4timedata(1)
           subpoint_lon_clo = subpoint_lon_sat(1)
           write(6,*)
           write(6,*)'Returning ',var_2d,' for ',csatid(1), 
     &                ' ',asc9_tim,' 4'
           return

        elseif(nsats.le.0)then
           write(6,*)'No lvd fields found. Returning  no data'

        endif


        istatus=0
        return
        end
c
c--------------------------------------------------------------------
c
        subroutine get_laps_lvd(c_sat_id,
     1      i4time_needed,i4tol,i4time_nearest
     1     ,i4time_sys
     1     ,var_2d,units_2d,comment_2d
     1     ,imax,jmax,field_2d,istatus)

!       Steve Albers            1996
!       This routine tries to read in the desired variable from all files
!       having the proper extension, picking the closest one within the
!       specified time window.
!
!       John Smart              1998
!       Modified original get_laps_2dvar for satellite lvd files in
!       lvd subdirectories. Subdirectories are those satellites
!       known to the laps system (see src/include/satdata_lvd.for,
!       and data/static/satellite_lvd.nl). As many as 'maxsat' 2d fields
!       can be returned depending on configuration specified in satellite_lvd.nl.
!       max_sat is defined in src/include/satellite_dims_lvd.inc.
!
!       John Smart              1999
!       Added multiple lvd file search for best domain cover for files that
!       are within i4tol. If all files have the same domain cover (like 100%)
!       then the first (and closest in time) file is returned.
!
!       John Smart              1999
!       return the time closest to requested time in the event that more
!       than 1 file with complete domain cover satisfies i4tol

        character*150 dir
        character*150 satdir
        character*31  EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer LVL_2d
        character*4 LVL_COORD_2d

        integer max_lvd_files
        parameter (max_lvd_files=20)

        real field_2d(imax,jmax)
        real field_2d_save(imax,jmax,max_lvd_files)

        integer max_files
        parameter (max_files = 600)

        character*9 asc9_tim
        character*9 asc9_time(max_files)

        real pctmiss(max_files)
        real pctmnmiss
        real r_missing_data
        integer jf,jr,imiss

        character*255 c_filespec
        character*125 comment_2d_save(max_lvd_files)
        character*120 c_fnames(max_files)
        integer i4times(max_files)
        integer i4times_data(max_lvd_files)
        integer i_selected(max_files)
        character*6 c_sat_id          !input satellite id's known to system
        logical lcont

        call make_fnam_lp(i4time_needed,asc9_tim,istatus)

        write(6,*)
        write(6,*)' Subroutine get_laps_lvd for ',c_sat_id,
     1            ' needed at ',asc9_tim,'...'

        ext = 'lvd'
        call get_directory(ext,dir,ldir)

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lext = i-1

        do j = 1,max_files
           i_selected(j) = 0
        enddo ! j

        satdir=dir(1:ldir)//trim(c_sat_id)//'/'
        lsdir=index(satdir,' ')-1
        c_filespec = satdir(1:lsdir)//'*.'//ext(1:lext)

        call get_r_missing_data(r_missing_data,istat)
        call get_file_times(c_filespec,max_files,c_fnames
     1                     ,i4times,i_nbr_files_ret,istatus)
        if(istatus .ne. 1)then
           write(6,*)'get_laps_lvd: Bad status returned '
     1              ,'from get_file_times'
           return
        endif

        lcont=.true.
        jf=0

        print*,'Number of files returned get_file_times',i_nbr_files_ret

!       Calculate minimum time difference of unselected files from desired time
!       with respect to both needed time and system time
50      i4_diff_min = 999999999
        i4_diff_min_sys = 999999999
        do j = 1,i_nbr_files_ret
           i4_diff = abs(i4times(j) - i4time_needed)
           i4_diff_sys=abs(i4times(j)-i4time_sys)      
           if(i_selected(j) .eq. 0)then
            i4_diff_min     = min(i4_diff    ,i4_diff_min)
            i4_diff_min_sys = min(i4_diff_sys,i4_diff_min_sys)
           endif
        enddo ! j

!       Check whether any unselected files lie within the time window
        if(i4_diff_min_sys .gt. i4tol)then
           write(6,*)' No remaining files found within ',i4tol
     1              ,' sec systime window ',ext(1:5),var_2d
           write(6,*)' Closest file (in seconds) is ',i4_diff_min_sys

           lcont=.false.

c          istatus = 0
c          return
        endif

        do j=1,i_nbr_files_ret

           i4_diff=    abs(i4times(j)-i4time_needed)
           i4_diff_sys=abs(i4times(j)-i4time_sys)      

!          Select file that has minimum time diff (of unselected)
           if(i4_diff.eq.i4_diff_min.and.
     1            i_selected(j).eq.0.and.
     1            i4_diff_sys.le.i4tol)               then

              i_selected(j) = 1
              lvl_2d = 0
              lvl_coord_2d = 'MSL'
              call make_fnam_lp(i4times(j),asc9_tim,istatus)

              write(6,11)trim(satdir),asc9_tim,ext(1:5),var_2d
11            format(' Reading 2d ',a,1x,a9,1x,a5,1x,a3)

              CALL READ_LAPS_DATA(i4times(j),satdir,EXT,imax
     1            ,jmax,1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D
     1            ,COMMENT_2D,field_2d,ISTATUS)

              if(istatus .ne. 1)then
                 write(6,*)' No field found at ',ext(1:10)
     1                       ,var_2d,' ',asc9_tim
                 go to 50
c
c we need to expect some missing data in the lvd fields.
c determine which lvd file within the time window has the best
c domain coverage.
c
              else   !  istatus = 1, check for domain cover maxima

                 jf=jf+1
                 if(jf .gt. max_lvd_files)then
                     write(6,*)' ERROR in get_laps_lvd: jf > '
     1                        ,max_lvd_files
                     stop
                 endif

                 imiss=0
                 do il = 1,imax
                 do jl = 1,jmax

!                   QC check
                    if(field_2d(il,jl) .gt. 1e10)then
                       field_2d(il,jl) = r_missing_data
                    endif

                    if(field_2d(il,jl) .eq. r_missing_data)then
                       imiss=imiss+1

c                           write(6,*)il,jl,
c    1                        ' Missing Data Value Detected in 2D Field'
c                           istatus = -1
c                           return

                    endif
                    field_2d_save(il,jl,jf)=field_2d(il,jl)
                 enddo ! j
                 enddo ! i

                 pctmiss(jf)=float(imiss)/float(imax*jmax)
                 asc9_time(jf)=asc9_tim
                 i4times_data(jf)=i4times(j)
                 comment_2d_save(jf)=comment_2d
                 write(6,*)'         pctmiss = ',pctmiss(jf)

              endif

c             return

              if(lcont)goto 50

           endif ! File is closest unread file to desired time

        enddo ! ith file
c
c full domain coverage takes precedence over time matching.
c
        pctmnmiss=1.0
        jr=0
        if(jf.gt.1)then
           do i=1,jf
              if(pctmiss(i).le.pctmnmiss)then
                 jr=jr+1
                 pctmnmiss=pctmiss(i)
              endif
           enddo
           if(jr.gt.1)then
              i4diffmn=1000.
              do i=1,jr
                 i4diff=abs(i4times_data(i)-i4time_needed)
                 if(i4diff.lt.i4diffmn)then
                    jf=i
                    i4diffmn=i4diff
                 endif
              enddo
           else
              jf=jr
           endif
        elseif(jf.eq.0)then
           istatus = 0
           return
        endif

        call move(field_2d_save(1,1,jf),field_2d,imax,jmax)
        print*,'Returning requested field from ',
     1' get_laps_lvd  ',c_sat_id,'   ',asc9_time(jf),'   ',var_2d
        i4time_nearest=i4times_data(jf)
        comment_2d=comment_2d_save(jf)
        istatus = 1

        return 
        end


        subroutine get_3dgrid_dname(directory_in
     1         ,i4time_needed,i4tol,i4time_nearest
     1         ,EXT,var_2d,units_2d
     1         ,comment_2d,imax,jmax,kmax,field_3d,istatus)

cdoc    Returns a 3-D analysis grid. Inputs include a directory, ext, and time window.

!       directory_in        Input      Slash at end is optional
!       i4time_needed       Input      Desired i4time
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       ext                 Input      3 character file extension
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

        character*9 asc9_tim

        character*(*) DIRECTORY_IN
        character*255 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(kmax),comment_2d
        character*10 units_3d(kmax),units_2d
        character*3 var_3d(kmax),var_2d
        integer LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real field_3d(imax,jmax,kmax)

        character*255 c_filespec

        logical ltest_vertical_grid

        write(6,*)' Subroutine get_3dgrid_dname...'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_3dgrid_dname: bad istatus, return'
            return
        endif

        call s_len(ext,len_ext)
        call s_len(directory_in,len_dir)

        if(directory_in(len_dir:len_dir) .ne. '/')then
            directory = directory_in(1:len_dir)//'/'
            len_dir = len_dir + 1
        else
            directory = directory_in
        endif

        c_filespec = directory(1:len_dir)//'*.'//ext(1:len_ext)
        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .le. i4tol)then
            call get_3d_dir_time(directory_in,i4time_nearest,EXT
     1                          ,var_2d,units_2d,comment_2d
     1                          ,imax,jmax,kmax,field_3d,istatus)

        else
            write(6,*)' No field found within window ',ext(1:10)
            istatus = 0
        endif

        return
        end


        subroutine get_3d_dir_time(directory_in,i4time
     1                            ,EXT,var_2d,units_2d
     1                            ,comment_2d
     1                            ,imax,jmax,kmax,field_3d,istatus)

cdoc    Returns a 3-D grid. Inputs include a directory, ext, and time.

!       directory_in        Input      Slash at end is optional
!       i4time              Input      Desired i4time
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       ext                 Input      file extension
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

        character*9 asc9_tim

        character*(*) DIRECTORY_IN
        character*255 DIRECTORY
        character*(*) EXT, var_2d

        character*125 comment_3d(kmax),comment_2d
        character*10 units_3d(kmax),units_2d
        character*3 var_3d(kmax)
        integer LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real field_3d(imax,jmax,kmax)

        character*255 c_filespec

        logical ltest_vertical_grid, l_is_vxx, l_compress_radar

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_3d_dir_time: bad istatus, return'
            return
        endif

        call s_len(ext,len_ext)
        call s_len(directory_in,len_dir)

        if(directory_in(len_dir:len_dir) .ne. '/')then
            directory = directory_in(1:len_dir)//'/'
            len_dir = len_dir + 1
        else
            directory = directory_in
        endif


        call make_fnam_lp(i4time,asc9_tim,istatus)

        len_high = max(45,len_dir)
        len_low = len_high - 44

        write(6,11)directory(len_low:len_high),asc9_tim,ext,var_2d
11      format('  get_3d_dir_time: ',a,1x,a,1x,a,1x,a)

        do k = 1,kmax
            units_3d(k)   = units_2d
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'MB'
            elseif(ltest_vertical_grid('SIGMA_HT'))then
                lvl_3d(k) = nint(zcoord_of_level(k))
                lvl_coord_3d(k) = 'M'
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(k) = var_2d

        enddo ! k

        call get_l_compress_radar(l_compress_radar,istatus)
        if(istatus .ne. 1)return

        if(.not. l_is_vxx(EXT))then
          CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1          kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1          COMMENT_3D,field_3d,ISTATUS)

        elseif(.not. l_compress_radar)then
          CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1          kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1          COMMENT_3D,field_3d,ISTATUS)
       
          if(ISTATUS .ne. 1)then
            write(6,*)' Attempting compressed radar data read instead'      
            call read_laps_compressed(i4time,directory,ext,imax,jmax
     1                               ,kmax,var_3d,lvl_3d,lvl_coord_3d       
     1                               ,units_3d,comment_3d,field_3d
     1                               ,istatus)
          endif

        else ! l_compress_radar is TRUE for vxx file
          call read_laps_compressed(i4time,directory,ext,imax,jmax
     1                               ,kmax,var_3d,lvl_3d,lvl_coord_3d       
     1                               ,units_3d,comment_3d,field_3d
     1                               ,istatus)

          if(ISTATUS .ne. 1)then
            write(6,*)' Attempting netCDF radar data read instead'      
            CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1          kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1          COMMENT_3D,field_3d,ISTATUS)
          endif

        endif

        comment_2d=comment_3d(1)
        units_2d=units_3d(1)

        return
        end

        function l_is_vxx(ext)

        character*(*) ext
        logical l_is_vxx

        l_is_vxx = .false.

        if(ext(1:1) .ne. 'v')return
        if(ext(2:2) .eq. '0')l_is_vxx = .true.
        if(ext(2:2) .eq. '1')l_is_vxx = .true.
        if(ext(2:2) .eq. '2')l_is_vxx = .true.
        if(ext(2:2) .eq. '3')l_is_vxx = .true.
        if(ext(2:2) .eq. '4')l_is_vxx = .true.
        if(ext(2:2) .eq. '5')l_is_vxx = .true.
        if(ext(2:2) .eq. '6')l_is_vxx = .true.
        if(ext(2:2) .eq. '7')l_is_vxx = .true.
        if(ext(2:2) .eq. '8')l_is_vxx = .true.
        if(ext(2:2) .eq. '9')l_is_vxx = .true.

        return
        end
c
c ------------------------------------------------------------------------
c This routine is used in get_maps_lapsgrid.f (subroutine get_modelfg_3d and
c other related routines). Returns bgmodel name corresponding to the
c setting of variable "cmodel" in lga namelist (background.nl).
c
      subroutine bgmodel_name(maxbgmodels,nbgm,bgmodelnames,istatus)

      implicit none

c     include 'bgdata.inc'

      integer       maxbgmodels
      character*256 bgpaths(maxbgmodels)
      character*132 cmodel(maxbgmodels)
      character*5   bgmodelnames(maxbgmodels)

      integer bgmodels(maxbgmodels)
      integer oldest_forecast
      integer max_forecast_delta
      integer forecast_length
      integer itime_inc,ntmin,ntmax
      integer nbgm
      logical sfc_bkgd
      logical use_analysis
      logical use_forecast
      logical smooth_fields
      logical lgb_only

      integer   istatus,i

      call get_background_info(bgpaths,bgmodels
     +,forecast_length
     +,use_analysis,use_forecast
     +,cmodel,itime_inc,smooth_fields
     +,sfc_bkgd
     +,ntmin,ntmax
     +,lgb_only)

      nbgm=0
      do i=1,maxbgmodels
         if(cmodel(i).ne.' ')then
            if(cmodel(i).eq.'RUC40_NATIVE'.or.
     +         cmodel(i).eq.'RUC20_NATIVE'.or.
     +         cmodel(i).eq.'RUC'              )then
               nbgm=nbgm+1
               bgmodelnames(nbgm)='ruc'
            endif
            if(cmodel(i).eq.'ETA48_CONUS'.or.
     +         cmodel(i).eq.'MesoEta_SBN')then
               nbgm=nbgm+1
               bgmodelnames(nbgm)='eta'
            endif
            if(cmodel(i).eq.'AVN_FSL_NETCDF'.or.
     +         cmodel(i).eq.'AVN_SBN_CYLEQ')then
               nbgm=nbgm+1
               bgmodelnames(nbgm)='avn'
            endif
            if(cmodel(i).eq.'MODEL_FUA'.or.
     +         cmodel(i).eq.'LAPS_FUA')then
               nbgm=nbgm+1
               bgmodelnames(nbgm)='fua'
            endif
            if(cmodel(i).eq.'CWB_20FA_LAMBERT_NF')then
               nbgm=nbgm+1
               bgmodelnames(nbgm)='nfs'
            endif
         endif
      enddo

      return
      end

cdis

        subroutine get_laps_multi_2d(i4time,EXT,var_2d,units_2d,
     1                  comment_2d,imax,jmax,nf,field_2d,istatus)

cdoc    Used to read in one or more surface grids with inputs of time and ext

        character*9 asc9_tim
        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d(nf)
        character*10 units_2d(nf)
        character*3 var_2d(nf)
        integer LVL_2d(nf)
        character*4 LVL_COORD_2d(nf)

        real field_2d(imax,jmax,nf)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_multi_2d: bad istatus, return'
            return
        endif

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)trim(directory),asc9_tim,ext,var_2d
11      format(' Reading 2d ',a,1x,a,1x,a,1x,a)

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  nf,nf,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!       Check for missing data
        do if = 1,nf
        do j = 1,jmax
        do i = 1,imax
            if(istatus .eq. 1)then
                if(field_2d(i,j,if) .eq. r_missing_data)then
                    write(6,*)' Missing Data Value Detected in 2D Field'
                    istatus = -1
                endif
            endif
        enddo ! i
        enddo ! j
        enddo ! if

        return
        end
