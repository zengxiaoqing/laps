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

!       This routine can be used to read in a surface grid of known time

        character*9 asc9_tim
        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_2d: bad istatus, return'
            return
        endif

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)directory(1:45),asc9_tim,ext,var_2d
11      format(' Reading 2d ',a,1x,a,1x,a,1x,a)

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!       Spot check for missing data
        if(istatus .eq. 1 .and.
     1            field_2d(imax/2,jmax/2) .eq. r_missing_data)then
            write(6,*)' Missing Data Value Detected in 2D Field'
            istatus = -1
        endif

        return
        end

        subroutine get_lapsdata_2d(i4time,i4_valid,EXT,var_2d,units_2d,
     1                  comment_2d,imax,jmax,field_2d,istatus)

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

!       This routine can be used to read in a surface grid of known time
!       by calling the new READ_LAPS routine

        character*9 asc9_tim
        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_lapsdata_2d: bad istatus, return'
            return
        endif

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)directory(1:45),asc9_tim,ext(1:5),var_2d
11      format(' Reading 2d ',a45,1x,a9,1x,a5,1x,a3)

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL READ_LAPS(I4TIME,i4_valid,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!       Spot check for missing data
        if(istatus .eq. 1 .and.
     1            field_2d(imax/2,jmax/2) .eq. r_missing_data)then
            write(6,*)' Missing Data Value Detected in 2D Field'
            istatus = -1
        endif

        return
        end

        subroutine get_laps_2dgrid(i4time_needed,i4tol,i4time_nearest
     1         ,EXT,var_2d,units_2d
     1         ,comment_2d,imax,jmax,field_2d,ilevel,istatus)

!       Steve Albers            1990

        character*9 asc9_tim

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        character*255 c_filespec

        logical ltest_vertical_grid

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_2dgrid: bad istatus, return'
            return
        endif

        call get_directory(ext,directory,len_dir)

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lenext = i-1

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

            write(6,11)directory(1:45),asc9_tim,ext(1:5),var_2d
11          format(' Reading 2d ',a45,1x,a9,1x,a5,1x,a3)

            CALL READ_LAPS_DATA(I4TIME_nearest,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!           Spot check for missing data
            if(istatus .eq. 1 .and.
     1                 field_2d(imax/2,jmax/2) .eq. r_missing_data)then
                write(6,*)' Missing Data Value Detected in 2D Field'
                istatus = -1
            endif

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

!       Steve Albers            1990
!           J Smart             1998

        character*9 asc9_tim

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

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

            write(6,11)directory(1:52),asc9_tim,ext(1:5),var_2d
11          format(' Reading 2d ',a52,1x,a9,1x,a5,1x,a3)

            CALL READ_LAPS_DATA(I4TIME_nearest,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

!           Spot check for missing data
            if(istatus .eq. 1 .and.
     1                 field_2d(imax/2,jmax/2) .eq. r_missing_data)then
                write(6,*)' Missing Data Value Detected in 2D Field'
c               istatus = -1
            endif

        else
            write(6,*)' No field found within window ',ext(1:10)
            istatus = 0

        endif

        return
        end




        subroutine get_laps_2dvar(i4time_needed,i4tol,i4time_nearest
     1         ,EXT,var_2d,units_2d
     1         ,comment_2d,imax,jmax,field_2d,ilevel,istatus)

!       Steve Albers            1996
!       This routine tries to read in the desired variable from all files
!       having the proper extension, picking the closest one within the
!       specified time window.
!
!       J Smart                 1998
!       added lvd subdirectory flexibility. Only one 2d satellite field returned.

        character*9 asc9_tim

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        integer max_files
        parameter (max_files = 600)
        character*255 c_filespec
        character*120 c_fnames(max_files)
        integer i4times(max_files)
        integer i_selected(max_files)

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
     1     ,field_2d,istatus)

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

                write(6,11)directory(1:45),asc9_tim,ext(1:5),var_2d
11              format(' Reading 2d ',a45,1x,a9,1x,a5,1x,a3)

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

!       Returns a 3D laps grid
!       i4time              Input      Desired i4time
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       ext                 Input      3 character file extension
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

!       Steve Albers            1990

        character*150 DIRECTORY
cc        character*31 EXT
        character*(*) EXT, var_2d

        character*125 comment_3d(kmax),comment_2d
        character*10 units_3d(kmax),units_2d
cc        character*3 var_3d(kmax),var_2d
        character*3 var_3d(kmax)
        integer*4 LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real*4 field_3d(imax,jmax,kmax)

        character*9 asc9_tim

        logical ltest_vertical_grid

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)directory(1:45),asc9_tim,ext,var_2d
11      format(' Reading 3D ',a,1x,a,1x,a,1x,a)

        do k = 1,kmax
            units_3d(k)   = units_2d
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'MB'
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(k) = var_2d

        enddo ! k

        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,field_3d,ISTATUS)

        comment_2d=comment_3d(1)
        units_2d=units_3d(1)

        return
        end

        subroutine get_lapsdata_3d(i4time,i4_valid,imax,jmax,kmax
     1  ,EXT,var_2d,units_2d,comment_2d,field_3d,istatus)

!       Returns a 3D laps grid
!       i4time              Input      Desired i4time initial
!       i4_valid            Input      i4time for valid data time
!       imax,jmax,kmax      Input      LAPS grid dimensions
!       ext                 Input      3 character file extension
!       var_2d              Input      Which Variable do you want?
!       units_2d            Output     Units of data
!       Comment_2d          Output     Comment block
!       field_3d            Output     3D grid

!       Steve Albers            1990

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(kmax),comment_2d
        character*10 units_3d(kmax),units_2d
        character*3 var_3d(kmax),var_2d
        integer*4 LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real*4 field_3d(imax,jmax,kmax)

        logical ltest_vertical_grid

        call get_directory(ext,directory,len_dir)

        call s_len(ext,len)
        write(6,11)directory(1:len_dir),ext(1:len),var_2d
11      format(' Reading 3d ',a,1x,a5,1x,a3)

        do k = 1,kmax
            units_3d(k)   = units_2d
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'MB'
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif


            var_3d(k) = var_2d

        enddo ! k

        CALL READ_LAPS(I4TIME,I4_VALID,DIRECTORY,EXT,imax,jmax,
     1  kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,field_3d,ISTATUS)

        comment_2d=comment_3d(1)
        units_2d=units_3d(1)

        return
        end

        subroutine get_laps_3dgrid(i4time_needed,i4tol,i4time_nearest,
     1          imax,jmax,kmax,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,istatus)

!       Returns a 3D laps grid
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

        real*4 field_3d(imax,jmax,kmax)

        character*255 c_filespec

        call get_directory(ext,directory,len_dir)

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lenext = i-1

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

        character*150 DIRECTORY
cc        character*31 EXT
        character*(*) EXT

        character*125 comment_2d
        character*10 units_2d
cc        character*3 var_2d
        character*(*) var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        call get_directory(ext,directory,len_dir)

        write(6,11)directory,ext,var_2d
11      format(' Writing 2d ',a50,1x,a5,1x,a3)

        lvl_2d = 0
        lvl_coord_2d = 'MSL'

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,field_2d,ISTATUS)

        return
        end

        subroutine put_laps_3d(i4time,EXT,var_2d,units_2d,
     1                          comment_2d,field_3d,ni,nj,nk)

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(nk),comment_2d
        character*10 units_3d(nk),units_2d
        character*3 var_3d(nk),var_2d
        integer*4 LVL_3d(nk)
        character*4 LVL_COORD_3d(nk)

        real*4 field_3d(ni,nj,nk)

        call get_directory(ext,directory,len_dir)

        write(6,11)directory,ext(1:5),var_2d
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

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(nk*nf),comment_2d(nf)
        character*10 units_3d(nk*nf),units_2d(nf)
        character*3 var_3d(nk*nf),var_2d(nf)
        integer*4 LVL_3d(nk*nf)
        character*4 LVL_COORD_3d(nk*nf)

        real*4 field_3d(ni,nj,nk,nf)

        istatus = 0

        call get_directory(ext,directory,len_dir)

        do l = 1,nf
            write(6,11)directory,ext(1:5),var_2d(l)
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

        subroutine put_laps_multi_2d(i4time,EXT,var_a,units_a,
     1                          comment_a,field_2d,ni,nj,nf,istatus)

        integer*4 MAX_FIELDS
        parameter (MAX_FIELDS = 10)

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_a(nf)
        character*10 units_a(nf)
        character*3 var_a(nf)
        integer*4 LVL_2d(MAX_FIELDS)
        character*4 LVL_COORD_2d(MAX_FIELDS)

        real*4 field_2d(ni,nj,nf)

        istatus = 0

        if(nf .gt. MAX_FIELDS)then
            write(6,*)' Too many fields in put_laps_multi_2d'
        endif

        call get_directory(ext,directory,len_dir)

        do l = 1,nf
            write(6,11)directory,ext(1:5),var_a(l)
11          format(' Writing 2d ',a50,1x,a5,1x,a3)
        enddo ! l

        do l = 1,nf
            lvl_2d (l)= 0
            lvl_coord_2d(l) = 'MSL'
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

        subroutine open_lapsprd_file(lun,i4time,ext,istatus)

!       1997   Steve Albers

        character*(*)    ext
        character*150    directory
        character*13 filename13

!       Test for proper length of extension
        call s_len(ext,len_ext)
        if(len_ext .ne. 3)then
           write(6,*)' Error in open_lapsprd_file: ext has length'
     1               ,len_ext
           istatus = 0
           return
        endif

        call get_directory(ext,directory,len_dir)

        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))       
     1          ,status='unknown',err=998)
        go to 999

 998    write(6,*)' Error in open_lapsprd_file: cannot open the product'
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
        if(len_ext .ne. 3)then
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
        if(len_ext .ne. 3)then
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
     2          'cannot open the file',
     3          directory(1:len_dir)//filename13(i4time,ext(1:3))       
        else
           istatus = 1
        endif

        return
        end ! open_lapsprd_file_append




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
     1     ,field_2d,istatus)
c
c J. Smart 2-98.
c This routine acquires satellite data from the lvd subdirectories
c and makes decisions about what the best data is to return as the
c 2d field (for var_2d).
c
        implicit none

        integer imax,jmax,maxsat

        real*4    field_2d_lvd(imax,jmax,maxsat)
        real*4    field_2d(imax,jmax)

        integer   isats(maxsat)
        integer   nsats

        integer   i4time_needed
        integer   i4time_nearest
        integer   i4tol

        integer   i
        integer   istatus 
        integer   jstatus

        character comment_2d*125
        character units_2d*10
        character var_2d*3
        character c_sat_id(maxsat)*6   !satellite id's known to system
        character csatid(maxsat)*6   !satellite id's returned from routine

        istatus=0   !No data found
        nsats=0
        do i=1,maxsat
           if(isats(i).eq.1)then
              call get_laps_lvd(c_sat_id(i),
     &                 i4time_needed,i4tol,i4time_nearest,
     &                 var_2d,units_2d,comment_2d,
     &                 imax,jmax,field_2d,jstatus)

              if(jstatus.ne.1)then
                 write(6,*)'No data returned from get_laps_lvd',
     &               ' for ',c_sat_id(i)
              else
                 nsats=nsats+1
                 csatid(nsats)=c_sat_id(i)
                 call move(field_2d,field_2d_lvd(1,1,nsats),imax,jmax)
              endif
           endif
        enddo
c
c this section can make decisions about which satellite data
c to return in the event there is more than 1 2d field.
c
        if(nsats.gt.1)then
           write(6,*)'Found data for ',nsats,' satellites'
           write(6,*)'Returning ',var_2d,' for ',csatid(1),' only'
        elseif(nsats.eq.1)then
           write(6,*)'Found ',var_2d,' for ',nsats,' satellite'
        elseif(nsats.le.0)then
           write(6,*)'No lvd fields found. Returning  no data'
           return
        endif

c       call move_3dto2d(field_2d_lvd,1,field_2d,imax,jmax,maxsat)
        call move(field_2d_lvd(1,1,1),field_2d,imax,jmax)

        istatus=1
        return
        end
c
c--------------------------------------------------------------------
c
        subroutine get_laps_lvd(c_sat_id,
     1      i4time_needed,i4tol,i4time_nearest
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

        character*9 asc9_tim

        character*150 dir
        character*150 satdir
        character*31  EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        integer max_files
        parameter (max_files = 600)
        character*255 c_filespec
        character*120 c_fnames(max_files)
        integer i4times(max_files)
        integer i_selected(max_files)
        character*6 c_sat_id          !input satellite id's known to system

        ext = 'lvd'
        call get_directory(ext,dir,ldir)

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lext = i-1

        do j = 1,max_files
           i_selected(j) = 0
        enddo ! j

        satdir=dir(1:ldir)//c_sat_id//'/'
        lsdir=index(satdir,' ')-1
        c_filespec = satdir(1:lsdir)//'*.'//ext(1:lext)

        call get_file_times(c_filespec,max_files,c_fnames
     1                     ,i4times,i_nbr_files_ret,istatus)
        if(istatus .ne. 1)then
           write(6,*)'get_laps_2dvar: Bad status returned '
     1              ,'from get_file_times'
           return
        endif

50      i4_diff_min = 999999999
        do j = 1,i_nbr_files_ret
           i4_diff = abs(i4times(j) - i4time_needed)
           if(i_selected(j) .eq. 0)then
            i4_diff_min = min(i4_diff,i4_diff_min)
           endif
        enddo ! j

        if(i4_diff_min .gt. i4tol)then
           write(6,*)' No remaining files found within ',i4tol
     1              ,' sec time window ',ext(1:5),var_2d

           istatus = 0
           return
        endif

        do j=1,i_nbr_files_ret

           i4_diff=abs(i4times(j)-i4time_needed)
           if(i4_diff .eq. i4_diff_min .and. i_selected(j) .eq. 0)then

              i_selected(j) = 1
              lvl_2d = 0
              lvl_coord_2d = 'MSL'
              call make_fnam_lp(i4times(j),asc9_tim,istatus)

              write(6,11)satdir(1:lsdir),asc9_tim,ext(1:5),var_2d
11            format(' Reading 2d ',a61,1x,a9,1x,a5,1x,a3)

              CALL READ_LAPS_DATA(i4times(j),satdir,EXT,imax
     1            ,jmax,1,1,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D
     1            ,COMMENT_2D,field_2d,ISTATUS)

              if(istatus .ne. 1)then
                 write(6,*)' No field found at ',ext(1:10)
     1                       ,var_2d,' ',asc9_tim
                 go to 50

c we need to expect some missing data in the lvd fields.
c
c             else   !  istatus = 1, check for missing data
c                do il = 1,imax
c                do jl = 1,jmax
c                   if(field_2d(il,jl) .eq. r_missing_data)then
c                           write(6,*)il,jl,
c    1                        ' Missing Data Value Detected in 2D Field'
c                           istatus = -1
c                           return
c                   endif
c                enddo ! j
c                enddo ! i

              endif

              return

           endif ! File is closest unread file to desired time

        enddo ! ith file

        return 
        end

