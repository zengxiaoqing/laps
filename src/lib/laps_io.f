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

        include 'lapsparms.inc'

        character*9 asc9_tim
        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)directory(1:45),asc9_tim,ext(1:5),var_2d
11      format(' Reading 2d ',a45,1x,a9,1x,a5,1x,a3)

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

        include 'lapsparms.inc'

        character*9 asc9_tim
        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

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

        include 'lapsparms.inc'

        character*9 asc9_tim

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        character*255 c_filespec

        call get_directory(ext,directory,len_dir)

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lenext = i-1

        c_filespec = directory(1:len_dir)//'*.'//ext(1:lenext)
        call get_file_time(c_filespec,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .le. i4tol)then
            if(ilevel .ne. 0)then
                if(vertical_grid .eq. 'HEIGHT')then
                    lvl_2d = zcoord_of_level(k)/10
                    lvl_coord_2d = 'MSL'
                elseif(vertical_grid .eq. 'PRESSURE')then
                    lvl_2d = ilevel
                    lvl_coord_2d = 'MB'
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



        subroutine get_laps_2dvar(i4time_needed,i4tol,i4time_nearest
     1         ,EXT,var_2d,units_2d
     1         ,comment_2d,imax,jmax,field_2d,ilevel,istatus)

!       Steve Albers            1996
!       This routine tries to read in the desired variable from all files
!       having the proper extension, picking the closest one within the
!       specified time window.

        include 'lapsparms.inc'

        character*9 asc9_tim

        character*150 DIRECTORY
        character*31 EXT

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

        call get_directory(ext,directory,len_dir)

        do i = 1,max_files
            i_selected(i) = 0
        enddo ! i

        do i = 1,31
            if(ext(i:i) .eq. ' ')goto20
        enddo
20      lenext = i-1

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
            write(6,*)' No remaining files found within time window '
     1                ,ext(1:5),var_2d
            istatus = 0
            return
        endif

        do i = 1,i_nbr_files_ret
            i4_diff = abs(i4times(i) - i4time_needed)
            if(i4_diff .eq. i4_diff_min .and. i_selected(i) .eq. 0)then
                i_selected(i) = 1

                if(ilevel .ne. 0)then
                    if(vertical_grid .eq. 'HEIGHT')then
                        lvl_2d = zcoord_of_level(k)/10
                        lvl_coord_2d = 'MSL'
                    elseif(vertical_grid .eq. 'PRESSURE')then
                        lvl_2d = ilevel
                        lvl_coord_2d = 'MB'
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

        include 'lapsparms.inc'

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_3d(NZ_L_MAX),comment_2d
        character*10 units_3d(NZ_L_MAX),units_2d
        character*3 var_3d(NZ_L_MAX),var_2d
        integer*4 LVL_3d(NZ_L_MAX)
        character*4 LVL_COORD_3d(NZ_L_MAX)

        real*4 field_3d(imax,jmax,kmax)

        character*9 asc9_tim

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)directory(1:45),asc9_tim,ext(1:5),var_2d
11      format(' Reading 3D ',a45,1x,a9,1x,a5,1x,a3)

        do k = 1,kmax
            units_3d(k)   = units_2d
            if(vertical_grid .eq. 'HEIGHT')then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(vertical_grid .eq. 'PRESSURE')then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'MB'
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

        include 'lapsparms.inc'

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_3d(NZ_L_MAX),comment_2d
        character*10 units_3d(NZ_L_MAX),units_2d
        character*3 var_3d(NZ_L_MAX),var_2d
        integer*4 LVL_3d(NZ_L_MAX)
        character*4 LVL_COORD_3d(NZ_L_MAX)

        real*4 field_3d(imax,jmax,kmax)

        call get_directory(ext,directory,len_dir)

        call s_len(ext,len)
        write(6,11)directory(1:len_dir),ext(1:len),var_2d
11      format(' Reading 3d ',a,1x,a5,1x,a3)

        do k = 1,kmax
            units_3d(k)   = units_2d
            if(vertical_grid .eq. 'HEIGHT')then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(vertical_grid .eq. 'PRESSURE')then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'MB'
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

        include 'lapsparms.inc'

        character*150 DIRECTORY
        character*31 EXT

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

!       include 'lapsparms.inc'

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_2d
        integer*4 LVL_2d
        character*4 LVL_COORD_2d

        real*4 field_2d(imax,jmax)

        call get_directory(ext,directory,len_dir)

        write(6,11)directory,ext(1:5),var_2d
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

        include 'lapsparms.inc'

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_3d(NZ_L_MAX),comment_2d
        character*10 units_3d(NZ_L_MAX),units_2d
        character*3 var_3d(NZ_L_MAX),var_2d
        integer*4 LVL_3d(NZ_L_MAX)
        character*4 LVL_COORD_3d(NZ_L_MAX)

        real*4 field_3d(ni,nj,nk)

        call get_directory(ext,directory,len_dir)

        write(6,11)directory,ext(1:5),var_2d
11      format(' Writing 3d ',a50,1x,a5,1x,a3)

        do k = 1,nk
            units_3d(k)   = units_2d
            comment_3d(k) = comment_2d
            if(vertical_grid .eq. 'HEIGHT')then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(vertical_grid .eq. 'PRESSURE')then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'HPA'
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

        include 'lapsparms.inc'

        integer*4 MAX_FIELDS,MAX_DIM_3D
        parameter (MAX_FIELDS = 10)

        parameter (MAX_DIM_3D = MAX_FIELDS * NZ_L_MAX)


        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_3d(MAX_DIM_3D),comment_2d(nf)
        character*10 units_3d(MAX_DIM_3D),units_2d(nf)
        character*3 var_3d(MAX_DIM_3D),var_2d(nf)
        integer*4 LVL_3d(MAX_DIM_3D)
        character*4 LVL_COORD_3d(MAX_DIM_3D)

        real*4 field_3d(ni,nj,nk,nf)

        istatus = 0

        if(nf .gt. MAX_FIELDS)then
            write(6,*)' Too many fields in put_laps_multi_3d'
        endif

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
            if(vertical_grid .eq. 'HEIGHT')then
                lvl_3d(iscript_3d) = zcoord_of_level(k)/10
                lvl_coord_3d(iscript_3d) = 'MSL'
            elseif(vertical_grid .eq. 'PRESSURE')then
                lvl_3d(iscript_3d) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(iscript_3d) = 'HPA'
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

        include 'lapsparms.inc'

        integer*4 MAX_FIELDS
        parameter (MAX_FIELDS = 10)

        character*150 DIRECTORY
        character*31 EXT

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

