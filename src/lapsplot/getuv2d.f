
        subroutine get_uv_2d(i4time,k_level,uv_2d,ext,imax,jmax
     1                      ,fcst_hhmm,istatus)

!       97-Aug-17     Ken Dritz     Used commenting to (temporarily) hardwire
!                                   VERTICAL_GRID to 'PRESSURE' (without
!                                   accessing VERTICAL_GRID)
!       97-Aug-17     Ken Dritz     Removed include of lapsparms.for

        character*150 DIRECTORY
        character*31 EXT
        character*20 c_model

        character*125 comment_2d(2)
        character*10 units_2d(2)
        character*3 var(2)
        integer*4 LVL_2d(2)
        character*4 LVL_COORD_2d(2), fcst_hhmm
        character*9 a9time
        character*13 a13_time

        real*4 uv_2d(imax,jmax,2)

        call get_directory(ext,directory,len_dir)
        call s_len(ext,lext)
        if(ext(1:lext).eq.'balance')then
           ext='lw3'
           directory=directory(1:len_dir)//'lw3'//'/'
        endif

        do k = 1,2
          units_2d(k)   = 'M/S'
          comment_2d(k) = '3DWIND'

          if(k_level .gt. 0)then
!           if(vertical_grid .eq. 'HEIGHT')then
!               lvl_2d(k) = zcoord_of_level(k_level)/10
!               lvl_coord_2d(k) = 'MSL'
!           elseif(vertical_grid .eq. 'PRESSURE')then
                lvl_2d(k) = zcoord_of_level(k_level)/100
                lvl_coord_2d(k) = 'HPA'
!           endif
            var(1) = 'U3' ! newvar 'U3', oldvar = 'U'
            var(2) = 'V3' ! newvar 'V3', oldvar = 'V'

          else
            var(1) = 'U'
            var(2) = 'V'
            lvl_2d(k) = 0
            lvl_coord_2d(k) = 'AGL'

          endif

        enddo ! k

        if(ext(1:3) .eq. 'lga' .or. ext(1:3) .eq. 'fua')then
            call get_laps_cycle_time(laps_cycle_time,istatus)

            call input_background_info(
     1                              ext                     ! I
     1                             ,directory,c_model       ! O
     1                             ,i4time                  ! I
     1                             ,laps_cycle_time         ! I
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O
            if(istatus.ne.1)return

            CALL READ_LAPS(i4_initial,i4_valid,DIRECTORY,EXT
     1                    ,imax,jmax,2,2,VAR,LVL_2d,LVL_COORD_2d
     1                    ,UNITS_2d,COMMENT_2d,uv_2d,ISTATUS)
            IF(ISTATUS .ne. 1)THEN
                write(6,*)
     1          ' Sorry, file has not yet been generated this cycle'
                stop
            else
                write(6,*)
     1          ' 2d - LAPS U and V analysis successfully read in'
     1          ,lvl_2d(1)        
            endif

            i4time = i4_valid

        else

!           Read in 2d W array
            CALL READ_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,2,2,
     1          VAR,LVL_2d,LVL_COORD_2d,UNITS_2d,COMMENT_2d,
     1          uv_2d,ISTATUS)
            IF(ISTATUS .ne. 1)THEN
                write(6,*)
     1      ' Sorry, file has not yet been generated this hour'
                stop
            else
                write(6,*)
     1    ' 2d - LAPS U and V analysis successfully read in',lvl_2d(1)        
            endif

        endif

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
        integer*4 LVL_2d(nf)
        character*4 LVL_COORD_2d(nf)

        real*4 field_2d(imax,jmax,nf)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' get_laps_multi_2d: bad istatus, return'
            return
        endif

        call get_directory(ext,directory,len_dir)

        call make_fnam_lp(i4time,asc9_tim,istatus)

        write(6,11)directory(1:45),asc9_tim,ext,var_2d
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
