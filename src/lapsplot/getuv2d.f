
        subroutine get_uv_2d(i4time,k_level,uv_2d,ext,imax,jmax
     1                      ,fcst_hhmm,c_model,istatus)

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
        integer LVL_2d(2)
        character*4 LVL_COORD_2d(2)
        character*5 fcst_hhmm
        character*9 a9time

        real uv_2d(imax,jmax,2)

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
     1          ' Sorry, unable to read in wind variables ',var          
                return
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

