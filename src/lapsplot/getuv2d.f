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

        subroutine get_uv_2d(i4time,k_level,uv_2d,ext,imax,jmax,istatus)

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
        character*4 LVL_COORD_2d(2)
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
            if(ext(1:3) .eq. 'fua')then
                write(6,205)ext(1:3)
 205            format(/'  Enter model [e.g. mm5] for ',a3,' file: ',$)       

                read(5,206)c_model
 206            format(a)

                call get_directory(ext,directory,len_dir)
                call s_len(c_model,len_model)
                directory = directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

            endif

            write(6,211)ext(1:3)
 211        format(/'  Enter yydddhhmmHHMM for ',a3,' file: ',$)

            read(5,1)a13_time
 1          format(a13)
            call get_fcst_times(a13_time,I4TIME,i4_valid,i4_fn)

            CALL READ_LAPS(I4TIME,i4_valid,DIRECTORY,EXT,imax,jmax,2,2,       
     1          VAR,LVL_2d,LVL_COORD_2d,UNITS_2d,COMMENT_2d,
     1          uv_2d,ISTATUS)
            IF(ISTATUS .ne. 1)THEN
                write(6,*)
     1          ' Sorry, file has not yet been generated this hour'
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
