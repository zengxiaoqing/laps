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
        subroutine get_w_3d(i4time,imax,jmax,kmax,w_3d,ext,istatus)

!       Returns vertical wind component from a 3D gridded LAPS dataset

!       Steve Albers            1990

        include 'lapsparms.inc'

        character*50 DIRECTORY
        character*31 EXT     ! Input extension of file (normally 3 characters)

        character*125 comment_3d(NZ_L_MAX)
        character*10 units_3d(NZ_L_MAX)
        character*3 var_w(NZ_L_MAX)
        integer*4 LVL_3d(NZ_L_MAX)
        character*4 LVL_COORD_3d(NZ_L_MAX)

        logical l_convert

        common /lapsplot_omega/l_convert

        real*4 w_3d(imax,jmax,kmax)

        call get_directory(ext,directory,len_dir)

        do k = 1,kmax
            units_3d(k)   = 'PA/S'
            comment_3d(k) = '3DWIND'
            if(vertical_grid .eq. 'HEIGHT')then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(vertical_grid .eq. 'PRESSURE')then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'HPA'
            endif

            if(i4time .gt. 954547200)then
                if(ext(1:3) .eq. 'lco')then
                    var_w(k) = 'COM' ! newvar = 'COM', oldvar = 'OM'
                else
                    var_w(k) = 'OM'
                endif
            else
                write(6,*)' Setting L_convert to .false.'
                l_convert = .false.
                var_w(k) = 'W'
            endif

        enddo ! k


!       Read in 3d W array
        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,
     1          EXT,imax,jmax,kmax,kmax,
     1          var_w,LVL_3d,LVL_COORD_3d,UNITS_3d,COMMENT_3d,
     1          w_3d,ISTATUS)
        IF(ISTATUS .ne. 1)THEN
            write(6,*)' Sorry, file has not yet been generated this hour
     1'
        else
            write(6,*)' 3D - LAPS W analysis successfully read in'
        endif


        return
        end
