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
        subroutine get_uv_3d(i4time,imax,jmax,kmax,u_3d,v_3d,ext,istatus
     1)

!       Returns U and V wind components from a 3D gridded LAPS dataset

!       Steve Albers            1990

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*31 EXT      ! Input extension of file (normally 3 characters)

        character*125 comment_3d(kmax)
        character*10 units_3d(kmax)
        character*3 var_u(kmax)
        character*3 var_v(kmax)
        integer LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        real u_3d(imax,jmax,kmax)
        real v_3d(imax,jmax,kmax)

        call get_directory(ext,directory,len_dir)

        do k = 1,kmax
            units_3d(k)   = 'M/S'
            comment_3d(k) = '3DWIND'
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

            var_u(k) = 'U3' ! newvar = 'U3', oldvar = 'U'
            var_v(k) = 'V3' ! newvar = 'V3', oldvar = 'V'

        enddo ! k


!       Read in 3d U array
        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,
     1          EXT,imax,jmax,kmax,kmax,
     1          VAR_u,LVL_3d,LVL_COORD_3d,UNITS_3d,COMMENT_3d,
     1          u_3d,ISTATUS)
        IF(ISTATUS .ne. 1)THEN
            write(6,*)
     1  ' Sorry, file has not yet been generated this hour'
        else
            write(6,*)
     1  ' 3D - U analysis successfully read in ',EXT(1:3)
        endif

!       Read in 3d V array
        CALL READ_LAPS_DATA(I4TIME,DIRECTORY,
     1          EXT,imax,jmax,kmax,kmax,
     1          VAR_v,LVL_3d,LVL_COORD_3d,UNITS_3d,COMMENT_3d,
     1          v_3d,ISTATUS)
        IF(ISTATUS .ne. 1)THEN
            write(6,*)
     1  ' Sorry, file has not yet been generated this hour'
!           stop
        else
            write(6,*)
     1  ' 3D - V analysis successfully read in ',EXT(1:3)
        endif


        return
        end
