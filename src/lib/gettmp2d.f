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
        subroutine get_temp_2d(i4time_needed,i4time_tol,i4time_nearest
     1                          ,ilevel_mb,imax,jmax,temp_2d,istatus)

!       Steve Albers            1990

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment_2d
        character*10 units_2d
        character*3 var_t
        integer LVL_2d
        character*4 LVL_COORD_2d

        real temp_2d(imax,jmax)

        character*255 c_filespec

        write(6,*)
        write(6,*)' Getting Temperature Analysis at ',ilevel_mb

        ext = 'lt1'
        call get_directory(ext,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

        call get_file_time(c_filespec
     1  ,i4time_needed,i4time_nearest)

        if(abs(i4time_needed - i4time_nearest) .gt. i4time_tol)then
            write(6,*)' No file available within requested time window'
            istatus = 0
            return
        endif

        units_2d  = 'K'
        if(ltest_vertical_grid('HEIGHT'))then
            lvl_2d = zcoord_of_level(k)/10
            lvl_coord_2d = 'MSL'
        elseif(ltest_vertical_grid('PRESSURE'))then
            lvl_2d = ilevel_mb
            lvl_coord_2d = 'MB'
        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return
        endif

        var_t = 'T3'  ! newvar = 'T3', oldvar = 'T'

!       Read in 2d T array
        CALL READ_LAPS_DATA(I4TIME_NEAREST,DIRECTORY,
     1          EXT,imax,jmax,1,1,
     1          var_t,LVL_2d,LVL_COORD_2d,UNITS_2d,COMMENT_2d,
     1          temp_2d,ISTATUS)

!       Test for bad temperatures
        if(temp_2d(29,29) .lt. 173.)istatus = 0

        return

        end

