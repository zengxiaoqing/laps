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

        subroutine put_wind_2d(i4time,DIRECTORY,EXT,var,units,
     1                  comment,wind_2d,imax,jmax,istatus)

        character*150 DIRECTORY
        character*31 EXT

        character*125 comment,comment_2d(2)
        character*10 units,units_2d(2)
        character*3 var,var_2d(2)
        integer*4 LVL,LVL_2d(2)
        character*4 LVL_COORD,LVL_COORD_2d(2)

        real*4 wind_2d(imax,jmax,2)

        write(6,11)directory,ext(1:5)
11      format(' Writing 2d Wind ',a50,1x,a5,1x,a3)

        lvl = 0
        lvl_coord = 'MSL'

        var_2d(1) = 'U'
        var_2d(2) = 'V'

        do k = 1,2
            comment_2d(k) = comment
            LVL_2d(k) = LVL
            LVL_Coord_2d(k) = LVL_Coord
            units_2d(k) = 'M/S'
        enddo

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1  2,2,VAR_2D,LVL_2D,LVL_COORD_2D,UNITS_2D,
     1                     COMMENT_2D,wind_2d,ISTATUS)

        return
        end
