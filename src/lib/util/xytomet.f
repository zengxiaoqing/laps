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
        SUBROUTINE   XY_TO_MET_XM( X,
     1                     Y,
     1                     RANGE,
     1                     DIR,
     1                     ISTATUS )

!  XY_TO_MET_XM  converts cartesian coordinates to meteorological (polar)
!                coordinates.

!       Bob Lipschutz   12-APR-1983     Original version
!       Windsor         15-Aug-1985     Fixed so differentiates btwn 0 and 180


        REAL*4   DEG_PER_RAD
        PARAMETER (DEG_PER_RAD  = 180. / 3.14159)


        REAL*4
     1  X,              ! X-coordinate.
     1  Y,              ! Y-coordinate (in same units as X).
     1  RANGE,          ! Output:  Range in same units as X and Y).
     1  DIR             ! Output:  Meteorological degrees.


        RANGE  = SQRT( X*X + Y*Y )

        IF ( X .ne. 0)  THEN

            DIR  = ATAN2( Y,X )
!       print *, 'after atan2 : ',dir

            DIR  = DIR * DEG_PER_RAD
            IF( DIR .LT. 0 )  DIR  = 360. + DIR

!         ... Convert to met. degrees.

            DIR  = 450. - DIR
            IF( DIR .GT. 360. )  DIR  = DIR - 360.

          ELSE

            If (y .lt. 0) then

                dir = 180.0

             else

               DIR  = 0.0

            Endif

        END IF

!       dir = float(nint(dir))

!       print *,' x,y : ',x,y,' deg : ',dir

        ISTATUS   = 1
        RETURN
        END
