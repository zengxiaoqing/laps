cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
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


        REAL   DEG_PER_RAD
        PARAMETER (DEG_PER_RAD  = 180. / 3.14159)


        REAL
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
