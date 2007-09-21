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

        subroutine trilinear_laps(ri,rj,rk,imax,jmax,kmax,array_3d
     1                           ,result)

cdoc    Interpolate 3-d array to find the field value at a fractional grid
cdoc    point.

        real array_3d(imax,jmax,kmax)

        k_low = int(rk)
        if(k_low .ge. kmax)k_low = kmax - 1 ! Extrapolate if on high end

        k_high = k_low + 1
        frac_k = rk - k_low

        call bilinear_laps(ri,rj,imax,jmax,array_3d(1,1,k_low )
     1                    ,result_low)
        call bilinear_laps(ri,rj,imax,jmax,array_3d(1,1,k_high)
     1                    ,result_high)

!       Interpolate in the vertical
        result = result_low  * (1.0 - frac_k)     +
     1           result_high *        frac_k

        return
        end

        subroutine trilinear_interp_extrap(ri,rj,rk,imax,jmax,kmax
     1                                    ,array_3d,result,istatus)

cdoc    Interpolate 3-d array to find the field value at a fractional grid
cdoc    point. This one allows you to extrapolate very slightly outside the 
cdoc    grid.

        real array_3d(imax,jmax,kmax)

        k_low = int(rk)
        if(k_low .ge. kmax)k_low = kmax - 1 ! Extrapolate if on high end

        k_high = k_low + 1
        frac_k = rk - k_low

        call bilinear_interp_extrap(ri,rj,imax,jmax,array_3d(1,1,k_low )
     1                             ,result_low,istatus)
        call bilinear_interp_extrap(ri,rj,imax,jmax,array_3d(1,1,k_high)
     1                             ,result_high,istatus)

!       Interpolate in the vertical
        result = result_low  * (1.0 - frac_k)     +
     1           result_high *        frac_k

        return
        end
