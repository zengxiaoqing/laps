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

        subroutine trilinear_laps(ri,rj,rk,imax,jmax,kmax,array_3d
     1                           ,result)

cdoc    Interpolate 3-d array to find the field value at a fractional grid
cdoc    point.

        real*4 array_3d(imax,jmax,kmax)

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

        real*4 array_3d(imax,jmax,kmax)

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
