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
	subroutine getdudv_lam(lov,lap,dx,dy,lat_orig,lon_orig,
     1du,dv,u_orig,v_orig)

c	this mid-level routine performs the same task as its polar
c	counterpart except it is written for lambert.

c	Written by Dan Birkenheuer January 1995

        implicit none

        real lov,dx,dy,lap,lat_orig,lon_orig

        real du,dv,v1,v2,dlat,u_orig,v_orig
	
        real polat,polon
        real pi,dummy


        pi = acos(-1.0)
        polon = lov * pi/180.
        polat = lap * pi/180.  ! by definition

c	increment in up-down (n-s) [lat] dir
c       this is preferred since this is immune to cos(lat) adjustments for 
c	distance that one finds in the lon direction.

        dlat =  dy/2. /6.3712e3 ! km to radians
        dlat = dlat * 180./pi ! convert to degrees
        dummy = (lap+dlat)*pi/180.
        call getuv_lam (polat,polon,dummy,polon,dummy,v1)

        dummy = (lap-dlat)*pi/180.
        call getuv_lam (polat,polon,dummy,polon,dummy,v2)

        dv = v1-v2

c ----------------- make the grid square -----
        du = dv
c --------------------------------------------

c	now compute the addition grid parameters

c	



        call getuv_lam (polat,polon,lat_orig*pi/180.,lon_orig*pi/180.,
     1	u_orig,v_orig)
 

        return
        end
