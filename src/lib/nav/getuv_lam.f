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
	subroutine getuv_lam (polat,polon,lat,lon,u,v)

c 	derived for Lambert conformal conic with one standard parallel
c 	routine assumes all angles in radians... 
c 	compare polat = theta_zero (ncar manual) 
c 	compare polon = phi_zero (ncar)

c	Written by Dan Birkenheuer February 1994
c       J Smart 6-96  modified the arguments for trig functions by computing
c 			them first and storing in variable term

	real*4 polat,polon,lat,lon,u,v
	real*8 n,pi
        real*8 term
        real*8 rtanterm,rexpterm

	pi = acos(-1.)
c
        term = pi/2.-polat
	n = cos(term)
        term = pi/4.-lat/2.
	rtanterm = tan(term)
        rexpterm = rtanterm**n
        term = n*(lon-polon)
	u = rexpterm* sin (term)
	v = (-1.)*rexpterm*cos(term)

	return
	end
