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
	subroutine uv_ij (ny,u1,v1,du,dv,u,v,ri,rj)

c	this routine takes inital u1 v1 (lower left) map projection coords,
c	the increment in u and v (du and dv) (parent grid)
c	then given ANY set of u and v it generates real
c	representations of i and j in the parent grid corresponding to the 
c	u and v sets.

c	Written by Dan Birkenheuer February 1994

	implicit none

	real u1,v1,du,dv,u,v,ri,rj
	integer ny

	ri = (u-u1)/du  +1.0
	rj = float(ny) - (v-v1)/dv

	return
	end
