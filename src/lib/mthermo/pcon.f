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
	function pcon(p,t,tc)
c
c   this function returns the pressure pcon (mb) at the lifted condensa-
c   tion level, given the initial pressure p (mb) and temperature t
c   (celsius) of the parcel and the temperature tc (celsius) at the lcl.
c
c	baker,schlatter	17-may-1982	original version
c
c   the algorithm is exact.  it makes use of the formula for the
c   potential temperatures corresponding to t at p and tc at pcon.
c   these two potential temperatures are equal.
c
	data akapi/3.5037/
c   akapi = (specific heat at constant pressure for dry air) /
c	    (gas constant for dry air)
c
c   convert t and tc to kelvin temperatures.
	tk = t+273.16
	tck = tc+273.16
	pcon = p*(tck/tk)**akapi
	return
	end
