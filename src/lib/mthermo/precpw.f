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
	function precpw(td,p,n)
c
c   this function computes total precipitable water precpw (cm) in a
c   vertical column of air based upon sounding data at n levels:
c	   td = dew point (celsius)
c	   p = pressure (millibars)
c
c	baker,schlatter	17-may-1982	original version
c
c   calculations are done in cgs units.
	dimension td(n),p(n)
c   g = acceleration due to the earth's gravity (cm/s**2)
	data g/980.616/
c   initialize value of precipitable water
	pw = 0.
	nl = n-1
c   calculate the mixing ratio at the lowest level.
	wbot = wmr(p(1),td(1))
	do 5 i=1,nl
	wtop = wmr(p(i+1),td(i+1))
c   calculate the layer-mean mixing ratio (g/kg).
	w = 0.5*(wtop+wbot)
c   make the mixing ratio dimensionless.
	wl = .001*w
c   calculate the specific humidity.
	ql = wl/(wl+1.)
c   the factor of 1000. below converts from millibars to dynes/cm**2.
	dp = 1000.*(p(i)-p(i+1))
	pw = pw+(ql/g)*dp
	wbot = wtop
    5	continue
	precpw = pw
	return
	end
