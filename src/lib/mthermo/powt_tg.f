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
	function powt(t,p,td)

c   this function yields wet-bulb potential temperature powt
c   (celsius), given the following input:
c	   t = temperature (celsius)
c	   p = pressure (millibars)
c	   td = dew point (celsius)

c	baker, schlatter  17-may-1982	  original version.

	data cta,akap/273.15,0.28541/

c   cta = difference between kelvin and celsius temperatures
c   akap = (gas constant for dry air) / (specific heat at
c	   constant pressure for dry air)
c   compute the potential temperature (celsius)

	pt = (t+cta)*(1000./p)**akap-cta

c   compute the lifting condensation level (lcl).

	tc = tcon(t,td)

c   for the origin of the following approximation, see the documen-
c   tation for the wobus function.

	powt = pt-wobf(pt)+wobf(tc)
	return
	end
