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
	function esat(t)
c
c   this function returns the saturation vapor pressure over
c   water (mb) given the temperature (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
c   tions of selected meteorlolgical parameters for cloud physics prob-
c   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
c   electronics command, white sands missile range, new mexico 88002.
c
	tk = t+273.16
	p1 = 11.344-0.0303998*tk
	p2 = 3.49149-1302.8844/tk
	c1 = 23.832241-5.02808*alog10(tk)
	esat = 10.**(c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/tk)
	return
	end
