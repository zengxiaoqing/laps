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
	function tv(t,td,p)
c
c   this function returns the virtual temperature tv (celsius) of
c   a parcel of air at temperature t (celsius), dew point td
c   (celsius), and pressure p (millibars). the equation appears
c   in most standard meteorological texts.
c
c	baker,schlatter	17-may-1982	original version
c
	data cta,eps/273.16,0.62197/
c   cta = difference between kelvin and celsius temperatures.
c   eps = ratio of the mean molecular weight of water (18.016 g/mole)
c	  to that of dry air (28.966 g/mole)
	tk = t+cta
c   calculate the dimensionless mixing ratio.
	w = .001*wmr(p,td)
	tv = tk*(1.+w/eps)/(1.+w)-cta
	return
	end
