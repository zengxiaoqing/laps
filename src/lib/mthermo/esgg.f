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
	function esgg(t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the saturation vapor pressure over liquid
c   water esgg (millibars) given the temperature t (celsius). the
c   formula used, due to goff and gratch, appears on p. 350 of the
c   smithsonian meteorological tables, sixth revised edition, 1963,
c   by roland list.

	data cta,ews,ts/273.15,1013.246,373.15/

c   cta = difference between kelvin and celsius temperatures
c   ews = saturation vapor pressure (mb) over liquid water at 100c
c   ts = boiling point of water (k)

	data c1,      c2,      c3,      c4,       c5,       c6
     1	/ 7.90298, 5.02808, 1.3816e-7, 11.344, 8.1328e-3, 3.49149 /
	tk = t+cta

c   goff-gratch formula

	rhs = -c1*(ts/tk-1.)+c2*alog10(ts/tk)-c3*(10.**(c4*(1.-tk/ts))
     1	      -1.)+c5*(10.**(-c6*(ts/tk-1.))-1.)+alog10(ews)
	esw = 10.**rhs
	if (esw.lt.0.) esw = 0.
	esgg = esw
	return
	end
