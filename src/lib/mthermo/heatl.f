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
	function heatl(key,t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the latent heat of
c		evaporation/condensation         for key=1
c		melting/freezing		 for key=2
c		sublimation/deposition		 for key=3
c   for water. the latent heat heatl (joules per kilogram) is a 
c   function of temperature t (celsius). the formulas are polynomial
c   approximations to the values in table 92, p. 343 of the smithsonian
c   meteorological tables, sixth revised edition, 1963 by roland list.
c   the approximations were developed by eric smith at colorado state
c   university.
c   polynomial coefficients

	data a0,a1,a2/ 3337118.5,-3642.8583, 2.1263947/
	data b0,b1,b2/-1161004.0, 9002.2648,-12.931292/
	data c0,c1,c2/ 2632536.8, 1726.9659,-3.6248111/
	hltnt = 0.
	tk = t+273.15
	if (key.eq.1) hltnt=a0+a1*tk+a2*tk*tk
	if (key.eq.2) hltnt=b0+b1*tk+b2*tk*tk
	if (key.eq.3) hltnt=c0+c1*tk+c2*tk*tk
	heatl = hltnt
	return
	end
