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
	function esice(t)
c
c   this function returns the saturation vapor pressure with respect to
c   ice esice (millibars) given the temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the formula used is based upon the integration of the clausius-
c   clapeyron equation by goff and gratch.  the formula appears on p.350
c   of the smithsonian meteorological tables, sixth revised edition,
c   1963.
c
	data cta,eis/273.16,6.1071/
c   cta = difference between kelvin and celsius temperature
c   eis = saturation vapor pressure (mb) over a water-ice mixture at 0c
	data c1,c2,c3/9.09718,3.56654,0.876793/
c   c1,c2,c3 = empirical coefficients in the goff-gratch formula
	if (t.le.0.) go to 5
	esice = 99999.
	write(6,3)esice
    3	format(' saturation vapor pressure for ice cannot be computed',
     1	       /' for temperature > 0c. esice =',f7.0)
	return
    5	continue
c   freezing point of water (k)
	tf = cta
	tk = t+cta
c   goff-gratch formula
	rhs = -c1*(tf/tk-1.)-c2*alog10(tf/tk)+c3*(1.-tk/tf)+alog10(eis)
	esi = 10.**rhs
	if (esi.lt.0.) esi = 0.
	esice = esi
	return
	end
