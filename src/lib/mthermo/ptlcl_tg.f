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
	subroutine ptlcl(p,t,td,pc,tc)

c   this subroutine estimates the pressure pc (mb) and the temperature
c   tc (celsius) at the lifted condensation level (lcl), given the
c   initial pressure p (mb), temperature t (celsius) and dew point
c   (celsius) of the parcel.  the approximation is that lines of 
c   constant potential temperature and constant mixing ratio are
c   straight on the skew t/log p chart.
c   teten's formula for saturation vapor pressure as a function of
c   pressure was used in the derivation of the formula below.  for
c   additional details, see math notes by t. schlatter dated 8 sep 81.
c   t. schlatter, noaa/erl/profs program office, boulder, colorado,
c   wrote this subroutine.

c	baker, schlatter  17-may-1982	  original version.

c   akap = (gas constant for dry air) / (specific heat at constant
c	   pressure for dry air)
c   cta = difference between kelvin and celsius temperatures

	data akap,cta/0.28541,273.15/
	c1 = 4098.026/(td+237.3)**2
	c2 = 1./(akap*(t+cta))
	pc = p*exp(c1*c2*(t-td)/(c2-c1))
	tc = t+c1*(t-td)/(c2-c1)
	return
	end
