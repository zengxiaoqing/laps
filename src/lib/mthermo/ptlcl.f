cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
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
