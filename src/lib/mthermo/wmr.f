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

	function wmr(p,t)
c
c   this function approximates the mixing ratio wmr (grams of water
c   vapor per kilogram of dry air) given the pressure p (mb) and the
c   temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c	the formula used is given on p. 302 of the
c   smithsonian meteorological tables by roland list (6th edition).
c
c   eps = ratio of the mean molecular weight of water (18.016 g/mole)
c         to that of dry air (28.966 g/mole)
	data eps/0.62197/
c   the next two lines contain a formula by herman wobus for the 
c   correction factor wfw for the departure of the mixture of air
c   and water vapor from the ideal gas law. the formula fits values
c   in table 89, p. 340 of the smithsonian meteorological tables,
c   but only for temperatures and pressures normally encountered in
c   in the atmosphere.
	x = 0.02*(t-12.5+7500./p)
	wfw = 1.+4.5e-06*p+1.4e-03*x*x
	fwesw = wfw*esw(t)
	r = eps*fwesw/(p-fwesw)
c   convert r from a dimensionless ratio to grams/kilogram.
	wmr = 1000.*r
	return
	end
