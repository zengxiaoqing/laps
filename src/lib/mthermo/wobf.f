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

	function wobf(t)
c
c   this function calculates the difference of the wet bulb potential
c   temperatures for saturated and dry air given the temperature.
c
c	baker,schlatter	17-may-1982	original version
c
c        let wbpts = wet-bulb potential temperature for saturated
c   air at temperature t (celsius). let wbptd = wet-bulb potential
c   temperature for completely dry air at the same temperature t.
c   the wobus function wobf (in degrees celsius) is defined by
c                      wobf(t) = wbpts-wbptd.
c   although wbpts and wbptd are functions of both pressure and
c   temperature, their difference is a function of temperature only.
c
c        to understand why, consider a parcel of dry air at tempera-
c   ture t and pressure p. the thermodynamic state of the parcel is
c   represented by a point on a pseudoadiabatic chart. the wet-bulb
c   potential temperature curve (moist adiabat) passing through this
c   point is wbpts. now t is the equivalent temperature for another
c   parcel saturated at some lower temperature tw, but at the same
c   pressure p.  to find tw, ascend along the dry adiabat through
c   (t,p). at a great height, the dry adiabat and some moist
c   adiabat will nearly coincide. descend along this moist adiabat
c   back to p. the parcel temperature is now tw. the wet-bulb
c   potential temperature curve (moist adiabat) through (tw,p) is wbptd.
c   the difference (wbpts-wbptd) is proportional to the heat imparted
c   to a parcel saturated at temperature tw if all its water vapor
c   were condensed. since the amount of water vapor a parcel can
c   hold depends upon temperature alone, (wbptd-wbpts) must depend
c   on temperature alone.
c
c        the wobus function is useful for evaluating several thermo-
c   dynamic quantities.  by definition:
c		    wobf(t) = wbpts-wbptd.               (1)
c   if t is at 1000 mb, then t is a potential temperature pt and
c   wbpts = pt. thus
c		    wobf(pt) = pt-wbptd.                 (2)
c   if t is at the condensation level, then t is the condensation
c   temperature tc and wbpts is the wet-bulb potential temperature
c   wbpt. thus
c		    wobf(tc) = wbpt-wbptd.               (3)
c   if wbptd is eliminated from (2) and (3), there results
c		    wbpt = pt-wobf(pt)+wobf(tc).
c   if wbptd is eliminated from (1) and (2), there results
c		    wbpts = pt-wobf(pt)+wobf(t).
c
c        if t is an equivalent potential temperature ept (implying
c   that the air at 1000 mb is completely dry), then wbpts = ept
c   and wbptd = wbpt. thus
c		    wobf(ept) = ept-wbpt.
c   this form is the basis for a polynomial approximation to wobf.
c   in table 78 on pp.319-322 of the smithsonian meteorological
c   tables by roland list (6th revised edition), one finds wet-bulb
c   potential temperatures and the corresponding equivalent potential
c   temperatures listed together. herman wobus, a mathematician for-
c   merly at the navy weather research facility, norfolk, virginia,
c   and now retired, computed the coefficients for the polynomial
c   approximation from numbers in this table.
c
c                                    notes by t.w. schlatter
c                                    noaa/erl/profs program office
c                                    august 1981
c
	x = t-20.
	if (x.gt.0.) go to 10
	pol = 1.                 +x*(-8.8416605e-03
     1	     +x*( 1.4714143e-04  +x*(-9.6719890e-07
     2       +x*(-3.2607217e-08  +x*(-3.8598073e-10)))))
	wobf = 15.130/pol**4
	return
   10	continue
	pol = 1.                 +x*( 3.6182989e-03
     1       +x*(-1.3603273e-05  +x*( 4.9618922e-07
     2       +x*(-6.1059365e-09  +x*( 3.9401551e-11
     3       +x*(-1.2588129e-13  +x*( 1.6688280e-16)))))))
	wobf = 29.930/pol**4+0.96*x-14.8
	return
	end
