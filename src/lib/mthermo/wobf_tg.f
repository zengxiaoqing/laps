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
	function wobf(t)

c   this function calculates the difference of the wet-bulb potential
c   temperatures for saturated and dry air given the temperature.

c	baker, schlatter  17-may-1982	  original version.

c        let wbpts = wet-bulb potential temperature for saturated
c   air at temperature t (celsius). let wbptd = wet-bulb potential
c   temperature for completely dry air at the same temperature t.
c   the wobus function wobf (in degrees celsius) is defined by
c                      wobf(t) = wbpts-wbptd.
c   although wbpts and wbptd are functions of both pressure and
c   temperature, their difference is a function of temperature only.

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
