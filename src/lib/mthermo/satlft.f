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

	function satlft(thw,p)

c	baker, schlatter  17-may-1982	  original version.

c   input:  thw = wet-bulb potential temperature (celsius).
c		  thw defines a moist adiabat.
c	    p = pressure (millibars)
c   output: satlft = temperature (celsius) where the moist adiabat
c		  crosses p

	data cta,akap/273.15,0.28541/

c   cta = difference between kelvin and celsius temperatures
c   akap = (gas constant for dry air) / (specific heat at constant
c	    pressure for dry air)

c        the algorithm below can best be understood by referring to a
c   skew-t/log p chart.  it was devised by herman wobus, a mathemati-
c   cian formerly at the navy weather research facility but now retired.
c   the value returned by satlft can be checked by referring to table
c   78, pp.319-322, smithsonian meteorological tables, by roland list
c   (6th revised edition).
c

	if (p.ne.1000.) go to 5
	satlft = thw
	return
    5	continue

c   compute tone, the temperature where the dry adiabat with value thw
c   (celsius) crosses p.

	pwrp = (p/1000.)**akap
	tone = (thw+cta)*pwrp-cta

c   consider the moist adiabat ew1 through tone at p.  using the defini-
c   tion of the wobus function (see documentation on wobf), it can be
c   shown that eone = ew1-thw.

	eone = wobf(tone)-wobf(thw)
	rate = 1.
	go to 15

c   in the loop below, the estimate of satlft is iteratively improved.

   10	continue

c   rate is the ratio of a change in t to the corresponding change in
c   e.  its initial value was set to 1 above.

	rate = (ttwo-tone)/(etwo-eone)
	tone = ttwo
	eone = etwo
   15	continue

c   ttwo is an improved estimate of satlft.

	ttwo = tone-eone*rate

c   pt is the potential temperature (celsius) corresponding to ttwo at p

	pt = (ttwo+cta)/pwrp-cta

c   consider the moist adiabat ew2 through ttwo at p. using the defini-
c   tion of the wobus function, it can be shown that etwo = ew2-thw.

	etwo = pt+wobf(ttwo)-wobf(pt)-thw

c   dlt is the correction to be subtracted from ttwo.

	dlt = etwo*rate
	if (abs(dlt).gt.0.1) go to 10
	satlft = ttwo-dlt
	return
	end
