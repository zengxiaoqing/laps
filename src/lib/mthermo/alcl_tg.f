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
	function alcl(t,td,p)

c	g.s. stipanuk     1973      	  original version.
c	reference stipanuk paper entitled:
c            "algorithms for generating a skew-t, log p
c	     diagram and computing selected meteorological
c	     quantities."
c	     atmospheric sciences laboratory
c	     u.s. army electronics command
c	     white sands missile range, new mexico 88002
c	     33 pages
c	baker, schlatter  17-may-1982	 

c   this function returns the pressure alcl (mb) of the lifting conden-
c   sation level (lcl) for a parcel initially at temperature t (celsius)
c   dew point td (celsius) and pressure p (millibars). alcl is computed
c   by an iterative procedure described by eqs. 8-12 in stipanuk (1973),
c   pp.13-14.
c   determine the mixing ratio line through td and p.

	aw = w(td,p)

c   determine the dry adiabat through t and p.

	ao = o(t,p)

c   iterate to locate pressure pi at the intersection of the two
c   curves. pi has been set to p for the initial guess.

 3	continue
	   pi = p
	   do 4 i= 1,10
	      x= .02*(tmr(aw,pi)-tda(ao,pi))
	      if (abs(x).lt.0.01) go to 5
 4	      pi= pi*(2.**(x))
 5	   alcl= pi
	return

	entry alclm(t,td,p)
c   for entry alclm only, t is the mean potential temperature (celsius)
c   and td is the mean mixing ratio (g/kg) of the layer containing the
c   parcel.

	aw = td
	ao = t
	go to 3
	end
