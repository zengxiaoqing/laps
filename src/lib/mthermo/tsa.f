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
	function tsa(os,p)

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

c   this function returns the temperature tsa (celsius) on a saturation
c   adiabat at pressure p (millibars). os is the equivalent potential
c   temperature of the parcel (celsius). sign(a,b) replaces the
c   algebraic sign of a with that of b.
c   b is an empirical constant approximately equal to 0.001 of the latent 
c   heat of vaporization for water divided by the specific heat at constant
c   pressure for dry air.

	data b/2.6518986/
	a= os+273.15

c   tq is the first guess for tsa.

	tq= 253.15

c   d is an initial value used in the iteration below.

	d= 120.

c   iterate to obtain sufficient accuracy....see table 1, p.8
c   of stipanuk (1973) for equation used in iteration.

	do 1 i= 1,12
	   tqk= tq-273.15
	   d= d/2.
	   x= a*exp(-b*w(tqk,p)/tq)-tq*((1000./p)**.286)
	   if (abs(x).lt.1e-7) go to 2
	   tq= tq+sign(d,x)
 1	continue
 2	tsa= tq-273.15
	return
	end
