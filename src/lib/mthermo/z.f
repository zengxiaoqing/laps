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
	function z(pt,p,t,td,n)

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

c   this function returns the thickness of a layer bounded by pressure
c   p(1) at the bottom and pressure pt at the top.
c   on input:
c	p = pressure (mb).  note that p(i).gt.p(i+1).
c	t = temperature (celsius)
c	td = dew point (celsius)
c	n = number of levels in the sounding and the dimension of
c	    p, t and td
c   on output:
c	z = geometric thickness of the layer (m)
c   the algorithm involves numerical integration of the hydrostatic
c   equation from p(1) to pt. it is described on p.15 of stipanuk
c   (1973).

	dimension t(1),p(1),td(1),tk(100)

c	c1 = .001*(1./eps-1.) where eps = .62197 is the ratio of the
c			      molecular weight of water to that of
c			      dry air. the factor 1000. converts the
c			      mixing ratio w from g/kg to a dimension-
c			      less ratio.
c	c2 = r/(2.*g) where r is the gas constant for dry air
c		      (287 kg/joule/deg k) and g is the acceleration
c		      due to the earth's gravity (9.8 m/s**2). the
c		      factor of 2 is used in averaging two virtual
c		      temperatures.

	data c1/.0006078/,c2/14.64285/
	do 5 i= 1,n
	   tk(i)= t(i)+273.15
    5	continue
	z= 0.0
	if (pt.lt.p(n)) go to 20
	i= 0
   10	i= i+1
	j= i+1
	if (pt.ge.p(j)) go to 15
	a1= tk(j)*(1.+c1*w(td(j),p(j)))
	a2= tk(i)*(1.+c1*w(td(i),p(i)))
	z= z+c2*(a1+a2)*(alog(p(i)/p(j)))
	go to 10
   15	continue
	a1= tk(j)*(1.+c1*w(td(j),p(j)))
	a2= tk(i)*(1.+c1*w(td(i),p(i)))
	z= z+c2*(a1+a2)*(alog(p(i)/pt))
	return
 20	z= -1.0
	return
	end
