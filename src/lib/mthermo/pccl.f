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

	function pccl(pm,p,t,td,wbar,n)
c
c   this function returns the pressure at the convective condensation
c   level given the appropriate sounding data.
c
c	baker,schlatter	17-may-1982	original version
c
c   on input:
c	p = pressure (millibars). note that p(i).gt.p(i+1).
c	t = temperature (celsius)
c	td = dew point (celsius)
c	n = number of levels in the sounding and the dimension of
c	    p, t and td
c	pm = pressure (millibars) at upper boundary of the layer for
c	     computing the mean mixing ratio. p(1) is the lower 
c	     boundary.
c   on output:
c	pccl = pressure (millibars) at the convective condensation level
c	wbar = mean mixing ratio (g/kg) in the layer bounded by
c	       pressures p(1) at the bottom and pm at the top
c   the algorithm is decribed on p.17 of stipanuk, g.s.,1973:
c   "algorithms for generating a skew-t log p diagram and computing
c   selected meteorological quantities," atmospheric sciences labora-
c   tory, u.s. army electronics command, white sands missile range, new
c   mexico 88002.
	dimension t(1),p(1),td(1)
	if (pm.ne.p(1)) go to 5
	wbar= w(td(1),p(1))
	pc= pm
	if (abs(t(1)-td(1)).lt.0.05) go to 45
	go to 25
    5	continue
	wbar= 0.
	k= 0
   10	continue
	k = k+1
	if (p(k).gt.pm) go to 10
	k= k-1
	j= k-1
	if(j.lt.1) go to 20
c   compute the average mixing ratio....alog = natural log
	do 15 i= 1,j
	   l= i+1
   15      wbar= (w(td(i),p(i))+w(td(l),p(l)))*alog(p(i)/p(l))
     *          +wbar
   20	continue
	l= k+1
c   estimate the dew point at pressure pm.
	tq= td(k)+(td(l)-td(k))*(alog(pm/p(k)))/(alog(p(l)/p(k)))
	wbar= wbar+(w(td(k),p(k))+w(tq,pm))*alog(p(k)/pm)
	wbar= wbar/(2.*alog(p(1)/pm))
c   find level at which the mixing ratio line wbar crosses the
c   environmental temperature profile.
   25	continue
	do 30 j= 1,n
	   i= n-j+1
	   if (p(i).lt.300.) go to 30
c   tmr = temperature (celsius) at pressure p (mb) along a mixing
c   	  ratio line given by wbar (g/kg)
	   x= tmr(wbar,p(i))-t(i)
	   if (x.le.0.) go to 35
   30	continue
	pccl= 0.0
	return
c  set up bisection routine
   35	l = i
	i= i+1
	del= p(l)-p(i)
	pc= p(i)+.5*del
	a= (t(i)-t(l))/alog(p(l)/p(i))
	do 40 j= 1,10
	   del= del/2.
	   x= tmr(wbar,pc)-t(l)-a*(alog(p(l)/pc))
c   the sign function replaces the sign of the first argument
c   with that of the second.
   40	pc= pc+sign(del,x)
   45	pccl = pc
	return
	end
