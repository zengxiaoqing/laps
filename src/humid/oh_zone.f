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
      
      subroutine oh_zone (i_p, ozone, nl, repeat, istatus)

c supplies ozone profile when there is no data to match a given input array
c this routine will interplate once and thenceforth return the once determined
c array.  Based in part on gimtau.f, this code is better since it actually
c can extrapolate if pressure ranges go outside of the standard atmospheric
c data range.  For LAPS this should provide a good OZONE profile even in
c situations where the profile changes from grid point to grid point.

c Note that the routine is devised to take this situation in hand.  If
c repeat is 0 it will only work for a onetime call. If repeat is switched to 
c 1, it will recompute the profile each call.

c one reason for repeat might be where the calling routine modifies the
c ozone profile every call.

      implicit none
      save

c external variables
      integer nl
      integer istatus
      integer repeat
      real i_p (nl), ozone (nl)

c internal variables
      real prer(40)
      real pref(40)
      real oref(40)
      real tref(40)
      real wref(40)
      integer i

c internal logical varaibles
      integer first_time
      data first_time/0/


c $ Transmittance-Model Reference = U.S. Standard Atmosphere, 1976


      data prer/ .1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,
     1 50.,60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,
     1 430.,475.,500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./
      data tref/
     1  231.70, 245.22, 263.35, 270.63, 264.07, 257.93, 249.51, 243.65,
     1  239.24, 232.64, 228.07, 225.00, 223.13, 221.72, 220.54, 217.28,
     1  216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.72,
     1  220.85, 228.58, 235.38, 241.45, 244.81, 249.48, 251.95, 258.32,
     1  262.48, 266.40, 268.61, 274.21, 278.74, 282.97, 284.71, 287.50/
      data wref/
     1   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,
     1   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,
     1   0.003,  0.003,  0.003,  0.003,  0.003,  0.004,  0.005,  0.014,
     1   0.036,  0.089,  0.212,  0.331,  0.427,  0.588,  0.699,  1.059,
     1   1.368,  1.752,  1.969,  2.741,  3.366,  3.976,  4.255,  4.701/
      data oref/
     1 0.65318,1.04797,2.13548,3.82386,5.26768,6.11313,7.35964,7.75004,
     1 7.82119,7.56126,6.92006,6.10266,5.55513,5.15298,4.59906,2.86792,
     1 2.29259,1.80627,1.28988,0.93973,0.72277,0.54848,0.46009,0.29116,
     1 0.16277,0.09861,0.06369,0.05193,0.04718,0.04097,0.03966,0.03614,
     1 0.03384,0.03342,0.03319,0.03249,0.03070,0.0287,0.02805,0.02689/

      if( repeat .eq. 1) first_time=0


      if (first_time.eq.0) then

c this type of interpolation only works when inteval spacing is non
c exponential  therefore must convert pressure coordinate to log space.

       do i = 1,40
        pref(i) = log(prer(i))
       enddo
 
       do i = 1,nl  !fill each level via polonomial interpolation log p.
         call interpolate_ozone(pref,oref,40,
     1                         log(i_p(i)),ozone(i),istatus)
          if (istatus .ne. 1) then
           write(6,*) 'using ozone next best ozone guess'
           if (i .ne. 1) ozone (i) = ozone (i-1)
          endif
         ozone(i) = max(0.0,ozone(i))
       enddo !i loops through all laps levels
       first_time = 1
      else
       return
      endif

      return
      
      end


      subroutine interpolate_ozone(xa,ya,n,x,y,istatus)

c based on program found in numerical recipes.
c this routine will not only interpolote but extrapolate using rational
c funcitons.

      integer istatus
      integer n,nmax
      real dy,x,y,xa(n),ya(n),tiny
      parameter (nmax=100,tiny=1.e-25)
      integer i,m,ns
      real dd,h,hh,t,w,c(nmax),d(nmax)
      istatus = 1
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
c the tiny code here helps prevent the zero-over-zero situation that may
c occur if not checked.
        c(i)=ya(i)
        d(i)=ya(i)+tiny
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
c note that variable h here will never be zero (and hence cause no problems
c being in the denominator in computing t.
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.) then
             write (6,*) 'failure in ozone interpolation ... pole'
             istatus = 0
             return
          endif
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
