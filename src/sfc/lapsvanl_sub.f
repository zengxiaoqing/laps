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

c
c=====  Here are John's subroutines...(abandon hope, ye who enter)
c
	subroutine vortdiv(u,v,vort,div,imax,jmax,dx,dy)
c this routine computes vorticity and divergence from u and v winds
c using centered differences.
	real*4 u(imax,jmax),v(imax,jmax),vort(imax,jmax)
	real*4 div(imax,jmax),dx(imax,jmax),dy(imax,jmax)
c
	do j=2,jmax-1
	do i=2,imax-1
	    div(i,j) = (u(i,j-1) - u(i-1,j-1)) / dx(i,j)
     &               + (v(i-1,j) - v(i-1,j-1)) / dy(i,j)
	    vort(i,j) = (v(i,j) - v(i-1,j)) / dx(i,j)
     &                - (u(i,j) - u(i,j-1)) / dy(i,j)
	enddo !i
	enddo !j
	call bounds(div,imax,jmax)
	call bounds(vort,imax,jmax)
c
	return
	end
c
c
	subroutine channel(u,v,topo,imax,jmax,top,pblht,dx,dy,
     &                     z,div)
c
c=====================================================================
c
c	Routine to channel winds around terrain features.
c       Includes option to conserve surface convergence from 
c       raw wind data since channeling routine acts to eliminate
c       convergence totally.
c
c       Original:   J. McGinley, 2nd half of 20th Century
c       Changes:    P. Stamus  26 Aug 1997  Changes for dynamic LAPS
c
c=====================================================================
c
	real*4 u(imax,jmax),v(imax,jmax),z(imax,jmax)
	real*4 dx(imax,jmax),dy(imax,jmax),top(imax,jmax),topo(imax,jmax)
	real*4 div(imax,jmax)
c
	real*4 phi(imax,jmax),ter(imax,jmax),du(imax,jmax),dv(imax,jmax)
	real*4 dpbl(imax,jmax),b(imax,jmax),c(imax,jmax)
c
	call zero(phi,imax,jmax)	! zero the work arrays
	call zero(ter,imax,jmax)
	call zero(du,imax,jmax)
	call zero(dv,imax,jmax)
	call zero(dpbl,imax,jmax)
	call zero(b,imax,jmax)
	call zero(c,imax,jmax)
c
	do j=1,jmax
	do i=1,imax
	  dpbl(i,j) = top(i,j) - topo(i,j)
	enddo !i
	enddo !j
c
	do j=2,jmax-1
	do i=2,imax-1
	  dzx2 = (dpbl(i+1,j)+dpbl(i+1,j+1))*.5
	  dzx1 = (dpbl(i-1,j)+dpbl(i-1,j-1))*.5
	  dzy2 = (dpbl(i-1,j)+dpbl(i,j))*.5
	  dzy1 = (dpbl(i-1,j-1)+dpbl(i,j-1))*.5
	  u2 = u(i,j-1)
	  u1 = u(i-1,j-1)
	  v2 = v(i-1,j)
	  v1 = v(i-1,j-1)
          zbar=(dzx2+dzy2+dzx1+dzy1)*.25
          zbars=zbar*zbar
          b(i,j)=2.*(dzx2-dzx1)/dx(i,j)/zbar
          c(i,j)=2.*(dzy2-dzy1)/dy(i,j)/zbar
	  ter(i,j) = -((u2*dzx2-u1*dzx1)/dx(i,j)
     1                +(v2*dzy2-v1*dzy1)/dy(i,j)) /zbars
     2                -div(i,j)/zbar
	enddo !i
	enddo !j
c
c	print *,' Calculating the solution for streamfunction'
	call zero(phi,imax,jmax)
	call leib(phi,ter,100,.1,imax,jmax,z,b,c,z,z,dx,dy,1.)
c
c.....	Adjust the winds.
c
	do j=1,jmax-1
	do i=1,imax-1
	  dzx2 = (dpbl(i,j)+dpbl(i,j+1))*.5
	  dzy2 = (dpbl(i,j)+dpbl(i+1,j+1))*.5
	  du(i,j) = (phi(i+1,j+1) - phi(i,j+1)) / dx(i,j) * dzx2
	  u(i,j) = du(i,j) + u(i,j)
	  dv(i,j) = (phi(i+1,j+1) - phi(i+1,j)) / dy(i,j) * dzy2
	  v(i,j) = dv(i,j) + v(i,j)
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine frict(fu,fv,u,v,uo,vo,imax,jmax,ak,akk)
c
	real*4 fu(imax,jmax),fv(imax,jmax),u(imax,jmax),v(imax,jmax)
	real*4 vo(imax,jmax),uo(imax,jmax),akk(imax,jmax)
c
	do j=1,jmax
	do i=1,imax
	    uu = u(i,j) * .75 + uo(i,j) * .25
	    vv = v(i,j) * .75 + vo(i,j) * .25
	    fu(i,j) = ak * uu * abs(uu) * akk(i,j)
	    fv(i,j) = ak * vv * abs(vv) * akk(i,j)
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine nonlin(nu,nv,u,v,uo,vo,imax,jmax,dx,dy)
c
	real*4 nu(imax,jmax),nv(imax,jmax),u(imax,jmax),v(imax,jmax)
	real*4 uo(imax,jmax),vo(imax,jmax),dx(imax,jmax),dy(imax,jmax)
c
	do j=2,jmax-1
	do i=2,imax-1
	    dudx = (u(i+1,j)-u(i-1,j)) / dx(i,j) * .375 +
     &           (uo(i+1,j)-uo(i-1,j)) / dx(i,j) * .125
	    dvdy = (v(i,j+1)-v(i,j-1)) / dy(i,j) * .375 +
     &           (vo(i,j+1)-vo(i,j-1)) / dy(i,j) * .125
	    dudy = (u(i,j+1)-u(i,j-1)) / dy(i,j) * .375 +
     &           (uo(i,j+1)-uo(i,j-1)) / dy(i,j) * .125
	    dvdx = (v(i+1,j)-v(i-1,j)) / dx(i,j) * .375 +
     &           (vo(i+1,j)-vo(i-1,j)) / dx(i,j) * .125
	    uu = (u(i,j)+u(i,j-1)+u(i+1,j)+u(i+1,j-1)) * .1875 +
     &         (uo(i,j)+uo(i,j-1)+uo(i+1,j)+uo(i+1,j-1)) * .0675
	    vv = (v(i,j)+v(i-1,j)+v(i,j+1)+v(i-1,j+1)) * .1875 +
     &         (vo(i,j)+vo(i-1,j)+vo(i,j+1)+vo(i-1,j+1)) * .0675
	    utt = u(i,j) * .75 + uo(i,j) * .25
	    vtt = v(i,j) * .75 + vo(i,j) * .25
	    nu(i,j) = utt * dudx + vv * dudy
	    nv(i,j) = uu * dvdx + vtt * dvdy
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine leib(sol,force,itmax,erf,imax,jmax,a,b,c,d,e,
     &                   dx,dy,dz)
c
	real*4 sol(imax,jmax),force(imax,jmax),a(imax,jmax)
	real*4 b(imax,jmax),c(imax,jmax),d(imax,jmax),e(imax,jmax)
	real*4 dx(imax,jmax), dy(imax,jmax)
c
	ovr = 1.
	reslmm = 0.
	erb = 0.
c  first guess here
	ittr = 0
	do 1 it=1,itmax
	  ertm = 0.
	  ermm = 0.
	  ia = 0
	  corlm = 0.
	  do 2 j=2,jmax-1
	  do 2 i=2,imax-1
	    dx2 = dx(i,j) * 2.
	    dxs = dx(i,j) * dx(i,j)
	    dy2 = dy(i,j) * 2.
	    dys = dy(i,j) * dy(i,j)
	    aa = a(i,j)
	    bb = b(i,j)
	    cc = c(i,j)
	    dd = d(i,j)
	    cortm = (-2. / dxs) - (2. / dys) + e(i,j)
20	    res = (sol(i+1,j) + sol(i-1,j)) / dxs +
     &            (sol(i,j+1) + sol(i,j-1)) / dys +
     &            (cortm * sol(i,j)) - force(i,j) + bb * 
     &            ((sol(i+1,j) - sol(i-1,j)) / dx2) + cc * 
     &            (sol(i,j+1) - sol(i,j-1)) / dy2
	    cor = res / cortm
	    if(abs(cor) .gt. erf) ia = 1
	    if(abs(cor) .gt. corlm) corlm = abs(cor)
	    sol(i,j) = sol(i,j) - cor * ovr
2	  continue
5	  ittr = ittr + 1
	  cor5 = corlm
	  if(ittr .ne. 5) go to 15
	  ittr = 0
	  rho = (cor5 / cor0) ** .2
	  if(rho .gt. 1) go to 16
	  ovr = 2. / (1. + sqrt(1. - rho))
16	  continue
	  cor0 = cor5
15	  if(ia .ne. 1) go to 4
	  if(it .ne. 1) go to 1
	  corlmm = corlm
	  cor0 = corlmm
1	continue
4	continue
	reslm = corlm * cortm
	write(6,1001) it,reslm,corlm,corlmm,erb
	write(6,1002) ovr
1002	format(1x,'OVR RLXTN CONST AT FNL ITTR = ',e10.4)
1001	format(1x,'ITERATIONS= ',i4,' MAX RESIDUAL= ',e10.3,
     & ' MAX CORRECTION= ',e10.3, ' FIRST ITER MAX COR= ',e10.3,
     & 'MAX BNDRY ERROR= ',e10.3)
c
	return
	end
c
c
	subroutine spline(t,to,tb,alf,alf2,beta,a,s_in,cormax,err,imax,
     &                    jmax,roi,bad_mult,imiss,mxstn,obs_error,name)
c
c*******************************************************************************
c	LAPS spline routine...based on one by J. McGinley.
c
c	Changes:  
c	  P. Stamus	10-18-90  Started to clean code. Made alf2/alf2o arrays.
c			11-11-91  Pass in dummy work arrays.
c			07-27-93  Changes for new barnes2 routine. 
c                       07-20-95  Put wt calcs here...call to dynamic_wts.
c                       08-26-97  Changes for dynamic LAPS. Pass in obs_error.
c         J. McGinley   09-22-98  Changes to fully exploit background info.
c          and P.Stamus           Routine modified to always use a background
c                                 field for QC and 1st guess.  Spline solves for
c                                 a solution difference from 1st guess.  Removed
c                                 2 Barnes calls. When no satellite data for HSM,
c                                 'a' weight set to zero.  Bkg added back into
c                                 t, to, and s arrays on exit.
c         P. Stamus     09-29-98  Calc std dev from just obs (not boundaries+obs)
c
c*******************************************************************************
c
	real*4 t(imax,jmax), to(imax,jmax), s(imax,jmax), s_in(imax,jmax)
	real*4 RESS(100), tb(imax,jmax) , alf2(imax,jmax)
c
	real*4 fnorm(0:imax-1,0:jmax-1)
c	real*4 alf2o(imax,jmax)  !work array
c
	character name*10
	logical iteration
c
!	write(9,910)
!910	format(' in spline routine')
	imiss = 0
	ovr = 1.4
	iflag = 0
	itmax = 100	! max number of iterations
	zeros = 1.e-30
	smsng = 1.e37
	cormax = 1.
	beta = 0.
	call move(s_in, s, imax, jmax)
	call zero(t, imax,jmax)
c
c.....	first guess use barnes
c
	npass = 1
c
c.....  Count the number of observations in the field (not counting the 
c.....  boundaries.
c
	n_obs_var = 0
	do j=3,jmax-2
	do i=3,imax-2
	  if(to(i,j) .ne. 0.) n_obs_var = n_obs_var + 1
	enddo !i
	enddo !j
	if(n_obs_var .eq. 0) then
	  print *,' +++ In SPLINE: All zeros in data array. +++'
	  imiss = 1
	  return
	else
	   print *,' Observations in data array: ', n_obs_var
	endif
c
	if(name.ne.'NOPLOT' .and. name(1:3).ne.'TB8') then
	  write(9,912)
 912	  format(' data passed into spline:')
	endif
c
c.....	Data check algorithm and computation of difference 
c.....  from background.
c
	sum = 0.
	cnt = 0.
	sum1 = 0.
	icnt = 0
c
c.....	Compute standard deviation of the obs
c
	do j=3,jmax-2
	do i=3,imax-2
	  if(to(i,j) .eq. 0.) go to 99
	  sum = sum + ((to(i,j) - tb(i,j)) ** 2)
	  cnt = cnt + 1.
99	continue
        enddo !i
	enddo !j
c
	if(cnt .eq. 0.) then
	  print *, ' +++ In SPLINE: all zeros in data array. +++'
	  return
	else
	  std = sqrt(sum / cnt)
	endif
	if(std .eq. 0.) then
	  write(6,927) t(21,21)
 927	  format(1x,'++ Uniform input to spline routine with value: ',
     &           e12.4,' ++')
	  std = zeros
	endif
c
c.....  Bad data defined as deviating 'bad_mult' sigmas from 1st guess
c
	bad = bad_mult * std
	iflag = 0
	print *,' std dev: ',std,', bad value: ',bad
c
c.....  eliminate bad data   
c
	do j=3,jmax-2
	do i=3,imax-2
	  if(to(i,j) .eq. 0.) go to 98
	  diff = to(i,j) - tb(i,j)
	  if(abs(diff) .lt. bad) then
	     to(i,j) = diff
	  else
	     iflag = 1
	     write(6,1099) i,j,to(i,j), diff
	     to(i,j) = 0.
	     if(i.ne.1 .and. i.ne.imax .and. j.ne.1 .and. j.ne.jmax)
     &         n_obs_var = n_obs_var - 1
	  endif
 98	  continue 
        enddo !i
	enddo !j
 1099	format(1x,'bad data at i,j ',2i5,': value ',e12.4,', diff ',e12.4)
c
	if(n_obs_var .gt. 0) then
	   print *,
     &       ' Observations in data array after spline QC: ',n_obs_var
	else
	   print *,' No obs in data array.  Returning.'
	   return
	endif
c
c.....  Ensure that HSM weights are zero if there is no HSM field
c
	isat_flag = 0
	do j=1,jmax
	do i=1,imax
	   if(s(i,j) .ne. 0.) then
	      isat_flag = 1
	      s(i,j) = s(i,j) - tb(i,j) !diff from background
	   endif
	enddo !i
	enddo !j
c
c.....  Set the weights for the spline.
c
	alf = 100.
	beta = 75.
	a = 50.
	if(isat_flag .eq. 0) a = 0
	alf2a = (1./9.) * alf
c
	write(6,9995) alf, beta, alf2a, a
 9995	format(5x,'Using spline wts: alf, beta, alf2a, a = ',4f10.4)
c
c.....  Now do the spline.
c
	iteration = .true.
	do it=1,itmax
	  cormax = 0.
	  if(iteration) then
	     do j=3,jmax-2
	     do i=3,imax-2
		alfo = alf
		alf2o = alf2a
		if(to(i,j) .eq. 0.) alfo = 0.
		dtxx = t(i+1,j) + t(i-1,j) - 2. * t(i,j)
		dtyy = t(i,j+1) + t(i,j-1) - 2. * t(i,j)
		d4t = 20. * t(i,j) 
     &            - 8. * (t(i+1,j) + t(i,j+1) + t(i-1,j) + t(i,j-1))
     &            + 2.*(t(i+1,j+1)+t(i+1,j-1) + t(i-1,j+1) + t(i-1,j-1))
     &            + (t(i+2,j) + t(i-2,j) + t(i,j+2) + t(i,j-2))
		d2t = dtxx + dtyy
		dtx = (t(i+1,j) - t(i-1,j)) * .5
		dty = (t(i,j+1) - t(i,j-1)) * .5
		sxx = (s(i+1,j) + s(i-1,j) - 2. * s(i,j))
		syy = (s(i,j+1) + s(i,j-1) - 2. * s(i,j))
c       
		res = d4t - a * (d2t - sxx - syy) / beta
     &               + alfo/beta * (t(i,j) - to(i,j)) ! stations
     &               + alf2o/beta * t(i,j)            ! background
		cortm = 20. + a*4./beta + alfo/beta + alf2o/beta
		tcor = abs(res / cortm)
		t(i,j) = t(i,j) - res / cortm * ovr
		if(tcor .le. cormax) go to 5
		cormax = tcor
		ress(it) = tcor
c	write(6,1010) i,j,res,cortm,tcor
c1010	format(1x,2i5,3e12.4)
c	write(6,1009)beta,d4t,d2t,dtxy,dtx,dty,gam,sxx,syy,sxy,sx,sy
 5		continue
	     enddo !i
	     enddo !j

c
c	write(6,1000) it,cormax
	     if(cormax .lt. err) iteration = .false.
	     corhold = cormax
	     ithold = it
	  endif
	enddo !it
c
c
c.....  Add backgrounds back to t, to, and s
c
	do j=1,jmax
	do i=1,imax
	   t(i,j) = t(i,j) + tb(i,j)
	   s(i,j) = s(i,j) + tb(i,j)
	   if(to(i,j) .ne. 0.) to(i,j) = to(i,j) + tb(i,j)
	enddo !i
	enddo !j
c
 6	write(6,1000) ithold ,corhold !it, cormax
 1000	format(1x,i4,e12.4)
	if(name.ne.'NOPLOT' .and. name(1:3).ne.'TB8') then
	   write(9,923)
 923	   format(1x,' solution after spline')
	endif
	if(cormax.lt.err .and. it.eq.1) return
 3	continue
c       
!	print *,' leaving spline'
	return
	end
c
c
	subroutine meso_anl(u,v,p,t,td,theta,dx,dy,q,qcon,qadv,
     &                      thadv,tadv,ni,nj)
c
c*******************************************************************************
c
c	Routine to calculate derived quantities from the LAPS surface 
c	analysis.  From the Meso program written by Mark Jackson...derived 
c	from AFOS MESO sometime during fall of 1988.....?????
c
c	Input units:  u, v  -- m/s          Output:  q           -- g/kg
c	              p     -- mb                    qcon, qadv  -- g/kg/s
c	              theta -- K                     tadv, thadv -- K/s
c	              t, td -- K
c
c	Changes:
c		P. A. Stamus	04-21-89  Changed for use in lapsvanl.
c				05-02-89  VORT calc in main program.
c				05-08-89  Working version...really.
c				05-11-89  Add Moisture advect. calc.
c				04-16-90  Added temp adv.
c				10-30-90  Added boundary routine.
c				04-10-91  Bugs, bugs, bugs....sign/unit errs.
c
c*****************************************************************************
c
	real*4 dx(ni,nj), dy(ni,nj)
c
	real*4 p(ni,nj), td(ni,nj), u(ni,nj), v(ni,nj), thadv(ni,nj)
	real*4 qcon(ni,nj), q(ni,nj), theta(ni,nj), qadv(ni,nj)
	real*4 t(ni,nj), tadv(ni,nj)
c
c
c.....	Calculate mixing ratio.
c.....	Units:  g / kg
c
	do 10 j = 1,nj
	do 10 i = 1,ni
	  tdp = td(i,j) - 273.15                  ! convert K to C
	  tl = (7.5 * tdp) / (237.3 + tdp)
	  e = 6.11 * 10. ** tl
	  drprs = 1. / (p(i,j) - e)   !invert to avoid further divisions.
	  if(p(i,j) .le. 0.0) then
!	    write(6,990) i,j,p(i,j)
	    drprs = 1. / (p(i-1,j-1) - e)
	  endif	
 	  q(i,j) = 622. * e * drprs  !mixing ratio using (0.622*1000) for g/kg.
10	continue
990	format(1x,'Point ',2i4,' has a pressure of ',e12.4)
c
c.....	Compute moisture flux convergence on the laps grid.
c.....	Units:  g / kg / sec
c
	do 30 j=2,nj-1
	do 30 i=2,ni-1
	  ddx1 = ((q(i,j-1) + q(i,j)) * .5) * u(i,j-1)
	  ddx2 = ((q(i-1,j) + q(i-1,j-1)) * .5) * u(i-1,j-1)
	  ddx = (ddx1 - ddx2) / dx(i,j)
	  ddy1 = ((q(i-1,j) + q(i,j)) * .5) * v(i-1,j)
	  ddy2 = ((q(i,j-1) + q(i-1,j-1)) * .5) * v(i-1,j-1)
	  ddy = (ddy1 - ddy2) / dy(i,j)
	  qcon(i,j) = - ddx - ddy
30	continue
	call bounds(qcon,ni,nj)
c
c.....	Compute Theta advection on the laps grid.
c.....	Units:  deg K / sec
c
	do 40 j=2,nj-1
	do 40 i=2,ni-1
	  dth1 = (theta(i,j) - theta(i-1,j)) / dx(i,j)
	  dth2 = (theta(i,j-1) - theta(i-1,j-1)) / dx(i,j)
	  dtdx = (u(i,j-1) + u(i-1,j-1)) * (dth1 + dth2) * .25
	  dth3 = (theta(i,j) - theta(i,j-1)) / dy(i,j)
	  dth4 = (theta(i-1,j) - theta(i-1,j-1)) / dy(i,j)
	  dtdy = (v(i-1,j) + v(i-1,j-1)) * (dth3 + dth4) * .25
	  thadv(i,j) = - dtdx - dtdy   ! deg K/sec
40	continue
	call bounds(thadv,ni,nj)
c
c.....	Compute temperature advection.
c.....	Units:  deg K / sec
c
	do 45 j=2,nj-1
	do 45 i=2,ni-1
	  dth1 = (t(i,j) - t(i-1,j)) / dx(i,j)
	  dth2 = (t(i,j-1) - t(i-1,j-1)) / dx(i,j)
	  dtdx = (u(i,j-1) + u(i-1,j-1)) * (dth1 + dth2) * .25
	  dth3 = (t(i,j) - t(i,j-1)) / dy(i,j)
	  dth4 = (t(i-1,j) - t(i-1,j-1)) / dy(i,j)
	  dtdy = (v(i-1,j) + v(i-1,j-1)) * (dth3 + dth4) * .25
	  tadv(i,j) = - dtdx - dtdy     ! deg K/sec
45	continue
	call bounds(tadv,ni,nj)
c
c.....	Compute Moisture advection on the laps grid.
c.....	Units:  g / kg / sec
c
	do 50 j=2,nj-1
	do 50 i=2,ni-1
	  dqa1 = (q(i,j) - q(i-1,j)) / dx(i,j)
	  dqa2 = (q(i,j-1) - q(i-1,j-1)) / dx(i,j)
	  dqdx = (u(i,j-1) + u(i-1,j-1)) * (dqa1 + dqa2) * .25
	  dqa3 = (q(i,j) - q(i,j-1)) / dy(i,j)
	  dqa4 = (q(i-1,j) - q(i-1,j-1)) / dy(i,j)
	  dqdy = (v(i-1,j) + v(i-1,j-1)) * (dqa3 + dqa4) * .25
	  qadv(i,j) = - dqdx - dqdy     ! g/kg/sec
50	continue
	call bounds(qadv,ni,nj)
c
c.....	Send the fields back to the main program.
c
	return
	end
c
c
	subroutine barnes2(t,imax,jmax,to,smsng,mxstn,npass,fnorm)
c
	real*4 to(imax,jmax), t(imax,jmax), val(imax*jmax)
	real*4 fnorm(0:imax-1,0:jmax-1)
	real*4 h1(imax,jmax), h2(imax,jmax)  !work arrays
c
	integer*4 iob(imax*jmax), job(imax*jmax), dx, dy
c
c
	call zero(h1,imax,jmax)
	call zero(h2,imax,jmax)
c
	badd = 1.e6 - 2.	! bad data value
c
c.....	loop over field npass times 
c
!	print *,' *** In BARNES2 ***'
        ncnt = 0
        do j=1,jmax
        do i=1,imax
          if (to(i,j).ne.0. .and. to(i,j).lt.badd) then
            ncnt = ncnt + 1
            iob(ncnt) = i
            job(ncnt) = j 
            val(ncnt) = to(i,j)
cc	     write(6,999) ncnt, i, j, to(i,j)
          endif
        enddo !i
        enddo !j
	if(ncnt .eq. 0) then
	  print *,' *** NCNT = 0 in BARNES2. ***'
	  return
	endif
 999	format('   ncnt: ',i4,' at ',2i5,f10.3)
c
	do ipass=1,npass
c
	  do j=1,jmax
	  do i=1,imax
	    sum = 0.
	    sumwt = 0.
	    sum2 = 0.
	    sumwt2 = 0.
	    do n=1,ncnt
	      dy = abs(j - job(n))
	      dx = abs(i - iob(n))
	      sum2 = fnorm(dx,dy) * val(n) + sum2
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
c
	    if(sumwt2 .eq. 0.) then
	      if(ipass .eq. 1) then
c		   print *,' got into wierd loop.............'
	            sum2 = 0.
	            sumwt2 = 0.
	            do n=1,ncnt
	              dx = abs(i - iob(n))
	              dy = abs(j - job(n))
	              sum2 = (fnorm(dx,dy) + .01) * val(n) + sum2
	              sumwt2 = sumwt2 + (fnorm(dx,dy) + .01)
	            enddo !n
	            if(sumwt2 .ne. 0.) go to 490
	      else
	        go to 500
	      endif
	    endif
c
490	    continue
	    t(i,j) = sum2 / sumwt2
c
500	  continue
          enddo !i
	  enddo !j
c
	  if(ipass .eq. 2) then
	    call move(h2,to,imax,jmax)
	    call diff(h1,t,t,imax,jmax)
!	    write(9,915)
!915	    format(' after 2nd barnes pass')
!	    write(9,909) rom2
!909	    format(' radm2 = ',f8.4)
	    go to 550
	  endif
!	  write(9,912)
912	  format(' after 1st barnes pass')
	  if(npass .eq. 1) go to 550
	  call move(t,h1,imax,jmax)
 	  call move(to,h2,imax,jmax)
	  do n=1,ncnt
	    val(n) = t(iob(n),job(n)) - val(n)  
	  enddo !n
550	continue
        enddo !ipass
c
!	print *,' *** BARNES2 Done. ***'
	return
	end
c
c
	subroutine barnes_wide(t,imax,jmax,ii,jj,t_ob,numsta,smsng,
     &                         mxstn,npass,fnorm)
c
c.....	Routine to do a Barnes analysis that will consider stations in
c.....	the 't_ob' array that are outside the boundaries of the 't' array.
c
	real*4 t(imax,jmax), t_ob(mxstn) 
	real*4 fnorm(0:imax-1,0:jmax-1)
	real*4 h1(imax,jmax), val(mxstn)
c	
	integer*4 iob(mxstn), job(mxstn), ii(mxstn), jj(mxstn)
	integer*4 dx, dy 
c
!	print *,' *** In BARNES_wide ***'
	call zero(h1,imax,jmax)
	im1 = imax - 1
	jm1 = jmax - 1
c
c.....	loop over field npass times 
c
	ncnt = 0
	do n=1,numsta
          if (t_ob(n).ne.0. .and. t_ob(n).ne.smsng) then
	    ncnt = ncnt + 1
	    iob(ncnt) = ii(n)
	    job(ncnt) = jj(n)
	    val(ncnt) = t_ob(n)
          endif
	enddo !n 
c
	if(ncnt .eq. 0) then
	  print *,' +++ NCNT = 0 in BARNES_WIDE. +++'
	  return
	endif
	write(6,900) ncnt, numsta
900	format('   Selected ',i4,' obs out of ',i4,' total.')
c
	do ipass=1,npass
c
	  do j=1,jmax
	  do i=1,imax
	    sum2 = 0.
	    sumwt2 = 0.
	    do n=1,ncnt
	      dy = min(abs(j - job(n)), jm1) 
	      dx = min(abs(i - iob(n)), im1) 
	      sum2 = fnorm(dx,dy) * val(n) + sum2
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
	    if(sumwt2 .eq. 0.) then
	      if(ipass .eq. 1) then
		print *,' barneswide wierd loop...........'
	        sum2 = 0.
	        sumwt2 = 0.
	        do n=1,ncnt
	          dx = min(abs(i - iob(n)), im1) 
	          dy = min(abs(j - job(n)), jm1) 
	          sum2 = fnorm(dx,dy) * val(n) + sum2
	          sumwt2 = sumwt2 + fnorm(dx,dy) 
	        enddo !n
	        if(sumwt2 .ne. 0.) go to 490 
	      else
	        go to 500
	      endif 
	    endif 
c
490	    continue 
c
	    t(i,j) = sum2 / sumwt2
c
500 	  continue
          enddo !i
	  enddo !j
c
	  if(ipass .eq. 2) then
	    call diff(h1,t,t,imax,jmax)
!	    write(9,915)
!915	    format(' after 2nd barnes pass')
!	    write(9,909) rom2
!909	    format(' radm2 = ',f8.4)
	    go to 550
	  endif
!	  write(9,912)
912	  format(' after 1st barnes pass')
	  if(npass .eq. 1) go to 550
	  call move(t,h1,imax,jmax)
c
	  do n=1,ncnt
	    if(iob(n).lt.1 .or. iob(n).gt.imax  .or.
     &         job(n).lt.1 .or. job(n).gt.jmax) then
	      val(n) = 0. ! which is 't_ob(n)-val(n)' at stns outside the grid
	    else
	      val(n) = t(iob(n),job(n)) - val(n)
	    endif
	  enddo !n
550	continue
        enddo !ipass
c
!	print *,'   leaving barnes_wide'
	return
	end
c
c
	subroutine bounds(x,imax,jmax)
c
c.....	Routine to fill in the boundaries of an array.  Just uses the
c.....	interior points for now.
c
	real*4 x(imax,jmax)
c
	do i=1,imax
	  x(i,1) = x(i,2)
	  x(i,jmax) = x(i,jmax-1)
	enddo !i
	do j=1,jmax
	  x(1,j) = x(2,j)
	  x(imax,j) = x(imax-1,j)
	enddo !j
c
	x(1,1) = x(2,2)
	x(1,jmax) = x(2,jmax-1)
	x(imax,1) = x(imax-1,2)
	x(imax,jmax) = x(imax-1,jmax-1)
c
	return
	end
c
c
	subroutine make_cssi(t,td,pmsl,u,v,cssi,ni,nj,badflag)
c
c======================================================================
c
c	Routine to calculate the CSSI (Rodgers and Maddox 81) at each 
c       LAPS gridpoint.  The temp and dewpt enter in deg F, the MSL 
c       pressure in mb, and the wind components in m/s, which have to be 
c       converted to speed and direction in kts.
c
c	Original version: 05-03-91  Peter A. Stamus NOAA/FSL
c	Changes:          11-11-91  Pass in dummy arrays.
c                         08-26-97  Changes for dynamic LAPS
c
c======================================================================
c
	real*4 t(ni,nj), td(ni,nj), pmsl(ni,nj), u(ni,nj), v(ni,nj)
	real*4 cssi(ni,nj)
c
	real*4 spd(ni,nj), dir(ni,nj)  !work arrays
c
c.....	Start.  Convert u,v in m/s to spd/dir in kts.
c
	call windconvert(u,v,dir,spd,ni,nj,badflag)
	call conv_ms2kt(spd,spd,ni,nj)	
c
c.....	Calculate each of the 4 terms involved, then combine.
c
	do j=1,nj
	do i=1,ni
	  term1 = t(i,j) - 60.		! temperature
	  term2 = 2. * (td(i,j) - 45.)	! moisture
	  term3 = abs(1010.0 - pmsl(i,j))  ! pressure: abs of diff 
	  if(dir(i,j).gt.180. .and. dir(i,j).lt.360.) then	! west wind
	    term4 = -2. * spd(i,j)
	  else							! east wind
	    if(td(i,j) .ge. 45.) then				! that's moist
	      term4 = 2. * spd(i,j)
	    else						! that's not...
	      term4 = spd(i,j)
	    endif
	  endif
c
	  cssi(i,j) = term1 + term2 - term3 + term4
c
	enddo !i
	enddo !j
c
	return
	end
c
c
	subroutine windconvert(uwind,vwind,direction,speed,
     &                         ni,nj,badflag)
c
c======================================================================
c
c       Given wind components, calculate the corresponding speed and 
c       direction.  Hacked up from the windcnvrt_gm program.
c
c
c       Argument     I/O   Type       Description
c      --------	     ---   ----   -----------------------------------
c       UWind         I    R*4A    U-component of wind
c       VWind         I	   R*4A    V-component of wind
c       Direction     O    R*4A    Wind direction (meteoro. degrees)
c       Speed         O    R*4A    Wind speed (same units as input)
c       ni,nj         I    I       Grid dimensions
c       badflag       I    R*4     Bad flag value
c
c       Notes:
c       1.  If magnitude of UWind or VWind > 500, set the speed and 
c           direction set to the badflag value.
c
c       2.  Units are not changed in this routine.
c
c======================================================================
c
	real*4  uwind(ni,nj), vwind(ni,nj)
	real*4  direction(ni,nj), speed(ni,nj)
c
	do j=1,nj
	do i=1,ni
	   if(abs(uwind(i,j)).gt.500. .or. 
     &                           abs(vwind(i,j)).gt.500.) then
	      speed(i,j) = badflag
	      direction(i,j) = badflag
c
	   elseif(uwind(i,j).eq.0.0 .and. vwind(i,j).eq.0.0) then
	      speed(i,j) = 0.0
	      direction(i,j) = 0.0			!Undefined
c
	   else
	      speed(i,j) = 
     &          sqrt(uwind(i,j)*uwind(i,j) + vwind(i,j)*vwind(i,j))  !speed
	      direction(i,j) = 
     &          57.2957795 * (atan2(uwind(i,j),vwind(i,j))) + 180.   !dir
	   endif
	enddo !i
	enddo !j
c
	return
	end
c
c
        subroutine enhance_vis(i4time,vis,hum,topo,ni,nj,kcloud)
c
c==============================================================================
c
c       Routine to call other routines to adjust the visibility analysis
c       based on other data (radar, cloud, etc.).
c            ** May want to put the spline call in here someday...
c
c       Original:  ??-??-93  Peter A. Stamus  NOAA/FSL
c       Changes:   02-03-94  Rewritten
c                  08-26-97  Changes for dynamic LAPS.
c       
c       Notes:
c          1.  The variables here are:
c                  i4time = Time for this analysis.
c                  vis    = Visibility (units are not changed) 
c                  hum    = Relative humidity (0 to 100 percent)
c                  topo   = LAPS topography (meters)
c                  ni,nj  = LAPS grid dimensions
c                  kcloud = Cloud grid dimension in vertical
c          2.  Units of visibility (miles, meters) are not changed
c              in this routine.
c
c==============================================================================
c
        real*4 vis(ni,nj), hum(ni,nj), topo(ni,nj)
        real*4 vismod(ni,nj)  !work array
c
        print *,' In enhance_vis routine...'
c
c..... Radar adjustment.            ! still disabled...2-3-94 pas (no data)
c
c       call constant(vis,-10.,ni,nj)
c       call get_radar_visibility(i4time,vis,istatus)
c
c
c..... Low cloud/humidity adjustment.
c.....................................
c..... Get the modification array using the cloud data from LC3 and the
c..... surface relative humidity.  The multiply the visibilities by the
c..... modification factor to get the adjusted visibility.
c
        call get_vismods(i4time,hum,topo,vismod,ni,nj,kcloud)
c
        do j=1,nj
        do i=1,ni
           vis(i,j) = vis(i,j) * vismod(i,j)
        enddo !i
        enddo !j
c
c..... that's it...
c
        print *,' Enhance_vis routine...done.'
        return
        end
c
c
        subroutine get_vismods(i4time,hum,topo,vismod,ni,nj,kcloud)
c
c==============================================================================
c
c       Routine to set a visibility adjustment based on low cloud data from
c       the LAPS 3-D cloud analysis and the LAPS surface humidity analysis.
c       The adjustment is the percentage (0-1) that the visibility is 
c       reduced for given cloud amounts/humidities.  The adjustment is put
c       into the vismod array, and is multiplied by the vis array in the
c       calling routine (enhance_vis).
c
c       Original:  ??-??-93  Peter A. Stamus
c       Changes:   02-03-94  Rewritten
c                  08-26-97  Changes for dynamic LAPS
c	           08-19-98  Initialize k_hold array.
c
c==============================================================================
c
        real*4 hum(ni,nj), vismod(ni,nj), topo(ni,nj)
	real*4 cld_hts(kcloud), cld_pres(kcloud)
        real*4 clouds_3d(ni,nj,kcloud)
c
        integer*4 k_hold(ni,nj), lvl(kcloud)
c
        character ext*31, var(kcloud)*3, comment(kcloud)*125
        character units(kcloud)*10, lvl_coord(kcloud)*4, dir*256
c
c..... Start by setting up the default values for vismod (1.0=no adjustment)
c
        call constant(vismod,1.0,ni,nj)
c
c..... Get the low cloud data from the nearest LC3 file (timewise).
c
        icnt = 0
        i4time_c = i4time
        do k=1,kcloud
           lvl(k) = k
           var(k) = 'lc3'
        enddo !k
        ext = 'lc3'
	call get_directory('lc3', dir, len)
 500    call read_laps_data(i4time_c,dir,ext,ni,nj,kcloud,kcloud,
     &       var,lvl,lvl_coord,units,comment,clouds_3d,istatus)
c
        if(istatus .ne. 1) then
           if(istatus .eq. 0) then  !no data
              if(icnt .lt. 1) then  !just try 1-hr for now.
           print *,' No data for given i4time...trying 1-hr earlier.'
               i4time_c = i4time_c - 3600
               icnt = icnt + 1
               go to 500
              else
               print *,' LC3 data too old.'
               print *,' No visibility modification done.'
               return
              endif 
           else
              print *,' Bad return from LC3 read: istatus = ',istatus
              print *,' No visibility modification done.'
              return
           endif
        endif
c
c..... Check time difference.  Don't use if cloud analysis is too old.
c
c        if((i4time - i4time_nearest) .gt. 5400) then  !1.5 hours
c           print *,' LC3 data too old.'
c           print *,' No visibility modification done.'
c           return
c        endif
c
c..... Decode the cloud heights and pressures.
c
        print *,' Got cloud data.'
        do k=1,kcloud
           read(comment(k),100,err=999) cld_hts(k), cld_pres(k)
100        format(2e20.7)
        enddo !k
c
c..... Find the vertical level from 'cld_hts' just below the surface.  Will
c..... then start cloud checks at the next level up.
c
        do j=1,nj
        do i=1,ni
	   k_hold(i,j) = 0
           do k=1,kcloud
              if(cld_hts(k) .lt. topo(i,j)) k_hold(i,j) = k
           enddo !k
        enddo !i
        enddo !j
c
c..... Now loop over the grid and check the lowest 3 levels above the 
c..... surface.  Check for fog first, then check for low clouds.
c
        do j=1,nj
        do i=1,ni
           k_start = k_hold(i,j) + 1     ! 1st level above the surface
           k_end = k_start + 3           ! this is usually within 300 m
c
c.....     Check for fog in the layer just above the surface.
c
           if(clouds_3d(i,j,k_start) .gt. 0.65) then
              if(hum(i,j) .gt. 70.) vismod(i,j) = 0.90
              if(hum(i,j) .gt. 80.) vismod(i,j) = 0.75
              if(hum(i,j) .gt. 90.) vismod(i,j) = 0.55
              if(hum(i,j) .gt. 95.) vismod(i,j) = 0.35
              go to 200
           endif
c
c.....     If no fog in lowest layer, find the maximum value in the 3 levels
c.....     above the surface.  Then adjust the vismod based on humidity.
c
           amax_layer = 0.
           do k=k_start,k_end
              if(clouds_3d(i,j,k) .gt. amax_layer) then
                 amax_layer = clouds_3d(i,j,k)
              endif
           enddo !k
c
           if(amax_layer .gt. 0.65) then
              if(hum(i,j) .gt. 70.) vismod(i,j) = 0.95
              if(hum(i,j) .gt. 80.) vismod(i,j) = 0.80
              if(hum(i,j) .gt. 90.) vismod(i,j) = 0.60
              if(hum(i,j) .gt. 95.) vismod(i,j) = 0.40
           endif
c
200     continue
        enddo !i
        enddo !j
c
c..... That is all.
c
        return ! normal return
c
 999    print *,' Error reading comment field from LC3.'
        print *,' No visibility modification done.'
        return
        end
c
c
	subroutine dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
c
c=====================================================================
c
c     Routine to calculate the weights to be used in the Barnes
c     analysis.  The data density is used to set the cutoff for
c     the response function.  Then that cutoff is used to calculate
c     the exp, based on differences so that no additional distance
c     calculations are required in the Barnes routine.  All of this
c     is done in gridpoint space.
c
c     Original:  07-14-95  P. Stamus, NOAA/FSL
c     Changes:   P. Stamus  08-28-97  Declare dx,dy integers.
c                           09-10-97  Bag include. Pass in fnorm.
c
c     Notes:
c
c       1.  If variable 'rom2' is passed in as zero, it is calculated
c           from the data density.  Otherwise, the value passed in is
c           used in the weight calculation.
c
c       2.  The response for 2d waves is hard-wired in this routine.
c           This is the 'con' variable, and comes from the eqn:
c                     D = exp -(pi**2 R**2)/lamba**2
c           If we set D (the response) to our desired cutoff, set 
c           lamba to the desired wavelength in gridpt space (2d),
c           then solve for R in terms of d, we get the 'con' value
c           (i.e.,  R = (con)*d).  Here are some values for different
c           cutoffs:
c                     D = 0.01     R = 1.36616d
c                         0.10     R = 0.96602d
c                         0.25     R = 0.74956d
c                         0.50     R = 0.53002d
c
c=====================================================================
c
	real*4 fnorm(0:imax-1,0:jmax-1)
	integer*4 dx, dy
c
c.... First, find the area that each ob covers in gridpt space (this
c.... of course assumes a uniform coverage).
c
cc	con = 0.96602     ! resp of 0.10
cc	con = 0.74956     ! resp of 0.25
	con = 0.53002     ! resp of 0.50
	if(rom2 .eq. 0.) then
	   area = float(imax * jmax) / n_obs_var
	   d = sqrt( area )
	   rom2 = 1. / ((con * con) * (d * d))
	   write(6,900) n_obs_var, area, d, rom2
 900	   format(1x,'Num obs: ',i5,'  Area: ',f8.2,'  d: ',f8.2,
     &       '  rom2: ',f8.5)
	else
	   d = sqrt(1./(con * con * rom2))
	   write(6,902) rom2, d
 902	   format(' Using preset rom2 of: ',f8.5,'  Calc d: ',f8.2)
	endif
c
c.... Now calculate the weights for all the possible distances.
c
	pi = 4. * atan(1.)
	fno = 1. / (sqrt(2. * pi))
c
	do dy=0,jmax-1
	do dx=0,imax-1
	   rr = float(dx*dx + dy*dy)
	   fnorm(dx,dy) = fno * (exp( -(rr * rom2)))
	enddo !dx
	enddo !dy
c
c.... That's it.
c
	return
	end
c
c
	subroutine calc_beta(d,obs_error,beta)
c
c=======================================================================
c
c       Routine to calculate the 'beta' weight for the spline.  'Beta'
c       is calculated based on the 'd' from the gridpt to data distance
c       and an expected observation error for the partictular variable.
c
c
c       Original:  07-19-95  P. Stamus, NOAA/FSL
c       Changes:
c
c=======================================================================
c
	pi = 4. * atan( 1. )
	pi4 = pi * pi * pi * pi
	d4 = d * d * d * d
	alpha = -99.9
	beta = -99.9
c
	if(obs_error .ne. 0.) then
	   alpha = 1. / (obs_error * obs_error)
	else
	   print *,' **ERROR. obs_error = 0 in CALC_BETA.**'
	   go to 100
	endif
c
	beta = alpha * d4 / pi4
c
 100	continue
	write(6,900) obs_error, d, alpha, beta
 900	format(1x,'obs error: ',f9.4,'  d: ',f9.4,
     &                      '  alpha: ',f9.4,'  beta: ',f15.4)
	if(beta .eq. 0.) then
	   print *,' **ERROR. beta = 0 in CALC_BETA.**'
	   beta = -99.9
	endif
c
	return
	end
