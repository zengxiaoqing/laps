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
c
	subroutine mean_lapse(num_sfc,elev,t,td,a_t,b_t,a_td,b_td,
     &                        hbar,badflag)
c
c*******************************************************************************
c
c	Routine to calculate a mean lapse rate from observed surface data.
c	Values returned are the coefficients for the regression equations
c	for the temperatures and dew points.  Also calculates the mean 
c	elevation for the surface stations.
c
c	Changes:
c		P.A. Stamus	12-01-88  Original (from J. McGinley)
c				12-22-88  Added consistency check on td.
c				05-24-90  Rig to return std lapse rate.
c                               08-25-97  Changes for dynamic LAPS.
c
c	Inputs/Outputs:
c
c	   Variable     Var Type    I/O   Description
c	  ----------   ----------  ----- -------------
c	   num_sfc         I         I    Number of surface stations.
c	   elev            RA        I    Station elevation.
c	   t               RA        I    Temperature.
c	   td              RA        I    Dew point temperature.
c	   a_t             R         O    'a' regression value for temp.(intrcp)
c          b_t             R         O    'b'      "       "    "    "  (lapse)
c	   a_td            R         O    'a'      "       "    "  dewpt.
c	   b_td            R         O    'b'      "       "    "    "
c	   hbar            R         O    Mean elevation of the stations.
c          badflag         R         I    Bad flag value.
c
c	User Notes:
c
c	1. Units are not changed in this routine.
c
c*******************************************************************************
c
	real*4 elev(num_sfc), t(num_sfc), td(num_sfc)
c
c.....	Set up storage variables.
c
	cnt = 0.
	cntd = 0.
	sumht = 0.
	sumh = 0.
	sumt = 0.
	sumh2 = 0.
	sumt2 = 0.
	sumtd = 0.
	sumhtd = 0.
c
c.....	Gather sums and then calculate the 'a' and 'b' for the regression
c.....	equation y = bz + a, for both the temperature and dew point.  The
c.....	'a' is the intercept with sea level, and the 'b' is the lapse rate.
c.....	Also calculate the mean elevation of the stations.
c
	do 10 i=1,num_sfc
	  if(elev(i).le.badflag .or. t(i).le.badflag) go to 10
	  sumht = (elev(i) * t(i)) + sumht
	  sumh = elev(i) + sumh
	  sumh2 = (elev(i) * elev(i)) + sumh2
	  sumt = t(i) + sumt
	  cnt = cnt + 1.
c
	  if(td(i) .le. badflag) go to 10
	  sumtd = td(i) + sumtd
	  sumhtd = (elev(i) * td(i)) + sumhtd
	  cntd = cntd + 1.
10	continue
c
	b_t = (cnt*sumht - sumh*sumt) / (cnt*sumh2 - sumh*sumh)
	a_t = (sumt - b_t * sumh) / cnt
c
	b_td = (cntd*sumhtd - sumh*sumtd) / (cntd*sumh2 - sumh*sumh)
	a_td = (sumtd - b_td * sumh) / cntd
c
	hbar = sumh / cnt
c
c.....	Do a consistency check on the dewpoint regression.  If msl intercept 
c.....	is below zero or if the dewpoint lapse rate is positive, set td slope 
c.....	to t slope and slide intercept over.
c
	if(a_td.lt.0. .or. b_td.gt.0. .or. a_td.gt.a_t) then
	  write(6,900)
900	  format(1x,'++ Suspect dewpoint regression in MEAN_LAPSE. ++')
	  write(6,901) a_td, b_td
901	  format(1x,'  MSL intercept = ',e12.4,'  Slope = ',e12.4)
	  b_td = b_t * .6
	  a_td = (sumtd / cntd) - b_td * hbar  !mean td - new td lapse @ mean h
	  write(6,902)
902	  format(3x,'  Setting values to:')
	  write(6,901) a_td, b_td
	endif
c
c.....	Now change that stuff to std lapse rate...in deg F...temp. fix
c
	b_t = -.01167
	a_t = 59.
	b_td = -.007
	a_td = 40.
c
c.....	End of routine
c
	return
	end
c
c
	real function alt_2_sfc_press(alt,elev)
c
c*****************************************************************************
c
c	Routine to reduce an SAO altimeter to surface pressure.  Have to 
c	use MSL pressure and Standard temperature since that's how they
c	got altimeter in the first place.
c
c	Original:	07-12-88 (from CNVRT_SFC_DATA_2CMM by M. McCoy)
c							Peter A. Stamus
c
c	Changes:	07-20-88	Changed from subroutine to function.
c
c*****************************************************************************
c
	data const/0.190284/, const1/5.2553026/		!const1 = 1/const
	data alapse/0.0065/, tstd/288.15/, amslp/1013.25/
c
	badflag = -99.9
c
	if(alt .eq. badflag) then
	  alt_2_sfc_press = badflag
	  return
	endif
c
	term1 = (alapse * elev) / tstd
	term2 = 1 - ( ( (amslp / alt) ** const) * term1)
	alt_2_sfc_press = alt * (term2 ** const1)
c
	return
	end
c
c
	subroutine reduce_p(temp,dewp,pres,elev,lapse_t,
     &                          lapse_td,redpres,ref_lvl,badflag)
c
c
c================================================================================
c   This routine is designed to reduce the mesonet plains stations' pressure
C   reports to the elevation of the Boulder mesonet station, namely 1612 m.  The
C   hydrostatic equation is used to perform the reduction, with the assumption
C   that the mean virtual temperature in the layer between the station in ques-
C   tion and Boulder can approximated by the station virtual temperature.  This
C   is a sufficient approximation for the mesonet plains stations below about
C   7000 ft.  For the mountain and higher foothill stations (Estes Park,
C   Rollinsville, Ward, Squaw Mountain, and Elbert), a different technique needs
C   to be utilized in order to decently approximate the mean virtual temperature
C   in the layer between the station and Boulder. Ideas to do this are 1) use
C   the free air temperature from a Denver sounding, or 2) use the data from
C   each higher station to construct a vertical profile of the data and iterate
C   downward to Boulder.

C	D. Baker	 2 Sep 83  Original version.
C	J. Wakefield	 8 Jun 87  Changed ZBOU from 1609 to 1612 m.
c	P. Stamus 	27 Jul 88  Change ranges for good data tests.
c	P. Stamus	05 Dec 88  Added lapse rate for better computation
c	 				of mean virtual temps.
c			19 Jan 89  Fixed error with lapse rates (sheeze!)
c			19 Dec 89  Change reduction to 1500 m.
c			20 Jan 93  Version with variable reference level.
c                       25 Aug 97  Changes for dynamic LAPS.
c
c================================================================================
c
	real lapse_t, lapse_td
!	DATA GOR,ZBOU,CTV/.03414158,1612.,.37803/
	data gor,ctv/.03414158,.37803/
		!GOR= acceleration due to gravity divided by the dry air gas
		!     constant (9.8/287.04)
		!F2M= conversion from feet to meters
		! *** zbou is now the standard (reduction) level 12-19-89 ***
		!CTV= 1-EPS where EPS is the ratio of the molecular weight of
		!     water to that of dry air.


C** Check input values......good T, Td, & P needed to perform the reduction.
	if(dewp.gt.temp .or. pres.le.620. .or. pres.gt.1080. .or.
     &      temp.lt.-30. .or. temp.gt.120. .or. dewp.lt.-35. .or.
     &      dewp.gt.90.) then
	 redpres = badflag	!FLAG VALUE RETURNED FOR BAD INPUT
	 return
	endif

	dz= elev - ref_lvl	!thickness (m) between station & reference lvl
	dz2 = 0.5 * dz		! midway point in thickness (m)
	t_mean = temp - (lapse_t * dz2)	! temp at midpoint (F)
	td_mean = dewp - (lapse_td * dz2)	! dewpt at midpoint (F)
	TD= 0.55556 * (td_mean - 32.)		! convert F to C
	e= esw(td)		!saturation vapor pressure
	tkel= 0.55556 * (t_mean - 32.) + 273.15	! convert F to K
	tv= tkel/(1.-ctv*e/pres)	!virtual temperature

	redpres= pres*exp(gor*(dz/tv))	!corrected pressure
c
	return
	end
c
c
	subroutine extract(a,imax,jmax,i,j,ix,jy)
c
c.....	Routine designed to zero out a grouping of points about a 
c.....  named point i,j.  All points from i-ix to i+ix, and 
c.....	j-jy to j+jy will be zeroed.  This is aimed at processing 
c.....	Band 8 temperatures to remove cloud edges.
c
	real*4 a(imax,jmax)
c
	jmx=j+jy
	imx=i+ix
	jmn=j-jy
	imn=i-ix
	if(jmx.gt.jmax) jmx=jmax
	if(imx.gt.imax) imx=imax
	if(jmn.lt.1)jmn=1
	if(imn.lt.1)imn=1
c
	do jj=jmn,jmx
	do ii=imn,imx
	  a(ii,jj) = 0.
	enddo !ii
	enddo !jj
c
	return
	end
c
c
	subroutine decompwind_gm(dd,ff,ucomp,vcomp,status)
C***Decompose vector wind into U and V

C	J. Wakefield	16 Sep 83	Original version
c	P. Stamus       29 Apr 93	Unix version

C Argument	I/O	Type			Description
C --------	---	----	-----------------------------------------------
C DD		 I	R*4	Wind direction (meteorological degrees)
C FF		 I	R*4	Wind speed
C UComp		 O	R*4	U-component of wind
C VComp		 O	R*4	V-component of wind
C Status	 O	I*4	Standard system status

	parameter	(flag=1e37)

	real*4		dd,ff,ucomp,vcomp
	integer*4	status

	Status = 1

	if(ff .eq. .0) then
	 ucomp = .0
	 vcomp = .0
	elseif(dd .ge. .0 .and. dd .le. 360) then
	 angle = .01745239 * dd			!Radians
	 ucomp = -ff * sin(angle)
	 vcomp = -ff * cos(angle)
	else
	 ucomp = flag
	 vcomp = flag
	 status = 0
	endif
c
	return
	end
c
c
      subroutine stats(x,ni,nj)
c
c========================================================================
c
c     Routine to calculate some statistical scores for a 2-d grid.
c     The following scores are calculated and printed:
c
c     Max, Min     Absolute Average     Std. Deviation     Mean
c     Z-score for Max, Min              Range              Skewness
c     Variance     Kurtosis
c
c     Original: 04-86 (from an old program)
c     LAPS version:  07-26-94  P. Stamus
c
c========================================================================
c
      real*4 x(ni,nj)
c
c..... Start by zeroing some counters.
c
      amean = 0.
      ave = 0.
      sum = 0.
      sum_a = 0.
      sum_v = 0.
      sum_sk = 0.
      sum_kr = 0.
      range = 0.
      st_dev = 0.
      var = 0.
      amax = -1.e25
      amin =  1.e25
      pts = float(ni * nj)
c
c..... Calculate the means and range.
c
      do j=1,nj
      do i=1,ni
         if(x(i,j) .lt. amin) then
            amin = x(i,j)
            imin = i
            jmin = j
         endif
         if(x(i,j) .gt. amax) then
            amax = x(i,j)
            imax = i
            jmax = j
         endif
         sum = sum + x(i,j)
         sum_a = sum_a + abs(x(i,j))
      enddo !i
      enddo !j
c
      amean = sum / pts
      ave = sum_a / pts
      range = amax - amin
c
c..... Now calculate the variance, stdev, etc.
c
      do j=1,nj
      do i=1,ni
         dif = x(i,j) - amean
         dif2 = dif * dif
         sum_v = sum_v + dif2
         sum_sk = sum_sk + (dif2 * dif)
         sum_kr = sum_kr + (dif2 * dif2)
      enddo !i
      enddo !j
c
      var = sum_v / (pts - 1)

c .... test divide by zero ... added by Dan Birkenheuer (11/15)

      if (var .ne. 0.0) then
        st_dev = sqrt( var )
        z_max = (amax - amean) / st_dev
        z_min = (amin - amean) / st_dev
        coef_sk = (sum_sk / pts) / sqrt( (sum_v / pts) ** 3 )
        coef_kr = (sum_kr / pts) / ( (sum_v / pts) ** 2 )
      endif
c
c..... Write out the stats.
c
	write(6,900) amean, ave, st_dev, coef_sk
 900  format(1x,'Mean:',g12.4,2x,'AbsAve:',g12.4,2x,'StDev:',g12.4,
     &       2x,'Skew:',g12.4)
c
	write(6,910) amax,imax,jmax,amin,imin,jmin
 910  format(1x,'Max:',g12.4,' @ ',i4,',',i4,2x,'Min:',g12.4,' @ ',
     &       i4,',',i4)
c
	write(6,920) z_max, z_min, range
 920	format(1x,'Z-Max:',g12.4,2x,'Z-Min:',g12.4,2x,'Range:',g12.4)
c
c.... That's it.  Let's go home.
c
      return
      end
c
c
	subroutine clouds(imax,jmax,topo,t,smsng,tb8,
     &                    i4time,laps_cycle_time,lat,lon,
     &                    r_missing_data)
c
c*************************************************************************
c
c	Routine to process band 8 brightness temps for clouds.
c	
c	Changes: 11-01-91  Changes for new LAPS grids.
c                07-20-94  New version.    S. Albers
c                08-25-97  Changes for dynamic LAPS.   P. Stamus
c
c*************************************************************************
c
	real*4 t(imax,jmax), tb8(imax,jmax), lapse_t
        real*4 lat(imax,jmax), lon(imax,jmax), topo(imax,jmax)
c
        real*4 cvr_snow(imax,jmax), t_gnd_k(imax,jmax)  !work arrays
	real*4 t_est(imax,jmax), dtb8(imax,jmax)        !work arrays
c
c
	call zero(t_est,imax,jmax)

	do 11 j=1,jmax
	do 11 i=1,imax
	  terr = topo(i,j)
          t_est(i,j) = t(i,j)
11	continue


        call get_ground_temperature(i4time,laps_cycle_time
     &                    ,imax,jmax,lat,lon,r_missing_data
     &                             ,cvr_snow,t_est,t_gnd_k)


!       Call the ground temperature routine

	do 21 j=1,jmax
	do 21 i=1,imax

	  if(tb8(i,j).ne.smsng .and. tb8(i,j).ne.0.) then
!           t_compare = t_est(i,j)
            t_compare = t_gnd_k(i,j)
	    dtb8(i,j) = tb8(i,j) - t_compare
	  endif

21	continue
c
c.....  Now use the std dev data to check the Band 8 data for clouds.
c
	icnt = 0

!       Check number of points thrown out as a function of offset to check
!       for IR navigation errors

        thresh1 = 10.

	do j=1,jmax
	do i=1,imax

!         Estimate whether tb8 - t < -1.5 * stdt
	  if(dtb8(i,j) .lt. -thresh1) then	             ! probably clds

c.....	set the point and surrounding points to zero...clouds.
	    tb8(i,j) = 0.                         ! Set to zero as a cloud flag
	  call extract(tb8,imax,jmax,i,j,7,7)
	    icnt = icnt + 1

          endif

	enddo	! on i
	enddo	! on j
	write(6,951) icnt, thresh1
951	format(1x,i6,' points removed for thresh = ',f8.1)
	return
	end
c
c
	subroutine lp_fire_danger (ni,nj,lp_10krh,lp_10kt,lp_10kws,
     &                             soil_moist,snow_cover,topo,
     &                        ismoist,isnow,lp_fire_index,i_status)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	Routine to caluclate fire danger from LAPS surface data.
c
c	Original version:  Matt Kelsch  05-18-93
c	Changes:           Pete Stamus  10-14-93  Set up for LAPS use.
c                                       07-27-94  Unix version.
c                                       07-29-94  Change units on tests.
c                                       02-24-95  Add snow, topo caps.
c                                       08-25-97  Changes for dynamic LAPS.
c
c	Notes:
c
c       The if-then structure of this routine will produce a 0 to 20 unit
c	scale of the fire danger based on the current conditions observed 
c	within LAPS.  Each unit represents an increase in fire danger ...
c
c	The four components to the fire danger index and the number of
c	points out of 20 each component may contribute are given below
c	(note the relative humidity and wind speed are most influential):
c		i.   relative humidity (7),
c		ii.  wind speed (7),
c		iii. soil moisture (3),
c		iv.  temperature (3).
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	implicit none
c
	real*4 lp_10krh(ni,nj), lp_10kt(ni,nj), lp_10kws(ni,nj)
	real*4 snow_cover(ni,nj), soil_moist(ni,nj), topo(ni,nj)
	real*4 lp_fire_index(ni,nj)
c
	integer*4 i_lp,	j_lp, i_rh, i_temp, i_wspeed, i_soil,
     &            i_status
c
c	...BEGIN...

c	... loop over all of the LAPS gridpoints
c
	do j_lp = 1,nj
	  do i_lp = 1,ni
c
c	... first, the relative humidity component, I_RH, of fire danger...
c       ... the RH enters this routine in a range of 0 to 100 Percent ...
c
	    i_rh = 0
c	
	    if (lp_10krh(i_lp,j_lp) .ge. 75.) i_rh = -1 
c
	    if (lp_10krh(i_lp,j_lp) .lt. 60.) then
	      i_rh = i_rh + 1
	      if (lp_10krh(i_lp,j_lp) .lt. 50.) then
	        i_rh = i_rh + 1
	        if (lp_10krh(i_lp,j_lp) .lt. 41.) then
	          i_rh = i_rh + 1
	          if (lp_10krh(i_lp,j_lp) .lt. 33.) then
	            i_rh = i_rh + 1
	            if (lp_10krh(i_lp,j_lp) .lt. 25.) then
	              i_rh = i_rh + 1
	              if (lp_10krh(i_lp,j_lp) .lt. 17.) then
	                i_rh = i_rh + 1
	                if (lp_10krh(i_lp,j_lp) .lt. 9.) then
	                  i_rh = i_rh + 1
	                endif
	              endif
		    endif
	          endif
	        endif
	      endif
	    endif

c	... second, the wind speed component, I_WIND, of fire danger ...
c       ... the wind speed enters this routine in m/s...
	    i_wspeed = 0
c
	    if (lp_10kws(i_lp,j_lp) .gt. 1.79) then  ! 4 mph
	      i_wspeed = i_wspeed + 1
	      if (lp_10kws(i_lp,j_lp) .gt. 3.58) then  ! 8 mph
	        i_wspeed = i_wspeed + 1
	        if (lp_10kws(i_lp,j_lp) .gt. 5.36) then  ! 12 mph
	          i_wspeed = i_wspeed + 1
	          if (lp_10kws(i_lp,j_lp) .gt. 6.71) then  ! 15 mph
	            i_wspeed = i_wspeed + 1
	            if (lp_10kws(i_lp,j_lp) .gt. 8.05) then  ! 18 mph
	              i_wspeed = i_wspeed + 1
	              if (lp_10kws(i_lp,j_lp) .gt. 9.39) then  ! 21 mph
	                i_wspeed = i_wspeed + 1
	                if (lp_10kws(i_lp,j_lp) .gt. 10.73) then  ! 24 mph
	                  i_wspeed = i_wspeed + 1
	                endif
	              endif
		    endif
	          endif
	        endif
	      endif
	    endif

c	... third, the soil moisture component, I_SOIL, of fire danger ...

	    i_soil = 0
c
	    if (ismoist .ne. 1) go to 600     ! no soil moisture analysis
	    if (soil_moist(i_lp,j_lp) .ge. 75.) i_soil = -1 
	    if (soil_moist(i_lp,j_lp) .lt. 50.) then
	      i_soil = i_soil + 1
	      if (soil_moist(i_lp,j_lp) .lt. 35.) then
	        i_soil = i_soil + 1
	        if (soil_moist(i_lp,j_lp) .lt. 20.) then
	          i_soil = i_soil + 1
	        endif
	      endif
	    endif

c	... fourth, the temperature component, I_TEMP, of fire danger ...
c       ... temp enters this routine in deg K
 600	    i_temp = 0
c
	    if (lp_10kt(i_lp,j_lp) .ge. 298.15) then  ! 77 F
	      i_temp = i_temp + 1
	      if (lp_10kt(i_lp,j_lp) .ge. 303.15) then  ! 86 F
	        i_temp = i_temp + 1
	        if (lp_10kt(i_lp,j_lp) .ge. 308.15) then  ! 95 F
	          i_temp = i_temp + 1
	        endif
	      endif
	    endif

c	... the LAPS fire danger product is the sum of the four components ...

	    lp_fire_index(i_lp,j_lp) = 
     &                         i_rh + i_wspeed + i_soil + i_temp

c.....  If there's no soil moisture analysis, increase index by 7.5% to
c.....  allow for at least some influence.

	    if(ismoist .ne. 1) then
	       lp_fire_index(i_lp,j_lp) = lp_fire_index(i_lp,j_lp) + 
     &                   0.075 * lp_fire_index(i_lp,j_lp)
	    endif

c.....  Cap the index under the following conditions:
c.....  1. The elevation is greater than treeline (about 11000 feet)

	    if(topo(i_lp,j_lp) .gt. 3350.) then
	       if(lp_fire_index(i_lp,j_lp) .gt. 10.)
     &                             lp_fire_index(i_lp,j_lp) = 10.
	    endif

c.....  2. The snow cover is greater than 25%.

	    if(isnow .ne. 1) go to 700        ! no snow data
	    if(snow_cover(i_lp,j_lp) .gt. .25) then
	       if(lp_fire_index(i_lp,j_lp) .gt. 5.)
     &                             lp_fire_index(i_lp,j_lp) = 5.
	    endif

c.....  Check for index below zero.

 700	    continue
	    if (lp_fire_index(i_lp,j_lp) .lt. 0.0) then
	        lp_fire_index(i_lp,j_lp) = 0.0
	    endif

	  enddo !i_lp
	enddo!j_lp

c	...END...

	i_status = 1
	return
c
97	i_status = -1
	print *,' *** Bad status in LP_FIRE_DANGER ***'
c
	return
	end
c
c
      subroutine heat_index(t,rh,hi,ni,nj,badflag)
c
c====================================================================
c
c     Routine to calculate a heat index.  Based on a formula 
c     by Lans Rothfusz, NWS.  Seems to provide valid HI numbers
c     for temperatures above 75 deg F.
c
c     Original:  07-18-95  P. Stamus, NOAA/FSL
c     Changes:  P. Stamus  08-25-97  Return badflag if Temp < 75F
c                                    Change units returned to K.
c                          01-20-98  T in as deg K.
c
c     Notes:
c
c       1.  Inputs:
c                    rh = Relative Humidity (0 to 100 %)
c                    t  = Temperature (deg K)
c	             ni, nj  = Grid dimensions
c                    badflag = Bad flag value
c
c           Output:
c                    hi = Heat Index (deg K)
c
c       2.  If the temperature is below 75 deg F, no heat index is
c           calculated and the point is set to "badflag".
c
c====================================================================
c
      real*4 t(ni,nj), rh(ni,nj), hi(ni,nj)
c
      do j=1,nj
      do i=1,ni
c
	 temp = ( 1.8 * (t(i,j) - 273.15) ) + 32.         ! K to F
c
	 if(temp .lt. 75.) then
	    hi(i,j) = badflag

	 else
	    rh1 = rh(i,j)                                  ! %
	    rh2 = rh1 * rh1
	    t1 = temp          
	    t2 = t1 * t1
c
	    heat = -42.379 + (2.04901523  * t1)
     &                     + (10.1433312  * rh1)
     &                     - (0.22475541  * t1  * rh1)
     &                     - (6.83783e-3  * t2)
     &                     - (5.481717e-2 * rh2)
     &                     + (1.22874e-3  * t2  * rh1)
     &                     + (8.52e-4     * rh2 * t1)
     &                     - (1.99e-6     * t2  * rh2)
c
	    hi(i,j) = ((heat - 32.) * 0.55555555) + 273.15 ! F to K
c
	 endif

      enddo !i
      enddo !j
c
c.... That's it.  Let's go home.
c
      return
      end
c
c
      subroutine verify(field,ob,stn,n_obs_b,title,iunit,
     &                  ni,nj,mxstn,x1a,x2a,y2a,ii,jj,
     &                  field_ea,badflag)
c
c======================================================================
c
c     Routine to interpolate a field back to station locations, to 
c     compare the analysis to the original obs.
c
c     Original: P.Stamus, NOAA/FSL  08-07-95
c     Changes:  
c               P.Stamus  08-14-95  Added mean.
c                         08-25-97  Changes for dynamic LAPS
c                         05-13-98  Added expected accuracy counts.
c
c     Notes:
c
c======================================================================
c
	real*4 field(ni,nj), ob(mxstn), interp_ob
	real*4 x1a(ni), x2a(nj), y2a(ni,nj)
	integer*4 ii(mxstn), jj(mxstn)
	character title*40, stn(mxstn)*3, stn_mx*3, stn_mn*3
c
c.... Start.
c
	num = 0
	num_ea1 = 0
	num_ea2 = 0
	num_ea3 = 0
	abs_diff = 0.
	sum = 0.
	amean = 0.
	diff_mx = -1.e30
	diff_mn = 1.e30
	write(iunit,900) title
 900	format(/,/,2x,a40,/)
c
	ea1 = field_ea
	ea2 = field_ea * 2.
	ea3 = field_ea * 3.
c
c....   Find the 2nd derivative table for use by the splines.
c
	call splie2(x1a,x2a,field,ni,nj,y2a)
c
c....   Now call the spline for each station in the grid.
c
	do i=1,n_obs_b
	   if(ii(i).lt.1 .or. ii(i).gt.ni) go to 500
	   if(jj(i).lt.1 .or. jj(i).gt.nj) go to 500
	   aii = float(ii(i))
	   ajj = float(jj(i))
	   call splin2(x1a,x2a,field,y2a,ni,nj,aii,ajj,interp_ob)
c
	   if(ob(i) .le. badflag) then
	      diff = badflag
	   else
	      diff = interp_ob - ob(i)
	      sum = diff + sum
	      adiff = abs(diff)
	      abs_diff = abs_diff + adiff
	      num = num + 1
c
	      if(adiff .gt. diff_mx) then
		 diff_mx = adiff
		 stn_mx = stn(i)
	      endif
	      if(adiff .lt. diff_mn) then
		 diff_mn = adiff
		 stn_mn = stn(i)
	      endif
c
c.....  Count how many stns are within the exp accuracy (and multiples)
c
	      if(adiff .le. ea1) num_ea1 = num_ea1 + 1
	      if(adiff .le. ea2) num_ea2 = num_ea2 + 1
	      if(adiff .le. ea3) num_ea3 = num_ea3 + 1
c
	   endif
c
	   write(iunit,905) i, stn(i), interp_ob, ob(i), diff
 905	   format(5x,i3,1x,a3,1x,3f10.2)
c
 500	enddo !i
c
c.... Get the average diff over the obs.
c     
	ave_diff = badflag
	amean = badflag
	if(num .ne. 0) amean = sum / float(num)
	if(num .ne. 0) ave_diff = abs_diff / float(num)
	write(iunit,909) amean, num
 909	format(/,'    Mean difference: ',f10.2,' over ',i4,' stations.')
	write(iunit,910) ave_diff, num
 910	format(' Average difference: ',f10.2,' over ',i4,' stations.')
	write(iunit,920) diff_mx, stn_mx
 920	format(' Maximum difference of ',f10.2,' at ',a3)
	write(iunit,925) diff_mn, stn_mn
 925	format(' Minimum difference of ',f10.2,' at ',a3)
	write(iunit, 930)
 930	format(' ')
c
	write(iunit, 950) field_ea
 950	format(' Number of obs within multiples of exp acc of ',f8.2)
	percent = -1.
	anum = float(num)
	if(num .ne. 0) percent = (float(num_ea1) / anum) * 100.
	write(iunit, 952) num_ea1, num, percent
 952	format(10x,'1x exp accuracy: ',i5,' of ',i5,' (',f5.1,'%)')
	if(num .ne. 0) percent = (float(num_ea2) / anum) * 100.
	write(iunit, 953) num_ea2, num, percent
 953	format(10x,'2x exp accuracy: ',i5,' of ',i5,' (',f5.1,'%)')	
	if(num .ne. 0) percent = (float(num_ea3) / anum) * 100.
	write(iunit, 954) num_ea3, num, percent
 954	format(10x,'3x exp accuracy: ',i5,' of ',i5,' (',f5.1,'%)')
	write(iunit, 930)
	write(iunit, 930)
	write(iunit, 931)
 931	format(1x,'===============================================')
c
	return
	end
c
c
	subroutine bkgwts(lat,lon,topo,numsfc,lat_s,lon_s,elev_s,
     &                    rii,rjj,wt,ni,nj,mxstn)
c
c***************************************************************************
c
c	This routine finds the distance from each gridpoint to each station
c	to find the station density and set up the assmilation weight array.
c
c	Changes:
c 	P.A. Stamus	12-13-96  Original (from build_sfc_static)
c                       08-25-97  Changes for dynamic LAPS
c
c       Note:  'topo' (the LAPS grid elevations) is passed in since we may
c              want to have an elevation-based weight in the future.
c***************************************************************************
c
c
c..... Arrays for the OBS file input data
c
	real*4 lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
c
c.....	Grids for the outputs, weights, and stuff 
c
	real*4 wt(ni,nj), gpmean(ni,nj)
c
c..... LAPS Lat/lon grids.
c
	real*4 lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real*4 rii(mxstn), rjj(mxstn)
c
c
c.....  Set up constants.  Open file for background info.
c
	print *,' '
	print *,' Calculating Background weights:'
c
	imax = ni
	jmax = nj
c
c.....	range of values for the background weights
c
cc	small = .1
cc	alarge = .5
	small = 100.
	alarge = 500.
c
c.....	Find the distance from each grid point to each station
c
!	write(9,899)
!899	format(/,10x,'Within 2',3x,'Within 5',3x,'Within 10',3x,'Within 20',
!     &         3x,' Mean Distance')
	dist_max = -1.e30
	dist_min = 1.e30
	do j=1,jmax
	do i=1,imax
	  num50 = 0
	  num40 = 0
	  num30 = 0
	  num20 = 0
	  num10 = 0
	  num5 = 0
	  num2 = 0
	  amean = 0.
	  numsta = 0
	do ista=1,numsfc
	  if(rii(ista).le.0. .or. rjj(ista).le.0.) go to 12
	  if(rii(ista).gt.ni .or. rjj(ista).gt.nj) go to 12
	  distx = float(i) - rii(ista)
	  disty = float(j) - rjj(ista)
	  dist = sqrt(distx*distx + disty*disty)
c
cc	  if(dist .le. 2.) num2 = num2 + 1
cc	  if(dist .le. 5.) num5 = num5 + 1
cc	  if(dist .le. 10.) num10 = num10 + 1
cc	  if(dist .le. 20.) num20 = num20 + 1
cc	  if(dist .le. 30.) num30 = num30 + 1
cc	  if(dist .le. 40.) num40 = num40 + 1
cc	  if(dist .le. 50.) num50 = num50 + 1
c
	  amean = amean + dist
	  numsta = numsta + 1
 12    enddo !ista
c
	gpmean(i,j) = amean / float(numsta)
	if(gpmean(i,j) .lt. dist_min) then
	  dist_min = gpmean(i,j)
	  min_i = i
	  min_j = j
	endif
	if(gpmean(i,j) .gt. dist_max) then
	  dist_max = gpmean(i,j)
	  max_i = i
	  max_j = j
	endif
c
!	write(9,900) i,j,num2,num5,num10,num20,gpmean(i,j)
!	write(9,900) i,j,num20,num30,num40,num50,gpmean(i,j)
!900	format(i3,i3,':',3x,i4,7x,i4,7x,i4,8x,i4,5x,f12.1)
c
        enddo !i
        enddo !j
c
	write(6,905)
905	format(/,' Distance:')
	write(6,902) 'MAX',dist_max,max_i,max_j
	write(6,902) 'MIN',dist_min,min_i,min_j
902	format(3x,'The ',a3,' of ',f9.1,' occured at gridpoint ',2i4)
c
c.....	Now calculate the weight array by scaling the mean distances
c.....	so that the larger the mean, the larger the weight.
c
	r = dist_min / dist_max
	con =  (small - alarge * r) / (1. - r)
	x = alarge - con
	wt_min = 1.e30
	wt_max = -1.e30
c
	do j=1,jmax
	do i=1,imax
	  wt(i,j) = (gpmean(i,j) / dist_max) * x + con
	  if(wt(i,j) .lt. wt_min) then
	    wt_min = wt(i,j)
	    min_i = i
	    min_j = j
	  endif
	  if(wt(i,j) .gt. wt_max) then
	    wt_max = wt(i,j)
	    max_i = i
	    max_j = j
	  endif
        enddo !i
        enddo !j
c
	write(6,906)
906	format(/,' Weights:')
	write(6,902) 'MAX',wt_max,max_i,max_j
	write(6,902) 'MIN',wt_min,min_i,min_j
c
	print *,' '
	print *,' Normal completion of BKGWTS.'
c
	return
	end
c
c
      subroutine aplot(field, ni, nj)
c
c======================================================================
c
c     Routine to do a quick, scaled ASCII plot of a 2-d field.
c     Based on a version by D. Birkenheuer, NOAA/FSL
c
c     The variables are:
c          field      the field to plot
c          ni, nj     the dimensions of the field
c          amx, amn   the max/min of the values of the field
c
c     Original: 03-06-98  P. Stamus, NOAA/FSL
c
c======================================================================
c
      real*4 field(ni,nj)
c
      character line(ni)*1
c
c
c.....  Start by clearing out the line.
c
      do i=1,ni
         line(i) = ' '
      enddo !i
c
c.....  Figure out the max and min in the field.
c
      amax = -9.e20
      amin =  9.e20
      do j=1,nj
      do i=1,ni
         if(field(i,j) .gt. amax) then
            amax = field(i,j)
            i_max = i
            j_max = j
         endif
         if(field(i,j) .lt. amin) then
            amin = field(i,j)
            i_min = i
            j_min = j
         endif
      enddo !i
      enddo !j
c
      print *,' '
      write(6,901) amax, i_max, j_max, amin, i_min, j_min
 901  format(1x,' Field Max: ',f15.6,' at ',i4,',',i4,
     &       '   Field Min: ',f15.6,' at ',i4,',',i4)
      print *,' '
c        
c.....  Now at each point along the row, scale the value and put
c.....  the proper character in the line at that point.  Then write
c.....  the line and move to the next one down (remember to go from
c.....  the top of the domain to the bottom).
c
      do j=nj,1,-1
      do i=1,ni
c
	 ifld = nint( field(i,j) )
	 fld = float( ifld )
	 ich = ifix( amod(fld, 10.) )
c
	 if(ich .eq. 0) then
	    ich = ifix( amod(fld, 100.) )
	    ich = ich / 10
	    write(line(i), 905) ich
 905        format(i1)
	 elseif(ich.ge.1 .and. ich.lt.5) then
	    line(i) = '.'
	 elseif(ich.ge.5 .and. ich.le.9) then
	    line(i) = ':'
	 else
	    print *,' Bad value in plot routine'
	 endif
c
      enddo !i
c
	 write(6,*) (line(k),k=1,ni)
c
      enddo !j
c
c.....  That's it...lets go home.
c
      print *,' '
      print *,' '
      return
      end
c
c
       subroutine check_field_2d(x,ni,nj,fill_val,istatus)
c
c========================================================================
c
c     Routine to check a 2-d field for NaN's, out of range values, and
c     other bad stuff.  The istatus flag returns what we found to the
c     calling routine:
c
c               istatus =  1  Field ok.
c                          0  Missing field (all 'fill_val')
c                         -1  Found NaN's.
c
c     Some stats are calculated and printed in the log file:
c
c            Max, Min  Absolute Average  Std. Deviation  Mean  Range  
c
c     Original: 06-10-98  (from "stats.f").  P.A. Stamus, NOAA/FSL
c     Changes:  
c
c========================================================================
c
      real*4 x(ni,nj)
c
      istatus = 1
c
c..... Start by zeroing some counters.
c
      amean = 0.
      ave = 0.
      sum = 0.
      sum_a = 0.
      sum_v = 0.
      range = 0.
      st_dev = 0.
      var = 0.
      amax = -1.e25
      amin =  1.e25
      pts = float(ni * nj)
c
c..... Check the field for NaN's
c
      call checknan_2d(x,ni,nj,nan_flag)
      if(nan_flag .ne. 1) then
         print *,'   *** ERROR. NaNs found in field. ***'
         istatus = -1
         return
      endif
c
c..... Check for an empty field.
c
      do j=1,nj
      do i=1,ni
         if(x(i,j) .ne. fill_val) go to 100
      enddo !i
      enddo !j
      print *,'   ** WARNING. Empty field. **'
      istatus = 0
      return
c
 100  continue
c
c..... Calculate the means and range.
c
      do j=1,nj
      do i=1,ni
         if(x(i,j) .lt. amin) then
            amin = x(i,j)
            imin = i
            jmin = j
         endif
         if(x(i,j) .gt. amax) then
            amax = x(i,j)
            imax = i
            jmax = j
         endif
         sum = sum + x(i,j)
         sum_a = sum_a + abs(x(i,j))
      enddo !i
      enddo !j
c
      amean = sum / pts
      ave = sum_a / pts
      range = amax - amin
c
c..... Now calculate the variance, stdev, etc.
c
      do j=1,nj
      do i=1,ni
         dif = x(i,j) - amean
         dif2 = dif * dif
         sum_v = sum_v + dif2
      enddo !i
      enddo !j
c
      var = 0.0
      st_dev = -99.9
      if(pts .ne. 1) var = sum_v / (pts - 1)
      if (var .gt. 0.0) st_dev = sqrt( var )
c
c..... Write out some stat info.
c
      write(6,900) amax,imax,jmax,amin,imin,jmin,range
 900  format(5x,'Max: ',g12.4,' at ',2i4,'  Min: ',
     &         g12.4,' at ',2i4,'  Range: ',g12.4)
c
      write(6,910) amean, ave, st_dev
 910  format(5x,'Mean:',g12.4,2x,'AbsAve:',g12.4,2x,'StDev:',g12.4)
c
c
c.... That's it.  Let's go home.
c
      return
      end
c
c
      subroutine checknan_2d(x,ni,nj,nan_flag)
c
c     Routine to check a real 2-d array for NaN's.
c
      real*4 x(ni,nj)
c
      nan_flag = 1
c
      do j=1,nj
      do i=1,ni
         if( nan( x(i,j) ) .eq. 1) then
            print *,' ** ERROR. Found a NaN at ', i, j
            nan_flag = -1
            return
         endif
      enddo !i
      enddo !j
c
      return
      end
c
c
      subroutine checknan_3d(x,ni,nj,nk,nan_flag)
c
c     Routine to check a real 3-d array for NaN's.
c
      real*4 x(ni,nj,nk)
c
      nan_flag = 1
c
      do k=1,nk
      do j=1,nj
      do i=1,ni
         if( nan( x(i,j,k) ) .eq. 1) then
            print *,' ** ERROR. Found a NaN at ', i, j, k
            nan_flag = -1
            return
         endif
      enddo !i
      enddo !j
      enddo !k
c
      return
      end
c
c
	subroutine get_background_sfc(i4time_in,var_in,bkg_ext,
     &                bkg_time,bkg_field,laps_cycle_time,ni,nj,
     &                istatus)

c
c
c*****************************************************************************
c
c       Routine to get a surface background field for a given variable.  
c       This routine checks several different files for an appropiate field
c       to use, and checks the field for NaN's and other problems before
c       returning.  In addition to returning the appropiate background
c       field, the routine also returns the extension and time of the
c       background file.
c
c	History:
c		P. Stamus, NOAA/FSL 
c                   06-16-98  Original version.
c                   09-10-98  Change SFM filename read to be like LGB.
c
c
c       Notes:
c               Inputs:    i4time_in  - time of analysis    (int)
c                             var_in  - requested variable  (ch*4)
c                    laps_cycle_time  - analysis interval   (int)
c                             ni, nj  - grid dimensions     (int)
c
c               Outputs:    bkg_field - 2d background field    (r*4)
c                             bkg_ext - ext of bkg file used  (ch*31)
c                            bkg_time - time of bkg file used  (int)
c                          bkg_status - Status:
c                                          1 - normal return
c                                          0 - No background files found
c
c
c*****************************************************************************
c
	real*4 bkg_field(ni,nj)
c
	integer*4 i4time_in, lvl_in, bkg_time, ni, nj, bkg_status
c
	character bkg_ext*31, var_in*4, var(3)*3
	character bkg_dir*256, filename*9
	character units*10, lvlc*4, comment*125
c
	integer*4 max_files
	parameter(max_files = 300)
	character fnames(max_files)*80
	character dir*150, ext*31, filespec*255, filename13*13
c
c
c.....	Start here.  
c
	bkg_status = 0
	lvl_in = 0      ! surface level
	fill_val = 1.e37
	i4time_prev = i4time_in - laps_cycle_time
	call constant(bkg_field, fill_val, ni,nj)
c
c.....  Set up var array for the different names used in the different
c.....  files, based on the 'var_in' variable.
c
	if(var_in .eq. 'TEMP') then
	   var(1) = 'T  '  ! SFM variable
	   var(2) = 'TSF'  ! LGB variable
	   var(3) = 'T  '  ! LSX variable
	elseif(var_in .eq. 'SFCP') then
	   var(1) = 'PS '  ! SFM variable
	   var(2) = 'PSF'  ! LGB variable
	   var(3) = 'PS '  ! LSX variable
	elseif(var_in .eq. 'MSLP') then
	   var(1) = 'MSL'  ! SFM variable
	   var(2) = 'SLP'  ! LGB variable
	   var(3) = 'MSL'  ! LSX variable
	elseif(var_in .eq. 'DEWP') then
	   var(1) = 'TD '  ! SFM variable
	   var(2) = 'DSF'  ! LGB variable
	   var(3) = 'TD '  ! LSX variable
	elseif(var_in .eq. 'REDP') then
	   var(1) = 'P  '  ! SFM variable
	   var(2) = 'P  '  ! LGB variable
	   var(3) = 'P  '  ! LSX variable
	elseif(var_in .eq. 'VISB') then
	   var(1) = 'VIS'  ! SFM variable
	   var(2) = 'VIS'  ! LGB variable
	   var(3) = 'VIS'  ! LSX variable
	else
	   print *,
     &       ' Invalid variable request sent to get_background_sfc'
	   return
	endif
c
c.....  Get the background data.  Try for a SFM forecast first.  If not
c.....  available, try the LGB file.  If that's not there either, use a
c.....  previous LAPS analysis.  If nothings available, print a warning.
c
	isfm_bk = 1
	print *,' Trying for SFM background '
c
c	bkg_dir = '../lapsprd/rsf/'
	bkg_ext = 'rsf'
	i4time_bk = i4time_in

	call get_directory('rsf',bkg_dir,len)
	filespec = bkg_dir(1:len) // '*.'// bkg_ext
	call get_file_names(filespec,numfiles,fnames,max_files,istatus)

	i_best_file = 0
	i4_ftime_min = 9999999

	do i=1,numfiles
	   call get_directory_length(fnames(i), lend)
	   call get_time_length(fnames(i), lenf)
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
	   if(i4valid .eq. i4time_in) then
	      i4_fcst_time = i4valid - i4init
	      if(i4_fcst_time .lt. i4_ftime_min) then
		 i4_fcst_time = i4_ftime_min
		 i_best_file = i
	      endif
	   endif
	enddo !i
c
	if(i_best_file .gt. 0) then !found one

	   i = i_best_file
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
c
 110	   call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &        var(1),lvl_in,lvlc,units,comment,bkg_field,istatus)
  	   if(istatus .ne. 1) then
	      print *,' ERROR reading SFM file at ', filename13
	      isfm_bk = 0
	      go to 200
	   endif
c
	else
c
	   print *,' No SFM file with proper valid time.'
	   isfm_bk = 0
	   go to 200

	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found SFM background at ',filename13,
     &          '...checking field.'
	call check_field_2d(bkg_field,ni,nj,fill_val,istatus)
	if(istatus .eq. 1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,
     &    '  Problem with SFM background, check status = ', istatus
	endif
c
c.....  Try LGB.
c
 200	continue
	ilgb_bk = 1
	print *,' Trying for LGB background '
c	bkg_dir = '../lapsprd/lgb/'
	bkg_ext = 'lgb'
	i4time_bk = i4time_in

	call get_directory('lgb',bkg_dir,len)
	filespec = bkg_dir(1:len) // '*.'// bkg_ext
	call get_file_names(filespec,numfiles,fnames,max_files,istatus)

	i_best_file = 0
	i4_ftime_min = 9999999

	do i=1,numfiles
	   call get_directory_length(fnames(i), lend)
	   call get_time_length(fnames(i), lenf)
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
	   if(i4valid .eq. i4time_in) then
	      i4_fcst_time = i4valid - i4init
	      if(i4_fcst_time .lt. i4_ftime_min) then
		 i4_fcst_time = i4_ftime_min
		 i_best_file = i
	      endif
	   endif
	enddo !i
c
	if(i_best_file .gt. 0) then !found one

	   i = i_best_file
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
c
 210	   call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &        var(2),lvl_in,lvlc,units,comment,bkg_field,istatus)
  	   if(istatus .ne. 1) then
	      print *,' ERROR reading LGB file at ', filename13
	      ilgb_bk = 0
	      go to 300
	   endif
c
	else
c
	   print *,' No LGB file with proper valid time.'
	   ilgb_bk = 0
	   go to 300

	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found LGB background at ',filename13,
     &          '...checking field.'
	call check_field_2d(bkg_field,ni,nj,fill_val,istatus)
	if(istatus .eq. 1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,
     &    '  Problem with LGB background, check status = ', istatus
	endif
c
c.....	Try the previous LSX.
c
 300	ilaps_bk = 1
	ilaps_loop = 1
	print *,' Trying for previous LSX background '
c	bkg_dir = '../lapsprd/lsx/'
	call get_directory('lsx',bkg_dir,len)
	bkg_ext = 'lsx'
	i4time_bk = i4time_prev
c
 310	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var(3),lvl_in,lvlc,units,comment,bkg_field,istatus)
	if(istatus .ne. 1) then
	   if(ilaps_loop .gt. 3) then     ! give up
	      print *,'  No LSX background available. '
	      ilaps_bk = 0
	      go to 500
	   else
	      ilaps_loop = ilaps_loop + 1
	      i4time_bk = i4time_bk - laps_cycle_time
	      go to 310
	   endif
	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
        call make_fnam_lp(i4time_bk,filename,istat)  ! make filename	
	print *,'  Found LSX background at ',filename,
     &          '...checking field.'
	call check_field_2d(bkg_field,ni,nj,fill_val,istatus)
	if(istatus .eq. 1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,'  Problem with LSX background, check status = ', istatus
	endif
c
 500	continue
c
c.....	That's about it...let's go home.
c
	return
	end
c
c
	subroutine get_bkgwind_sfc(i4time_in,bkg_ext,bkg_time,bkg_u,bkg_v,
     &                laps_cycle_time,ni,nj,istatus)

c
c
c*****************************************************************************
c
c       Routine to get a surface background wind field.  This routine checks 
c       several different files for an appropiate fields to use, and checks 
c       the fields for NaN's and other problems before returning.  The u and
c       v components must be from the same file at the same time to be vaild.
c       In addition to returning the appropiate background u and v, this 
c       routine also returns the extension and time of the background file.
c
c	History:
c		P. Stamus, NOAA/FSL 
c                   07-10-98  Original version.
c                   09-10-98  Change SFM filename read to be like LGB.
c
c
c       Notes:
c               Inputs:    i4time_in  - time of analysis    (int)
c                    laps_cycle_time  - analysis interval   (int)
c                             ni, nj  - grid dimensions     (int)
c
c               Outputs: bkg_u, bkg_v - 2d background component fields (r*4)
c                             bkg_ext - ext of bkg file used  (ch*31)
c                            bkg_time - time of bkg file used  (int)
c                          bkg_status - Status:
c                                          1 - normal return
c                                          0 - No background files found
c
c
c*****************************************************************************
c
	real*4 bkg_u(ni,nj), bkg_v(ni,nj)
c
	integer*4 i4time_in, lvl_in, bkg_time, ni, nj, bkg_status
c
	character bkg_ext*31, var_u*3, var_v*3
	character bkg_dir*256, filename*9
	character units*10, lvlc*4, comment*125
c
	integer*4 max_files
	parameter(max_files = 300)
	character fnames(max_files)*80
	character dir*150, ext*31, filespec*255, filename13*13
c
c.....	Start here.  
c
	bkg_status = 0
	lvl_in = 0      ! surface level
	fill_val = 1.e37
	i4time_prev = i4time_in - laps_cycle_time
	call constant(bkg_u, fill_val, ni,nj)
	call constant(bkg_v, fill_val, ni,nj)
	var_u = 'U  '
	var_v = 'V  '
	print *,' '
	print *,' Getting background wind data....'
c
c
c.....  Get the background data.  Try for a SFM forecast first.  If not
c.....  available, try the LGB file.  If that's not there either, use a
c.....  previous LAPS analysis.  If nothings available, print a warning.
c
	isfm_bk = 1
	print *,' Trying for SFM background '
c	bkg_dir = '../lapsprd/rsf/'
	bkg_ext = 'rsf'
	i4time_bk = i4time_in

	call get_directory('rsf',bkg_dir,len)
	filespec = bkg_dir(1:len) // '*.'// bkg_ext
	call get_file_names(filespec,numfiles,fnames,max_files,istatus)

	i_best_file = 0
	i4_ftime_min = 9999999

	do i=1,numfiles
	   call get_directory_length(fnames(i), lend)
	   call get_time_length(fnames(i), lenf)
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
	   if(i4valid .eq. i4time_in) then
	      i4_fcst_time = i4valid - i4init
	      if(i4_fcst_time .lt. i4_ftime_min) then
		 i4_fcst_time = i4_ftime_min
		 i_best_file = i
	      endif
	   endif
	enddo !i
c
	if(i_best_file .gt. 0) then !found one

	   i = i_best_file
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
c
 110	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)

  	   if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	      print *,' ERROR reading SFM file at ', filename13
	      isfm_bk = 0
	      go to 200
	   endif
c
	else
c
	   print *,' No SFM file with proper valid time.'
	   isfm_bk = 0
	   go to 200

	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found SFM backgrounds...checking fields.'
	call check_field_2d(bkg_u,ni,nj,fill_val,istatus_u)
	call check_field_2d(bkg_v,ni,nj,fill_val,istatus_v)
	if(istatus_u.eq.1 .and. istatus_v.eq.1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,
     &   '   Problem with SFM backgrounds, check status = ',
     &   istatus_u, istatus_v
	endif
c
c.....  Try LGB.
c
 200	continue
	var_u = 'USF'
	var_v = 'VSF'
	ilgb_bk = 1
	print *,' Trying for LGB background '
c	bkg_dir = '../lapsprd/lgb/'
	bkg_ext = 'lgb'
	i4time_bk = i4time_in

	call get_directory('lgb',bkg_dir,len)
	filespec = bkg_dir(1:len) // '*.'// bkg_ext
	call get_file_names(filespec,numfiles,fnames,max_files,istatus)

	i_best_file = 0
	i4_ftime_min = 9999999

	do i=1,numfiles
	   call get_directory_length(fnames(i), lend)
	   call get_time_length(fnames(i), lenf)
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
	   if(i4valid .eq. i4time_in) then
	      i4_fcst_time = i4valid - i4init
	      if(i4_fcst_time .lt. i4_ftime_min) then
		 i4_fcst_time = i4_ftime_min
		 i_best_file = i
	      endif
	   endif
	enddo !i
c
	if(i_best_file .gt. 0) then !found one

	   i = i_best_file
	   filename13 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename13,i4init,i4valid,i4fn)
c
 210	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)

  	   if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	      print *,' ERROR reading LGB file at ', filename13
	      ilgb_bk = 0
	      go to 300
	   endif
c
	else
c
	   print *,' No LGB file with proper valid time.'
	   ilgb_bk = 0
	   go to 300

	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found LGB backgrounds at ',filename13,
     &          '...checking fields.'
	call check_field_2d(bkg_u,ni,nj,fill_val,istatus_u)
	call check_field_2d(bkg_v,ni,nj,fill_val,istatus_v)
	if(istatus_u.eq.1 .and. istatus_v.eq.1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,
     &    '  Problem with LGB backgrounds, check status = ', 
     &    istatus_u, istatus_v
	endif
c
c.....	Try LWM (sfc interpolated from 3d)
c
 300	ilwm_bk = 1
	ilwm_loop = 1
	print *,' Trying for LWM wind background '
c	bkg_dir = '../lapsprd/lwm/'
	call get_directory('lwm',bkg_dir,len)
	bkg_ext = 'lwm'
	i4time_bk = i4time_prev
	var_u = 'SU '
	var_v = 'SV '
c
 310	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)
c
	if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	   if(ilwm_loop .gt. 3) then     ! give up
	      print *,'  No LWM background available. '
	      ilwm_bk = 0
	      go to 400
	   else
	      ilwm_loop = ilwm_loop + 1
	      i4time_bk = i4time_bk - laps_cycle_time
	      go to 310
	   endif
	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found LWM backgrounds...checking fields.'
	call check_field_2d(bkg_u,ni,nj,fill_val,istatus_u)
	call check_field_2d(bkg_v,ni,nj,fill_val,istatus_v)
	if(istatus_u.eq.1 .and. istatus_v.eq.1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,'  Problem with LWM background, check status = ', 
     &            istatus_u, istatus_v
	endif
c
c.....	Try the previous LSX.
c
 400	ilaps_bk = 1
	ilaps_loop = 1
	print *,' Trying for previous LSX background '
c	bkg_dir = '../lapsprd/lsx/'
	call get_directory('lsx',bkg_dir,len)
	bkg_ext = 'lsx'
	i4time_bk = i4time_prev
	var_u = 'U  '
	var_v = 'V  '
c
 410	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
 	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)
c
	if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	   if(ilaps_loop .gt. 3) then     ! give up
	      print *,'  No LSX background available. '
	      ilaps_bk = 0
	      go to 500
	   else
	      ilaps_loop = ilaps_loop + 1
	      i4time_bk = i4time_bk - laps_cycle_time
	      go to 410
	   endif
	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found LSX backgrounds...checking fields.'
	call check_field_2d(bkg_u,ni,nj,fill_val,istatus_u)
	call check_field_2d(bkg_v,ni,nj,fill_val,istatus_v)
	if(istatus_u.eq.1 .and. istatus_v.eq.1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,'  Problem with LSX background, check status = ', 
     &             istatus_u, istatus_v
	endif
c
 500	continue
c
c.....	That's about it...let's go home.
c
	return
	end

