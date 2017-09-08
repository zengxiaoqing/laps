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
	implicit none
	integer num_sfc
	real a_t, b_t, a_td, b_td, hbar, badflag
	real elev(num_sfc), t(num_sfc), td(num_sfc), elev_i
	real cnt, cntd, sumht, sumh, sumt, sumh2, sumt2, sumtd, sumhtd
	integer i
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

          if(abs(elev(i)) .lt. 1e-10)then ! prevents underflow if elev is zeros
              elev_i = 0.
          else
              elev_i = elev(i)
          endif

	  sumht = (elev_i * t(i)) + sumht
	  sumh = elev_i + sumh
	  sumh2 = (elev_i * elev_i) + sumh2
	  sumt = t(i) + sumt
	  cnt = cnt + 1.
c
	  if(td(i) .le. badflag) go to 10
	  sumtd = td(i) + sumtd
	  sumhtd = (elev_i * td(i)) + sumhtd
	  cntd = cntd + 1.
10	continue
c
        if(cnt .gt. 0.0)then
            hbar = sumh / cnt
        else
            hbar = 0.
        endif

        if(hbar .eq. 0.)then
            write(6,*)' WARNING in mean_lapse: station elevations'
     1               ,' and/or hbar = 0.'
            write(6,*)' Skipping lapse rate regression'
            goto990
        endif
c
        if(cntd .eq. 0.)then
            write(6,*)' WARNING in mean_lapse: cntd = 0.!'
     1               ,' dewpoints all badflag/missing'
            write(6,*)' Skipping lapse rate regression'
            goto990
        endif
c
	b_t = (cnt*sumht - sumh*sumt) / (cnt*sumh2 - sumh*sumh)
	a_t = (sumt - b_t * sumh) / cnt
c
	b_td = (cntd*sumhtd - sumh*sumtd) / (cntd*sumh2 - sumh*sumh)
	a_td = (sumtd - b_td * sumh) / cntd
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
990     write(6,*)' mean_lapse: Using std lapse rate...'
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
        subroutine mean_pressure(p_s,n_sfc,p_bk,ni,nj,badflag,pbar)
c
c*******************************************************************************
c
c       Routine to calculate the mean surface or reduced pressure.
c       Based on 'mean_press' routine.
c
c       Changes:
c               P.A.Stamus    11-23-99       Original
c
c       Inputs/Outputs:
c
c          Variable     Var Type     I/O     Description
c         ----------   ----------   -----   -------------
c          p_s             RA         I      Pressure observations.
c          n_sfc           I          I      Number of surface observations.
c          p_bk            RA         I      Pressure background array.
c          ni, nj          I          I      Domain dimensions
c          badflag         R          I      Bad flag value.
c          pbar            R          O      Mean pressure of all obs or bkg.
c
c       User Notes:
c
c       1.  Units are not changed in this routine.
c
c*******************************************************************************
c
        real p_s(n_sfc), p_bk(ni,nj)
c
c.....  Find the mean pressure of the observations.
c
        sump = 0.
        cntp = 0.
c
        do i=1,n_sfc
          if(p_s(i) .le. 0.) go to 1
          sump = sump + p_s(i)
          cntp = cntp + 1.
 1	  continue
	enddo !i
c
	if(cntp .eq. 0.) then
	   print *,
     &     '  WARNING. No pressure observations. Trying background.'
	else
	   pbar = sump / cntp
	   return
	endif
c
c.....  If no pressure obs, try the background.  If bkg no good,
c.....  set the pbar to badflag.
c
	sump = 0.
	cntp = 0.
c
	do j=1,nj
	do i=1,ni
	   if(p_bk(i,j) .le. 0.) go to 2
	   sump = sump + p_bk(i,j)
	   cntp = cntp + 1.
 2	   continue
	enddo !i
	enddo !j
c
	if(cntp .eq. 0.) then
	   print *,
     &     '  WARNING. No pressure background either.'
	   pbar = badflag
	else
	   pbar = sump / cntp
	endif
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
	implicit none
	real alt, elev, badflag
	real term1, term2
	real const, const1, alapse, tstd, amslp
	parameter(const=0.190284,const1=1./const, alapse=0.0065,
     +            tstd=288.15,amslp=1013.25)
cc	data const/0.190284/, const1/5.2553026/		!const1 = 1/const
cc	data alapse/0.0065/, tstd/288.15/, amslp/1013.25/
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
	subroutine find_ij(lat_s,lon_s,lat,lon,numsta,mxsta,
     &                     ni,nj,ii,jj,rii,rjj)
c
c======================================================================
c
c       Routine to find the i,j locations for each station.  Do not "round"
c       the ii,jj's "up"...straight truncation puts the ob at the proper
c       grid point on the major grid.
c
c       Orginal:  P. Stamus NOAA/FSL c.1990
c       Changes:
c
c       Changed by Steve Albers in 2003 as it seems to me that the major
c       (non-staggered) grid can best be emplaced by each thermo ob by 
c       allowing a rounding up to occur.
c
c======================================================================
c
	real lat_s(mxsta), lon_s(mxsta)
        real lat(ni,nj), lon(ni,nj)
        real rii(mxsta), rjj(mxsta)
	integer ii(mxsta), jj(mxsta)
c
	do ista=1,numsta
          call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat,lon,
     &       ni,nj,rii(ista),rjj(ista),istatus)
	  ii(ista) = nint(rii(ista))
	  jj(ista) = nint(rjj(ista))
	enddo !ista
c
	return
	end
c
c
c
	subroutine extract(a,imax,jmax,i,j,ix,jy)
c
c.....	Routine designed to zero out a grouping of points about a 
c.....  named point i,j.  All points from i-ix to i+ix, and 
c.....	j-jy to j+jy will be zeroed.  This is aimed at processing 
c.....	Band 8 temperatures to remove cloud edges.
c
	implicit none
	integer imax, jmax, i,j,ix,jy
	real a(imax,jmax)
	integer jmx, imx, jmn, imn, ii, jj
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

	implicit none
	real flag, angle
	parameter	(flag=1e37)

	real		dd,ff,ucomp,vcomp
	integer	status

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
      implicit none

      integer ni,nj
      real x(ni,nj)
      real amean, ave, sum, sum_a, sum_v, sum_sk, sum_kr, range, st_dev
      real var, amax, amin, pts, dif, dif2, z_max, z_min, coef_sk
      real coef_kr
      integer i,j, imax, imin, jmax, jmin
      
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

      if (var .gt. 0. .and. sum_v .gt. 0. .and. pts .gt. 0.) then
        st_dev = sqrt( var )
        z_max = (amax - amean) / st_dev
        z_min = (amin - amean) / st_dev
        coef_sk = (sum_sk / pts) / ( (sum_v / pts) ** 1.5 )
        coef_kr = (sum_kr / pts) / ( (sum_v / pts) ** 2   )
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
!       write(6,920) z_max, z_min, range
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
c                09-30-98  Housekeeping.
c                04-12-99  More housekeeping.
c
c*************************************************************************
c
	integer imax, jmax, i, j
	real t(imax,jmax), tb8(imax,jmax)
        real lat(imax,jmax), lon(imax,jmax), topo(imax,jmax)
c
        real cvr_snow(imax,jmax), t_gnd_k(imax,jmax)  !work arrays
	real t_est(imax,jmax), dtb8(imax,jmax)        !work arrays
c
c
	call zero(t_est,imax,jmax)
	call constant(dtb8,-99.,imax,jmax)

	do j=1,jmax
	do i=1,imax
	   terr = topo(i,j)
	   t_est(i,j) = t(i,j)
	enddo !i
	enddo !j

        call get_ground_temperature(i4time,laps_cycle_time
     &                    ,imax,jmax,lat,lon,r_missing_data
     &                             ,cvr_snow,t_est,t_gnd_k)


!       Call the ground temperature routine

	do j=1,jmax
	do i=1,imax

	   if(tb8(i,j).ne.smsng .and. tb8(i,j).ne.0.) then
c             t_compare = t_est(i,j)
	      t_compare = t_gnd_k(i,j)
	      dtb8(i,j) = tb8(i,j) - t_compare
	   endif

	enddo !i
	enddo !j
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
	   if(dtb8(i,j) .lt. -thresh1) then ! probably clds
	      
c.....	set the point and surrounding points to zero...clouds.
	      tb8(i,j) = 0.	! Set to zero as a cloud flag
	      call extract(tb8,imax,jmax,i,j,7,7)
	      icnt = icnt + 1
          endif
	enddo	! on i
	enddo	! on j
c
	write(6,951) icnt, thresh1
951	format(1x,i6,' points removed for thresh = ',f8.1)
c
	return
	end
c
c
	subroutine lp_fire_danger (ni,nj,lp_10krh,lp_10kt,lp_10kws,
     &                             soil_moist,snow_cover,topo,ldf,
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
	integer ni,nj
	real lp_10krh(ni,nj), lp_10kt(ni,nj), lp_10kws(ni,nj)
	real snow_cover(ni,nj), soil_moist(ni,nj), topo(ni,nj)
	real lp_fire_index(ni,nj), ldf(ni,nj)
c
	integer i_lp,	j_lp, i_rh, i_temp, i_wspeed, i_soil,
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

c.....  Apply land fraction as a mask
            lp_fire_index(i_lp,j_lp) 
     1    = lp_fire_index(i_lp,j_lp) * ldf(i_lp,j_lp)

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
      subroutine heat_index(t,rh,hi,ni,nj,r_missing_data)
c
c====================================================================
c
c     Routine to calculate a heat index.  Based on a formula 
c     by Lans Rothfusz, NWS.  Seems to provide valid HI numbers
c     for temperatures above 75 deg F.
c
c     Original:  07-18-95  P. Stamus, NOAA/FSL
c     Changes:  P. Stamus  08-25-97  Return r_missing_data if Temp < 75F
c                                    Change units returned to K.
c                          01-20-98  T in as deg K.
c
c     Notes:
c
c       1.  Inputs:
c                    rh = Relative Humidity (0 to 100 %)
c                    t  = Temperature (deg K)
c	             ni, nj  = Grid dimensions
c                    r_missing_data = Bad flag value
c
c           Output:
c                    hi = Heat Index (deg K)
c
c       2.  If the temperature is below 75 deg F, no heat index is
c           calculated and the point is set to "r_missing_data".
c
c====================================================================
c
      integer ni,nj
      real t(ni,nj), rh(ni,nj), hi(ni,nj)
c
      do j=1,nj
      do i=1,ni
c
	 temp = ( 1.8 * (t(i,j) - 273.15) ) + 32.         ! K to F
c
	 if(temp .lt. 75.) then
	    hi(i,j) = r_missing_data

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
	subroutine bkgwts(lat,lon,topo,numsfc,lat_s,lon_s,elev_s,     ! I
     &                    rii,rjj,                                    ! I
     &                    wt,                                         ! O
     &                    ni,nj,mxstn,                                ! I
     &                    istatus)                                    ! O
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
        integer mxstn, ni, nj
	real lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
c
c.....	Grids for the outputs, weights, and stuff 
c
	real wt(ni,nj), gpmean(ni,nj)
c
c..... LAPS Lat/lon grids.
c
	real lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real rii(mxstn), rjj(mxstn)
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
	dist_max = -1.e15
	dist_min = +1.e15
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

          if(.false.)then
	      do ista=1,numsfc
                  riii = rii(ista)
                  rjjj = rjj(ista)

  	          if(riii.le.0. .or. rjjj.le.0.) go to 12
	          if(riii.gt.ni .or. rjjj.gt.nj) go to 12

                  distx = float(i) - riii
                  disty = float(j) - rjjj
	          dist = sqrt(distx*distx + disty*disty)
c
                  amean = amean + dist
	          numsta = numsta + 1
 12           enddo !ista

          else ! approximate as final output is hopefully unaffected
              numsta = 1
              dist = ni/2
              amean = dist

          endif
c
          if(numsta .eq. 0)then
            write(6,*)' Aborting subroutine bkgwts, no obs in domain...'
            istatus = 0
            return
          else
   	    gpmean(i,j) = amean / float(numsta)
          endif

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
!	  write(9,900) i,j,num2,num5,num10,num20,gpmean(i,j)
!	  write(9,900) i,j,num20,num30,num40,num50,gpmean(i,j)
!900	  format(i3,i3,':',3x,i4,7x,i4,7x,i4,8x,i4,5x,f12.1)
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

        if(r .ne. 1.)then
	    con =  (small - alarge * r) / (1. - r)
        endif

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
        istatus = 1
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
      integer ni,nj
      real field(ni,nj)
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
c     Changes:  09-24-98  P.A. Stamus, NOAA/FSL
c                  Added zero to missing field check.
c               10-02-98  P.A. Stamus, NOAA/FSL
c                  Added check for fill_val to max/min/ave calcs.
c
c
c========================================================================
c
      integer ni,nj
      real x(ni,nj)
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
      pts = 0
c
c..... Check the field for NaN's
c
      call check_nan2(x,ni,nj,nan_flag)
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
         if(x(i,j).ne.fill_val .and. x(i,j).ne.0.) go to 100
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
	 if(x(i,j) .ne. fill_val) then
	    pts = pts + 1
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
	 endif
      enddo !i
      enddo !j
c
      amean = sum / pts
      ave = sum_a / pts
      range = amax - amin
c
c..... Now calculate the variance, stdev, etc.
c
      pts = 0
      do j=1,nj
      do i=1,ni
	 if(x(i,j) .ne. fill_val) then
	    pts = pts + 1
	    dif = x(i,j) - amean
	    dif2 = dif * dif
	    sum_v = sum_v + dif2
	 endif
      enddo !i
      enddo !j
c
      var = 0.0
      st_dev = -99.9
      if(pts .gt. 1) var = sum_v / (pts - 1)
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
	subroutine get_background_sfc(i4time_in,var_in,bkg_ext,
     &                bkg_time,bkg_field,laps_cycle_time,ni,nj,
     &                bkg_status)
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
c                   09-30-98  Housekeeping.
c                   12-02-98  Change 'istatus' to 'bkg_status' in call list.
c                   01-28-99  Skip LGB until can figure out how to deal 
c                               with bias in its fields.
c                   06-11-99  Turn LGB back on.
c                   07-25-99  Set LGB MSL variable to MSL so won't use for now.
c                   09-16-99  Try LGB MSL again.
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
        include 'lapsparms.for' ! MAX_BACKGROUND_FILES

	integer ni,nj
	real bkg_field(ni,nj)
	real bkg_field_dum(ni,nj)
c
	integer i4time_in, lvl_in, bkg_time, bkg_status
c
	character bkg_ext*31, var_in*4, var(3)*3
	character bkg_dir*256, filename*9
	character units*10, lvlc*4, comment*125
c
	integer max_files
	parameter(max_files = MAX_BACKGROUND_FILES)
	character fnames(max_files)*256
	character filespec*255
c
c
c.....	Start here.  
c
	call tagit('get_background_sfc', 19990725)
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
cc	   var(2) = 'MSL'  ! LGB variable
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


        if(.true.)then ! new way
           call get_modelfg_2d(i4time_in,var(2),ni,nj,bkg_field,istatus)
           if(istatus .ne. 1)then
	       print *,' No LGB/FSF file with proper valid time.'
               bkg_status = 0
	       go to 300
           endif

           i4time_bk = i4time_in

        endif

c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'...checking field.'
	call check_field_2d(bkg_field,ni,nj,fill_val,istatus)

	if(istatus .eq. 1) then
	   bkg_status = 1
	   bkg_time = i4time_bk
	   return
	else
	   print *,
     &    '  Problem with LGB/RSF background, check status = ', istatus       
	endif
c
c.....	Try the previous LSX.
c
 300	if(.true.)return ! Using the previous LSX analysis as a background can
                         ! promote gradient overshoots so we'll turn this off
                         ! for now. One future option would be a climo field or
                         ! to analyze only on larger scales with this type of
                         ! cycling. The latter could be accomplished by setting
                         ! an artifically high value to the rms iteration 
                         ! thresholds.

        ilaps_bk = 1
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
     &                laps_cycle_time,ni,nj,bkg_status)
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
c                   09-30-98  Housekeeping.
c                   12-02-98  Change 'istatus' to 'bkg_status' in call list.
c                   01-28-99  Skip LGB until can figure out how to deal 
c                               with bias in its fields.
c                   06-11-99  Turn LGB back on.
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
        include 'lapsparms.for' ! MAX_BACKGROUND_FILES

	integer ni,nj
	real bkg_u(ni,nj), bkg_v(ni,nj)
c
	integer i4time_in, lvl_in, bkg_time, bkg_status
c
	character bkg_ext*31, var_u*3, var_v*3
	character bkg_dir*256
	character units*10, lvlc*4, comment*125
c
	integer max_files
	parameter(max_files = MAX_BACKGROUND_FILES)
	character fnames(max_files)*256
	character filespec*255, filename14*14
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
	call tagit('get_bkgwind_sfc', 19990611)
c
c
c.....  Get the background data.  Try for LWM file first. Next, try for an 
c.....  SFM forecast.  If not available, try the LGB file.  If that's not there
c.....  either, use a previous LAPS analysis.  If nothings available, print a 
c.....  warning.
c
c.....	Try LWM (sfc interpolated from 3d)
c
	ilwm_loop = 1
	print *,' Trying for LWM wind background '
c	bkg_dir = '../lapsprd/lwm/'
	call get_directory('lwm',bkg_dir,len)
	bkg_ext = 'lwm'
	i4time_bk = i4time_in
	var_u = 'SU '
	var_v = 'SV '
c
 310	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
	call read_laps_data(i4time_bk,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)
c
	if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	   if(ilwm_loop .gt. 0) then     ! give up
	      print *,'  No LWM background available. '
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
 400	continue
c
c.....  Try SFM
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
	   filename14 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename14,i4init,i4valid,i4fn)
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
	   filename14 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename14,i4init,i4valid,i4fn)
c
 110	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)

  	   if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	      print *,' ERROR reading SFM file at ', filename14
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
	print *,'  Found SFM background at ',filename14,
     &          '...checking field.'
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
 200	continue
c
c.....  Try LGB.
c
	var_u = 'USF'
	var_v = 'VSF'
	imodel_bk = 1
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
	   filename14 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename14,i4init,i4valid,i4fn)
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
	   filename14 = fnames(i)(lend+1:lenf)
	   call get_fcst_times(filename14,i4init,i4valid,i4fn)
c
 210	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_u,lvl_in,lvlc,units,comment,bkg_u,istatus_u)
	call read_laps(i4init,i4valid,bkg_dir,bkg_ext,ni,nj,1,1,
     &   var_v,lvl_in,lvlc,units,comment,bkg_v,istatus_v)

  	   if(istatus_u.ne.1 .or. istatus_v.ne.1) then
	      print *,' ERROR reading LGB file at ', filename14
	      imodel_bk = 0
	      go to 300
	   endif
c
	else
c
	   print *,' No LGB file with proper valid time.'
	   imodel_bk = 0
	   go to 300

	endif
c
c.....  Check the field for NaN's and other bad stuff.
c
	print *,'  Found LGB backgrounds at ',filename14,
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

 300	if(.true.)return ! Using the previous LSX analysis as a background can
                         ! promote gradient overshoots so we'll turn this off
                         ! for now. One future option would be a climo field or
                         ! to analyze only on larger scales with this type of
                         ! cycling. The latter could be accomplished by setting
                         ! an artifically high value to the rms iteration 
                         ! thresholds.
c
c.....	Try the previous LSX.
c
        ilaps_bk = 1
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
c
c
	subroutine tagit(name, code)
c
c*****************************************************************************
c
c       Routine to print a tracking code for the 'name' routine into the log.
c
c       Original:   P. Stamus, NOAA/FSL   04 Aug 1999
c
c       Notes:
c
c*****************************************************************************
c
	integer len
	integer code
	character name*(*)
c
	call s_len(name, len)
c
	write(6,*) ' Welcome to ',name(1:len),';  tag: ',code
c
	return
	end

        subroutine pstn_anal(back_mp,back_sp,prin_bk,prin,ni,nj
     1                      ,stnp_bk,stnp)

!       This routine solves for the 'stnp' analysis by assuming that ratio
!       stnp/stnp_bk is equal to prin/prin_bk. This provides an indirect
!       way to get the observations to modify the STNP background to
!       come up with a STNP analysis. We are thus feeding off of any PRIN obs
!       that contributed to the PRIN analysis.

        integer back_mp,back_sp

!       PRIN is input pressure array, either REDP field (mb)
        real prin_bk(ni,nj)                                    ! I
        real prin(ni,nj)                                       ! I
        real prin_diff(ni,nj)                                  ! L
        
        real stnp_bk(ni,nj)                                    ! I
        real stnp(ni,nj)                                       ! I/O


        write(6,*)' Subroutine pstn_anal (using pressure in mb)'

        call get_r_missing_data(r_missing_data,istatus)

!       Check bounds of prin/stnp fields
        if(back_sp .eq. 1)then
            call array_minmax(stnp_bk,ni,nj,rmin,rmax,r_missing_data)
            if(rmin .lt. 300. .or. rmax .gt. 1100)then
                write(6,*)' ERROR: STNP bkgnd range out of bounds'      
     1                   ,rmin,rmax
                stop
            endif
        endif

        call array_minmax(prin,ni,nj,rmin,rmax,r_missing_data)
        if(rmin .lt. 850. .or. rmax .gt. 1100)then
            write(6,*)' Warning: input P analysis range out of bounds'      
     1               ,rmin,rmax
            write(6,*)' Check if expected for reduced pressure height'      
        endif

        if(back_mp .eq. 1)then
            call array_minmax(prin_bk,ni,nj,rmin,rmax,r_missing_data)
            if(rmin .lt. 850. .or. rmax .gt. 1100)then
                write(6,*)' Warning: input P bkgnd range out of bounds'       
     1                   ,rmin,rmax
                write(6,*)
     1                 ' Check if expected for reduced pressure height'
            endif

            call diff(prin,prin_bk,prin_diff,ni,nj)
            call array_minmax(prin_diff,ni,nj,rmin,rmax,r_missing_data)       
            if(rmin .lt. -50. .or. rmax .gt. +50.)then
                write(6,*)
     1                  ' Warning: input P bkg diff range out of bounds'       
     1                   ,rmin,rmax
            endif
        endif

        if(back_mp .ne. 1 .or. back_sp .ne. 1)then
            write(6,*)' Skipping stnp adjustment'
            return
        endif

!       Adjust stnp field
        write(6,*)' Performing stnp adjustment using prin/prin_bk'

        do i = 1,ni
        do j = 1,nj
            if(       prin_bk(ni,nj) .ne. r_missing_data       
     1          .and. prin(ni,nj)    .ne. r_missing_data              
     1          .and. stnp_bk(ni,nj) .ne. r_missing_data   )then

                stnp(i,j) = ( prin(i,j) / prin_bk(i,j) ) * stnp_bk(i,j)       

            else
                stnp(i,j) = stnp_bk(i,j)

            endif
      
        enddo ! j
        enddo ! i

        return
        end

c
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
c                               11-15-99  Clean up mix ratio calc (QC check).
c
c*****************************************************************************
c
	real dx(ni,nj), dy(ni,nj)
c
	real p(ni,nj), td(ni,nj), u(ni,nj), v(ni,nj), thadv(ni,nj)
	real qcon(ni,nj), q(ni,nj), theta(ni,nj), qadv(ni,nj)
	real t(ni,nj), tadv(ni,nj)
c
	integer qbad
c
c
c.....	Calculate mixing ratio.
c.....	Units:  g / kg
c
	qbad = 0  !q ok
	do j=1,nj
	do i=1,ni
	   tdp = td(i,j) - 273.15          ! convert K to C
	   tl = (7.5 * tdp) / (237.3 + tdp)
	   e = 6.11 * 10. ** tl
	   if(p(i,j).le.0.0 .or. p(i,j).eq.e) then
	      write(6,990) i,j,p(i,j)
	      q(i,j) = badflag
	      qbad = 1                     !have a bad q field
	   else
	      drprs = 1. / (p(i,j) - e)    !invert to avoid further divisions.
	      q(i,j) = 622. * e * drprs	   !mixing ratio using (0.622*1000) for g/kg.
	   endif	
	enddo !i
	enddo !j
990	format(1x,' ERROR. Bad pressure in mixing ratio calc at point ',
     &         2i5,'-- pressure: ',f12.4,' calculated e: ',f12.4)
c
c.....	Compute moisture flux convergence on the laps grid.
c.....	Units:  g / kg / sec
c
	if(qbad .eq. 1) then
	   call constant(qcon, badflag, ni,nj)
	   go to 30
	endif
	do j=2,nj-1
	do i=2,ni-1
	  ddx1 = ((q(i,j-1) + q(i,j)) * .5) * u(i,j-1)
	  ddx2 = ((q(i-1,j) + q(i-1,j-1)) * .5) * u(i-1,j-1)
	  ddx = (ddx1 - ddx2) / dx(i,j)
	  ddy1 = ((q(i-1,j) + q(i,j)) * .5) * v(i-1,j)
	  ddy2 = ((q(i,j-1) + q(i-1,j-1)) * .5) * v(i-1,j-1)
	  ddy = (ddy1 - ddy2) / dy(i,j)
	  qcon(i,j) = - ddx - ddy
	enddo !i
	enddo !j
	call bounds(qcon,ni,nj)
 30	continue
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
	if(qbad .eq. 1) then
	   call constant(qadv, badflag, ni,nj)
	   go to 50
	endif
	do j=2,nj-1
	do i=2,ni-1
	  dqa1 = (q(i,j) - q(i-1,j)) / dx(i,j)
	  dqa2 = (q(i,j-1) - q(i-1,j-1)) / dx(i,j)
	  dqdx = (u(i,j-1) + u(i-1,j-1)) * (dqa1 + dqa2) * .25
	  dqa3 = (q(i,j) - q(i,j-1)) / dy(i,j)
	  dqa4 = (q(i-1,j) - q(i-1,j-1)) / dy(i,j)
	  dqdy = (v(i-1,j) + v(i-1,j-1)) * (dqa3 + dqa4) * .25
	  qadv(i,j) = - dqdx - dqdy     ! g/kg/sec
	enddo !i
	enddo !j
	call bounds(qadv,ni,nj)
 50	continue
c
c.....	Send the fields back to the main program.
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
c                         07-13-99  Change stn character arrays.
c                                     Rm *4 from declarations.
c                         11-15-99  Add writes to log file.
c
c     Notes:
c
c======================================================================
c
	integer ni,nj,mxstn
	real field(ni,nj), ob(mxstn), interp_ob
	real x1a(ni), x2a(nj), y2a(ni,nj)
	integer ii(mxstn), jj(mxstn)
	character title*(*), stn(mxstn)*20, stn_mx*5, stn_mn*5
        logical l_bilinear_interp 
c
c.... Start.
c
!       This will interpolate the grid to the stations using a significantly
!       faster bilinear interpolation instead of the splines.
        l_bilinear_interp = .true.

	num = 0
	num_ea1 = 0
	num_ea2 = 0
	num_ea3 = 0
	abs_diff = 0.
	sum = 0.
        sumsq = 0.
	amean = 0.
	rms = 0.
	diff_mx = -1.e30
	diff_mn = 1.e30
	print *,' '
	write(6,900) title
	write(iunit,900) title
 900	format(/,2x,a,/'          sta      i   j       '
     1        ,'grid      obs       diff')
c
	ea1 = field_ea
	ea2 = field_ea * 2.
	ea3 = field_ea * 3.

        if(.not. l_bilinear_interp)then
c
c....       Find the 2nd derivative table for use by the splines.
	    call splie2(x1a,x2a,field,ni,nj,y2a)

        endif
c
c....   Now call the spline for each station in the grid.
c
	do i=1,n_obs_b
	   if(ii(i).lt.1 .or. ii(i).gt.ni) go to 500
	   if(jj(i).lt.1 .or. jj(i).gt.nj) go to 500
	   aii = float(ii(i))
	   ajj = float(jj(i))

           if(l_bilinear_interp)then
               call bilinear_laps(aii,ajj,ni,nj,field,interp_ob)
           else
	       call splin2(x1a,x2a,field,y2a,ni,nj,aii,ajj,interp_ob)
           endif
c
	   if(ob(i) .le. badflag) then
	      diff = badflag
	   else
	      diff = ob(i) - interp_ob 
	      sum = diff + sum
              sumsq = sumsq + diff**2
	      adiff = abs(diff)
	      abs_diff = abs_diff + adiff
	      num = num + 1
c
	      if(adiff .gt. diff_mx) then
		 diff_mx = adiff
		 stn_mx = stn(i)(1:5)
	      endif
	      if(adiff .lt. diff_mn) then
		 diff_mn = adiff
		 stn_mn = stn(i)(1:5)
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
	   write(iunit,905) i, stn(i)(1:5), ii(i), jj(i), 
     1                      interp_ob, ob(i), diff
 905	   format(4x,i5,1x,a5,1x,i4,1x,i4,1x,3f10.2)
c
 500	enddo !i
c
c.... Get the average diff over the obs.
c     
	ave_diff = badflag
	amean = badflag
	if(num .gt. 0) amean = sum / float(num)
	if(num .gt. 0) ave_diff = abs_diff / float(num)
	if(num .gt. 0 .and. sumsq .ge. 0.)rms = sqrt(sumsq / float(num))       

	write(6,909) amean, num
	write(iunit,909) amean, num
 909	format(/,' Mean difference:    ',f10.2,' over ',i6,' stations.')       

	write(6,910) ave_diff, num
	write(iunit,910) ave_diff, num
 910	format(' Average difference: ',f10.2,' over ',i6,' stations.')

	write(6,915) rms, num
	write(iunit,915) rms, num
 915	format(' RMS difference:     ',f10.2,' over ',i6,' stations.')

	write(iunit,920) diff_mx, stn_mx
 920	format(' Maximum difference of ',f10.2,' at ',a5)

	write(iunit,925) diff_mn, stn_mn
 925	format(' Minimum difference of ',f10.2,' at ',a5)

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
	write(6, 931)
	write(iunit, 931)
 931	format(1x,'===============================================')
c
	return
	end
c
c
c
c
      subroutine splie2(x1a,x2a,ya,m,n,y2a)

c	15 May 1991  birkenheuer

      dimension x1a(m),x2a(n),ya(m,n),y2a(m,n),ytmp(n*50),y2tmp(n*50)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline_db(x2a,ytmp,n,1.e30,1.e30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      end
c
c
      subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)

c	15 May 1991 birkenheuer

       dimension x1a(m),x2a(n),ya(m,n),y2a(m,n),ytmp(n*50),
     &          y2tmp(n*50),yytmp(n*50)
       do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12     continue
       call spline_db(x1a,yytmp,m,1.e30,1.e30,y2tmp)
       call splint(x1a,yytmp,y2tmp,m,x1,y)
       return
       end
c
c
      subroutine spline_db(x,y,n,yp1,ypn,y2)

c	15 may 1991  birkenheuer

      dimension x(n),y(n*50),y2(n*50),u(n*50)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then   ! test for overflow condition
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
c
c
      subroutine splint(xa,ya,y2a,n,x,y)


c	15 May 1991 Birkenheuer

      dimension xa(n),ya(n*50),y2a(n*50)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     1      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
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
c       P. Stamus       15 Nov 99  Change checks of incoming values so the
c                                    AF can run over Mt. Everest and Antarctica.
c
c       Notes:  This routine may or may not be giving reasonable results over
c               extreme areas (Tibet, Antarctica).  As noted above, there are
c               questions about use of the std lapse rate to do these reductions.
c               15 Nov 99
c
c================================================================================
c
	implicit none
	real lapse_t, lapse_td, temp, dewp, pres, elev, redpres, ref_lvl
	real badflag
	real gor, ctv
	parameter(gor=0.03414158,ctv=0.37803)
	real dz, dz2, t_mean, td_mean, td, e, tkel, tv, esw
!	DATA GOR,ZBOU,CTV/.03414158,1612.,.37803/
!	data gor,ctv/.03414158,.37803/
		!GOR= acceleration due to gravity divided by the dry air gas
		!     constant (9.8/287.04)
		!F2M= conversion from feet to meters
		! *** zbou is now the standard (reduction) level 12-19-89 ***
		!CTV= 1-EPS where EPS is the ratio of the molecular weight of
		!     water to that of dry air.


C** Check input values......good T, Td, & P needed to perform the reduction.
cc	if(dewp.gt.temp .or. pres.le.620. .or. pres.gt.1080. .or.
cc     &      temp.lt.-30. .or. temp.gt.120. .or. dewp.lt.-35. .or.
cc     &      dewp.gt.90.) then
	if(dewp.gt.temp .or. pres.le.275. .or. pres.gt.1150. .or.
     &      temp.lt.-130. .or. temp.gt.150. .or. dewp.lt.-135. .or.
     &      dewp.gt.100.) then
	   print *,' Warning. Bad input to reduce_p routine.',dewp,temp
     &                                                       ,pres   
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

	subroutine read_sfc_verif_history(lun,r_missing_data             ! I
     1                                   ,mo,mf,mt                       ! I
     1                                   ,stn_a,bkg_a,obs_a,diff_a,nsta  ! O
     1                                   ,istatus)

        character*150 dirname,filename
        character*17 basename
        character*100  line
        character*5    stn

        integer mo,mf,mt

        character*5 stn_a(mo)
        real bkg_a(mo,mf,mt)
        real obs_a(mo,mf,mt)
        real diff_a(mo,mf,mt)

        integer mapcount(mo)

        write(6,*)' Subroutine read_sfc_verif_history'

        ierr = 0

        call get_sfc_badflag(badflag,istatus)

        bkg_a = r_missing_data
        obs_a = r_missing_data
        diff_a = r_missing_data

        nsta = 0
        ifield = 0

!       dirname = '/data/fab/parallel/laps/data/log/qc'
!       call s_len(dirname,lend)
        call get_directory('log',dirname,lend)
        dirname = dirname(1:lend)//'/qc'
        lend = lend + 3

        do itime = 1,mt

          write(basename,11)itime-1
 11	  format('laps_sfc.ver.',i2.2,'00')

          filename = dirname(1:lend)//'/'//basename

          write(6,*)' Open ',filename

          open(lun,file=filename,status='old',err=950)

          ifound_last = 0
          iblock = 0
          mapcount = 0 ! array

 100	  read(lun,101,err=990,end=900)line
 101	  format(a)

          call s_len2(line,lenl)

          if(itime .le. 1 .and. ifield .le. 2 .and. nsta .le. 50)then
              write(6,102)lenl,line(1:65)
 102	      format(24x,' newline is ',i3,a)
          endif

          if(lenl .eq. 56)then
 	    read(line,905,err=906) i, stn(1:5), ilaps, jlaps, 
     1                      bkg, ob, diff
 905	    format(4x,i5,1x,a5,1x,i4,1x,i4,1x,3f10.2)
            ifound = 1

            goto 910

 906	    ierr = ierr + 1
            if(ierr .le. 100)then
                write(6,*)' Error reading line with length 56:'
                write(6,*)line
            endif

            goto 100 ! try another line

          else
            ifound = 0

          endif

 910	  if(ifound_last .eq. 0 .and. ifound .eq. 1)then
            I4_elapsed = ishow_timer()
            iblock = iblock + 1
            write(6,*)' Start text block',iblock
            nmap = 0
            nsearch = 0
            nnew = 0

          elseif(ifound_last .eq. 1 .and. ifound .eq. 0)then
            write(6,104)iblock,nmap,nsearch,nsta,nnew
 104        format(' End text block',i4,3x     
     1            ,' nmap/nsearch/nsta/nnew=',4i7)                     

          endif
      
          if(ifound .eq. 1 .and. iblock .ne. (iblock/2)*2 )then
            ifield = iblock/2 + 1

            if(itime .eq. 1 .and. ifield .eq. 1)then
              nsta = nsta + 1
              ista = nsta
              stn_a(ista) = stn(1:5)
              nnew = nnew + 1
              if(nnew .le. 50)then
                  write(6,*)' new station = ',nsta,stn_a(ista)
              endif
              mapcount(i) = nsta

            else ! repeating loop   
              if(mapcount(i) .eq. 0)then ! search for station in list
                do ii = 1,nsta
                  if(stn_a(ii) .eq. stn)then ! match                
                     ista = ii
                     goto120
                  endif
                enddo

                nsta = nsta + 1          ! no match, add station to list
                ista = nsta
                stn_a(ista) = stn(1:5)
                nnew = nnew + 1
                if(nnew .le. 50)then
                  write(6,*)' new station = ',nsta,stn_a(ista)
                endif                                             

 120            mapcount(i) = ista
    	        nsearch = nsearch + 1

              else ! obtain station number from mapcount
                ista = mapcount(i) 
                nmap = nmap + 1

              endif

            endif

            if(ifield .le. mf)then
                if(itime .le. 1 .and. ifield .le. 2 
     1                          .and. nsta .le. 50)then
                  write(6,103)iblock,ifield,ista,line(1:56)
 103	          format(1x,' place text in odd block ',i3,i3,i6,1x,a)
                endif

                if(diff .ne. badflag)then
                  bkg_a(ista,ifield,itime)  = bkg
                  obs_a(ista,ifield,itime)  = ob
                  diff_a(ista,ifield,itime) = diff
                endif
            else ! skip rest of fields
                write(6,*)' Skipping rest of fields for this time'
                goto 900
                
            endif

          endif

          ifound_last = ifound

          goto 100

 900	  continue ! end of file

          close(lun)

          goto 960

 950	  write(6,*)' Error opening itime ',itime

 960	  write(6,*)' end time: itime/nsta = ',itime,nsta
          write(6,*)

        enddo ! itime

 990	continue

        return
        end

	subroutine sfc_verif_qc(r_missing_data                  ! I
     1                         ,mo,mf,mt                        ! I
     1                         ,stn_a,bkg_a,obs_a,diff_a,nsta   ! I
     1                         ,bias_a,obs_mean,obs_std)        ! O

!       Calculate obs minus background statistics, and obs statistics

!       Note that wind directions are inserted into the 'obs_a' array for output

        character*5 stn_a(mo)
        real bkg_a(mo,mf,mt)
        real obs_a(mo,mf,mt)
        real diff_a(mo,mf,mt)

        real bias_a(mo,mf)
        real obs_mean(mo,mf)
        real obs_std(mo,mf)

        write(6,*)' Subroutine sfc_verif_qc'

        bias_a = r_missing_data
        obs_mean = r_missing_data
        obs_std = r_missing_data

!       Field array indices
        iv_t = 1
        iv_td = 2
        iv_u = 3
        iv_v = 4
        iv_spd = 5
        iv_mslp = 6
        iv_redp = 7
        iv_tgd = 8
        iv_dir = 9

        do ifield = 1,2 ! t and td
          do ista = 1,nsta
            icount = 0
            sum = 0.
            sumsq = 0.
            diff_sum   = 0.
            diff_sumsq = 0.
            do itime = 1,mt
              if(diff_a(ista,ifield,itime) .ne. r_missing_data)then
                icount = icount + 1             
                diff_sum   = diff_sum   + diff_a(ista,ifield,itime)             
                diff_sumsq = diff_sumsq + diff_a(ista,ifield,itime)**2
              endif

              sum   = sum   + obs_a(ista,ifield,itime)             
              sumsq = sumsq + obs_a(ista,ifield,itime)**2

            enddo

            frac = float(icount) / float(mt)

            if(frac .ge. 0.5)then ! enough reports to process station stats
              bias = diff_sum / float(icount)
              bias_a(ista,ifield) = bias
              if(ista .le. 50)then
                write(6,*)stn_a(ista),frac
     1                   ,' bias = ',bias                        
              endif

              rmean = sum / float(icount)
              obs_mean(ista,ifield) = rmean

              obs_std(ista,ifield) = 
     1          sqrt( (sumsq / float(icount)) - rmean**2)

            else 
              if(ista .le. 50)then
                write(6,*)stn_a(ista),frac
              endif

            endif

          enddo ! ista
        enddo ! ifield

!       Obtain wind direction
        write(6,*)' Obtaining wind directions'
        do ista = 1,nsta
          do itime = 1,mt
            u = obs_a(ista,iv_u,itime)
            v = obs_a(ista,iv_v,itime)
            if(u .ne. r_missing_data .and. v .ne. r_missing_data)then
               dir = atan3d(-u,-v)
               obs_a(ista,iv_dir,itime) = dir
            endif 
          enddo 
        enddo

!       Compute standard deviation of wind direction/speed
        do ifield = 5,9,4
          do ista = 1,nsta
            icount = 0
            sum   = 0.
            sumsq = 0.
            do itime = 1,mt
              if(obs_a(ista,ifield,itime) .ne. r_missing_data)then
                icount = icount + 1             
                sum   = sum   + obs_a(ista,ifield,itime)             
                sumsq = sumsq + obs_a(ista,ifield,itime)**2
              endif
            enddo ! itime

            frac = float(icount) / float(mt)

!           (stddev ala the IDL code)
!           MEAN[*,*] = SUM[*,*] / nmembers
            if(icount .ge. 2 .and. frac .ge. 0.5)then
              rmean = sum / float(icount)
              obs_mean(ista,ifield) = rmean

!             SPREAD[*,*] = SQRT( (SUMSQ[*,*] / nmembers ) - MEAN[*,*]^2.0 )
              obs_std(ista,ifield) = 
     1          sqrt( (sumsq / float(icount)) - rmean**2)

            endif

          enddo ! ista
        enddo ! ifield

        return
        end
