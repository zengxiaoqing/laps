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
        subroutine regress_precip(num_sfc,radar,gauge,a_t,b_t,rbar,gbar
     1                           ,istatus)
c
c*******************************************************************************
c
c       Routine to calculate regression coefficients from gauge/radar data
c       (formerly t and elev data).

c       Values returned are the coefficients for the regression equations
c       for the temperatures and dew points.  Also calculates the mean
c       elevation for the surface stations.
c
c       Changes:
c               P.A. Stamus     12-01-88        Original (from J. McGinley)
c
c       Inputs/Outputs:
c
c          Variable     Var Type    I/O   Description
c         ----------   ----------  ----- -------------
c          num_sfc         I         I    Number of surface stations.
c          radar           RA        I    Station elevation.
c          gauge           RA        I    Temperature.
c          a_t             R         O    'a'      "       "    "    "  (slope)
c          b_t             R         O    'b' regression value for temp.(intrcp)
c          rbar            R         O    Mean elevation of the stations.
c
c       User Notes:
c
c       1. Units are not changed in this routine.
c
c*******************************************************************************
c
        real radar(num_sfc), gauge(num_sfc)                           

c
        badflag = -99.9
c
c.....  Set up storage variables.
c
        cnt = 0.
        sumht = 0.
        sumh = 0.
        sumt = 0.
        sumh2 = 0.
        sumt2 = 0.
c
c.....  Gather sums and then calculate the 'a' and 'b' for the regression
c.....  equation y = az + b, for both the temperature and dew point.  The
c.....  'b' is the intercept with sea level, and the 'a' is the lapse rate.
c.....  'y' represents the gauge value and 'z' is radar/analyzed
c.....  Also calculate the mean elevation of the stations.
c
        istatus = 0

        gmax = 0.
        rmax = 0.
        gmin = 999.
        rmin = 999.

        write(6,*)'     gauge analyzed'
        do 10 i=1,num_sfc
!         if(radar(i).le.badflag .or. gauge(i).le.badflag) go to 10
          write(6,5)i,gauge(i),radar(i)
5         format(i3,2f8.3)
          sumht = (radar(i) * gauge(i)) + sumht
          sumh = radar(i) + sumh
          sumh2 = (radar(i) * radar(i)) + sumh2
          sumt = gauge(i) + sumt
          cnt = cnt + 1.

          rmin = min(radar(i),rmin)
          rmax = max(radar(i),rmax)
          gmin = min(gauge(i),gmin)
          gmax = max(gauge(i),gmax)

          istatus = 1
10      continue

        if(cnt .eq. 0.)then
            write(6,*)' Count = 0'
            istatus = 0
            return
        endif

!       Slope    
        denominator = (cnt*sumh2 - sumh*sumh)
        if(denominator .ne. 0.)then
            a_t = (cnt*sumht - sumh*sumt) / denominator                
        else
            a_t = 0.
            istatus = 0
        endif

!       Intercept
        b_t = (sumt - a_t * sumh) / cnt
c
        rbar = sumh / cnt
        gbar = sumt / cnt

        write(6,*)'num_sfc,rbar,gbar,a_t,b_t'
        write(6,*)num_sfc,rbar,gbar,a_t,b_t

        write(6,11)gmin,gmax
11      format('  Gauge range = ',2f9.3)    
        write(6,12)rmin,rmax
12      format('  Radar range = ',2f9.3)    

        if(a_t .lt. 0.1 .OR. a_t .gt. 10.)then
           write(6,*)' Warning, slope is ill conditioned'
           istatus = 0
        endif
c
c.....  End of routine
c
        return
        end
