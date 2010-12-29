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
        subroutine mean_lapse_old(num_sfc,elev,t,td,a_t,b_t,a_td,b_td
     1                           ,hbar)
c
c*******************************************************************************
c
c       Routine to calculate a mean lapse rate from observed surface data.
c       Values returned are the coefficients for the regression equations
c       for the temperatures and dew points.  Also calculates the mean
c       elevation for the surface stations.
c
c       Changes:
c               P.A. Stamus     12-01-88        Original (from J. McGinley)
c                               12-22-88        Added consistency check on td.
c
c       Inputs/Outputs:
c
c          Variable     Var Type    I/O   Description
c         ----------   ----------  ----- -------------
c          num_sfc         I         I    Number of surface stations.
c          elev            RA        I    Station elevation.
c          t               RA        I    Temperature.
c          td              RA        I    Dew point temperature.
c          a_t             R         O    'a' regression value for temp.(intrcp)
c          b_t             R         O    'b'      "       "    "    "  (lapse)
c          a_td            R         O    'a'      "       "    "  dewpt.
c          b_td            R         O    'b'      "       "    "    "
c          hbar            R         O    Mean elevation of the stations.
c
c       User Notes:
c
c       1. Units are not changed in this routine.
c
c*******************************************************************************
c
        real elev(num_sfc), t(num_sfc), td(num_sfc)
c
        badflag = -99.9
c
c.....  Set up storage variables.
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
c.....  Gather sums and then calculate the 'a' and 'b' for the regression
c.....  equation y = bz + a, for both the temperature and dew point.  The
c.....  'a' is the intercept with sea level, and the 'b' is the lapse rate.
c.....  Also calculate the mean elevation of the stations.
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
10      continue
c
        b_t = (cnt*sumht - sumh*sumt) / (cnt*sumh2 - sumh*sumh)
        a_t = (sumt - b_t * sumh) / cnt
c
        b_td = (cntd*sumhtd - sumh*sumtd) / (cntd*sumh2 - sumh*sumh)
        a_td = (sumtd - b_td * sumh) / cntd
c
        hbar = sumh / cnt
c
c.....  Do a consistency check on the dewpoint regression.  If msl intercept
c.....  is below zero or if the dewpoint lapse rate is positive, set td slope
c.....  to t slope and slide intercept over.
c
        if(a_td.lt.0. .or. b_td.gt.0.) then
          write(6,900)
900       format(1x,'++ Suspect dewpoint regression in MEAN_LAPSE. ++')
          write(6,901) a_td, b_td
901       format(1x,'  MSL intercept = ',e12.4,'  Slope = ',e12.4)
          b_td = b_t
          diff = (sumt / cnt) - (sumtd / cntd)
          a_td = a_t - diff
          write(6,902)
902       format(1x,'  Setting values to:')
          write(6,901) a_td, b_td
        endif
c
c.....  End of routine
c
        return
        end
