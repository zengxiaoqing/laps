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
        subroutine windspeed(u,v,ff,ni,nj)
c
c...............................................................................
c
c       Routine to calculate the wind speed from the wind components.
c
c       Changes:
c               P.A. Stamus     01-11-89        Original (based on Baker's)
c                               07-20-89        Removed conversion to kts.
c                               08-21-89        Bad data check.
c                               09-20-89        Add implicit none.
c                               11-15-90        Bag bad data ck. Vectorize.
c
c       Inputs/Outputs:
c
c          Variable     Var Type     I/O   Description
c         ----------   ----------   ----- -------------
c          u, v            RA         I    Wind components (m/sec).
c          ni, nj          I          I    Grid dimensions in x, y.
c          ff              RA         O    Wind speed (m/sec).
c
c       User Notes:
c
c          1.  Units are not changed (output units = input units).
c
c...............................................................................
c
        implicit none
        integer ni, nj, i, j
        real u(ni,nj), v(ni,nj), ff(ni,nj)
c
c.....  Compute the wind speed at each point.
c
        do j=1,nj
        do i=1,ni
          ff(i,j) = sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j))
        enddo !i
        enddo !j
c
        return
        end
