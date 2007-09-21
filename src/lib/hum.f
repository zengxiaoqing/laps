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
        subroutine hum(t,td,rh,ni,nj,sat_t,sat_td)
c
C       G.S. Stipanuk     1973            Original version.
C       Reference Stipanuk paper entitled:
C            "ALGORITHMS FOR GENERATING A SKEW-T, LOG P
C            DIAGRAM AND COMPUTING SELECTED METEOROLOGICAL
C            QUANTITIES."
C            ATMOSPHERIC SCIENCES LABORATORY
C            U.S. ARMY ELECTRONICS COMMAND
C            WHITE SANDS MISSILE RANGE, NEW MEXICO 88002
C            33 PAGES
C       Baker, Schlatter  17-MAY-1982
c       P. Stamus       02-08-89        Changed to subroutine.
c                       06-14-89        Pass in the sat_t/sat_td work arrays.
c                       08-21-89        Add bad data check.
c                       09-20-89        Add implicit none.
c                       04-13-90        Drop scaling.
c                       11-15-90        Bag bad data ck.  Vectorize.
c
C   THIS FUNCTION RETURNS RELATIVE HUMIDITY (0.00 - 1.00) GIVEN THE
C   TEMPERATURE T AND DEW POINT TD (K).  AS CALCULATED HERE,
C   RELATIVE HUMIDITY IS THE RATIO OF THE ACTUAL VAPOR PRESSURE TO
C   THE SATURATION VAPOR PRESSURE.
c
c.....  NOTE:  The temperature and dewpoint must be in degrees K.
c
        implicit none
        integer ni, nj, i, j
        real t(ni,nj), td(ni,nj), rh(ni,nj)
        real sat_t(ni,nj), sat_td(ni,nj)
c
c.....  Calculate the actual and saturation vapor pressure, then the rh.
c
        call esat1(t,sat_t,ni,nj)
        call esat1(td,sat_td,ni,nj)
c
        do j=1,nj
        do i=1,ni
            rh(i,j) = (sat_td(i,j) / sat_t(i,j))
        enddo !i
        enddo !j
c
        return
        end
