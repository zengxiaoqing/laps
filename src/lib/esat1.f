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
        subroutine esat1(t,es,ni,nj)
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
c       P. Stamus         01-26-89      Changed to subroutine.
c                         08-21-89      Add bad data check.
c                         09-20-89      Add implicit none.
c                         11-14-90      Bad bad data ck. Vectorize.
c
C   THIS FUNCTION RETURNS THE SATURATION VAPOR PRESSURE OVER
C   WATER (MB) GIVEN THE TEMPERATURE (K).
C   THE ALGORITHM IS DUE TO NORDQUIST, W.S.,1973: "NUMERICAL APPROXIMA-
C   TIONS OF SELECTED METEORLOLGICAL PARAMETERS FOR CLOUD PHYSICS PROB-
C   LEMS," ECOM-5475, ATMOSPHERIC SCIENCES LABORATORY, U.S. ARMY
C   ELECTRONICS COMMAND, WHITE SANDS MISSILE RANGE, NEW MEXICO 88002.
c
c.....  NOTE: T must be in degrees K coming into this routine.
c
        implicit none
        integer ni, nj, i, j
        real p1, p2, c1, term
        real t(ni,nj), es(ni,nj)
c
        do j=1,nj
        do i=1,ni
          p1 = 11.344 - 0.0303998 * t(i,j)
          p2 = 3.49149 - 1302.8844 / t(i,j)
          c1 = 23.832241 - 5.02808 * alog10(t(i,j))
          term = c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/t(i,j)
          es(i,j) = 10. ** term
        enddo !i
        enddo !j
c
        return
        end
