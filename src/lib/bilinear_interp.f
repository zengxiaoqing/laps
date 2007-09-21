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
        subroutine bilinear_interp (n1,m1,g1,n2,m2,g2,dn,dm)

c       author: d. birkenheuer 18 oct 90
c       checked vectorizing      8 November 1990    (D. Birkenheuer)

c       routine is designed to interpolate from a less dense grid (g1)
c       to a grid of higher density (g2)  where g1 is NESTED in g2.

c       for current purposes, g1 is lma array (8,8), n1=8, m1=8
c       g2 is laps array (57,57) n2=57,m2=57, and dn=dm=8 the spacing
c       between the nested gridpoints in the g2 array.

c       dn or dm can be computed by:
c       (n2-1)/(n1-1) = dn.....and similarly for dm

        implicit none

        integer n1,n2,m1,m2,dn,dm
        real g1(n1,m1),g2(n2,m2)

        real factorI, factorJ, fracI, fracJ, sum
        integer intI,intJ,i,j


        do j = 1,m2
        factorJ = float(j-1)/float(dm) +1.
        intJ = int(factorJ)
        fracJ = factorJ-intJ

        do i = 1,n2
        factorI = float(i-1)/float(dn) +1.
        intI = int(factorI)
        fracI = factorI-intI

        sum = 0.

        g2(i,j) = g1(intI,intJ) * (1.-fracI) * (1.-fracJ)
        sum = sum + (1.-fracI) * (1.-fracJ)

        if (intI+1 .le. n1) then
        g2(i,j) = g2(i,j) + g1(intI+1,intJ) * (fracI) * (1.-fracJ)
        sum = sum + (fracI) * (1.-fracJ)
        endif


        if (intJ+1 .le. m1) then
        g2(i,j) = g2(i,j) + g1(intI,intJ+1) * (fracJ) * (1.-fracI)
        sum = sum + (fracJ) * (1.-fracI)
        endif

        if (intI+1 .le. n1  .and. intJ+1  .le. m1)then
        g2(i,j) = g2(i,j) + g1(intI+1,intJ+1) * (fracI) * (fracJ)
        sum = sum + (fracI) * (fracJ)
        endif

        g2(i,j) = g2(i,j)/sum

        enddo
        enddo


        return
        end

