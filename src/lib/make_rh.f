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
        function make_rh (p,t,ssh,t_ref)

c       this function is designed to compute (rh) from basic variables
c       p (mb), t(c) and q (g/kg) to give rh in fraction.  The reference
c       temperature t_ref (c) is used to describe the temperature at which
c       the liquid and ice phase change occurs.

c       Dan Birkenheuer    14 May 1993

        implicit none

        real p  !ambient pressure (mb)
        real t   !ambient temp (c)
        real ssh  !specific humidity in g/kg
        real t_ref  !phase reference temp (c)
        real make_rh  !routine output (fraction)
        external eslo
        external esice
        external es
        integer istatus
        real rmd
        real eslo  !function type
        real esice !function type
        real es    !function type
        real esat   !saturation vapor pressure of the environment

        real ew ! vapor pressure of the environment

        real mw_air !molecular weight of dry air
        real mw_vap  !molecular weight of water

        data mw_vap /18.0152/
        data mw_air /28.966/

        if (ssh .eq. 0.0) then   !the obvious limitation
                make_rh = 0.0
                return
        endif


c       first compute the saturation  vapor pressure of water

        if(t_ref .lt.-132.) t_ref = -132.  ! limit t_ref

        if (t .gt. t_ref .and. t .ge. -132. ) then ! liquid phase
c                                               eslo approx
                esat = eslo (t)

        elseif (t .le. t_ref) then ! assume ice phase
                esat = esice (t)

        else ! situation not covered (impossible)

           write(6,*) 'Warning.. t and t_ref conflict, missing data
     1 returned'
           call get_r_missing_data(rmd, istatus)

           esat = rmd

        endif

c       now derive ambient vap pres from ssh from ideal gas law and
c       definition of sh = mv/(mv+ma)

        ew = p/ ( ( (mw_vap/ssh*1000. - mw_vap)/ mw_air) + 1. )

c       now compute the relative humidity using the partial vapor
c       by simple ratio.

        make_rh = ew /esat  !  units are fraction (dimensionless)

        return
        end

