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
        function make_ssh (p,t,rh, t_ref)

c       this function is designed to compute (q) from basic variables
c       p (mb), t(c) and rh(fraction) to give (q) in (g/kg).  The reference
c       temperature t_ref (c) is used to describe the temperature at which
c       the liquid and ice phase change occurs.

        implicit none

        real p
        real t
        real rh
        real t_ref
        real make_ssh
        external eslo
        external esice
        external es
        real eslo  !function type
        real esice !function type
        real es    !function type

        real ew ! vapor pressure of the environment

        real mw_air !molecular weight of dry air
        real mw_vap  !molecular weight of water

        data mw_vap /18.0152/
        data mw_air /28.966/

c       stop the obvious

        if (t.lt.-200.) then
                make_ssh = 0.0
                return
        endif



c       first compute the ambient vapor pressure of water

        if (t_ref .lt. -132.) t_ref = -132.

        if (t .ge. t_ref .and. t .ge. -132. ) then ! liquid phase
c                                               eslo approx

                ew = eslo (t)

        else ! ice phase

                ew = esice (t)

        endif

c       now sat vap pres obtained compute local vapor pressure

        ew = ew * rh

c       now compute the specific humidity using the partial vapor
c       pressures of water vapor (ew) and dry air (p-ew).

        make_ssh = mw_vap *ew

        make_ssh = make_ssh / (make_ssh + mw_air*(p-ew) ) * 1000. ! g(kg)

        return
        end

