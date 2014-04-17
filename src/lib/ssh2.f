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
        function ssh2 (p,t,td,t_ref)

c       This module computes the specific humidity
c       based on the same function calls to the m-thermo library (those
c       deemed the most accurate and fastest).  Unlike the conventional
c       mthermo routines, this specific humidity routine differentiates
c       between the liquid and ice phases for computing the specific
c       humidity.

c       this module is designed to take similar inputs as the thermo routine
c       ssh () module -- except that it returns more than just saturation sh
c       it returns the sh at the given dewpoint.  If you desire saturation sh
c       simply supply the ambient temperature as the dewpoint.  The other
c       added parameter denotes the ice/liquid phase change temperature.

c       dan birkenheuer
c       may 11 1993

        implicit none

        real p     !atmospheric pressure in (mb)
        real t     !atmospheric temperature (c)
        real td    !atmospheric dewpoint temperature (c)
        real ssh2  !returned specific humidity (g/kg)
        external esice
        external eslo
        external es
        real esice !function type
        real eslo  !function type
        real es    !function type
        real t_ref !reference temperature for the phase type (c)
        real ew    !computed vapor pressure of water based on the phase and td
c                   (mb)
        real mw_vap !molecular weight of water vapor
        real mw_air !molecular weight of dry air

        data mw_vap /18.0152/
        data mw_air /28.966/


c       trap the obvious error conditions

        if (td .gt. t) then
           print *, 'error condition, td greater than t '
           ssh2 = 0.0
           return
        endif

        if (td .lt. -200. ) then
           write (6,*) 'td is less than -199.C, zero vapor press assngd'
           ssh2 = 0.0
           return
        endif

        if (t_ref .lt. -132.) t_ref = -132.

        if (t .gt. t_ref  .and. td.ge.-132.) then !assume liquid phase

           ew = eslo (td)

        else                    !assume ice phase

           ew = esice (td)

        endif

c       now that we know the vapor pressure of the water substance, now
c       compute the specific humidity.  Specific humidity can be computed
c       effectively from the ideal gas law, which is the assumption that this
c       module makes concerning the atmospheric density.  The temperature,
c       gas constant, and volume all cancel out, leaving the computation
c       dependent only on the partial pressures of air and water vapor and
c       their molecular weights.

        ssh2 = mw_vap * ew

        ssh2 = ssh2 /(ssh2 + mw_air * (p -ew) ) * 1000.  ! conver to g/kg

c modify routine to work at very low pressures.  11/9/94 db

        if(ew .gt. p ) ssh2 = 1000.

        return
        end

