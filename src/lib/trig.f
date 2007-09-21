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


        function sind(x)

        real x,sind,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        sind = sin(x*rpd)

        return
        end

        function cosd(x)

        real x,cosd,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        cosd = cos(x*rpd)

        return
        end

        function tand(x)

        real x,tand,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        tand = tan(x*rpd)

        return
        end

        function asind(x)

        real x,asind,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        asind = asin(x) / rpd

        return
        end

        function acosd(x)

        real x,acosd,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        acosd = acos(x) / rpd

        return
        end


        function atan2d(x,y)

        real x,y,atan2d,pi,rpd

        parameter (pi = 3.1415926535897932)
        parameter (rpd = pi / 180.)

        atan2_deg = atan2(x,y) / rpd

        return
        end








