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

        subroutine latlon_to_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj,ist
     1atus)

!       1991            Steve Albers
!       1994            Steve Albers - partially added lambert option
!       This routine assumes polar stereographic, or tangent lambert projection

        real*4 rlat                         ! Input Lat
        real*4 rlon                         ! Input Lon
        real*4 lat(ni,nj),lon(ni,nj)        ! Input (Arrays of LAT/LON)
        integer ni,nj                       ! Input (LAPS Dimensions)
        real*4 ri,rj                        ! Output (I,J on LAPS Grid)
        integer istatus                     ! Output

        real*4 polat,polon,n

        save init,umin,umax,vmin,vmax,polat,polon
        data init/0/

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(init .eq. 0)then
            if(c6_maproj .eq. 'plrstr')then ! polar stereo
                polat = 90.
                polon = standard_longitude

                a=polat-lat(1,1)
                b=180.-lon(1,1)+polon
                r = tand(a/2.)
                umin=r*sind(b)
                vmin=r*cosd(b)

                a=polat-lat(ni,nj)
                b=180.-lon(ni,nj)+polon
                r = tand(a/2.)
                umax=r*sind(b)
                vmax=r*cosd(b)

            else          ! lambert
                polat = standard_latitude
                polon = standard_longitude

                n = abs(cosd(90.-polat))
                r = (tand(45.-lat(1,1)/2.))**n
                umin = r* sind(n*(lon(1,1)-polon))
                vmin = (-1.)*r*cosd(n*(lon(1,1)-polon))

                n = cosd(90.-polat)
                r = (tand(45.-lat(ni,nj)/2.))**n
                umax = r* sind(n*(lon(ni,nj)-polon))
                vmax = (-1.)*r*cosd(n*(lon(ni,nj)-polon))


            endif



            write(6,101)umin,umax,vmin,vmax
101         format(1x,' Initializing latlon_to_rlaps',4f10.5)
            init = 1
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale

        if(c6_maproj .eq. 'plrstr')then ! polar stereo

c           compute ulaps and vlaps (Polar stereographic grid)

            a=polat-rlat
            b=180.-rlon+polon

            r = tand(a/2.)

            ulaps=r*sind(b)
            vlaps=r*cosd(b)

        else          ! lambert

            n = cosd(90.-polat)
            r = (tand(45.-rlat/2.))**n
            ulaps = r* sind(n*(rlon-polon))
            vlaps = (-1.)*r*cosd(n*(rlon-polon))

        endif

        ri = (ulaps - u0) / uscale
        rj = (vlaps - v0) / vscale

!       Set status if location of point rounded off is on the LAPS grid
        if(nint(ri) .ge. 1 .and. nint(ri) .le. ni .and.
     1     nint(rj) .ge. 1 .and. nint(rj) .le. nj         )then
            istatus = 1
        else
            istatus = 0
        endif

        return
        end


        function projrot_laps(rlon)

        real*4 n

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic
            projrot_laps = standard_longitude - rlon

        else       ! lambert conformal
            n = cosd(90.-standard_latitude)
            projrot_laps = n * (standard_longitude - rlon)

        endif

        return
        end



