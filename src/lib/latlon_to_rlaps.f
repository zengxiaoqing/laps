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

        subroutine latlon_to_rlapsgrid(rlat,rlon,lat,lon,ni,nj,ri,rj
     1                                ,istatus)

!       1991            Steve Albers
!       1994            Steve Albers - partially added lambert option
!       1997            Steve Albers - Added both kinds of lambert projections
!                                    - as well as mercator projection

!       This routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        real*4 rlat                         ! Input Lat
        real*4 rlon                         ! Input Lon
        real*4 lat(ni,nj),lon(ni,nj)        ! Input (Arrays of LAT/LON)
        integer ni,nj                       ! Input (LAPS Dimensions)
        real*4 ri,rj                        ! Output (I,J on LAPS Grid)
        integer istatus                     ! Output

        save init,umin,umax,vmin,vmax
        data init/0/

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(init .eq. 0)then
            call latlon_to_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_to_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' Initializing latlon_to_rlapsgrid',4f10.5)
            init = 1
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale

!       Compute ulaps and vlaps

        call latlon_to_uv(rlat,rlon,ulaps,vlaps,istatus)

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

        subroutine rlapsgrid_to_latlon(ri,rj,lat,lon,ni,nj,rlat,rlon
     1                                ,istatus)

!       1997            Steve Albers 

!       This routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        real*4 ri,rj                        ! Input (I,J on LAPS Grid)
        real*4 lat(ni,nj),lon(ni,nj)        ! Input (Arrays of LAT/LON)
        integer ni,nj                       ! Input (LAPS Dimensions)
        real*4 rlat                         ! Output Lat
        real*4 rlon                         ! Output Lon
        integer istatus                     ! Output

        save init,umin,umax,vmin,vmax
        data init/0/

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(init .eq. 0)then
            call latlon_to_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_to_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' Initializing rlapsgrid_to_latlon',4f10.5)
            init = 1
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale

!       Compute lat,lon
        ulaps = u0 + uscale * ri
        vlaps = v0 + vscale * rj

        call uv_to_latlon(ulaps,vlaps,rlat,rlon,istatus)

        istatus = 1

        return
        end

        
        subroutine latlon_to_uv(rlat,rlon,u,v,istatus)

!       1997            Steve Albers 

!       This routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
            slat1 = standard_latitude
            polat = standard_latitude2
            slon = standard_longitude

            call latlon_to_uv_ps(rlat,rlon,slat1,polat,slon,u,v)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
            slat1 = standard_latitude
            slat2 = standard_latitude2
            slon = standard_longitude

            call latlon_to_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            slat1  = standard_latitude
            cenlon = grid_cen_lon_cmn

            call latlon_to_uv_mc(rlat,rlon,slat1,cenlon,u,v)

        else
            write(6,*)'latlon_to_uv: unrecognized projection '
     1                ,c6_maproj       
            stop

        endif

        istatus = 1

        return
        end


        subroutine latlon_to_uv_ps(rlat,rlon,slat,polat,slon,u,v)

        if(polat .ne. +90.)then
            write(6,*)' latlon_to_uv_ps: ERROR, polat .ne. +90 ',polat      
        endif

        a=polat-rlat
        b=180.-rlon+slon

        r = tand(a/2.)

        u=r*sind(b)
        v=r*cosd(b)

        return
        end

        subroutine latlon_to_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        real*4 n

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        if(slat1 .ge. 0)then
            s = -1.
        else
            s = +1.
        endif

        if(slat1 .eq. slat2)then ! tangent lambert
            n = cosd(90.-slat1)
        else                       ! two standard latitudes
            n = alog(cosd(slat1)/cosd(slat2))/
     1          alog(tand(45.-s*slat1/2.)/tand(45.-s*slat2/2.))
        endif

        r = (tand(45.-rlat/2.))**n
        u = r* sind(n*angdif(rlon,slon))
        v = s*r*cosd(n*angdif(rlon,slon))

        return
        end

        subroutine latlon_to_uv_mc(rlat,rlon,slat,cenlon,u,v)

        real*4 pi, rpd

        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.
 
        a = 90.-rlat
        b = cosd(slat)

        u = angdif(rlon,cenlon) * rpd * b
        v = alog(1./tand(a/2.))       * b

        return
        end


        
        subroutine uv_to_latlon(u,v,rlat,rlon,istatus)

!       1997            Steve Albers 

!       This routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
            slat1 = standard_latitude
            polat = standard_latitude2
            slon = standard_longitude

            call uv_to_latlon_ps(u,v,slat1,polat,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
            slat1 = standard_latitude
            slat2 = standard_latitude2
            slon = standard_longitude

            call uv_to_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            slat1  = standard_latitude
            cenlon = grid_cen_lon_cmn

            call uv_to_latlon_mc(u,v,slat1,cenlon,rlat,rlon)

        else
            write(6,*)'uv_to_latlon: unrecognized projection '
     1                ,c6_maproj       
            stop

        endif

        istatus = 1

        return
        end

        subroutine uv_to_latlon_ps(u,v,slat,polat,slon,rlat,rlon)

        if(polat .ne. +90.)then
            write(6,*)' uv_to_latlon_ps: ERROR, polat .ne. +90 ',polat      
        endif

        dist=sqrt(u**2+v**2)

        if (dist .eq. 0) then
            rlat=90.
            rlon=-90.

        else
            rlat=atand(dist/2.)
            rlat=90.-2.*rlat
c
            if (u .eq. 0.) then
                rlon=90.

            else
                if (u .gt. 0.) then
                    rlon=atand(v/u)
                else
                    rlon=atand(v/u)+180.
                endif

            endif

        endif
c
        rlon=amod(rlon+630.,360.) - 180.
c
        return
        end

        subroutine uv_to_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        real*4 n

        if(slat1 .ge. 0)then
            s = -1.
        else
            s = +1.
        endif

        if(slat1 .eq. slat2)then   ! tangent lambert
            n = cosd(90.-slat1)
        else                       ! two standard latitudes
            n = alog(cosd(slat1)/cosd(slat2))/
     1          alog(tand(45.-s*slat1/2.)/tand(45.-s*slat2/2.))
        endif

        rlon=slon+atand(-s*u/v)/n
        rlat=(90.- 2.*atand((-v/cosd(n*(rlon-slon)))**(1./n)))/s      

        return
        end

        subroutine uv_to_latlon_mc(u,v,slat,cenlon,rlat,rlon)

        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

        b = cosd(slat)

        rlat_abs = 90. - atand(exp(-abs(v)/b)) * 2.

        if(v .gt. 0)then
            rlat =  rlat_abs
        else
            rlat = -rlat_abs
        endif

        rlon = u/b/rpd + cenlon
        rlon = mod(rlon+540.,360.) - 180.

        return
        end


        function projrot_laps(rlon)

        real*4 n

        include 'lapsparms.cmn'

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(    c6_maproj .eq. 'plrstr')then ! polar stereographic
            projrot_laps = standard_longitude - rlon

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            slat1 = standard_latitude
            slat2 = standard_latitude2

            if(slat1 .ge. 0)then
                s = -1.
            else
                s = +1.
            endif

            if(slat1 .eq. slat2)then ! tangent lambert
                n = cosd(90.-slat1)
            else                       ! two standard latitudes
                n = alog(cosd(slat1)/cosd(slat2))/
     1              alog(tand(45.-s*slat1/2.)/tand(45.-s*slat2/2.))
            endif

            n = cosd(90.-standard_latitude)


            projrot_laps = n * angdif(standard_longitude,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        else
            write(6,*)'projrot_laps: unrecognized projection ',c6_maproj       
            stop

        endif

        return
        end
