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

        real*4 slat1,slat2,slon,n,pi

        parameter (pi=3.1415926535897932)

        save init,umin,umax,vmin,vmax,slat1,slat2,slon,cenlon
        data init/0/

        include 'lapsparms.cmn'

        if(iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' ERROR, get_laps_config not called'
            stop
        endif

        if(init .eq. 0)then
            if(c6_maproj .eq. 'plrstr')then ! polar stereo
                slat1 = standard_latitude
                slon = standard_longitude

                call latlon_to_uv_ps(slat1,slon
     1                              ,lat(1,1),lon(1,1)
     1                              ,umin,vmin)

                call latlon_to_uv_ps(slat1,slon
     1                              ,lat(ni,nj),lon(ni,nj)
     1                              ,umax,vmax)

            elseif(c6_maproj .eq. 'lambrt')then ! lambert
                slat1 = standard_latitude
                slat2 = standard_latitude2
                slon = standard_longitude

                call latlon_to_uv_lc(slat1,slat2,slon
     1                              ,lat(1,1),lon(1,1)
     1                              ,umin,vmin)

                call latlon_to_uv_lc(slat1,slat2,slon
     1                              ,lat(ni,nj),lon(ni,nj)
     1                              ,umax,vmax)

            elseif(c6_maproj .eq. 'merctr')then ! mercator
                slat1  = standard_latitude
                cenlon = grid_cen_lon_cmn
                call latlon_to_uv_mc
     1                  (slat1,cenlon,lat(1,1),lon(1,1),umin,vmin)

                call latlon_to_uv_mc
     1                  (slat1,cenlon,lat(ni,nj),lon(ni,nj),umax,vmax)        

            else
                write(6,*)'latlon_to_rlaps: unrecognized projection '
     1                    ,c6_maproj       
                stop

            endif


            write(6,101)umin,umax,vmin,vmax
101         format(1x,' Initializing latlon_to_rlaps',4f10.5)
            init = 1
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale


!       Compute ulaps and vlaps

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
            call latlon_to_uv_ps(slat1,slon,rlat,rlon
     1                          ,ulaps,vlaps)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
            call latlon_to_uv_lc(slat1,slat2,slon,rlat,rlon
     1                          ,ulaps,vlaps)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call latlon_to_uv_mc(slat1,cenlon,rlat,rlon,ulaps,vlaps)       

        else
            write(6,*)'latlon_to_rlaps: unrecognized projection '
     1                ,c6_maproj       
            stop

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

        subroutine latlon_to_uv_ps(slat,slon,rlat,rlon,u,v)

        if(slat .gt. 0.)then
            polat = +90.
        elseif(slat .lt. 0.)then
            polat = -90.
        else
            write(6,*)' subroutine latlon_to_uv_ps: error, slat = 0.'
        endif

        a=polat-rlat
        b=180.-rlon+slon

        r = tand(a/2.)

        u=r*sind(b)
        v=r*cosd(b)

        return
        end

        subroutine latlon_to_uv_lc(slat1,slat2,slon,rlat,rlon,u,v)

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

        subroutine latlon_to_uv_mc(slat,cenlon,rlat,rlon,u,v)

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

