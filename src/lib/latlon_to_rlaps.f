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
!       1997            Steve Albers - Added local stereographic
!       2002            Dan Birkenheuer - added input test for valid lats

cdoc    This routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.
cdoc
cdoc    Values returned are on the LAPS grid where I ranges from 1,ni and
cdoc    J ranges from 1,nj. Istatus is set to 1 if a valid I,J was determined
cdoc    to be in the LAPS domain, except that a buffer of 0.5 grid points 
cdoc    around the perimenter is allowed. If the point is farther outside
cdoc    the domain, values of I,J will still be returned, while istatus is set
cdoc    to 0. I increased from grid west to grid east, and J increases from
cdoc    grid south to grid north.

        real rlat                         ! Input Lat
        real rlon                         ! Input Lon
        real lat(ni,nj),lon(ni,nj)        ! Input (Arrays of LAT/LON)
        integer ni,nj                     ! Input (LAPS Dimensions)
        real ri,rj                        ! Output (I,J on LAPS Grid)
        integer istatus                   ! Input / Output

        save init,umin,umax,vmin,vmax
        data init/0/

        include 'grid_fname.cmn'

        if (abs(rlat) > 90.000) then
           if(istatus .ne. 100)then
              write(6,*) 'rejecting invalid latitude ',rlat
           endif
           istatus = -1
           return
        endif

        if(init.ne.nest)then
            call latlon_to_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_to_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' Initializing latlon_to_rlapsgrid',4f10.5)
            init = nest
        endif

        uscale = (umax - umin) / (float(ni) - 1.)
        vscale = (vmax - vmin) / (float(nj) - 1.)

        u0 = umin - uscale
        v0 = vmin - vscale

!       Compute ulaps and vlaps

        call latlon_to_uv(rlat,rlon,ulaps,vlaps,istatus)

        if(uscale .eq. 0. .or. vscale .eq. 0.)then
            write(6,*)
     1      ' SEVERE ERROR: (u|v)scale = 0 in latlon_to_rlapsgrid'
            stop
        else
            ri = (ulaps - u0) / uscale
            rj = (vlaps - v0) / vscale
        endif

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

cdoc    This routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.

        real ri,rj                        ! Input (I,J on LAPS Grid)
        real lat(ni,nj),lon(ni,nj)        ! Input (Arrays of LAT/LON)
        integer ni,nj                       ! Input (LAPS Dimensions)
        real rlat                         ! Output Lat
        real rlon                         ! Output Lon
        integer istatus                     ! Output

        save init,umin,umax,vmin,vmax
        data init/0/

        include 'grid_fname.cmn'

        if(init .ne. nest)then
            call latlon_to_uv(lat(1,1),lon(1,1),umin,vmin,istatus)
            call latlon_to_uv(lat(ni,nj),lon(ni,nj),umax,vmax,istatus)

            write(6,101)umin,umax,vmin,vmax
101         format(1x,' Initializing rlapsgrid_to_latlon',4f10.5)
            init = nest 
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

cdoc    This routine assumes a polar stereographic, lambert conformal,
cdoc    or mercator projection.

        use mem_namelist, ONLY: c6_maproj
     1                  ,slat1=>standard_latitude
     1                  ,slat2=>standard_latitude2
     1                  ,polat=>standard_latitude2
     1                  ,slon=>standard_longitude

        include 'trigd.inc'

        include 'grid_fname.cmn'

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
            call latlon_to_uv_ps(rlat,rlon,slat1,polat,slon,u,v)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
            if(abs(slat2).eq.90. .and. slat2.ne.slat1)then
               print*,'Error: lambert slat2 = 90.'
               istatus=1
               return
            endif
c           slat1 = standard_latitude
c           slat2 = standard_latitude2
c           slon = standard_longitude

            call latlon_to_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif
c           slat1  = standard_latitude
c           cenlon = grid_cen_lon_cmn

            call latlon_to_uv_mc(rlat,rlon,slat1,cenlon,u,v)

        elseif(c6_maproj .eq. 'latlon')then ! latlon (cylindrical equidistant)
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif

            call latlon_to_uv_ll(rlat,rlon,cenlon,u,v)

        else
            write(6,*)'ERROR in latlon_to_uv: unrecognized projection '
     1                ,c6_maproj       
            stop

        endif

        istatus = 1

        return
        end


        subroutine latlon_to_uv_ps(rlat_in,rlon_in,slat,polat,slon,u,v)

        include 'trigd.inc'
        if(abs(polat) .eq. 90.)then ! pole at N/S geographic pole
            if(.true.)then
                polon = slon
                call GETOPS(rlat,rlon,rlat_in,rlon_in,polat,polon)
                rlon = rlon - 270.  ! Counteract rotation done in GETOPS

            else ! .false. (older simple method)
                rlat = rlat_in
                rlon = rlon_in

            endif

            b = rlon - slon         ! rotate relative to standard longitude
            s = polat / abs(polat)

        else                        ! local stereographic
            polon = slon
            call GETOPS(rlat,rlon,rlat_in,rlon_in,polat,polon)
            b = rlon - 270.         ! rlon has zero angle pointing east
                                    ! b has zero angle pointing south
            s = 1.0

        endif

        a=90.-rlat
        r = tand(a/2.)      ! Consistent with Haltiner & Williams 1-21

!       b = angle measured counterclockwise from -v axis (zero angle south)
        u =   r * sind(b)
        v = (-r * cosd(b)) * s

        return
        end

        subroutine latlon_to_uv_lc(rlat,rlon,slat1,slat2,slon,u,v)

        include 'trigd.inc'
        real n

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        call lambert_parms(slat1,slat2,n,s,rconst)

        r = (tand(45.-s*rlat/2.))**n
        u =    r*sind(n*angdif(rlon,slon))
        v = -s*r*cosd(n*angdif(rlon,slon))

        return
        end

        subroutine latlon_to_uv_mc(rlat,rlon,slat,cenlon,u,v)

        include 'trigd.inc'
        real pi, rpd

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

        subroutine latlon_to_uv_ll(rlat,rlon,cenlon,u,v)

        include 'trigd.inc'
        real pi, rpd

        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.
 
        u = angdif(rlon,cenlon) * rpd 
        v = rlat                * rpd

        return
        end


        
        subroutine uv_to_latlon(u,v,rlat,rlon,istatus)

        use mem_namelist, ONLY: c6_maproj
     1                  ,slat1=>standard_latitude
     1                  ,slat2=>standard_latitude2
     1                  ,polat=>standard_latitude2
     1                  ,slon=>standard_longitude

!       1997            Steve Albers 

!       This routine assumes a polar stereographic, lambert conformal,
!       or mercator projection.

        include 'trigd.inc'

        include 'grid_fname.cmn'

        if(c6_maproj .eq. 'plrstr')then ! polar stereo
c           slat1 = standard_latitude
c           polat = standard_latitude2
c           slon = standard_longitude

            call uv_to_latlon_ps(u,v,slat1,polat,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'lambrt')then ! lambert
c           slat1 = standard_latitude
c           slat2 = standard_latitude2
c           slon = standard_longitude

            call uv_to_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif
c           slat1  = standard_latitude
c           cenlon = grid_cen_lon_cmn

            call uv_to_latlon_mc(u,v,slat1,cenlon,rlat,rlon)

        elseif(c6_maproj .eq. 'latlon')then ! latlon (cylindrical equidistant)
            call get_grid_center(cenlat_dum,cenlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif

            call uv_to_latlon_ll(u,v,cenlon,rlat,rlon)

        else
            write(6,*)'ERROR in uv_to_latlon: unrecognized projection '
     1                ,c6_maproj
            stop

        endif

        istatus = 1

        return
        end

        subroutine uv_to_latlon_ps(u,v,slat,polat,slon
     1                                         ,rlat_out,rlon_out)

        include 'trigd.inc'

        if(abs(polat) .eq. 90.)then
            s = polat / abs(polat)
        else
            s = 1.0
        endif

        r=sqrt(u**2+v**2)

        if (r .eq. 0) then
            rlat=90.
            rlon=0.

        else                           
            a=2.* atand(r)               ! From Haltiner & Williams 1-21
            rlat=90.- a
            rlon = atan2d(s*v,u)
            rlon = rlon + 90.
        endif

        if(.true.)then ! Rotate considering where the projection pole is
            polon = slon
!           This routine will rotate the longitude even if polat = +90.
            call PSTOGE(rlat,rlon,rlat_out,rlon_out,polat,polon)
        else
            rlat_out = rlat
            rlon_out = rlon + slon
        endif


        rlon_out = amod(rlon_out+540.,360.) - 180. ! Convert to -180/+180 range


        return
        end

        subroutine uv_to_latlon_lc(u,v,slat1,slat2,slon,rlat,rlon)

        include 'trigd.inc'
        real n

        call lambert_parms(slat1,slat2,n,s,rconst)

!       rlon=slon + atand(-s*u/v) /n
!       rlat=(90.- 2.*atand((-  v/cosd(n*(rlon-slon)))**(1./n)))/s      

        angle  = atan2d(u,-s*v)
        rlat = (90.- 2.*atand((-s*v/cosd(angle))**(1./n))) / s      
        rlon = slon + angle / n

        rlon = mod(rlon+540.,360.) - 180.          ! Convert to -180/+180 range

        return
        end

        subroutine uv_to_latlon_mc(u,v,slat,cenlon,rlat,rlon)

        include 'trigd.inc'
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
        rlon = mod(rlon+540.,360.) - 180.          ! Convert to -180/+180 range

        return
        end

        subroutine uv_to_latlon_ll(u,v,cenlon,rlat,rlon)

        include 'trigd.inc'
        parameter (pi=3.1415926535897932)
        parameter (rpd=pi/180.)

        rlat = v/rpd

        rlon = u/rpd + cenlon
        rlon = mod(rlon+540.,360.) - 180.          ! Convert to -180/+180 range

        return
        end


        function projrot_latlon(rlat,rlon,istatus)
        include 'trigd.inc'

cdoc    1997 Steve Albers    Calculate map projection rotation, this is the
cdoc                         angle between the y-axis (grid north) and
cdoc                         true north. Units are degrees.
!
!                            projrot_laps = (true north value of wind direction
!                                          - grid north value of wind direction)
       
!       Added 8/4/2000 to make sure these are declared even if not passed in
        integer istatus
        real rlat, rlon

        real n

        save init
        data init/0/

        character*6  c6_maproj

        include 'grid_fname.cmn'

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus.ne.1)then
           print*,'Error returned from get_c6_maproj'
           return
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic

            call get_standard_longitude(polon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_longitude'
               return
            endif
            call get_standard_latitudes(stdlat,polat,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_latitudes'
               return
            endif

            call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif


            if(polat .eq. +90.)then
                projrot_laps = polon - rlon

            elseif(polat .eq. -90.)then
                projrot_laps = rlon - polon 

            else ! abs(polat) .ne. 90.
                if(grid_cen_lat .eq. polat .and. 
     1             grid_cen_lon .eq. polon)then ! grid centered on proj pole

                    if(init .eq. 0)then
                        write(6,*)
     1                   ' NOTE: local stereographic projection.'
                        write(6,*)
     1                   ' Using approximation for "projrot_laps",'
     1                  ,' accurate calculation not yet in place.'
                        init = 1
                    endif

                    rn = cosd(90.-polat)
                    projrot_laps = rn * angdif(polon,rlon)      

                elseif(.true.)then
                    if(init .eq. 0)then
                        write(6,*)' ERROR in projrot_laps: '
                        write(6,*)' This type of local'
     1                  ,' stereographic projection not yet supported.'
                        write(6,*)' Grid should be centered on'
     1                  ,' projection pole.'
                        init = 1
                    endif

                    projrot_laps = 0.
         
                else ! .false.
!                   Find dx/lat and dy/lat, then determine projrot_laps

                endif

            endif ! polat

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_latitudes'
               return
            endif
            call get_standard_longitude(slon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_longitude'
               return
            endif

            call lambert_parms(slat1,slat2,n,s,rconst)

            call get_standard_longitude(stdlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_c6_maproj'
               return
            endif

            projrot_laps = n * s * angdif(stdlon,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        elseif(c6_maproj .eq. 'latlon')then ! latlon
            projrot_laps = 0.

        else
            write(6,*)
     1            'ERROR in projrot_latlon: unrecognized projection '      
     1               ,c6_maproj       
            stop

        endif

        projrot_latlon = projrot_laps

        istatus = 1 
        return
        end

        subroutine projrot_latlon_2d(rlat,rlon,ni,nj,projrot_laps
     1                                              ,istatus)
        include 'trigd.inc'

cdoc    1997 Steve Albers    Calculate map projection rotation, this is the
cdoc                         angle between the y-axis (grid north) and
cdoc                         true north. Units are degrees.
!
!                            projrot_laps = (true north value of wind direction
!                                          - grid north value of wind direction)
       
!       Added 8/4/2000 to make sure these are declared even if not passed in
        integer istatus
        real rlat(ni,nj), rlon(ni,nj), projrot_laps(ni,nj)

        real n

        save init
        data init/0/

        character*6  c6_maproj

        include 'grid_fname.cmn'

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus.ne.1)then
           print*,'Error returned from get_c6_maproj'
           return
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic

            call get_standard_longitude(polon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_longitude'
               return
            endif
            call get_standard_latitudes(stdlat,polat,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_latitudes'
               return
            endif

            call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif


            if(polat .eq. +90.)then
                projrot_laps(:,:) = polon - rlon(:,:)

            elseif(polat .eq. -90.)then
                projrot_laps(:,:) = rlon(:,:) - polon 

            else ! abs(polat) .ne. 90.
                if(grid_cen_lat .eq. polat .and. 
     1             grid_cen_lon .eq. polon)then ! grid centered on proj pole

                    if(init .eq. 0)then
                        write(6,*)
     1                   ' NOTE: local stereographic projection.'
                        write(6,*)
     1                   ' Using approximation for "projrot_laps",'
     1                  ,' accurate calculation not yet in place.'
                        init = 1
                    endif

                    rn = cosd(90.-polat)
                    do j = 1,nj
                    do i = 1,ni
                        projrot_laps(i,j) = rn * angdif(polon,rlon(i,j))
                    enddo ! i
                    enddo ! j

                elseif(.true.)then
                    if(init .eq. 0)then
                        write(6,*)' ERROR in projrot_laps: '
                        write(6,*)' This type of local'
     1                  ,' stereographic projection not yet supported.'
                        write(6,*)' Grid should be centered on'
     1                  ,' projection pole.'
                        init = 1
                    endif

                    projrot_laps = 0.
         
                else ! .false.
!                   Find dx/lat and dy/lat, then determine projrot_laps

                endif

            endif ! polat

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_latitudes'
               return
            endif
            call get_standard_longitude(slon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_longitude'
               return
            endif

            call lambert_parms(slat1,slat2,n,s,rconst)

            call get_standard_longitude(stdlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_c6_maproj'
               return
            endif

            do j = 1,nj
            do i = 1,ni
                projrot_laps(i,j) = n * s * angdif(stdlon,rlon(i,j))
            enddo ! i
            enddo ! j

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        elseif(c6_maproj .eq. 'latlon')then ! latlon
            projrot_laps = 0.

        else
            write(6,*)
     1         'ERROR in projrot_latlon_2d: unrecognized projection '
     1         ,c6_maproj       
            stop

        endif

        istatus = 1 
        return
        end



        function projrot_laps(rlon)
        include 'trigd.inc'

cdoc    This routine is being phased out. Please try to use 'projrot_latlon'.

!       1997 Steve Albers    Calculate map projection rotation, this is the
!                            angle between the y-axis (grid north) and
!                            true north. Units are degrees.
!
!                            projrot_laps = (true north value of wind direction
!                                          - grid north value of wind direction)
       
!       Added 8/4/2000 to make sure these are declared even if not passed in
        integer istatus
        real rlat, rlon

        real n

        save init
        data init/0/

        character*6  c6_maproj

        include 'grid_fname.cmn'

!       Difference between two angles, result is between -180. and +180.
        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus.ne.1)then
           print*,'Error returned from get_c6_maproj'
           return
        endif

        if(c6_maproj .eq. 'plrstr')then ! polar stereographic

            call get_standard_longitude(polon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_longitude'
               return
            endif
            call get_standard_latitudes(stdlat,polat,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_latitudes'
               return
            endif

            call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_grid_center'
               return
            endif


            if(polat .eq. +90.)then
                projrot_laps = polon - rlon

            elseif(polat .eq. -90.)then
                projrot_laps = rlon - polon 

            else ! abs(polat) .ne. 90.
                if(grid_cen_lat .eq. polat .and. 
     1             grid_cen_lon .eq. polon)then ! grid centered on proj pole

                    if(init .eq. 0)then
                        write(6,*)
     1                   ' NOTE: local stereographic projection.'
                        write(6,*)
     1                   ' Using approximation for "projrot_laps",'
     1                  ,' accurate calculation not yet in place.'
                        init = 1
                    endif

                    rn = cosd(90.-polat)
                    projrot_laps = rn * angdif(polon,rlon)      

                elseif(.true.)then
                    if(init .eq. 0)then
                        write(6,*)' ERROR in projrot_laps: '
                        write(6,*)' This type of local'
     1                  ,' stereographic projection not yet supported.'
                        write(6,*)' Grid should be centered on'
     1                  ,' projection pole.'
                        init = 1
                    endif

                    projrot_laps = 0.
         
                else ! .false.
!                   Find dx/lat and dy/lat, then determine projrot_laps

                endif

            endif ! polat

        elseif(c6_maproj .eq. 'lambrt')then ! lambert conformal

            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_latitudes'
               return
            endif
            call get_standard_longitude(slon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_standard_longitude'
               return
            endif

            call lambert_parms(slat1,slat2,n,s,rconst)

            call get_standard_longitude(stdlon,istatus)
            if(istatus.ne.1)then
               print*,'Error returned from get_c6_maproj'
               return
            endif

            projrot_laps = n * s * angdif(stdlon,rlon)

        elseif(c6_maproj .eq. 'merctr')then ! mercator
            projrot_laps = 0.

        elseif(c6_maproj .eq. 'latlon')then ! latlon
            projrot_laps = 0.

        else
            write(6,*)'ERROR in projrot_laps: unrecognized projection '
     1               ,c6_maproj       
            stop

        endif

        projrot_latlon = projrot_laps

!       istatus = 1 ! Not yet usable because of the entry
        return
        end


      subroutine check_domain(lat,lon,ni,nj,grid_spacing_m,intvl
     1                                                       ,istatus)

cdoc  This routine checks whether the lat/lon grid is consistent with
cdoc  map projection parameters as processed by latlon_to_rlapsgrid,
cdoc  and rlapsgrid_to_latlon. The grid size is also checked.
cdoc  This is a good sanity check of the NetCDF static file, namelist files,
cdoc  as well as various grid conversion routines.

!     1997 Steve Albers

      include 'trigd.inc'
      real pi, rpd
      parameter (pi=3.1415926535897932)
      parameter (rpd=pi/180.)

      real lat(ni,nj),lon(ni,nj)

      character*6  c6_maproj

      istatus = 1
      tolerance_m = 1000.

      diff_grid_max = 0.

      write(6,*)
      write(6,*)' subroutine check_domain: checking latlon_to_rlapsgrid'

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus.ne.1)then
         print*,'Error from get_c6_maproj'
         return
      endif

      do i = 1,ni,intvl
      do j = 1,nj,intvl
          call latlon_to_rlapsgrid(lat(i,j),lon(i,j),lat,lon,ni,nj
     1                                              ,ri,rj,istat)

          if(istat .ne. 1)then
              write(6,*)' Bad status from latlon_to_rlapsgrid'
              istatus = 0
              return
          endif

          diff_gridi = ri - float(i)
          diff_gridj = rj - float(j)
          diff_grid = sqrt(diff_gridi**2 + diff_gridj**2)

          if(diff_grid .gt. diff_grid_max)then
              diff_grid_max = diff_grid
              idmax = i
              jdmax = j
          endif

          diff_grid_max_m = diff_grid_max * grid_spacing_m

      enddo
      enddo

      write(6,*)' check_domain: max_diff (gridpoints) = ',diff_grid_max
     1         ,' at i/j',idmax,jdmax
      write(6,*)' check_domain: max_diff (approx m)   = '
     1                                                 ,diff_grid_max_m      

      if(diff_grid_max_m .gt. tolerance_m)then
          write(6,*)' WARNING: exceeded tolerance in check_domain'
     1               ,tolerance_m
          istatus = 0
      endif

!...........................................................................

      diff_ll_max = 0.

      write(6,*)' Checking rlapsgrid_to_latlon'

      do i = 1,ni,intvl
      do j = 1,nj,intvl
          ri = i
          rj = j
          call rlapsgrid_to_latlon(ri,rj,lat,lon,ni,nj
     1                                  ,rlat,rlon,istat)

          if(istat .ne. 1)then
              write(6,*)' Bad status from rlapsgrid_to_latlon'
              istatus = 0
              return
          endif

          diff_lli =  rlat - lat(i,j)
          diff_llj = (rlon - lon(i,j)) * cosd(lat(i,j))
          diff_ll = sqrt(diff_lli**2 + diff_llj**2)

          if(diff_ll .gt. diff_ll_max)then
              diff_ll_max = diff_ll
              idmax = i
              jdmax = j
          endif

          diff_ll_max_m = diff_ll_max * 110000. ! meters per degree

      enddo
      enddo

      write(6,*)' check_domain: max_diff (degrees) = ',diff_ll_max
     1         ,' at i/j',idmax,jdmax
      write(6,*)' check_domain: max_diff (approx m)   = ',diff_ll_max_m

      if(diff_ll_max_m .gt. tolerance_m)then
          write(6,*)' WARNING: exceeded tolerance in check_domain'
     1               ,tolerance_m
          istatus = 0
      endif

!...........................................................................

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1    ' Error calling get_earth_radius from check_domain'       
          return
      endif

      icen = ni/2+1
      jcen = nj/2+1

      diff_lat =  lat(icen,jcen+1) - lat(icen,jcen-1)
      diff_lon = (lon(icen,jcen+1) - lon(icen,jcen-1)) 
     1                        * cosd(lat(icen,jcen))

      dist = sqrt((diff_lat)**2 + (diff_lon)**2) / 2. * rpd * erad   

      if(abs(lat(icen,jcen)) .lt. 89.)then ! should be reasonably accurate
          write(6,*)
     1   ' measured grid spacing on earths surface at domain center is:'       
     1     ,dist       
      endif

!...........................................................................

      if(c6_maproj .ne. 'latlon')then
          call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen)
     1                                      ,dist_calc,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' Error calling get_grid_spacing_actual from check_domain'       
              return
          endif

          write(6,*)
     1 ' calculated grid spacing on earths surface at domain center is:'       
     1     ,dist_calc

      endif

!...........................................................................

      call latlon_to_xy(lat(1,1),lon(1,1),erad,x1,y1)
      call latlon_to_xy(lat(1,2),lon(1,2),erad,x2,y2)

      dist = sqrt((x2-x1)**2 + (y2-y1)**2)

      write(6,*)
     1 ' grid spacing on projection plane using "latlon_to_xy" is:'      
     1 ,dist
      
!...........................................................................

      if(c6_maproj .eq. 'lambrt')then
          call get_standard_latitudes(std_lat1,std_lat2,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' Error calling get_standard_latitudes from check_domain'
              return
          endif

          ave_lat = (std_lat1 + std_lat2) / 2.

          call get_grid_spacing_actual(ave_lat,lon(icen,jcen)
     1                                      ,dist_calc,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' Error calling get_grid_spacing_actual from check_domain'       
              return
          endif

          write(6,*)' calculated grid spacing on '
     1      ,'earths surface at average lambert latitude:'       
     1     ,ave_lat,dist_calc

      endif

!...........................................................................

      write(6,*)

      return
      end

      subroutine lambert_parms(slat1,slat2,n_out,s_out,rconst_out)

      include 'trigd.inc'
      real n,n_out

!     We only have to do the calculations once since the inputs are constants
      data init/0/
      save init,n,s,rconst 

      if(init .eq. 0)then ! Calculate saved variables
          if(slat1 .ge. 0)then
              s = +1.
          else
              s = -1.
          endif

          colat1 = 90. - s * slat1
          colat2 = 90. - s * slat2

          if(slat1 .eq. slat2)then ! tangent lambert
              n = cosd(90.-s*slat1)
              rconst =       tand(colat1)    / tand(colat1/2.)**n

          else                     ! two standard latitudes
              n = alog(cosd(slat1)/cosd(slat2))/
     1            alog(tand(45.-s*slat1/2.)/tand(45.-s*slat2/2.))
              rconst =      (sind(colat1)/n) / tand(colat1/2.)**n

          endif

          init = 1

      endif

      n_out = n
      s_out = s
      rconst_out = rconst

      return
      end

      subroutine get_grid_spacing_actual(rlat,rlon
     1                                  ,grid_spacing_actual_m,istatus)

cdoc  Calculate actual grid spacing at any given lat/lon location

      character*6 c6_maproj

      call get_standard_latitudes(slat1,slat2,istatus)
      if(istatus .ne. 1)then
          return
      endif

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1    ' Error calling get_grid_spacing from get_grid_spacing_actual'       
          return
      endif

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus .ne. 1)then
          return
      endif

      if(c6_maproj .eq. 'plrstr')then
          call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
      else
          grid_spacing_proj_m = grid_spacing_m
      endif

      call get_sigma(rlat,rlon,sigma,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1    ' Error calling get_sigma from get_grid_spacing_actual'       
          return
      endif

      grid_spacing_actual_m = grid_spacing_proj_m / sigma

      return
      end

      subroutine get_grid_spacing_actual_xy(rlat,rlon
     1                        ,grid_spacing_actual_mx
     1                        ,grid_spacing_actual_my
     1                        ,istatus)

cdoc  Calculate actual grid spacing (x,y directions) at any given lat/lon 
cdoc  location. This works for conformal or 'latlon' grids

      include 'trigd.inc'

      character*6 c6_maproj

      call get_standard_latitudes(slat1,slat2,istatus)
      if(istatus .ne. 1)then
          return
      endif

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1 ' Error calling get_grid_spacing from get_grid_spacing_actual_xy'       
          return
      endif

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus .ne. 1)then
          return
      endif

      if(c6_maproj .ne. 'latlon')then
          if(c6_maproj .eq. 'plrstr')then
              call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                     ,grid_spacing_proj_m)
          else
              grid_spacing_proj_m = grid_spacing_m
          endif

          call get_sigma(rlat,rlon,sigma,istatus)
          if(istatus .ne. 1)then
              write(6,*)
     1        ' Error calling get_sigma from get_grid_spacing_actual'       
              return
          endif

          grid_spacing_actual_mx = grid_spacing_proj_m / sigma
          grid_spacing_actual_my = grid_spacing_proj_m / sigma

      else
          grid_spacing_actual_mx = grid_spacing_m * cosd(rlat)
          grid_spacing_actual_my = grid_spacing_m 

      endif

      return
      end


        subroutine get_grid_spacing_cen(grid_spacing_cen_m,istatus)

cdoc    Calculate actual grid spacing at the center of the domain
cdoc    If we have a lat/lon domain the Y direction spacing will be used

        call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)
        if(istatus .ne. 1)return

        call get_grid_spacing_actual_xy(grid_cen_lat,grid_cen_lon
     1                                 ,grid_spacing_actual_mx
     1                                 ,grid_spacing_actual_my
     1                                 ,istatus)

        grid_spacing_cen_m = grid_spacing_actual_my
 
        return
        end


      subroutine get_ps_parms(slat1,slat2,grid_spacing_m                ! I
     1                       ,phi0,grid_spacing_proj_m)                 ! O

!     1998 Steve Albers

      include 'trigd.inc'

      logical l_secant
      data l_secant /.true./

!     Secant projections are described in "Principles of Meteorological 
!     Analysis", Saucier, p. 33. 'phi_std' is the value of phi at the 
!     "standard latitude" as specified by input parameter. 'phi0' is the value 
!     of phi on the actual projection plane utilized for internal map 
!     projection calculations. 'phi0' represents the standard latitude in the
!     more generic sense.

!     We will eventually use the secant projection assumption (unless we run 
!     into software problems)

!     Projection is tangent to earth's surface only if phi0 = 90.

      if(slat2 .eq. +90.)then     ! Projection pole is at geographic north pole
          phi_std = slat1       

      elseif(slat2 .eq. -90.)then ! Projection pole is at geographic south pole
          phi_std = -slat1

      else
          phi_std = +90.        ! We ignore standard lat for local sterographic
                                ! No need for this to be not equal to +90.
      endif

!     Calculate grid_spacing_proj_m: grid spacing in the projection plane
      if(l_secant)then
          grid_spacing_proj_m = grid_spacing_m ! equal to parameter value 
                                               ! in namelist files
          phi0 = phi_std

      else ! tangent projection
          grid_spacing_proj_m = grid_spacing_m * 
     1                          2. / (1. + sind(phi_std))
          phi0 = 90.

      endif

      return
      end
