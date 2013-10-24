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
	subroutine solar_normal(ni,nj,topo,dx,dy,lat,lon
     1                         ,sol_alt,sol_azi,alt_norm)

c       Compute solar altitude normal to the terrain

        include 'trigd.inc'

        angleunitvectors(a1,a2,a3,b1,b2,b3) = acosd(a1*b1+a2*b2+a3*b3)

	real topo(ni,nj)                 ! I Terrain elevation (m)
        real lat(ni,nj)                  ! I Lat (deg)
        real lon(ni,nj)                  ! I Lon (deg)
        real sol_alt(ni,nj)              ! I Solar Altitude (deg)
        real sol_azi(ni,nj)              ! I Solar Azimuth (deg)
        real alt_norm(ni,nj)             ! O Solar Alt w.r.t. terrain normal

	real dx(ni,nj)                   ! I Grid spacing in X direction (m)
	real dy(ni,nj)                   ! I Grid spacing in Y direction (m)

        real rot(ni,nj)                  ! L Rotation Angle (deg)
    
        write(6,*)' Subroutine solar_normal'

!       call get_grid_spacing_array(lat,lon,ni,nj,dx,dy)

!       Default value
        alt_norm = sol_alt

        call projrot_latlon_2d(lat,lon,ni,nj,rot,istatus)

	do j=2,nj-1
	do i=2,ni-1

!           Determine terrain slope
            dterdx = (topo(i,j)+topo(i,j-1)-topo(i-1,j)-topo(i-1,j-1)
     1                ) * .5 / dx(i,j)
            dterdy = (topo(i,j)+topo(i-1,j)-topo(i-1,j-1)-topo(i,j-1)
     1                ) * .5 / dy(i,j)

            terrain_slope = sqrt(dterdx**2 + dterdy**2)

            if(terrain_slope .gt. .001)then ! machine/terrain epsilon threshold

!             Direction cosines of terrain normal
              dircos_tx = dterdx / (sqrt(dterdx**2 + 1.))
              dircos_ty = dterdy / (sqrt(dterdy**2 + 1.))
              dircos_tz = 1.0 / sqrt(1.0 + terrain_slope**2)

              sol_azi_grid = sol_azi(i,j) - rot(i,j) 

              dircos_sx = sind(sol_azi_grid)
              dircos_sy = cosd(sol_azi_grid)
              dircos_sz = sind(sol_alt(i,j))

              result = angleunitvectors(dircos_tx,dircos_ty,dircos_tz
     1                                 ,dircos_sx,dircos_sy,dircos_sz)

              alt_norm(i,j) = 90. - result 
            
              if(i .eq. ni/2 .and. j .eq. nj/2)then
                write(6,*)'solar alt/az, dterdx, dterdy, alt_norm',
     1             sol_alt(i,j),sol_azi(i,j),dterdx,dterdy,alt_norm(i,j)        
                write(6,*)' dircos_t ',dircos_tx,dircos_ty,dircos_tz
                write(6,*)' dircos_s ',dircos_sx,dircos_sy,dircos_sz
                write(6,*)' rot = ',rot(i,j)
              endif

            else 
              alt_norm(i,j) = sol_alt(i,j) ! terrain virtually flat

              if(i .eq. ni/2 .and. j .eq. nj/2)then
                write(6,*)'solar alt/az, alt_norm',
     1             sol_alt(i,j),sol_azi(i,j),alt_norm(i,j)        
              endif

            endif

	enddo !i
	enddo !j

	return
	end


