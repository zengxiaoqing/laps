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
        subroutine comp_laps_vr(grid_ra_vel,u,v,ni,nj,nk,r_missing_data
     1  ,rms,lat,lon,rlat_radar,rlon_radar,rheight_radar)


        real*4 grid_ra_vel(ni,nj,nk)
        dimension u(ni,nj,nk),v(ni,nj,nk)
        real*4 lat(ni,nj),lon(ni,nj)
        real lat_grid,lon_grid

        nobs = 0
        residual = 0.

        write(6,2)
2       format(/'      Comparing Radial Velocities to LAPS'/
     1  1x,'   i   j   k  radar  laps  diff   ')

        do k = 1,nk
        height_grid = 0.
        do j = 1,nj
        do i = 1,ni

            if(         u(i,j,k) .ne. r_missing_data
     1  .and. grid_ra_vel(i,j,k) .ne. r_missing_data)then
                nobs = nobs + 1
                lat_grid = lat(i,j)
                lon_grid = lon(i,j)
                call latlon_to_radar(lat_grid,lon_grid,height_grid
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)
                call uvtrue_to_radar(u(i,j,k),
     1                       v(i,j,k),
     1                       t_radar,
     1                       r_radar,
     1                       azimuth)

                diff = r_radar - grid_ra_vel(i,j,k)
                residual = residual + diff ** 2
                if(float(nobs)/40. .eq. nobs/40)then
                    write(6,101)i,j,k,grid_ra_vel(i,j,k),r_radar,
     1                  diff
101                 format(1x,3i4,3f7.1,3i6)
                endif

            endif

        enddo ! i
        enddo ! j
        enddo ! k

        if(nobs .gt. 0)then
            rms = sqrt(residual/nobs)
        else
            rms = 0.
        endif

        write(6,*)' RMS between radial velocities & LAPS = ',nobs,rms
        write(15,*)' RMS between radial velocities & LAPS = ',nobs,rms

        return

        end

