cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
        subroutine comp_laps_vr(grid_ra_vel,u,v,ni,nj,nk,r_missing_data
     1  ,cgrid,rms
     1  ,nx_r,ny_r,ioffset,joffset                                  ! I
     1  ,lat,lon,rlat_radar,rlon_radar,rheight_radar)


        real grid_ra_vel(nx_r,ny_r,nk)
        dimension u(ni,nj,nk),v(ni,nj,nk)
        real lat(ni,nj),lon(ni,nj)
        real lat_grid,lon_grid
        integer ioffset,joffset                 

        character*(*) cgrid

        nobs = 0
        residual = 0.
        bias_sum = 0.

        write(6,2)
2       format(/'      Comparing Radial Velocities to LAPS'/
     1         6x,'   i   j   k   radar   laps  diff   ')

        do k = 1,nk
        height_grid = 0.
        do jo = 1,ny_r
        do io = 1,nx_r

          i = io + ioffset
          j = jo + joffset

          if(i .ge. 1 .and. i .le. ni .and. 
     1       j .ge. 1 .and. j .le. nj)then     

            if(         u(i,j,k) .ne. r_missing_data
     1  .and. grid_ra_vel(io,jo,k) .ne. r_missing_data)then
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

                diff = r_radar - grid_ra_vel(io,jo,k)
                residual = residual + diff ** 2
                bias_sum = bias_sum + diff
                if(float(nobs)/40. .eq. nobs/40 .or. nobs .le. 25)then       
                    write(6,101)nobs,i,j,k,grid_ra_vel(io,jo,k),r_radar
     1                         ,diff
101                 format(1x,i5,3i4,3f7.1,3i6)
                endif

            endif

          endif

        enddo ! i
        enddo ! j
        enddo ! k

        if(nobs .gt. 0)then
            rms = sqrt(residual/nobs)
            bias = bias_sum / nobs
        else
            rms = 0.
            bias = 0.
        endif

        write(6,102)cgrid,nobs,bias,rms
 102    format(' BIAS/RMS between radial velocities & ',a,' = '
     1        ,i5,' . . ',f6.1,' . . ',f6.1)

        return

        end

