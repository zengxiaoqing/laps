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
        subroutine plot_grid_2d(imax,jmax,lat,lon)

!       1997 Steve Albers

        real lat(imax,jmax),lon(imax,jmax)

        common /supmp6/ umin_n,umax_n,vmin_n,vmax_n

        size = 61. / 200.

        write(6,*)' subroutine plot_grid_2d...'
        write(6,*)' plotting with LAPS lat/lons converted to NCAR U/V'       

        do j = 1,jmax
        do i = 1,imax

!           call latlon_to_uv(lat(i,j),lon(i,j),u_l,v_l,istatus)
            call MAPTRN(lat(i,j),lon(i,j),u_n,v_n)
            call plot_gridpt(u_n,v_n,imax,jmax,size)

        enddo ! i
        enddo ! j


        write(6,11)umin_n,umax_n,vmin_n,vmax_n
11      format('  NCARG UMIN/UMAX/VMIN/VMAX',4f10.5)

        call latlon_to_uv(lat(1,1),lon(1,1),umin_l,vmin_l,istatus)
        call latlon_to_uv(lat(imax,jmax),lon(imax,jmax),umax_l,vmax_l
     1                                                   ,istatus)       

        write(6,21)umin_l,umax_l,vmin_l,vmax_l
21      format('  LAPS  UMIN/UMAX/VMIN/VMAX',4f10.5)

        write(6,31)umin_l/umin_n,umax_l/umax_n
     1            ,vmin_l/vmin_n,vmax_l/vmax_n
31      format('  ratio UMIN/UMAX/VMIN/VMAX',4f10.5)

        aspect_n = (vmax_n - vmin_n) / (umax_n - umin_n)
        aspect_l = (vmax_l - vmin_l) / (umax_l - umin_l)

        write(6,41)aspect_n,aspect_l,aspect_n/aspect_l
 41     format('  aspect NCARG/LAPS/ratio: ',3f10.5)

!       Call Check_domain for good measure
        call get_grid_spacing(grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' grid_spacing = ',grid_spacing_m

        call check_domain(lat,lon,imax,jmax,grid_spacing_m,1,istatus)
        write(6,*)' check_domain: istatus = ',istatus

        return
        end

        subroutine plot_gridpt(u_n,v_n,imax,jmax,relsize)

!       1997 Steve Albers
!       Note that the umin/umax come from the NCAR graphics subroutines while
!       the u/v values for the individual grid points come from the laps library.
!       This plot is therefore a good test of the consistency of these two
!       methods of calculating u and v.

        include 'lapsparms.cmn'

        common /supmp6/ umin_n,umax_n,vmin_n,vmax_n

!       This tries to keep the same size of barbs relative to the grid points

        call get_border(imax,jmax,x_1,x_2,y_1,y_2)
!       call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))
        call set(x_1,x_2,y_1,y_2,umin_n,umax_n,vmin_n,vmax_n)

        relsize = 0.3 * (umax_n - umin_n) / float(imax)

        du = relsize

!       Plot a '+'
        u = u_n
        u1 = u - du
        u2 = u + du

        v = v_n
        v1 = v - du
        v2 = v + du
        
        CALL LINE(U1,V,U2,V)
        CALL LINE(U,V1,U,V2)

    1 continue
      return
      end


