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
        subroutine plot_grid_2d(interval,size,imax,jmax,lat,lon)

!       1997 Steve Albers

        real*4 lat(imax,jmax),lon(imax,jmax)

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


