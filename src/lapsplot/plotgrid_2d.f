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

        real*4 lat(imax,jmax),lon(imax,jmax)

        common /supmp6/ umin,umax,vmin,vmax

        size = 61. / 200.

        write(6,*)' plot_grid_2d: UMIN/UMAX/VMIN/VMAX'
     1                           ,umin,umax,vmin,vmax

        do j = 1,jmax
        do i = 1,imax

!           ri = i
!           rj = j

            call latlon_to_uv(lat(i,j),lon(i,j),ri,rj,istatus)
           
            call plot_gridpt(ri,rj,imax,jmax,size)

        enddo ! i
        enddo ! j

        return
        end

      subroutine plot_gridpt(ri,rj,imax,jmax,relsize)

!     1997 Steve Albers
!     Note that the umin/umax come from the NCAR graphics subroutines while
!     the u/v values for the individual grid points come from the laps library.
!     This plot is therefore a good test of the consistency of these two
!     methods of calculating u and v.

      include 'lapsparms.cmn'

      common /supmp6/ umin,umax,vmin,vmax

!     call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!     write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
 1234 format(1x,4i5,4e12.4,i4)

!     This tries to keep the same size of barbs relative to the grid points

      call get_border(imax,jmax,x_1,x_2,y_1,y_2)
!     call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))
      call set(x_1,x_2,y_1,y_2,umin,umax,vmin,vmax)

      relsize = 0.3 * (umax - umin) / float(imax)

      du = relsize

!     Plot a '+'
      u = ri
      u1 = ri - du
      u2 = ri + du
      v = rj
      v1 = rj - du
      v2 = rj + du
        
      CALL LINE(U1,V,U2,V)
      CALL LINE(U,V1,U,V2)

    1 continue
      return
      end


