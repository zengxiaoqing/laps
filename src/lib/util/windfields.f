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
      subroutine get_grid_spacing_array(rlat,rlon,ni,nj
     1                                 ,grid_spacing_actual_mx
     1                                 ,grid_spacing_actual_my)

cdoc  Calculate actual grid spacing (x,y directions) at any given lat/lon 
cdoc  location. This works for conformal or 'latlon' grids

      include 'trigd.inc'

      real rlat(ni,nj),rlon(ni,nj)                                     ! I
      real grid_spacing_actual_mx(ni,nj),grid_spacing_actual_my(ni,nj) ! O

      character*6 c6_maproj

      write(6,*)' Subroutine get_grid_spacing_array'

      call get_standard_latitudes(slat1,slat2,istatus)
      if(istatus .ne. 1)then
          stop
      endif

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)
     1 ' Error calling get_grid_spacing from get_grid_spacing_actual_xy'       
          stop
      endif

      call get_c6_maproj(c6_maproj,istatus)
      if(istatus .ne. 1)then
          stop
      endif

      if(c6_maproj .ne. 'latlon')then
          if(c6_maproj .eq. 'plrstr')then
              call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                     ,grid_spacing_proj_m)
          else
              grid_spacing_proj_m = grid_spacing_m
          endif

          do i = 1,ni
          do j = 1,nj
              call get_sigma(rlat(i,j),rlon(i,j),sigma,istatus)
              if(istatus .ne. 1)then
                  write(6,*)
     1        ' Error calling get_sigma from get_grid_spacing_actual'       
                  stop
              endif

              grid_spacing_actual_mx(i,j) = grid_spacing_proj_m / sigma       
              grid_spacing_actual_my(i,j) = grid_spacing_proj_m / sigma
          enddo ! j
          enddo ! i

      else
          do i = 1,ni
          do j = 1,nj
              grid_spacing_actual_mx(i,j) = grid_spacing_m 
     1                                    * cosd(rlat(i,j))
              grid_spacing_actual_my(i,j) = grid_spacing_m 
          enddo ! j
          enddo ! i

      endif

      return
      end

c
c=====  Here are John's subroutines...(abandon hope, ye who enter)
c
	subroutine vortdiv(u,v,vort,div,imax,jmax,dx,dy)
c this routine computes vorticity and divergence from u and v winds
c using centered differences.
	real u(imax,jmax),v(imax,jmax),vort(imax,jmax)
	real div(imax,jmax),dx(imax,jmax),dy(imax,jmax)
c
	do j=2,jmax-1
	do i=2,imax-1
	    div(i,j) = (u(i,j-1) - u(i-1,j-1)) / dx(i,j)
     &               + (v(i-1,j) - v(i-1,j-1)) / dy(i,j)
	    vort(i,j) = (v(i,j) - v(i-1,j)) / dx(i,j)
     &                - (u(i,j) - u(i,j-1)) / dy(i,j)
	enddo !i
	enddo !j
	call bounds(div,imax,jmax)
	call bounds(vort,imax,jmax)
c
	return
	end
c
c
	subroutine bounds(x,imax,jmax)
c
c.....	Routine to fill in the boundaries of an array.  Just uses the
c.....	interior points for now.
c
	real x(imax,jmax)
c
	do i=1,imax
	  x(i,1) = x(i,2)
	  x(i,jmax) = x(i,jmax-1)
	enddo !i
	do j=1,jmax
	  x(1,j) = x(2,j)
	  x(imax,j) = x(imax-1,j)
	enddo !j
c
	x(1,1) = x(2,2)
	x(1,jmax) = x(2,jmax-1)
	x(imax,1) = x(imax-1,2)
	x(imax,jmax) = x(imax-1,jmax-1)
c
	return
	end
c

