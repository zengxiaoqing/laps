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
	subroutine up_mflux(ni,nj,nk,topo,u_3d,v_3d,q_3d,upslope_flux
     1                     ,pres_3d,r_missing_data)

c       Compute upslope moisture flux 

	real topo(ni,nj)                 ! I
        real u_3d(ni,nj,nk)              ! I
        real v_3d(ni,nj,nk)              ! I
        real q_3d(ni,nj,nk)              ! I
        real pres_3d(ni,nj,nk)           ! I
        real upslope_flux(ni,nj)         ! O

        real u_2d(ni,nj)
        real v_2d(ni,nj)
        real q_2d(ni,nj)

        value = 900. * 100.

!       Interpolate 3D fields to 900mb pressure level
	do j=1,nj
	do i=1,ni
            rk = rlevel_of_field(value,pres_3d,ni,nj,nk,i,j,istatus)
            level = nint(rk)
            u_2d = u_3d(i,j,level)
            v_2d = v_3d(i,j,level)
            q_2d = q_3d(i,j,level)
        enddo ! i
        enddo ! j

	do j=2,nj-1
	do i=2,ni-1
	   dterdx = (topo(i,j)+topo(i,j-1)-topo(i-1,j)-topo(i-1,j-1)
     1                ) * .5 / dx(i,j)
	   dterdy = (topo(i,j)+topo(i-1,j)-topo(i-1,j-1)-topo(i,j-1)
     1                ) * .5 / dy(i,j)

	   ubar = u_2d(i,j)
	   vbar = v_2d(i,j)

           if(ubar .ne. r_missing_data .and. 
     1        vbar .ne. r_missing_data      )then     

!              Calculate upslope wind component (m/s)
	       dvh = ubar * dterdx + vbar * dterdy

!              Calculate upslope moisture flux (m/s)
	       upslope_flux(i,j) = dvh * q(i,j)

           else
	       upslope_flux(i,j) = r_missing_data

           endif

	enddo !i
	enddo !j

	call bounds(upslope_flux,ni,nj)

	return
	end

