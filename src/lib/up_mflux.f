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
	subroutine up_mflux(ni,nj,nk,topo,dx,dy
     1                     ,u_3d,v_3d,tpw_2d,upslope_flux
     1                     ,ht_3d,r_missing_data)

c       Compute upslope moisture flux (using conventions in the PSD flux tool) 

	real topo(ni,nj)                 ! I Terrain elevation (m)
	real dx(ni,nj)                   ! I Grid spacing in X direction (m)
	real dy(ni,nj)                   ! I Grid spacing in Y direction (m)
        real u_3d(ni,nj,nk)              ! I U wind component in Grid X direction
        real v_3d(ni,nj,nk)              ! I V wind component in Grid Y direction
        real ht_3d(ni,nj,nk)             ! I
        real tpw_2d(ni,nj)               ! I
        real upslope_flux(ni,nj)         ! O

        real u_2d(ni,nj)
        real v_2d(ni,nj)
        
	do j=2,nj-1
	do i=2,ni-1

!         Controlling layer (defined relative to topography)
          ht_lo  = 1200. + topo(i,j)
          ht_hi  = 3200. + topo(i,j)

          if(topo(i,j) .le. ht_lo)then

!           Calculate mean wind over the layer       
            rk_lo = rlevel_of_field(ht_lo,ht_3d,ni,nj,nk,i,j,istatus) 
            rk_hi = rlevel_of_field(ht_hi,ht_3d,ni,nj,nk,i,j,istatus)  

            k_lo = int(rk_lo)
            k_hi = int(rk_hi)

!           Lower part
            frac_lo = rk_lo - float(k_lo)
            u_lo = u_3d(i,j,k_lo)*(1.-frac_lo)+u_3d(i,j,k_lo+1)*frac_lo
            v_lo = v_3d(i,j,k_lo)*(1.-frac_lo)+v_3d(i,j,k_lo+1)*frac_lo     

            ubar_lyr = (u_lo + u_3d(i,j,k_lo+1)) / 2.
            vbar_lyr = (v_lo + v_3d(i,j,k_lo+1)) / 2.
 
            ubar_sum = ubar_lyr
            vbar_sum = vbar_lyr

            sumk = 1. - frac_lo

!           Middle part
            do k = k_lo+1,k_hi
                ubar_lyr = (u_3d(i,j,k) + u_3d(i,j,k+1)) / 2.
                vbar_lyr = (v_3d(i,j,k) + v_3d(i,j,k+1)) / 2. 
                ubar_sum = ubar_sum + ubar_lyr
                vbar_sum = vbar_sum + vbar_lyr
                sumk = sumk + 1.0
            enddo ! k

!           Upper part  
            frac_hi = rk_hi - float(k_hi)
            u_hi = u_3d(i,j,k_hi)*(1.-frac_hi)+u_3d(i,j,k_hi+1)*frac_hi
            v_hi = v_3d(i,j,k_hi)*(1.-frac_hi)+v_3d(i,j,k_hi+1)*frac_hi     

            ubar_lyr = (u_3d(i,j,k_hi) + u_hi) / 2.
            vbar_lyr = (v_3d(i,j,k_hi) + v_hi) / 2.

            ubar_sum = ubar_sum + ubar_lyr
            vbar_sum = vbar_sum + vbar_lyr

            sumk = sumk + frac_hi

!           Divide to get the means
            ubar = ubar_sum / sumk
            vbar = vbar_sum / sumk            

!           Determine terrain slope
	    dterdx = (topo(i,j)+topo(i,j-1)-topo(i-1,j)-topo(i-1,j-1)
     1                ) * .5 / dx(i,j)
	    dterdy = (topo(i,j)+topo(i-1,j)-topo(i-1,j-1)-topo(i,j-1)
     1                ) * .5 / dy(i,j)

            terrain_slope = sqrt(dterdx**2 + dterdy**2)

            if(terrain_slope .gt. .001)then ! machine/terrain epsilon threshold
!             Calculate upslope wind component (m/s)
!             This is normalized by the terrain slope
              dvh = (ubar * dterdx + vbar * dterdy) / terrain_slope

!             Calculate upslope moisture flux (m**2/s)
	      upslope_flux(i,j) = dvh * tpw_2d(i,j)

            else ! assume cos(theta) = 1, so we just want moisture flux
              dvh = sqrt(ubar**2 + vbar**2)

	      upslope_flux(i,j) = dvh * tpw_2d(i,j)

            endif

          else
            upslope_flux(i,j) = r_missing_data

          endif
	enddo !i
	enddo !j

	call bounds(upslope_flux,ni,nj)

	return
	end

