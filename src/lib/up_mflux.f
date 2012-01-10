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
	subroutine up_mflux(ni,nj,nk,topo,ldf,dx,dy
     1                     ,u_3d,v_3d,tpw_2d,upslope_flux
     1                     ,ht_3d,r_missing_data)

c       Compute upslope moisture flux (using conventions in the PSD flux tool) 

	real topo(ni,nj)                 ! I Terrain elevation (m)
	real ldf(ni,nj)                  ! I Land Fraction
	real dx(ni,nj)                   ! I Grid spacing in X direction (m)
	real dy(ni,nj)                   ! I Grid spacing in Y direction (m)
        real u_3d(ni,nj,nk)              ! I U wind component in Grid X direction
        real v_3d(ni,nj,nk)              ! I V wind component in Grid Y direction
        real ht_3d(ni,nj,nk)             ! I
        real tpw_2d(ni,nj)               ! I
        real upslope_flux(ni,nj)         ! O

        real ht_lower(ni,nj)             ! L
        real ht_upper(ni,nj)             ! L
        real umean_2d(ni,nj)             ! L
        real vmean_2d(ni,nj)             ! L

        write(6,*)' Subroutine up_mflux'

        ht_lower(:,:)  =  750. + topo(:,:)
        ht_upper(:,:)  = 1250. + topo(:,:)

	call mean_wind_hlyr(ni,nj,nk,ht_lower,ht_upper
     1                     ,u_3d,v_3d,ht_3d,umean_2d,vmean_2d
     1                     ,r_missing_data)

	do j=2,nj-1
	do i=2,ni-1

          if(topo(i,j) .lt. ht_lower(i,j))then

            ubar = umean_2d(i,j)
            vbar = vmean_2d(i,j)

!           Determine terrain slope
	    dterdx = (topo(i,j)+topo(i,j-1)-topo(i-1,j)-topo(i-1,j-1)
     1                ) * .5 / dx(i,j)
	    dterdy = (topo(i,j)+topo(i-1,j)-topo(i-1,j-1)-topo(i,j-1)
     1                ) * .5 / dy(i,j)

            terrain_slope = sqrt(dterdx**2 + dterdy**2)

!           if(terrain_slope .gt. .001)then ! machine/terrain epsilon threshold

            if(ldf(i,j) .lt. .01 .and. abs(topo(i,j)) .lt. 10.)then ! ocean

!             Assume cos(theta) = 1, so we just want moisture flux

              dvh = sqrt(ubar**2 + vbar**2)

	      upslope_flux(i,j) = dvh * tpw_2d(i,j)

            else 

!             Calculate upslope wind component (m/s)
!             This is normalized by the terrain slope
              if(terrain_slope .gt. 0.)then
                  dvh = (ubar * dterdx + vbar * dterdy) / terrain_slope
              else
                  dvh = 0.
              endif

!             Calculate upslope moisture flux (m**2/s)
	      upslope_flux(i,j) = dvh * tpw_2d(i,j)

            endif

            if(i .eq. 10 .or. abs(upslope_flux(i,j)) .gt. 1e10)then ! write debugging info
              if(abs(upslope_flux(i,j)) .lt. 1e10)then
                  write(6,2)j,ubar,vbar,dvh,terrain_slope,tpw_2d(i,j)
     1                       ,upslope_flux(i,j)
2  	          format(i5,6f10.4)
              else
                  write(6,3)i,j,ubar,vbar,dvh,terrain_slope,tpw_2d(i,j)
     1                       ,upslope_flux(i,j)
3 		  format(' WARNING: large UMF ',2i5,6e13.4)
              endif
            endif 

          else
            upslope_flux(i,j) = r_missing_data

          endif
	enddo !i
	enddo !j

	call bounds(upslope_flux,ni,nj)

	return
	end


c
c
	subroutine mean_wind_hlyr(ni,nj,nk,ht_lower,ht_upper
     1                           ,u_3d,v_3d,ht_3d,umean_2d,vmean_2d
     1                           ,r_missing_data)       

c       Compute mean wind over a 2D height layer 

        real ht_lower(ni,nj)         ! I
        real ht_upper(ni,nj)         ! I
        real u_3d(ni,nj,nk)          ! I U wind component in Grid X direction
        real v_3d(ni,nj,nk)          ! I V wind component in Grid Y direction
        real ht_3d(ni,nj,nk)         ! I

        real umean_2d(ni,nj)         ! O
        real vmean_2d(ni,nj)         ! O

        umean_2d = r_missing_data
        vmean_2d = r_missing_data
        
	do j=1,nj
	do i=1,ni

!         Controlling layer (defined relative to topography)
          ht_lo  = ht_lower(i,j)
          ht_hi  = ht_upper(i,j)

          if(.true.)then

!           Calculate mass weighted mean wind over the height layer       
            rk_lo = rlevel_of_field(ht_lo,ht_3d(i,j,:),1,1,nk,1,1
     1                             ,istatus)        
            if(istatus .ne. 1)then
                write(6,*)' Bad status in rlevel_of_field for lo ',i,j
                goto 900
            endif

            rk_hi = rlevel_of_field(ht_hi,ht_3d(i,j,:),1,1,nk,1,1  
     1                             ,istatus)        
            if(istatus .ne. 1)then
                write(6,*)' Bad status in rlevel_of_field for hi ',i,j  
                goto 900
            endif

!           rk_lo = rlevel_of_logfield(ht_lo,ht_3d,ni,nj,nk,i,j,istatus)        
!           rk_hi = rlevel_of_logfield(ht_hi,ht_3d,ni,nj,nk,i,j,istatus)  

            k_lo = int(rk_lo)
            k_hi = int(rk_hi)

!           Lower part
            frac_lo = rk_lo - float(k_lo)
            u_lo = u_3d(i,j,k_lo)*(1.-frac_lo)+u_3d(i,j,k_lo+1)*frac_lo
            v_lo = v_3d(i,j,k_lo)*(1.-frac_lo)+v_3d(i,j,k_lo+1)*frac_lo     

            ubar_llyr = (u_lo + u_3d(i,j,k_lo+1)) / 2.
            vbar_llyr = (v_lo + v_3d(i,j,k_lo+1)) / 2.
 
            ubar_sum = ubar_llyr * (1. - frac_lo)
            vbar_sum = vbar_llyr * (1. - frac_lo)

            sumk = 1. - frac_lo

!           Middle part
            do k = k_lo+1,k_hi-1
                ubar_lyr = (u_3d(i,j,k) + u_3d(i,j,k+1)) / 2.
                vbar_lyr = (v_3d(i,j,k) + v_3d(i,j,k+1)) / 2. 
                ubar_sum = ubar_sum + ubar_lyr
                vbar_sum = vbar_sum + vbar_lyr
                sumk = sumk + 1.0
!               if(i .eq. 10)then ! write debugging info
!                   write(6,1)k,u_3d(i,j,k),u_3d(i,j,k+1),ubar_lyr
!    1                                                   ,ubar_sum    
!1		    format(30x,i5,4f10.4)
!               endif
            enddo ! k

!           Upper part  
            frac_hi = rk_hi - float(k_hi)
            u_hi = u_3d(i,j,k_hi)*(1.-frac_hi)+u_3d(i,j,k_hi+1)*frac_hi
            v_hi = v_3d(i,j,k_hi)*(1.-frac_hi)+v_3d(i,j,k_hi+1)*frac_hi     

            ubar_hlyr = (u_3d(i,j,k_hi) + u_hi) / 2.
            vbar_hlyr = (v_3d(i,j,k_hi) + v_hi) / 2.

            ubar_sum = ubar_sum + (ubar_hlyr * frac_hi)
            vbar_sum = vbar_sum + (vbar_hlyr * frac_hi)

            sumk = sumk + frac_hi

!           Divide to get the means
            umean_2d(i,j) = ubar_sum / sumk
            vmean_2d(i,j) = vbar_sum / sumk            

            if(abs(umean_2d(i,j)) .gt. 1e10)then ! write debugging info
              write(6,11)i,j,rk_lo,rk_hi,u_lo,u_hi,ubar_llyr,ubar_hlyr
     1                                           ,ubar_sum,ubar,sumk
 11	      format(' WARNING: large umean - ',2i5,9e13.4)
            endif 

          endif

900       continue

	enddo !i
	enddo !j

        write(6,*)' range of umean_2d = ',minval(umean_2d)
     1                                   ,maxval(vmean_2d)
        write(6,*)' range of vmean_2d = ',minval(vmean_2d)
     1                                   ,maxval(vmean_2d)

	return
	end
