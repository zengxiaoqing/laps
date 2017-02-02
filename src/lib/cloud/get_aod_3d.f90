
       subroutine get_aod_3d(pres_3d,heights_3d,topo_2d,ni,nj,nk,aod,aod_3d)

       use mem_namelist, ONLY: redp_lvl,aero_scaleht
       use mem_allsky, ONLY: mode_aero_cld

       include 'rad.inc'

       real aod_3d(ni,nj,nk) ! aerosol extinction coefficient
       real pres_3d(ni,nj,nk)
       real heights_3d(ni,nj,nk)
       real topo_2d(ni,nj)

       gas_scale_height = 8000.

       write(6,*)' subroutine get_aod_3d: aod/redp_lvl/aero_scaleht = ' &
                                         ,aod,redp_lvl,aero_scaleht

       do i = 1,ni
       do j = 1,nj
         sum_aod = 0.
         do k = 1,nk
!          pratio = pres_3d(i,j,k) / pref
           h_agl = heights_3d(i,j,k) - redp_lvl 
           alpha_low = (aod/aero_scaleht) * exp(-h_agl/aero_scaleht)
           if(heights_3d(i,j,k) .gt. h1_ha .and. & 
              heights_3d(i,j,k) .le. h2_ha .and. aod_ha .gt. 0.)then
               alpha_high = alpha_ha
           else
               alpha_high = 0.
           endif
           aod_3d(i,j,k) = alpha_low + alpha_high
         
           if(i .eq. ni/2 .and. j .eq. nj/2)then
               if(h_agl .gt. 0.)then
                   ave_aod = 0.5 * (aod_3d(i,j,k)+aod_3d(i,j,k-1))
                   sum_aod = sum_aod + ave_aod       & 
                           * (heights_3d(i,j,k) - heights_3d(i,j,k-1))
               endif
               write(6,101)k,h_agl,aod_3d(i,j,k),sum_aod
101            format('k,agl,aod_3d,sum',i5,f9.1,e15.8,f9.5)
           endif
         enddo ! k
       enddo ! j
       enddo ! i

       if(mode_aero_cld .gt. 1)then ! synthetic aerosol plume
         aod_3d(ni/2+12,nj/2,1:nk-8) = .001
       endif

       return
       end
