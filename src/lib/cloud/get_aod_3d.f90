
       subroutine get_aod_3d(pres_3d,heights_3d,topo_2d,ni,nj,nk,aod_3d)

       use mem_namelist, ONLY: aod

       real aod_3d(ni,nj,nk) ! aerosol extinction coefficient
       real pres_3d(ni,nj,nk)
       real heights_3d(ni,nj,nk)
       real topo_2d(ni,nj)

       gas_scale_height = 8000.
       aero_scale_height = 1500.

       pref = 85000.

       write(6,*)' subroutine get_aod_3d: aod = ',aod

       do i = 1,ni
       do j = 1,nj
         sum_aod = 0.
         do k = 1,nk
           pratio = pres_3d(i,j,k) / pref
           h_agl = heights_3d(i,j,k) - topo_2d(i,j)
           aod_3d(i,j,k) = aod/aero_scale_height * exp(-h_agl/aero_scale_height)

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

       return
       end
