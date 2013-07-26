
       subroutine get_aod_3d(pres_3d,ni,nj,nk,aod_3d)

       use mem_namelist, ONLY: aod

       real aod_3d(ni,nj,nk)
       real pres_3d(ni,nj,nk)

       do k = 1,nk
       do i = 1,ni
       do j = 1,nj
           if(pres_3d(i,j,k) .gt. 45000.)then
               aod_3d(i,j,k) = aod * 2.0
           else
               aod_3d(i,j,k) = aod * 0.0
           endif
       enddo ! j
       enddo ! i
       enddo ! k

       return
       end
