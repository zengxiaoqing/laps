c
c---------------------------------------------------------------------------
c
      subroutine tdcheck(nx_laps,ny_laps,td_sfc,tp_sfc,
     &icnt,i_mx,j_mx,dmax,dmin)

      integer icnt,i_mx,j_mx

      real td_sfc(nx_laps,ny_laps)
      real tp_sfc(nx_laps,ny_laps)
      real dmax,dmin

      call get_r_missing_data(r_missing_data,istatus)

      icnt=0
      diff_mx = -1.e30
      diff_mn =  1.e30
      i_mx = 0
      j_mx = 0
      i_mn = 0
      j_mn = 0
      do j=1,ny_laps
      do i=1,nx_laps
         if(td_sfc(i,j) .gt. tp_sfc(i,j) .AND. 
     1      td_sfc(i,j) .ne. r_missing_data)then
            diff = td_sfc(i,j) - tp_sfc(i,j)
            td_sfc(i,j)=tp_sfc(i,j)
            icnt=icnt+1
            if(diff .gt. diff_mx) then
               diff_mx = diff
               i_mx = i
               j_mx = j
            endif
            if(diff .lt. diff_mn) then
               diff_mn = diff
               i_mn = i
               j_mn = j
            endif
         endif
      enddo
      enddo

      return
      end
