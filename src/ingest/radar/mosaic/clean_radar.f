

      subroutine clean_radar(ref,ni,nj,nk)

      use mem_namelist, ONLY: grid_spacing_m, ref_base, r_missing_data

      real ref(ni,nj,nk)
      real ref_buf(ni,nj,nk)

      ref_buf(:,:,:) = ref(:,:,:)

      if(grid_spacing_m .gt. 3500.)then
          write(6,*)' Skip radar cleaning step'
          return
      else
          write(6,*)' Clean radar data from small-scale noise'
          if(grid_spacing_m .gt. 750.)then
              ib = 1
          else
              ib = 2
          endif
      endif

      nclean = 0
      rnumbox = float((2*ib+1)**2 - 1)

      do k = 1,nk
      do j = 1+ib,nj-ib
      do i = 1+ib,ni-ib
          sum9 = sum(ref(i-ib:i+ib,j-ib:j+ib,k))
          if(abs(sum9) .le. 1e6)then ! screen out summed missing data values
              sum8 = sum9 - ref(i,j,k)
              ave8 = sum8 / rnumbox
              if(ref(i,j,k) .gt. (ave8 + 20.) .OR.
     1           (ave8 .eq. ref_base .and. grid_spacing_m .le. 750.) 
     1                                                             )then        
                  ref_buf(i,j,k) = ave8 ! ref_base
                  nclean = nclean + 1
              endif
          endif
      enddo ! i
      enddo ! j
      enddo ! k

      ref(:,:,:) = ref_buf(:,:,:)

      write(6,*)' Number of grid points cleaned is ',nclean

      return
      end
