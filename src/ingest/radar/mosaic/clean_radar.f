

      subroutine clean_radar(ref,ni,nj,nk)

      use mem_namelist, ONLY: grid_spacing_m, ref_base

      real ref(ni,nj,nk)
      real ref_buf(ni,nj,nk)

      ref_buf(:,:,:) = ref(:,:,:)

      if(grid_spacing_m .gt. 3500.)then
          write(6,*)' Skip radar cleaning step'
          return
      else
          write(6,*)' Clean radar data from small-scale noise'
      endif

      nclean = 0

      do k = 1,nk
      do j = 2,nj-1
      do i = 2,ni-1
          sum9 = sum(ref(i-1:i+1,j-1:j+1,k))
          sum8 = sum9 - ref(i,j,k)
          ave8 = sum8 / 8.
          if(ref(i,j,k) .gt. ave8 + 20.)then
              ref_buf(i,j,k) = ave8 ! ref_base
              nclean = nclean + 1
          endif
      enddo ! i
      enddo ! j
      enddo ! k

      ref(:,:,:) = ref_buf(:,:,:)

      write(6,*)' Number of grid points cleaned is ',nclean

      return
      end
