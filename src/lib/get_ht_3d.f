 
      subroutine get_ht_1d(nk,ht_1d_out,istatus)

cdoc  Returns a 1-D grid of heights.

      include 'grid_fname.cmn'                          ! grid_fnam_common

      integer nk       
      real ht_1d_out(nk)      

      integer max_ht
      parameter (max_ht=150) 

      real heights(max_ht)

      integer init
      data init /0/

      save init, heights

      namelist /heights_nl/ heights

      character*150 static_dir,filename

      if(init .eq. 0)then
          call get_directory(grid_fnam_common,static_dir,len_dir)

          filename = static_dir(1:len_dir)//'/heights.nl'
 
          open(1,file=filename,status='old',err=900)
          read(1,heights_nl,err=901)
          close(1)

          init = 1

      endif ! init = 0

      do k = 1,nk
        ht_1d_out(k) = heights(k)
      enddo                  ! k

      do k = 2,nk
        if(heights(k) .le. heights(k-1))goto902 
        if(heights(k) .ne. int(heights(k)/1.) * 1.)goto903
      enddo                  ! k

!     write(6,*)' Success in get_ht_1d'
      istatus = 1
      return

  900 print*,'error opening file ',filename
      istatus = 0
      return

  901 print*,'error reading heights_nl in ',filename
      write(*,heights_nl)
      istatus = 0
      return

  902 print*,'error in heights_nl values in ',filename
      write(*,heights_nl)
      istatus = 0
      return

  903 print*,'error in heights_nl values in ',filename
      print*,'values should be an integral multiple of 1 meter'
      write(*,heights_nl)
      istatus = 0
      return

      end      

      subroutine get_ht_3d(ni,nj,nk,topo,ht_3d,istatus)

cdoc  Returns a 3-D grid of heights. This is useful if we have a non-uniform
cdoc  height grid or other type of arbitrary vertical grid.

      use mem_namelist, ONLY: vertical_grid

      integer ni,nj,nk       
      real ht_1d(nk)      
      real ht_3d(ni,nj,nk)      
      real topo(ni,nj)

      call get_ht_1d(nk,ht_1d,istatus)

      sigma_htop = 20000.
      sigma_hbot =     0.

      if(istatus .eq. 1)then
         if(vertical_grid .eq. 'SIGMA_HT')then
             
            do j = 1,nj
            do i = 1,ni
               frac = (sigma_htop - topo(i,j)) / 
     1                (sigma_htop - sigma_hbot)
               do k = 1,nk
                  ht_3d(i,j,k) = topo(i,j) + (ht_1d(k) * frac)
               enddo ! k
            enddo               ! i
            enddo               ! j

            write(6,*)' Success in get_ht_3d'
            istatus = 1
            return

         elseif(vertical_grid .eq. 'HEIGHT')then
            do k = 1,nk
               do j = 1,nj
               do i = 1,ni
                 ht_3d(i,j,k) = ht_1d(k)
               enddo               ! i
               enddo               ! j
            enddo                  ! k

            write(6,*)' Success in get_ht_3d'
            istatus = 1
            return

         else
            write(6,*)' ERROR in get_ht_3d, invalid vertical grid '
     1                ,vertical_grid
            istatus = 0
            return

         endif

         write(6,*)' Success in get_ht_3d'
         istatus = 1
         return

      else
         write(6,*)' No Success in get_ht_3d'
         istatus = 0
         return

      endif

      end
      

 
      subroutine get_sigma_1d(nk,sigma_1d_out,istatus)

cdoc  Returns a 1-D grid of sigmas.

      use mem_namelist, ONLY: vertical_grid

      include 'grid_fname.cmn'                          ! grid_fnam_common

      integer nk       
      real sigma_1d_out(nk)
      real ht_1d(nk)      

      integer max_sigma
      parameter (max_sigma=150) 

      real sigmas(max_sigma)

      integer init
      data init /0/

      save init, sigmas

      namelist /sigmas_nl/ sigmas

      character*150 static_dir,filename

      if(vertical_grid .eq. 'SIGMA_HT')then ! get heights and convert to sigma
          sigma_htop = 20000.
          sigma_hbot =     0.

          call get_ht_1d(nk,ht_1d,istatus)
     1                    
          do k = 1,nk
              sigma_1d_out(k) = ht_1d(k) / (sigma_htop - sigma_hbot)
          enddo ! k

          return
      endif

      if(init .eq. 0)then
          call get_directory(grid_fnam_common,static_dir,len_dir)

          filename = static_dir(1:len_dir)//'/sigmas.nl'
 
          open(1,file=filename,status='old',err=900)
          read(1,sigmas_nl,err=901)
          close(1)

          init = 1

      endif ! init = 0

      do k = 1,nk
        sigma_1d_out(k) = sigmas(k)
      enddo                  ! k

!     QC check
      do k = 2,nk
        if(sigmas(k) .ge. sigmas(k-1))goto902 
      enddo                  ! k

!     write(6,*)' Success in get_sigma_1d'
      istatus = 1
      return

  900 print*,'error opening file ',filename
      istatus = 0
      return

  901 print*,'error reading sigmas_nl in ',filename
      write(*,sigmas_nl)
      istatus = 0
      return

  902 print*,'error in sigmas_nl values in ',filename
      write(*,sigmas_nl)
      istatus = 0
      return

      end      
