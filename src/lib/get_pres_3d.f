 
      subroutine get_pres_3d_old(i4time,ni,nj,nk,pres_3d,istatus)

cdoc  Returns a 3-D grid of pressures. This is useful if we have a non-uniform
cdoc  pressure grid or other type of arbitrary vertical grid.

      integer ni,nj,nk       
      real*4 pres_3d(ni,nj,nk)      
      include 'lapsparms.cmn'

      write(6,*)' get_pres_3d: calling get_laps_config'

      call get_config(istatus)

c     call get_laps_config('nest7grid',istatus)

      if(istatus .ne. 1 .or. iflag_lapsparms_cmn .ne. 1)then
         write(6,*)' Error detected in calling get_config'
         istatus = 0
         return
      else
         write(6,*)' Success in calling get_config'
      endif

      call upcase(vertical_grid,vertical_grid)
      if(vertical_grid.eq.'PRESSURE') then
         do k = 1,nk
            pressure = PRESSURE_BOTTOM_L-(k-1)*PRESSURE_INTERVAL_L
            do j = 1,nj
               do i = 1,ni
                  pres_3d(i,j,k) = pressure
               enddo            ! i
            enddo               ! j
         enddo                  ! k
         istatus = 1
      else
         print*, 'ERROR: vertical scheme not supported in get_pres_3d '
     +        ,vertical_grid
         istatus = 0
      endif
      return
      end

      subroutine get_pres_1d(i4time,nk,pres_1d_out,istatus)

cdoc  Returns a 1-D grid of pressures. This is useful if we have a non-uniform
cdoc  pressure grid. This does not support an arbitrary vertical grid.

      include 'grid_fname.cmn'                          ! grid_fnam_common

      integer nk       
      real*4 pres_1d_out(nk)      

      integer max_p
      parameter (max_p=150) 

      real*4 pressures(max_p)

      integer*4 init
      data init /0/

      save init, pressures

      namelist /pressures_nl/ pressures

      character*150 static_dir,filename

      if(init .eq. 0)then
          call get_directory(grid_fnam_common,static_dir,len_dir)

          filename = static_dir(1:len_dir)//'/pressures.nl'
 
          open(1,file=filename,status='old',err=900)
          read(1,pressures_nl,err=901)
          close(1)

          init = 1

      endif ! init = 0

      do k = 1,nk
        if(pressures(k) .le. 0. .or. pressures(k) .gt. 150000.)goto902       
        pres_1d_out(k) = pressures(k)
      enddo                  ! k

      do k = 2,nk
        if(pressures(k) .ge. pressures(k-1))goto902       
      enddo                  ! k

!     write(6,*)' Success in get_pres_1d'
      istatus = 1
      return

  900 print*,'error opening file ',filename
      istatus = 0
      return

  901 print*,'error reading pressures_nl in ',filename
      write(*,pressures_nl)
      istatus = 0
      return

  902 print*,'error in pressures_nl values in ',filename
      write(*,pressures_nl)
      istatus = 0
      return

      end      

      subroutine get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)

cdoc  Returns a 3-D grid of pressures. This is useful if we have a non-uniform
cdoc  pressure grid or other type of arbitrary vertical grid.

      integer ni,nj,nk       
      real*4 pres_1d(nk)      
      real*4 pres_3d(ni,nj,nk)      

      call get_pres_1d(i4time,nk,pres_1d,istatus)

      if(istatus .eq. 1)then
         do k = 1,nk
            do j = 1,nj
            do i = 1,ni
              pres_3d(i,j,k) = pres_1d(k)
            enddo               ! i
            enddo               ! j
         enddo                  ! k

         write(6,*)' Success in get_pres_3d'
         return

      else
         write(6,*)' No Success in get_pres_3d'
         return

      endif

      end
      

      subroutine get_rep_pres_intvl(pres_3d,ni,nj,nk,rep_pres_intvl
     1                             ,istatus)
     
cdoc  This routine returns a representative pressure interval that is used
cdoc  in the free atmosphere. This is useful if we have a variable interval
cdoc  pressure grid (even perhaps an arbitrary vertical grid), especially if 
cdoc  there are more levels in the boundary layer.

      integer ni,nj,nk       
      real*4 pres_3d(ni,nj,nk)     

      i = ni/2
      j = nj/2

      free_atmos_pres = max(70000.,pres_3d(i,j,nk-1))

      free_atmos_zcoord = rlevel_of_field(free_atmos_pres,pres_3d
     1                                   ,ni,nj,nk,i,j,istatus)
      if(istatus .ne. 1)return

      klow = nint(free_atmos_zcoord)

      pres_diff = pres_3d(i,j,klow) - pres_3d(i,j,nk) 

      rep_pres_intvl = nint( pres_diff / float(nk-klow) )

      if(rep_pres_intvl .le. 0.)then
          write(6,*)' ERROR in get_rep_pres_intvl: rep_pres_intvl < 0 '
     1              ,rep_pres_intvl
          istatus = 0
          return
      endif

      istatus = 1
      return
      end
