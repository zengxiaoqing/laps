      subroutine get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)

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
      k = nk

      rep_pres_intvl = nint( pres_3d(i,j,k-1) - pres_3d(i,j,k) )

      if(rep_pres_intvl .le. 0.)then
          write(6,*)' ERROR in get_rep_pres_intvl: rep_pres_intvl < 0 '
     1              ,rep_pres_intvl
          istatus = 0
          return
      endif

      istatus = 1
      return
      end
