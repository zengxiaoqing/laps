      subroutine get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
      include 'lapsparms.cmn'
       
      real*4 pres_3d(ni,nj,nk)      

      if(iflag_lapsparms_cmn .ne. 1)then
         write(6,*)' get_pres_3d: calling get_laps_config'

         call get_laps_config('nest7grid',istatus)
         if(istatus .ne. 1 .or. iflag_lapsparms_cmn .ne. 1)then
            write(6,*)' Error detected in calling get_laps_config'
            istatus = 0
            return
         else
            write(6,*)' Success in calling get_laps_config'
         endif

      endif

      call upcase(vertical_grid,vertical_grid)
      if(vertical_grid.eq.'PRESSURE') then
         do k = 1,nk
            pressure = PRESSURE_BOTTOM_L+(k-1)*PRESSURE_INTERVAL_L
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

