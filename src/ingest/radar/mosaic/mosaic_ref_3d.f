
      subroutine mosaic_ref_multi(i_ra_count,maxradars,l_low_level,       ! I
     &                            radar_name,lat,lon,nx,ny,nz,            ! I
     &                            rlat_radar,rlon_radar,rheight_radar,    ! I
     &                            topo,rheight_laps,grid_ra_ref,          ! I
     &                            imosaic_3d,                             ! I
     &                            grid_mosaic_3dref,istatus)              ! O
c
c routine mosaics vxx radar files. Uses radar closest to the
c laps grid point. Depending on the value of l_low_level, uses 
c either the maximum reflectivity in the column (l_low_level = false)
c or uses the reflectivity on the first level above the laps elevation
c (l_low_level = true).
c
      Real*4    lat(nx,ny)
      Real*4    lon(nx,ny)
      Real*4    grid_ra_ref(nx,ny,nz,maxradars)
      Real*4    grid_mosaic_2dref(nx,ny)
      Real*4    grid_mosaic_3dref(nx,ny,nz)
      Real*4    topo(nx,ny)
      Real*4    rheight_laps(nx,ny,nz)
      Real*4    rlat_radar(maxradars)
      Real*4    rlon_radar(maxradars)
      Real*4    rheight_radar(maxradars)
      Real*4    ri(maxradars)
      Real*4    rj(maxradars)

      Logical   l_low_level
      Logical   found_height
      Logical   l_valid

      integer   lr_2d(nx,ny)

      Character*4 radar_name(maxradars)

      write(6,*)
      write(6,*)' Subroutine mosaic_ref_multi: imosaic_3d =', imosaic_3d       

      call get_ref_base(ref_base, istatus)
      call get_r_missing_data(r_missing_data, istatus)

!     Initialize
      grid_mosaic_2dref = r_missing_data
      grid_mosaic_3dref = r_missing_data
      lr_2d = 0
c
c first find the radar location in ri/rj laps-space.
c
      istatus = 1
      if(i_ra_count .ge. 1)then ! essentially all the time
         do k = 1,i_ra_count
            call latlon_to_rlapsgrid(rlat_radar(k),
     &                            rlon_radar(k),
     &                            lat,lon,
     &                            nx,ny,
     &                            ri(k),rj(k),
     &                            jstatus)
            if(jstatus.ne.1)then
               write(6,*)
     1               'Error computing ri/rj for radar (outside domain)?'       
               write(6,*)'Name: ',radar_name(k),ri(k),rj(k)
            endif
         enddo

         icntn=0
         icntp=0
         icntz=0
         do j = 1,ny
         do i = 1,nx

            r_min_dist = sqrt(float(nx*nx+ny*ny))
c
c find the valid radar with the minimum distance to the grid point in question.
c
            lr = 0 

            do l = 1,i_ra_count
               l_valid = .false.

               do k = 1,nz ! Look for non-missing reflectivity in column
                  if(grid_ra_ref(i,j,k,l) .ne. r_missing_data)then        
                     l_valid = .true.
                  endif
               enddo ! k

               ridist = float(i)-ri(l)
               rjdist = float(j)-rj(l)

               rijdist=sqrt(ridist*ridist + rjdist*rjdist)

               if(rijdist .lt. r_min_dist .and. l_valid)then
                  lr=l
                  r_min_dist = rijdist
               endif

            enddo ! l

            if(l_low_level)then

c use the altitude of the radar compared to the altitude of the laps
c level to find the appropriate level within grid_ra_ref as the data
c to mosaic
               found_height=.false.
               k=0
               do while (.not.found_height)
                  k=k+1
                  if(k.le.nz)then
                     if(rheight_laps(i,j,k).gt.topo(i,j))then
                        found_height=.true.

                        if(grid_ra_ref(i,j,k,lr).ne.ref_base .and.
     1                     grid_ra_ref(i,j,k,lr).ne.r_missing_data)then       
                           grid_mosaic_3dref(i,j,k)
     1                                           =grid_ra_ref(i,j,k,lr)

!                          Increment stats
                           if(grid_mosaic_3dref(i,j,k).lt.0.0)then
                              icntn=icntn+1
                           elseif(grid_mosaic_3dref(i,j,k).eq.0.0)then
                              icntz=icntz+1
                           else
                              icntp=icntp+1
                           endif

                        endif
                     endif
                  else
                    found_height=.true.
                  endif
               enddo

            else  ! by default this switch means to use max dBZ in vertical?

               r_dbzmax=ref_base

               if(lr .gt. 0)then
!                 Get max ref in column
                  do k=1,nz
                     if(grid_ra_ref(i,j,k,lr).ne.ref_base .and.
     1                  grid_ra_ref(i,j,k,lr).ne.r_missing_data)then
                        r_dbzmax=max(r_dbzmax,grid_ra_ref(i,j,k,lr))
                     endif
                  enddo

                  if(imosaic_3d .eq. 0)then      ! vrc output only
                     do k=1,nz 
                        grid_mosaic_2dref(i,j)=r_dbzmax 
                        grid_mosaic_3dref(i,j,k)=r_dbzmax 
                     enddo ! k 

                  elseif(imosaic_3d .eq. 1)then  ! vrz output only
                     do k=1,nz 
                        grid_mosaic_3dref(i,j,k)=grid_ra_ref(i,j,k,lr)       
                     enddo ! k 

                  elseif(imosaic_3d .eq. 2)then  ! both vrc & vrz
                     do k=1,nz 
                        grid_mosaic_2dref(i,j)=r_dbzmax 
                        grid_mosaic_3dref(i,j,k)=grid_ra_ref(i,j,k,lr)       
                     enddo ! k 

                  endif ! imosaic_3d

               endif ! lr

!              Increment stats
               if(r_dbzmax.ne.ref_base)then
                  if(r_dbzmax.lt.0.0)then
                     icntn=icntn+1
                  elseif(r_dbzmax.eq.0.0)then
                     icntz=icntz+1
                  else
                     icntp=icntp+1
                  endif
               endif

            endif ! low_level switch

            lr_2d(i,j) = lr

         enddo ! i
         enddo ! j

      endif ! i_ra_count > 1

      print*,'Statistics for this mosaic'
      print*,'--------------------------'
      print*,'Num points > 0.0 ',icntp
      print*,'Num points = 0.0 ',icntz
      print*,'Num points < 0.0 ',icntn

      intvl = int(nx/68) + 1

      write(6,*)
      write(6,*)' Map of radars used: intvl = ',intvl

      do j = ny,1,-intvl
          write(6,101)(lr_2d(i,j),i=1,nx,intvl) 
 101      format(75i2)
      enddo ! j

      return
      end
