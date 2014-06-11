
      subroutine mosaic_ref_multi(i_ra_count,maxradars,l_low_level,       ! I
     &                            radar_name,lat,lon,nx,ny,nz,            ! I
     &                            rlat_radar,rlon_radar,rheight_radar,    ! I
     &                            topo,rheight_laps,grid_ra_ref,          ! I
     &                            grid_ra_ref_offset,ioffset,joffset,     ! I
     &                            nx_r,ny_r,                              ! I
     &                            imosaic_3d,                             ! I
     &                            dist_multiradar_2d,                     ! I  
     &                            l_offset,                               ! I
     &                            grid_mosaic_2dref,grid_mosaic_3dref,    ! I/O
     &                            closest_radar_m,istatus)                ! O
c
c routine mosaics vxx radar files. Uses radar closest to the
c laps grid point. Depending on the value of l_low_level, uses 
c either the maximum reflectivity in the column (l_low_level = false)
c or uses the reflectivity on the first level above the laps elevation
c (l_low_level = true).
c
      Real    lat(nx,ny)
      Real    lon(nx,ny)
      Real    grid_ra_ref(nx,ny,nz,maxradars)
      Real    grid_ra_ref_offset(nx_r,ny_r,nz,maxradars)
      Real    grid_mosaic_2dref(nx,ny)
      Real    dist_multiradar_2d(nx,ny,maxradars)
      Real    closest_radar_m(nx,ny)
      Real    grid_mosaic_3dref(nx,ny,nz)
      Real    topo(nx,ny)
      Real    rheight_laps(nx,ny,nz)
      Real    rlat_radar(maxradars)
      Real    rlon_radar(maxradars)
      Real    rheight_radar(maxradars)
      Real    ri(maxradars)
      Real    rj(maxradars)
      real    dist_radar_m(maxradars)

      Logical   l_low_level
      Logical   found_height
      Logical   l_valid
      Logical   l_valid_latlon(maxradars)
      Logical   l_offset 
!     Parameter (l_offset = .true.)

      integer   lr_2d(nx,ny)                     ! closest radar
      Integer   ioffset(maxradars),joffset(maxradars)

      Character*4 radar_name(maxradars)

      write(6,*)
      write(6,*)' Subroutine mosaic_ref_multi: imosaic_3d =', imosaic_3d       

!     if(l_offset)then
!         write(6,*)' grid_ra_ref_offset range: '
!    1              ,minval(grid_ra_ref_offset)
!    1              ,maxval(grid_ra_ref_offset)
!     else
!         write(6,*)' grid_ra_ref range: '
!    1              ,minval(grid_ra_ref)
!    1              ,maxval(grid_ra_ref)
!     endif

      call get_ref_base(ref_base, istatus)
      call get_r_missing_data(r_missing_data, istatus)
      call get_grid_spacing_cen(grid_spacing_cen_m,istatus)

!     Initialize
      grid_mosaic_2dref = r_missing_data
      grid_mosaic_3dref = r_missing_data
      lr_2d = 0
c
c first find the radar location in ri/rj laps-space.
c
      istatus = 1
      if(i_ra_count .ge. 1)then ! essentially all the time

!        Determine which radars have valid lat/lon
         do k = 1,i_ra_count 
            if(rlat_radar(k) .eq. r_missing_data .or.
     1         rlon_radar(k) .eq. r_missing_data      )then
                write(6,*)' No valid or single lat/lon for radar ',k
     1                   ,' ',radar_name(k)
                l_valid_latlon(k) = .false.

            else
                call latlon_to_rlapsgrid(rlat_radar(k),
     &                                   rlon_radar(k),
     &                                   lat,lon,
     &                                   nx,ny,
     &                                   ri(k),rj(k),
     &                                   jstatus)
                if(jstatus.ne.1)then
                    write(6,*)
     1               'Error computing ri/rj for radar (outside domain)?'    
                endif
                write(6,*)radar_name(k),k,ri(k),rj(k)
     1                   ,rlat_radar(k),rlon_radar(k)
51              format('Valid lat/lon: ',a,i5,2f8.1,2f8.2)
                l_valid_latlon(k) = .true.

                if(l_offset)then
                    write(6,*)'      offsets ',ioffset(k),joffset(k)               
                endif

            endif

         enddo ! radars

!        Loop through all horizontal gridpoints to define best radar array
         icntn=0
         icntp=0
         icntz=0
         icntb=0
         do j = 1,ny
         do i = 1,nx

            r_min_dist_m = sqrt(float(nx*nx+ny*ny))*grid_spacing_cen_m       
c
c find the valid radar with the minimum distance to the grid point in question.
c
            lr = 0 

!           Loop through all radars
            do l = 1,i_ra_count
               l_valid = .false.

               if(.not. l_offset)then
                 do k = 1,nz ! Look for non-missing reflectivity in column
                   if(grid_ra_ref(i,j,k,l) .ne. r_missing_data)then     
                     l_valid = .true.
                   endif
                 enddo ! k

               else ! use offset array
                 io = i - ioffset(l)
                 jo = j - joffset(l)
                 if(io .lt. 1 .or. io .gt. nx_r .or. 
     1              jo .lt. 1 .or. jo .gt. ny_r)then
!                  write(6,*)' offset out of domain ',i,j,io,jo
!    1                      ,ioffset(l),joffset(l),radar_name(l)
!                  stop
                   continue

                 else   
                   do k = 1,nz ! Look for non-missing reflectivity in column
                     if(grid_ra_ref_offset(io,jo,k,l) .ne. 
     1                                  r_missing_data)then     
                       l_valid = .true.
                       goto 50
                     endif
                   enddo ! k

                 endif ! in bounds of offset array                  

               endif ! use full array (not offset)

 50            if(l_valid_latlon(l))then
                   ridist = float(i)-ri(l)
                   rjdist = float(j)-rj(l)
                   rijdist=sqrt(ridist*ridist + rjdist*rjdist)
                   dist_radar_m(l) = rijdist * grid_spacing_cen_m

               else ! No valid or single lat/lon, use distance array
                   dist_radar_m(l) = dist_multiradar_2d(i,j,l)

               endif

               if(dist_radar_m(l) .lt. r_min_dist_m .and. l_valid)then       
                  lr=l
                  r_min_dist_m = dist_radar_m(l)
                  closest_radar_m(i,j) = r_min_dist_m       
               endif

            enddo ! l

            lr_2d(i,j) = lr

         enddo ! i
         enddo ! j

         do l = 1,i_ra_count ! mosaic in one radar at a time
            do j = 1,ny
            do i = 1,nx

            if(l .eq. lr_2d(i,j))then ! closest radar at this grid point 
               r_dbzmax=ref_base

               if(.not. l_offset)then
!                 Get max ref in column
                  do k=1,nz
                     if(grid_ra_ref(i,j,k,l).ne.ref_base .and.
     1                  grid_ra_ref(i,j,k,l).ne.r_missing_data)then
                        r_dbzmax=max(r_dbzmax,grid_ra_ref(i,j,k,l))
                     endif
                  enddo

                  if(imosaic_3d .eq. 0)then      ! vrc output only
                     do k=1,nz 
                        grid_mosaic_2dref(i,j)=r_dbzmax 
                        grid_mosaic_3dref(i,j,k)=r_dbzmax 
                     enddo ! k 

                  elseif(imosaic_3d .eq. 1)then  ! vrz output only
                     grid_mosaic_3dref(i,j,:)=grid_ra_ref(i,j,:,l)   

                  elseif(imosaic_3d .eq. 2)then  ! both vrc & vrz
                     do k=1,nz 
                        grid_mosaic_2dref(i,j)=r_dbzmax 
                        grid_mosaic_3dref(i,j,k)=grid_ra_ref(i,j,k,l)   
                     enddo ! k 

                  endif ! imosaic_3d

               else ! l_offset
                 io = i - ioffset(l)
                 jo = j - joffset(l)
                 
                 if(io .ge. 1 .AND. io .le. nx_r .AND.
     1              jo .ge. 1 .AND. jo .le. ny_r)then

!                  Get max ref in column
                   do k=1,nz
                     if(grid_ra_ref_offset(io,jo,k,l).ne.ref_base .and.
     1                  grid_ra_ref_offset(io,jo,k,l).ne.r_missing_data
     1                                                             )then
                        r_dbzmax=max(r_dbzmax,
     1                               grid_ra_ref_offset(io,jo,k,l))
                     endif
                   enddo

                   if(imosaic_3d .eq. 0)then      ! vrc output only
                     grid_mosaic_2dref(i,j)=r_dbzmax 
                     grid_mosaic_3dref(i,j,:)=r_dbzmax 

                   elseif(imosaic_3d .eq. 1)then  ! vrz output only
                     grid_mosaic_3dref(i,j,:)=
     1               grid_ra_ref_offset(io,jo,:,l)       

                   elseif(imosaic_3d .eq. 2)then  ! both vrc & vrz
                     grid_mosaic_2dref(i,j)=r_dbzmax 
                     grid_mosaic_3dref(i,j,:)=
     1               grid_ra_ref_offset(io,jo,:,l)       

                   endif ! imosaic_3d

                 endif ! in bounds of offset array

               endif ! .true.

!              Increment stats
               if(r_dbzmax.ne.ref_base)then
                  if(r_dbzmax.lt.0.0)then
                     icntn=icntn+1
                  elseif(r_dbzmax.eq.0.0)then
                     icntz=icntz+1
                  else
                     icntp=icntp+1
                  endif
               else
                  icntb = icntb + 1
               endif

            endif ! .true.

            enddo ! i
            enddo ! j
         enddo ! l

      endif ! i_ra_count > 1

      print*,'Statistics for this mosaic'
      print*,'--------------------------'
      print*,'Num points > 0.0  ',icntp
      print*,'Num points = 0.0  ',icntz
      print*,'Num points < 0.0  ',icntn
      print*,'Num points = base ',icntb

      intvl = int(nx/80) + 1

      write(6,*)
      write(6,*)' Map of ',i_ra_count,' valid radars used:'

      do j = ny,1,-intvl
          write(6,101)(lr_2d(i,j),i=1,nx,intvl) 
 101      format(75i2)
      enddo ! j

      return
      end
