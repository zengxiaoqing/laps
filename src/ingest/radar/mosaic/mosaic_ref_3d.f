      subroutine mosaic_ref_multi(i_ra_count,maxradars,l_low_level,
     &radar_name,lat,lon,nx,ny,nz,rlat_radar,rlon_radar,rheight_radar,
     &topo,rheight_laps,grid_ra_ref,grid_mosaic_3dref,istatus)
c
c routine mosaics vxx radar files. Uses radar closest to the
c laps grid point. Depending on the value of l_low_level, uses 
c either the maximum reflectivity in the column (l_low_level = false)
c or uses the reflectivity on the first level above the laps elevation
c (l_low_level = true).
c
      implicit none

      Integer maxradars
      Integer i,j,k,kr
      Integer i_ra_count
      Integer nx,ny,nz
      Integer jstatus
      Integer istatus
      Integer icntn
      Integer icntp
      Integer icntz

      Real*4    lat(nx,ny)
      Real*4    lon(nx,ny)
      Real*4    grid_ra_ref(nx,ny,nz,maxradars)
      Real*4    grid_mosaic_3dref(nx,ny,nz)
      Real*4    topo(nx,ny)
      Real*4    rheight_laps(nx,ny,nz)
      Real*4    rlat_radar(maxradars)
      Real*4    rlon_radar(maxradars)
      Real*4    rheight_radar(maxradars)
      Real*4    ri(maxradars)
      Real*4    rj(maxradars)

      Real*4    r_dbzmax
      Real*4    rijdist
      Real*4    rjdist
      Real*4    ridist
      Real*4    r_min_dist
      Real*4    ref_base,r_missing_data

      Logical   l_low_level
      Logical   found_height

      Character*4 radar_name(maxradars)

      call get_ref_base(ref_base, istatus)
      call get_r_missing_data(r_missing_data, istatus)
c
c first find the radar location in ri/rj laps-space.
c
      istatus = 1
      if(i_ra_count .gt. 1)then
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
c find the radar with the minimum distance to the grid point in question.
c
            do k = 1,i_ra_count

               ridist=abs(float(i)-ri(k))
               rjdist=abs(float(j)-rj(k))

               rijdist=sqrt(ridist*ridist+rjdist+rjdist)

               if(rijdist.lt.r_min_dist)then
                  kr=k
                  r_min_dist=min(rijdist,r_min_dist)
               endif

            enddo

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
                        if(grid_ra_ref(i,j,k,kr).ne.ref_base .and.
     1                     grid_ra_ref(i,j,k,kr).ne.r_missing_data)then       
                         grid_mosaic_3dref(i,j,k)=grid_ra_ref(i,j,k,kr)
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

            else  ! by default this switch means to use max dBZ in vertical

               r_dbzmax=ref_base
               do k=1,nz
                  if(grid_ra_ref(i,j,k,kr).ne.ref_base .and.
     1               grid_ra_ref(i,j,k,kr).ne.r_missing_data)then
                     r_dbzmax=max(r_dbzmax,grid_ra_ref(i,j,k,kr))
                  endif
               enddo
               if(r_dbzmax.ne.ref_base)then
                  if(r_dbzmax.lt.0.0)then
                     icntn=icntn+1
                  elseif(r_dbzmax.eq.0.0)then
                     icntz=icntz+1
                  else
                     icntp=icntp+1
                  endif
               endif
               do k=1,nz 
                  grid_mosaic_3dref(i,j,k)=r_dbzmax 
               enddo ! k (SA added this k loop on 11/29/01)
            endif
         enddo
         enddo

      else
c
c no need to mosaic since only one radar!
c
        write(6,*)'Only 1 radar - no mosaic'

        icntn=0
        icntp=0
        icntz=0
        do j = 1,ny
        do i = 1,nx
          found_height=.false.
          k=0
          do while (.not.found_height)
            k=k+1
            if(k.le.nz)then
              if(rheight_laps(i,j,k).gt.topo(i,j))then
                 found_height=.true.
                 if(grid_ra_ref(i,j,k,kr).ne.ref_base .and.
     1              grid_ra_ref(i,j,k,kr).ne.r_missing_data)then
                   grid_mosaic_3dref(i,j,k)=grid_ra_ref(i,j,k
     +,i_ra_count)
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
        enddo
        enddo

      endif

      print*,'Statistics for this mosaic'
      print*,'--------------------------'
      print*,'Num points > 0.0 ',icntp
      print*,'Num points = 0.0 ',icntz
      print*,'Num points < 0.0 ',icntn
 
      return
      end
