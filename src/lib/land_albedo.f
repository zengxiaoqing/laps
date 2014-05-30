
      subroutine land_albedo(lu,ni,nj,albedo)

!     This routine takes land use categories and broadband albedo relationship
!     to derive a 3 color albedo

      real lu(ni,nj)          ! Land Use (USGS 24 category)
      real albedo(3,ni,nj)    ! Albedo (Red, Green, Blue)

!     See http://www.mmm.ucar.edu/mm5/mm5v2/landuse-usgs-tbl.html

      do i = 1,ni
      do j = 1,nj
          if(lu(i,j) .eq. 1.)then      ! Urban and Built-up Land   (.18 Gray)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.11
          elseif(lu(i,j) .eq. 2.)then  ! Dryland Cropland/Pasture  (.17 Brown)
              albedo(1,i,j) = 0.17
              albedo(2,i,j) = 0.17
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 3.)then  ! Irrigtd Cropland/Pasture  (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 4.)then  ! Mixed                     (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 5.)then  ! Cropland Grassland Mosaic (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 6.)then  ! Cropland Woodland Mosaic  (.16 Green)
              albedo(1,i,j) = 0.12
              albedo(2,i,j) = 0.16
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 7.)then  ! Grassland                 (.19 Brown)
              albedo(1,i,j) = 0.19
              albedo(2,i,j) = 0.19
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 8.)then  ! Shrubland                 (.22 Brown)
              albedo(1,i,j) = 0.22
              albedo(2,i,j) = 0.22
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 9.)then  ! Mixed Shrubland/Grassland (.20 Brown)
              albedo(1,i,j) = 0.20
              albedo(2,i,j) = 0.20
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 10.)then ! Savanna                   (.20 Brown)
              albedo(1,i,j) = 0.20
              albedo(2,i,j) = 0.20
              albedo(3,i,j) = 0.01
          elseif(lu(i,j) .eq. 11.)then ! Deciduous Broadleaf       (.16 Green)
              albedo(1,i,j) = 0.12
              albedo(2,i,j) = 0.16
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 12.)then ! Deciduous Needleleaf      (.14 Green)
              albedo(1,i,j) = 0.10
              albedo(2,i,j) = 0.14
              albedo(3,i,j) = 0.07
          elseif(lu(i,j) .eq. 13.)then ! Evergreen Broadleaf       (.12 Green)
              albedo(1,i,j) = 0.07
              albedo(2,i,j) = 0.09
              albedo(3,i,j) = 0.045
          elseif(lu(i,j) .eq. 14.)then ! Evergreen Needleleaf      (.12 Green)
              albedo(1,i,j) = 0.07
              albedo(2,i,j) = 0.09
              albedo(3,i,j) = 0.045
          elseif(lu(i,j) .eq. 15.)then ! Mixed Forest              (.13 Green)
              albedo(1,i,j) = 0.075
              albedo(2,i,j) = 0.10
              albedo(3,i,j) = 0.05
          elseif(lu(i,j) .eq. 16.)then ! water (.08)
              albedo(1,i,j) = 0.04
              albedo(2,i,j) = 0.04
              albedo(3,i,j) = 0.12
          elseif(lu(i,j) .eq. 19.)then ! Barren or Sparsely Vegetated (.25 Brown)
              albedo(1,i,j) = 0.25
              albedo(2,i,j) = 0.15
              albedo(3,i,j) = 0.07
          else ! default
              albedo(1,i,j) = 0.19
              albedo(2,i,j) = 0.19
              albedo(3,i,j) = 0.01
          endif                        
      enddo ! j
      enddo ! i

      return
      end

      subroutine land_albedo_bm(rlat_laps,rlon_laps,ni,nj,albedo
     1                         ,istatus)

      use ppm

!     Read and interpolate from sector of Blue Marble Image
!     See code in /scratch/staging/fab/albers/nasa

      real albedo(3,ni,nj)    ! Albedo (Red, Green, Blue)
      real result(ni,nj)      ! Interpolated image channel

!     convert -crop 2000x2000+12000+9000 -compress none
!     world.200404.3x21600x21600.A1.png world.200404.3x21600x21600.A1.crop2.ppm

      character*255 directory
      character*255 file
      integer u

      integer ncol,iwidth,iheight
      parameter (ncol=3)
      parameter (iwidth=2000)
      parameter (iheight=2000)

      real rlat_img(iwidth,iheight)
      real rlon_img(iwidth,iheight)
      real ri_img(ni,nj)            
      real rj_img(ni,nj)             
      real rlat_laps(ni,nj)
      real rlon_laps(ni,nj)

      integer img(ncol,iwidth,iheight)                                               
      real array_2d(iwidth,iheight)
      u = 11

      I4_elapsed = ishow_timer()

      write(6,*)' Subroutine land_albedo_bm...'

      call get_directory('static',directory,len_dir)
      file=trim(directory)//'world.200405.3x21600x21600.A1.crop.ppm'

      write(6,*)' Open ',trim(file)
      write(6,*)' dims ',ncol,iwidth,iheight
      
!     Read section of NASA Blue Marble Image in PPM format
      open(u,file=trim(file),status='old',err=999)
      call read_ppm (u,img,ncol,iwidth,iheight)
      close(u)

      pix_latlon = 1. / 240
      rlat_start =   90. - 11000. * pix_latlon
      rlon_start = -180. + 17000. * pix_latlon
      rlat_end = rlat_start - float(iheight)  * pix_latlon
      rlon_end = rlon_start + float(iwidth)  * pix_latlon

      write(6,*)' BM PPM lat range: ',rlat_start,rlat_end
      write(6,*)' BM PPM lon range: ',rlon_start,rlon_end

!     Fill arrays of image lat/lons and LAPS ri,rj
      istatus_img = 1
      do i = 1,ni     
      do j = 1,nj     
        rj_img(i,j) = (rlat_start - rlat_laps(i,j)) / pix_latlon
        ri_img(i,j) = (rlon_laps(i,j) - rlon_start) / pix_latlon
        if(rlat_laps(i,j) .gt. rlat_start .OR.
     1     rlat_laps(i,j) .lt. rlat_end   .OR.
     1     rlon_laps(i,j) .lt. rlon_start .OR. 
     1     rlon_laps(i,j) .gt. rlon_end)then
            istatus_img = 0
        endif
      enddo ! j
      enddo ! i

      if(istatus_img .eq. 0)then
          write(6,*)' WARNING: LAPS grid extends outside image'
      else
          write(6,*)' LAPS grid is contained within image'      
      endif

      do i = 1,iwidth,200
      do j = 1,iheight,200

        write(6,11)img(:,i,j),i,j
11      format(5i6,2f8.3)

      enddo ! j
      enddo ! i

      I4_elapsed = ishow_timer()

      write(6,*)' Bilinearly Interpolate'

!     Interpolate to LAPS grid using bilinear_laps_2d
      do ic = 1,3
        array_2d(:,:) = img(ic,:,:)
        call bilinear_laps_2d(ri_img,rj_img,iwidth,iheight,ni,nj 
     1                       ,array_2d,result)
        do i = 1,ni
        do j = 1,nj
          if(result(i,j) .le. 179.)then ! Based on NASA spline
            albedo(ic,i,j) = result(i,j) / 716.                  
          else
            albedo(ic,i,j) = 0.25 + (result(i,j) - 179.) / 122. 
          endif
        enddo ! j
        enddo ! i
      enddo ! ic

cdoc  Interpolate 2-d array to find the field values at fractional grid
cdoc  points.

!     real array_2d(imax,jmax)                         ! I
!     real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
!     real result(nx_laps,ny_laps)                     ! O

      write(6,*)' Interpolated albedo: ',albedo(:,ni/2,nj/2)

      istatus = 1
      I4_elapsed = ishow_timer()

      return

!     Error condition
999   istatus = 0
      return

      end
