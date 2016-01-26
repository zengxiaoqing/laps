
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
              albedo(3,i,j) = 0.10
          elseif(lu(i,j) .eq. 3.)then  ! Irrigtd Cropland/Pasture  (.18 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 4.)then  ! Mixed                     (.18 Green)
              albedo(1,i,j) = 0.14
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 5.)then  ! Cropland Grassland Mosaic (.18 Green)
              albedo(1,i,j) = 0.14
              albedo(2,i,j) = 0.18
              albedo(3,i,j) = 0.09
          elseif(lu(i,j) .eq. 6.)then  ! Cropland Woodland Mosaic  (.16 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.16
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 7.)then  ! Grassland                 (.19 Brown)
              albedo(1,i,j) = 0.19
              albedo(2,i,j) = 0.19
              albedo(3,i,j) = 0.10
          elseif(lu(i,j) .eq. 8.)then  ! Shrubland                 (.22 Brown)
              albedo(1,i,j) = 0.22
              albedo(2,i,j) = 0.22
              albedo(3,i,j) = 0.10
          elseif(lu(i,j) .eq. 9.)then  ! Mixed Shrubland/Grassland (.20 Brown)
              albedo(1,i,j) = 0.20
              albedo(2,i,j) = 0.20
              albedo(3,i,j) = 0.10
          elseif(lu(i,j) .eq. 10.)then ! Savanna                   (.20 Brown)
              albedo(1,i,j) = 0.20
              albedo(2,i,j) = 0.20
              albedo(3,i,j) = 0.10
          elseif(lu(i,j) .eq. 11.)then ! Deciduous Broadleaf       (.16 Green)
              albedo(1,i,j) = 0.13
              albedo(2,i,j) = 0.16
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 12.)then ! Deciduous Needleleaf      (.14 Green)
              albedo(1,i,j) = 0.11
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
              albedo(1,i,j) = 0.08
              albedo(2,i,j) = 0.08
              albedo(3,i,j) = 0.08
          elseif(lu(i,j) .eq. 19.)then ! Barren or Sparsely Vegetated (.25 Brown)
              albedo(1,i,j) = 0.25
              albedo(2,i,j) = 0.15
              albedo(3,i,j) = 0.07
          else ! default
              albedo(1,i,j) = 0.19
              albedo(2,i,j) = 0.19
              albedo(3,i,j) = 0.10
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
      character*255 file,file_bm 
      character*10  c10_fname /'nest7grid'/
      integer u,u_out
      logical l_there

      integer ncol,iwidth,iheight
      parameter (ncol=3)
!     parameter (iwidth=2000)
!     parameter (iheight=2000)

!     real rlat_img(iwidth,iheight)
!     real rlon_img(iwidth,iheight)
      real, allocatable :: rlat_img(:,:),rlon_img(:,:) 
      real ri_img(ni,nj)            
      real rj_img(ni,nj)             
      real rlat_laps(ni,nj)
      real rlon_laps(ni,nj)
      real topo_laps_dum(ni,nj)

!     integer img(ncol,iwidth,iheight)                                               
      integer, allocatable :: img(:,:,:)                   
!     real array_2d(iwidth,iheight)
      real, allocatable :: array_2d(:,:)                   
      u = 11

      write(6,*)' Subroutine land_albedo_bm...'

      call get_directory('static',directory,len_dir)
      file_bm=trim(directory)//'albedo_multispectral.dat'
      inquire(file=trim(file_bm),exist=l_there)

      if(l_there)then
        write(6,*)' Blue Marble binary remapped albedo file exists'
        open(u,file=trim(file_bm),form='unformatted' 
     1      ,status='old',err=999)
        write(6,*)' Successful open, now reading'
        read(u,err=999)albedo
        close(u)
        write(6,*)' first point is ',albedo(1,1,1)
        if(albedo(1,1,1) .gt. 1.0)then
          write(6,*)' ERROR bad first point'
          goto 999
        endif
      else
        write(6,*)
     1' Blue Marble binary file is absent, generate and create it'       

      file=trim(directory)//'world.200405.3x21600x21600.crop.ppm'

      write(6,*)' Open for reading ',trim(file)

!     Read section of NASA Blue Marble Image in PPM format
      open(u,file=trim(file),status='old',err=999)
      read(u,*)   
      read(u,*)iwidth,iheight
      rewind(u)
      write(6,*)' dynamic dims ',ncol,iwidth,iheight
      allocate(img(ncol,iwidth,iheight))
      call read_ppm (u,img,ncol,iwidth,iheight)
      close(u)

      allocate(rlat_img(iwidth,iheight))
      allocate(rlon_img(iwidth,iheight))
      allocate(array_2d(iwidth,iheight))

!     Consider dynamic means to get rlat_start and rlon_start
!     Use domain lat/lon bounds with a 0.2 deg cushion
      call get_domain_perimeter(ni,nj,c10_fname  
     1                  ,rlat_laps,rlon_laps,topo_laps_dum
     1                  ,0.2,rnorth,south,east,west,istatus)

      write(6,5)rnorth,south,east,west
5     format('  NSEW Domain 0.2deg perimeter',4f9.3)                 

      pix_latlon = 1. / 240.
!     rlat_start =   90. - 11000. * pix_latlon
!     rlon_start = -180. + 17000. * pix_latlon
      rlat_start = rnorth                         
      rlon_start = west                          
      rlat_end = rlat_start - float(iheight) * pix_latlon
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
          elseif(result(i,j) .le. 255.)then
            albedo(ic,i,j) = 0.25 + (result(i,j) - 179.) / 122. 
          else ! bad / default value
            albedo(ic,i,j) = 0.2
          endif
        enddo ! j
        enddo ! i
      enddo ! ic

      write(6,*)' Interpolated albedo: ',albedo(:,ni/2,nj/2)

      if(albedo(2,ni/2,nj/2) .gt. 1.0)then
        write(6,*)' ERROR bad interpolated albedo point'
        goto 999
      endif

!     Write albedo file in binary format
      u_out = 12
      file_bm=trim(directory)//'albedo_multispectral.dat'
      write(6,*)' Writing albedo file to '//file_bm
      open(u_out,file=trim(file_bm),form='unformatted'
     1    ,status='new',err=999)
      write(u_out)albedo
      close(u_out)

      endif ! binary file exists

!     Normal end
      istatus = 1
      goto 9999

!     Error condition
999   istatus = 0
      write(6,*)' error in land_albedo_bm'

9999  if(allocated(rlat_img))deallocate(rlat_img)
      if(allocated(rlon_img))deallocate(rlon_img)
      if(allocated(img))deallocate(img)
      if(allocated(array_2d))deallocate(array_2d)

      return
      end

      subroutine get_nlights(ni,nj,grid_spacing_m,lat,lon,gnd_glow)

!     Read in VIIRS imagery for night lights.
!     Merging of code in 'get_sfc_glow' and 'land_albedo_bm'

      use ppm

      real lat(ni,nj)
      real lon(ni,nj)
      real sfc_glow_c(3,ni,nj) ! Sfc Glow (Red, Green, Blue)
      real result(ni,nj)       ! Interpolated image channel
      real sfc_glow(ni,nj)     ! surface lighting intensity of clouds (nl)
      real gnd_glow(ni,nj)     ! ground lighting intensity (nl)      
      character*255 directory
      character*255 file,file_bm 
      character*10  c10_fname /'nest7grid'/
      integer u,u_out
      logical l_there

      integer ncol,iwidth,iheight
      parameter (ncol=3)

      real, allocatable :: rlat_img(:,:),rlon_img(:,:) 
      real ri_img(ni,nj)            
      real rj_img(ni,nj)             
      real rlat_laps(ni,nj)
      real rlon_laps(ni,nj)
      real topo_laps_dum(ni,nj)

!     integer img(ncol,iwidth,iheight)                                               
      integer, allocatable :: img(:,:,:)                   
!     real array_2d(iwidth,iheight)
      real, allocatable :: array_2d(:,:)                   
      u = 11

      write(6,*)' Subroutine get_nlights...'

      call get_directory('static',directory,len_dir)
      file_bm=trim(directory)//'nlights_multispectral.dat'
      inquire(file=trim(file_bm),exist=l_there)

      if(l_there)then
        write(6,*)' Nlights binary remapped sfc_glow_c file exists'
        open(u,file=trim(file_bm),form='unformatted' 
     1      ,status='old',err=999)
        write(6,*)' Successful open, now reading'
        read(u,err=999)sfc_glow_c
        close(u)
        write(6,*)' first point is ',sfc_glow_c(1,1,1)
        if(sfc_glow_c(1,1,1) .lt. 0.0)then
          write(6,*)' ERROR bad first point'
          goto 999
        endif
        sfc_glow = sfc_glow_c(2,:,:)

      else ! binary file is absent
        write(6,*)
     1' Nlights binary file is absent, generate and create it'       

       file=trim(directory)//'viirs_crop.ppm'

       write(6,*)' Open for reading ',trim(file)

!      Read section of VIIRS Image in PPM format
       open(u,file=trim(file),status='old',err=999)
       read(u,*)   
       read(u,*)iwidth,iheight
       read(u,*)ncol2
       rewind(u)
       write(6,*)' dynamic dims ',ncol,iwidth,iheight
       allocate(img(ncol,iwidth,iheight))
       call read_ppm (u,img,ncol,iwidth,iheight)
       close(u)

       allocate(rlat_img(iwidth,iheight))
       allocate(rlon_img(iwidth,iheight))
       allocate(array_2d(iwidth,iheight))

!      Consider dynamic means to get rlat_start and rlon_start
!      Use domain lat/lon bounds with a 0.2 deg cushion
       call get_domain_perimeter(ni,nj,c10_fname  
     1                  ,rlat_laps,rlon_laps,topo_laps_dum
     1                  ,0.2,rnorth,south,east,west,istatus)

       write(6,*)' NSEW',rnorth,south,east,west

       pix_latlon = 1. / 240.
!      rlat_start =   90. - 11000. * pix_latlon
!      rlon_start = -180. + 17000. * pix_latlon
       rlat_start = rnorth                         
       rlon_start = west                          
       rlat_end = rlat_start - float(iheight)  * pix_latlon
       rlon_end = rlon_start + float(iwidth)  * pix_latlon

       write(6,*)' NL PPM lat range: ',rlat_start,rlat_end
       write(6,*)' NL PPM lon range: ',rlon_start,rlon_end

!      Fill arrays of image lat/lons and LAPS ri,rj
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

       do i = 1,iwidth,400
       do j = 1,iheight,400

        write(6,11)img(:,i,j),i,j
11      format(5i6,2f8.3)

       enddo ! j
       enddo ! i

       write(6,*)' img max = ',maxval(img)

       I4_elapsed = ishow_timer()

       write(6,*)' Bilinearly Interpolate'

!      Interpolate to LAPS grid using bilinear_laps_2d
       do ic = 1,3
        array_2d(:,:) = img(1,:,:)
        call bilinear_laps_2d(ri_img,rj_img,iwidth,iheight,ni,nj 
     1                       ,array_2d,result)
        do i = 1,ni
        do j = 1,nj
!         Convert units and such
          sfc_glow(i,j) = result(i,j)
          sfc_glow_c(ic,i,j) = result(i,j)
        enddo ! j
        enddo ! i
       enddo ! ic

       write(6,*)' range of sfc_glow is',minval(sfc_glow)
     1                                  ,maxval(sfc_glow)      

cdoc   Interpolate 2-d array to find the field values at fractional grid
cdoc   points.

!      real array_2d(imax,jmax)                         ! I
!      real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
!      real result(nx_laps,ny_laps)                     ! O

       write(6,*)' Interpolated sfc_glow_c: ',sfc_glow_c(:,ni/2,nj/2)

       if(sfc_glow_c(2,ni/2,nj/2) .lt. 0.0)then
        write(6,*)' ERROR bad interpolated sfc_glow_c point'
        goto 999
       endif

!      Write sfc_glow_c file in binary format
       u_out = 12
       file_bm=trim(directory)//'nlights_multispectral.dat'
       write(6,*)' Writing sfc_glow_c file to '//file_bm
       open(u_out,file=trim(file_bm),form='unformatted'
     1     ,status='new',err=999)
       write(u_out)sfc_glow_c
       close(u_out)

      endif ! binary file exists

      gnd_glow = sfc_glow

!     Normal end
      istatus = 1
      goto 9999

!     Error condition
999   istatus = 0
      write(6,*)' error in get_nlights'

9999  if(allocated(rlat_img))deallocate(rlat_img)
      if(allocated(rlon_img))deallocate(rlon_img)
      if(allocated(img))deallocate(img)
      if(allocated(array_2d))deallocate(array_2d)

      return
      end

      subroutine compare_land_albedo(lu,ni,nj,albedo_usgs,albedo_bm
     1                                       ,albedo_static)

!     This routine compares USGS, Blue Marble, and static albedo fields

      real lu(ni,nj)               ! Land Use (USGS 24 category)
      real albedo_usgs(3,ni,nj)    ! Albedo (Red, Green, Blue)
      real albedo_bm(3,ni,nj)      ! Albedo (Red, Green, Blue)
      real albedo_static(ni,nj)    ! Albedo (broadband)


      write(6,*)' subroutine compare_land_albedo...'

      rpts = ni*nj
      albedo_mean_bm = sum(albedo_bm(2,:,:))/rpts
      albedo_mean_usgs = sum(albedo_usgs(2,:,:))/rpts
      albedo_mean_static = sum(albedo_static(:,:))/rpts

      write(6,11)albedo_mean_bm,albedo_mean_usgs,albedo_mean_static
11    format('  Mean albedo bm/usgs/static ',3f9.4)

      return
      end 
