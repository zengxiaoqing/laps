
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

      subroutine land_albedo_bm(rlat_laps,rlon_laps,ni,nj,i4time
     1                         ,albedo,bm_counts,istatus)

      use ppm
      use mem_namelist, ONLY: c6_maproj, grid_spacing_m
      use mem_allsky, ONLY: nc

      include 'wa.inc'

!     Read and interpolate from sector of Blue Marble Image
!     See code in /scratch/staging/fab/albers/nasa

      real albedo(3,ni,nj)    ! Albedo (Red, Green, Blue)
      real bm_counts(3,ni,nj) ! Counts (Red, Green, Blue)
      real albedo_buff(3,ni,nj) 
      real result(ni,nj)      ! Interpolated image channel
      character*13 cvt_i4time_wfo_fname13,c13_time
      character*2 c2_mn

!     convert -crop 2000x2000+12000+9000 -compress none
!     world.200404.3x21600x21600.A1.png world.200404.3x21600x21600.A1.crop2.ppm
!     convert -compress none
!     world.200408.3x5400x2700.png world.200408.3x5400x2700.ppm

      character*255 directory
      character*255 file            ! Blue Marble image data
      character*255 file_bm,file_dc ! remapped to model grid
      character*10  c10_fname /'nest7grid'/
      integer u,u_out
      logical l_there, l_global_bm, l_there_dc

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

      integer, allocatable :: img(:,:,:)                   
      real, allocatable :: array_2d(:,:)                   
      u = 11

      idb = ni/2 ! min(183,ni)
      jdb = nj/2 ! min(123,nj)

      c13_time = cvt_i4time_wfo_fname13(i4time)
      c2_mn = c13_time(5:6)

      write(6,*)' Subroutine land_albedo_bm...',c2_mn

      call get_directory('static',directory,len_dir)

      file_bm=trim(directory)//'albedo_multispectral_'//c2_mn//'.dat'
      inquire(file=trim(file_bm),exist=l_there)
      write(6,*)' File being inquired is ',trim(file_bm),' ',l_there

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
     1' Blue Marble binary file is absent, generate and create from ppm'

        if(grid_spacing_m .ge. 2500. .and.
     1     grid_spacing_m .le. 7400.      )then ! mid-range of grid spacing
          file=trim(directory)//'world.2004'
     1                      //c2_mn//'.global_montage.20.ppm'
          inquire(file=trim(file),exist=l_there)
          write(6,*)' File being inquired is ',trim(file),' ',l_there
        else ! outside mid-range of grid spacing
          l_there = .false.
        endif

        file_dc=trim(directory)//'vhires_dc_crop.ppm'       
        inquire(file=trim(file_dc),exist=l_there_dc)
        write(6,*)' File being inquired is ',trim(file_dc),' '
     1                                      ,l_there_dc      

        if(l_there_dc)then                    ! Descartes data
          pix_latlon_we = 1. / 865.954              
          pix_latlon_sn = 1. / 1130.26
          l_global_bm = .false.
          file = trim(file_dc)
          perimeter = 0.05

        elseif(l_there)then                   ! local domain (2.5km pixels)
          pix_latlon_we = 1. / 48.            ! 17280x8640 image
          pix_latlon_sn = pix_latlon_we       ! for 2.5km to 7.4km grid
          l_global_bm = .true.
          perimeter = 0.2

        elseif(c6_maproj .ne. 'latlon' .and.
     1         grid_spacing_m .le. 2500.)then ! local domain (~500m pixels)
                                              ! for 500m to 2.5km grid
          write(6,*)' Looking for 500m BM tile'
          file=trim(directory)//'world.2004' 
     1                        //c2_mn//'.3x21600x21600.crop.ppm'
          pix_latlon_we = 1. / 240. !                 (90x180 degree tile)
          pix_latlon_sn = pix_latlon_we 
          l_global_bm = .false.
          perimeter = 0.2

        elseif(grid_spacing_m .ge. 7400.)then ! 7.4km pixels for >7.4km grid
                                              ! (or latlon / global projection)
          file=trim(directory)//'world.2004'//c2_mn//'.3x5400x2700.ppm'
          pix_latlon_we = 1. / 15.                             
          pix_latlon_sn = pix_latlon_we 
          l_global_bm = .true.
          perimeter = 0.2

        else
          write(6,*)' conditions unsatisifed for BMNG land albedo data'
          goto 999

        endif

        write(6,*)' Open for reading ',trim(file)

!       Read section of NASA Blue Marble Image in PPM format
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

!       Consider dynamic means to get rlat_start and rlon_start
!       Use domain lat/lon bounds with a 0.2 deg cushion
        call get_domain_perimeter(ni,nj,c10_fname  
     1                  ,rlat_laps,rlon_laps,topo_laps_dum
     1                  ,perimeter,rnorth,south,east,west,istatus)

        write(6,5)perimeter,rnorth,south,east,west
5       format('  NSEW Domain perimeter',5f9.3)                 

        if(l_global_bm)then
          rlat_start = +90.                         
          rlon_start = -180.
          rlat_end = rlat_start - float(iheight) * pix_latlon_sn
          rlon_end = rlon_start + float(iwidth)  * pix_latlon_we
        else
          rlat_start = rnorth                         
          rlon_start = west                          
          rlat_end = rlat_start - float(iheight) * pix_latlon_sn
          rlon_end = rlon_start + float(iwidth)  * pix_latlon_we
        endif

        write(6,*)' BM PPM lat range: ',rlat_start,rlat_end
        write(6,*)' BM PPM lon range: ',rlon_start,rlon_end

!       Fill arrays of image lat/lons and LAPS ri,rj
        istatus_img = 1
        do i = 1,ni     
        do j = 1,nj     
          rj_img(i,j) = (rlat_start - rlat_laps(i,j)) / pix_latlon_sn
          ri_img(i,j) = (rlon_laps(i,j) - rlon_start) / pix_latlon_we
          if(rlat_laps(i,j) .gt. rlat_start .OR.
     1       rlat_laps(i,j) .lt. rlat_end   .OR.
     1       rlon_laps(i,j) .lt. rlon_start .OR. 
     1       rlon_laps(i,j) .gt. rlon_end)then
            istatus_img = 0
          endif
        enddo ! j
        enddo ! i

        write(6,*)' Image location of domain point',idb,jdb
     1                                ,ri_img(idb,jdb),rj_img(idb,jdb)

        write(6,*)' Image location of SW domain corner',1,1
     1                                ,ri_img(1,1),rj_img(1,1)        

        write(6,*)' Image location of NW domain corner',1,nj
     1                                ,ri_img(1,nj),rj_img(1,nj)        

        if(istatus_img .eq. 0)then
          write(6,*)' WARNING: LAPS grid extends outside image'
        else
          write(6,*)' LAPS grid is contained within image'      
        endif

        write(6,*)' Sampling of image ',rlat_start,rlon_start
        do i = 1,iwidth,120
!       do j = 1,iheight,200
        do j = iheight/2,iheight/2,1

          rlat_img_e = rlat_start - float(j) * pix_latlon_sn
          rlon_img_e = rlon_start + float(i) * pix_latlon_we

          write(6,11)i,j,rlat_img_e,rlon_img_e,img(:,i,j)
11        format(2i6,2x,2f9.3,2x,3i6)

        enddo ! j
        enddo ! i

        I4_elapsed = ishow_timer()

        write(6,*)' Bilinearly Interpolate'

!       Interpolate to LAPS grid using bilinear_laps_2d
        do ic = 1,3
          array_2d(:,:) = img(ic,:,:)
          call bilinear_laps_2d(ri_img,rj_img,iwidth,iheight,ni,nj 
     1                         ,array_2d,result)
          do i = 1,ni
          do j = 1,nj

!           http://earthobservatory.nasa.gov/Features/BlueMarble/bmng.pdf

            if(result(i,j) .le. 179.)then ! Based on NASA spline
!             albedo(ic,i,j) = result(i,j) / 716.                  
              albedo(ic,i,j) = .0010 * result(i,j)
     1                       + 6.92e-11 * result(i,j)**4.0
            elseif(result(i,j) .le. 255.)then
!             albedo(ic,i,j) = 0.25 + (result(i,j) - 179.) / 122. 
              albedo(ic,i,j) = 0.25 + (result(i,j) - 179.)    * .00258
     1                              + (result(i,j) - 179.)**2 * .000092
            else ! bad / default value
              albedo(ic,i,j) = 0.2
            endif

            bm_counts(ic,i,j) = result(i,j)

          enddo ! j
          enddo ! i
          write(6,*)' Color / Count / Albedo (observer)',ic
     1                          ,result(idb,jdb),albedo(ic,idb,jdb)          
        enddo ! ic

        write(6,*)' Initial albedo: ',albedo(:,idb,jdb)

!       Adjust colors based on wavelengths
!       BM data is .645, .555, .470 microns
        albedo_buff(:,:,:) = albedo(:,:,:)
        red_f = (wa(1) - wa(2)) / (.645 - wa(2))
        grn_f = 1.0 - red_f
        write(6,*)' Adjusting R color based on wavelengths',red_f,grn_f
        albedo(1,:,:) = albedo_buff(1,:,:) * red_f
     1                + albedo_buff(2,:,:) * grn_f    

        blu_f = (wa(3) - wa(2)) / (.470 - wa(2))
        grn_f = 1.0 - blu_f
        write(6,*)' Adjusting B color based on wavelengths',blu_f,grn_f
        albedo(3,:,:) = albedo_buff(3,:,:) * blu_f
     1                + albedo_buff(2,:,:) * grn_f    
        
        write(6,*)' Adjusted albedo: ',albedo(:,idb,jdb)

        if(albedo(2,idb,jdb) .gt. 1.0)then
          write(6,*)' ERROR bad interpolated albedo point'
          goto 999
        endif

!       Write albedo file in binary format
        u_out = 12
        write(6,*)' Writing albedo file to '//file_bm
        open(u_out,file=trim(file_bm),form='unformatted'
     1      ,status='new',err=999)
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

      subroutine get_nlights(ni,nj,grid_spacing_m,rlat_laps,rlon_laps
     1                      ,gnd_glow)

!     Read in VIIRS imagery for night lights.
!     Merging of code in 'get_sfc_glow' and 'land_albedo_bm'

      use ppm
      use mem_namelist, ONLY: c6_maproj 

      real sfc_glow_c(3,ni,nj) ! Sfc Glow (Red, Green, Blue)
      real result(ni,nj)       ! Interpolated image channel
      real sfc_glow(ni,nj)     ! surface lighting intensity of clouds (nl)
      real gnd_glow(ni,nj)     ! zenithal ground lighting intensity (wm2sr)
      character*255 directory
      character*255 file,file_bm 
      character*10  c10_fname /'nest7grid'/
      integer u,u_out
      logical l_there, l_global_nl

      integer ncol,iwidth,iheight
      parameter (ncol=3)

!     ri_img/rj_img tells us the image pixel of each laps grid point
      real ri_img(ni,nj)            
      real rj_img(ni,nj)             
      real rlat_laps(ni,nj)
      real rlon_laps(ni,nj)
      real topo_laps_dum(ni,nj)

      real nwcm2sr_to_wm2sr
      nwcm2sr_to_wm2sr(x) = x * 1e-5

!     integer img(ncol,iwidth,iheight)                                               
      integer, allocatable :: img(:,:,:)                   
!     real array_2d(iwidth,iheight)
      real, allocatable :: array_2d(:,:)                   
!     ri_laps/rj_laps tells us the laps grid point of each image pixel
      real, allocatable :: rlat_img(:,:),rlon_img(:,:) 
      real, allocatable :: ri_laps(:,:)                   
      real, allocatable :: rj_laps(:,:)                   

      write(6,*)' Subroutine get_nlights...'

      u = 11
      idb = min(893,ni) ! min(878,ni) or ni/2
      jdb = min(583,nj) ! min(383,nj) or nj/2

      call get_directory('static',directory,len_dir)
      file_bm=trim(directory)//'nlights_multispectral.dat'
      inquire(file=trim(file_bm),exist=l_there)
      write(6,*)' File being inquired is ',trim(file_bm),' ',l_there

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
     1' Nlights binary file is absent, generate and create it from ppm'       

       file=trim(directory)//'viirs_global_montage_20.ppm' 
       inquire(file=trim(file),exist=l_there)
       if(l_there)then                     ! local domain (2.5km pixels)
         pix_latlon = 1. / 48.             ! 17280x8640 image
         offset_lat = 0.     ! positional error in remapping
         offset_lon = 0.     !             "
         l_global_nl = .true.

       elseif(c6_maproj .ne. 'latlon')then ! local domain (500m pixels)
         file=trim(directory)//'viirs_crop.ppm' !         (120 degree tile)
         pix_latlon = 1. / 240.
         offset_lat = -.008  ! positional error in remapping
         offset_lon = +.0035 !             "
         l_global_nl = .false.

       else                     ! latlon global projection (9km pixels)
         file=trim(directory)//'viirs_global_montage.ppm' 
         pix_latlon = 1. / 12.
         offset_lat = 0.     ! positional error in remapping
         offset_lon = 0.     !             "
         l_global_nl = .true.

       endif

       write(6,*)' Open for reading ',trim(file)
       write(6,*)' pix_latlon,l_global_nl = ',pix_latlon,l_global_nl

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
       allocate(ri_laps(iwidth,iheight))
       allocate(rj_laps(iwidth,iheight))

!      Consider dynamic means to get rlat_start and rlon_start
!      Use domain lat/lon bounds with a 0.2 deg cushion
       if(l_global_nl)then
          rnorth = +90.
          south = -90.
          west = -180.
          east = +180.
       else
          call get_domain_perimeter(ni,nj,c10_fname  
     1                  ,rlat_laps,rlon_laps,topo_laps_dum
     1                  ,0.2,rnorth,south,east,west,istatus)
       endif

       write(6,*)' NSEW',rnorth,south,east,west

       rlat_start = rnorth - offset_lat                         
       rlon_start = west   - offset_lon                       
       rlat_end = rlat_start - float(iheight) * pix_latlon
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

       write(6,*)' img max before interpolation (nwcm2sr) = '
     1          ,maxval(img)

       I4_elapsed = ishow_timer()

       reslights_m = 10000.

!      Interpolate to LAPS grid using bilinear_laps_2d
       do ic = 1,3
        array_2d(:,:) = img(1,:,:)

        if(grid_spacing_m / reslights_m .lt. 1.5)then
!       if(.true.)then
          write(6,*)' Bilinearly Interpolate ',ic,grid_spacing_m
     1                                        ,reslights_m
          write(6,*)' lat/lon is ',rlat_laps(idb,jdb),rlon_laps(idb,jdb)
          write(6,*)' rij_img is ',idb,jdb,ri_img(idb,jdb)
     1                                    ,rj_img(idb,jdb)
          iimg = ri_img(idb,jdb)
          jimg = rj_img(idb,jdb)
          write(6,*)' array_2d is ',array_2d(iimg:iimg+1,jimg:jimg+1)

          call bilinear_laps_2d(ri_img,rj_img,iwidth,iheight,ni,nj 
     1                         ,array_2d,result)

          write(6,*)' result is ',result(idb,jdb)
     1                           ,nwcm2sr_to_wm2sr(result(idb,jdb))
          
        else ! pixel average interpolation

          if(ic .eq. 1)then

!           Loop through image array
            do i = 1,iwidth
            do j = 1,iheight

!             Fill lat/lon of each image pixel
              rlon_img(i,j) = rlon_start + float(i) * pix_latlon
              rlat_img(i,j) = rlat_start - float(j) * pix_latlon

!             Fill LAPS ri/rj of each image pixel
              call latlon_to_rlapsgrid(rlat_img(i,j),rlon_img(i,j)
     1                                ,rlat_laps,rlon_laps,ni,nj
     1                                ,ri_laps(i,j),rj_laps(i,j)
     1                                ,istatus)

            enddo ! j
            enddo ! i

          endif ! 1st color in loop

          write(6,*)' Pixel Ave Interpolate ',ic,grid_spacing_m
     1                                       ,reslights_m
          call pixelave_interp(ri_laps,rj_laps,iwidth,iheight,ni,nj 
     1                        ,0.,array_2d,result)
          write(6,*)' Max result value is (nanoWatts/cm2/sr)'
     1               ,maxval(result)
        endif ! interpolation method

!       This conversion should reduce values by 10^5
        do i = 1,ni
        do j = 1,nj
!         Input units are nanoWatts/cm2/sr
          arg1 = result(i,j)                     ! nanoWatts/cm2/sr
          sfc_glow(i,j) = nwcm2sr_to_wm2sr(arg1) ! wm2sr
          sfc_glow_c(ic,i,j) = sfc_glow(i,j)
        enddo ! j
        enddo ! i

       enddo ! ic

cdoc   Interpolate 2-d array to find the field values at fractional grid
cdoc   points.

!      real array_2d(imax,jmax)                         ! I
!      real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
!      real result(nx_laps,ny_laps)                     ! O

       write(6,*)' Interpolated sfc_glow_c (wm2sr): '
     1          ,sfc_glow_c(:,idb,jdb)

       if(sfc_glow_c(2,idb,jdb) .lt. 0.0)then
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

      write(6,*)' range of sfc_glow (wm2sr) is',minval(sfc_glow)
     1                                         ,maxval(sfc_glow)      

      gnd_glow = sfc_glow

!     Normal end
      istatus = 1
      write(6,*)' normal return from get_nlights'
      goto 9999

!     Error condition
999   istatus = 0
      gnd_glow = 0.
      write(6,*)' error in get_nlights'

9999  if(allocated(rlat_img))deallocate(rlat_img)
      if(allocated(rlon_img))deallocate(rlon_img)
      if(allocated(img))deallocate(img)
      if(allocated(array_2d))deallocate(array_2d)
      if(allocated(ri_laps))deallocate(ri_laps)
      if(allocated(rj_laps))deallocate(rj_laps)

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


      subroutine pixelave_interp(ri_img,rj_img,ni1,nj1,ni2,nj2
     1                          ,r_missing_data,array_2d,result)

      real ri_img(ni1,nj1)
      real rj_img(ni1,nj1)
      real pixsum(ni2,nj2)
      integer npix(ni2,nj2)
      real array_2d(ni1,nj1)
      real result(ni2,nj2)

      pixsum = 0.
      npix = 0

      rimin = 0.5
      rimax = float(ni2) + 0.5
      rjmin = 0.5
      rjmax = float(nj2) + 0.5

      do i = 1,ni1
      do j = 1,nj1

        ri = ri_img(i,j)
        rj = rj_img(i,j)

        if(ri .gt. rimin .and. ri .lt. rimax .and.
     1     rj .gt. rjmin .and. rj .lt. rjmax      )then

          iout = nint(ri)
          jout = nint(rj)
       
          npix(iout,jout) = npix(iout,jout) + 1
          pixsum(iout,jout) = pixsum(iout,jout) + array_2d(i,j)
 
        endif
      enddo ! j
      enddo ! i

      do iout = 1,ni2
      do jout = 1,nj2
        if(npix(iout,jout) .gt. 0)then
          result(iout,jout) = pixsum(iout,jout) 
     1                      / float(npix(iout,jout))
        else
          result(iout,jout) = r_missing_data
        endif
      enddo ! jout
      enddo ! iout      

      return
      end 
