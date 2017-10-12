

      subroutine get_topo_1s(ni,nj,grid_spacing_m,r8lat,r8lon,topo
     1                      ,rnorth,south,east,west                      ! I
     1                      ,path_to_topt1s,istatus)

!     Read in topo tiles
!     Presently called from 'gridgen_model' software
      use ppm
      use mem_namelist, ONLY: c6_maproj 

      parameter (maxtiles = 100)

      real*8 r8lat(ni,nj), r8lon(ni,nj)
      real sfc_glow_c(3,ni,nj) ! Sfc Glow (Red, Green, Blue)
      real result(ni,nj)       ! Interpolated image channel
      real result_tile(ni,nj)  ! Interpolated image channel
      real sfc_glow(ni,nj)     ! surface lighting intensity of clouds (nl)
      real topo(ni,nj)     ! zenithal ground lighting intensity (wm2sr)
      real ri_img(ni,nj)
      real rj_img(ni,nj)
      character*255 directory
      character*255 file,file_bm,file_a(maxtiles)
      character*(*) path_to_topt1s
      character*10  c10_fname /'nest7grid'/
      character*20 cdum
      character*7 basename
      integer u,u_out
      logical l_there, l_global_nl

      integer ncol,iwidth,iheight
      parameter (ncol=3)

      real rlat_laps(ni,nj)
      real rlon_laps(ni,nj)
      real topo_laps_dum(ni,nj)
      real pixsum(ni,nj)
      integer npix(ni,nj)

      real nwcm2sr_to_wm2sr
      nwcm2sr_to_wm2sr(x) = x * 1e-5

!     integer img(ncol,iwidth,iheight)                                               
      integer, allocatable :: img(:,:,:)                   
!     real array_2d(iwidth,iheight)
      real, allocatable :: array_2d(:,:)                   

!     ri_mdl/rj_mdl tells us the model grid point of each image pixel
!     they are allocated on the image grid
      real, allocatable :: ri_mdl(:,:),rj_mdl(:,:)                   
      real, allocatable :: rlat_img(:,:),rlon_img(:,:) 
      u = 11

      rlat_laps = r8lat; rlon_laps = r8lon

      write(6,*)' Subroutine get_topo_1s...',ni,nj
      if(nj .le. 0)then
          write(6,*)' returning early'
          istatus = 0
          return
      endif

      call get_r_missing_data(r_missing_data,istatus)
      result = r_missing_data ! initialize

      call get_directory('static',directory,len_dir)
      file_bm=trim(directory)//'topo1s_tile.dat'
      inquire(file=trim(file_bm),exist=l_there)
      write(6,*)' File being inquired is ',trim(file_bm),' ',l_there

      if(l_there)then ! this block not presently used
        write(6,*)' topo_1s tile remapped topo file exists'
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
       write(6,*)' setup array of topo_1s asc files'

       if(.true.)then

         pix_latlon = 1. / 240.
         offset_lat = -.000  ! positional error in remapping
         offset_lon = +.000  !             "
         l_global_nl = .false.

         write(6,*)' perimeter NSEW:',rnorth,south,east,west
         ilat_start = int(south) + 1
         ilat_end   = int(rnorth) + 1

         ilon_start = -floor(east)
         ilon_end   = -floor(west)

         write(6,*)' ilat range dynamic',ilat_start,ilat_end
         write(6,*)' ilon range dynamic',ilon_start,ilon_end

         itile = 0
         do ilat = ilat_start,ilat_end
         do ilon = ilon_start,ilon_end
             itile = itile + 1
             write(basename,11)ilat,ilon
11           format('n',i2,'w',i3)          
             file_a(itile) =
     1           trim(path_to_topt1s)//'/'//basename//'.asc'
             ntiles = itile
         enddo ! ilon
         enddo ! ilat

       endif

!      Process asc tiles derived from ArcGIS 30m topo data
!      For example, gdal_translate can convert 'adf' files to 'asc'        
       do itile=1,ntiles

        file=file_a(itile)

        write(6,*)' Open for reading ',trim(file)
        write(6,*)' pix_latlon,l_global_nl = ',pix_latlon,l_global_nl

!       Read terrain tile in PPM format
        open(u,file=trim(file),status='old',err=999)
        read(u,*)cdum,iwidth
        read(u,*)cdum,iheight
        read(u,*)cdum,xllcorner
        read(u,*)cdum,yllcorner
        read(u,*)cdum,cellsize
        read(u,*)cdum,value_nodata
        write(6,*)' dynamic dims ',iwidth,iheight
        if(itile .eq. 1)allocate(array_2d(iwidth,iheight))
        write(6,*)' reading array_2d'
        read(u,*)array_2d
        close(u)
        write(6,*)' 1s tile read successful ',itile,ntiles

        if(itile .eq. 1)then
          allocate(rlat_img(iwidth,iheight))
          allocate(rlon_img(iwidth,iheight))
          allocate(ri_mdl(iwidth,iheight))
          allocate(rj_mdl(iwidth,iheight))
        endif

        rlat_start = yllcorner ! rsouth 
        rlon_start = xllcorner ! west   
        rlat_end = rlat_start + float(iheight) * cellsize
        rlon_end = rlon_start + float(iwidth)  * cellsize

        write(6,*)' TILE ASC lat range: ',rlat_start,rlat_end
        write(6,*)' TILE ASC lon range: ',rlon_start,rlon_end

        pix_latlon = cellsize

!       Fill arrays of image ri,rj for each LAPS grid point
        istatus_img = 1
        do i = 1,ni ! Loop through LAPS grid points    
        do j = 1,nj     

!        The dimensions of this are questionable
         rj_mdl(i,j) = (rlat_laps(i,j) - rlat_start) / pix_latlon
         ri_mdl(i,j) = (rlon_laps(i,j) - rlon_start) / pix_latlon

         rj_img(i,j) = (rlat_laps(i,j) - rlat_start) / pix_latlon
         rj_img(i,j) = 3600. - rj_img(i,j) + 1.0
         ri_img(i,j) = (rlon_laps(i,j) - rlon_start) / pix_latlon

         if(rlat_laps(i,j) .gt. rlat_start .OR.
     1      rlat_laps(i,j) .lt. rlat_end   .OR.
     1      rlon_laps(i,j) .lt. rlon_start .OR. 
     1      rlon_laps(i,j) .gt. rlon_end)then
             istatus_img = 0
         endif
        enddo ! j
        enddo ! i

        if(istatus_img .eq. 0)then
          write(6,*)' NOTE: LAPS grid extends outside image tile'
        else
          write(6,*)' LAPS grid is contained within image tile'      
        endif

        write(6,*)' selected image topo points for tile',itile
        write(6,*)'    i     j   j1    topo       lat      lon'
        do j = 1,iheight,400
        do i = 1,iwidth,400

          j1 = iheight+1-j  
          rlon_img_tmp = rlon_start + float(i) * pix_latlon
          rlat_img_tmp = rlat_start + float(j) * pix_latlon

          write(6,21)i,j,j1,array_2d(i,j),rlat_img_tmp,rlon_img_tmp
21        format(3i6,3f9.2)

        enddo ! i
        enddo ! j

        write(6,*)' array_2d range before interpolation (topo) = '
     1          ,minval(array_2d),maxval(array_2d)

        I4_elapsed = ishow_timer()

        reslights_m = 30.

!       Interpolate to LAPS grid using bilinear_laps_2d
!       array_2d(:,:) = img(1,:,:)

        if(grid_spacing_m / reslights_m .lt. 1.5)then
!       if(.false.)then
          write(6,*)' Bilinearly Interpolate Tile ',itile,grid_spacing_m
     1                                             ,reslights_m
          result_tile = r_missing_data

!         Find interpolated image values on the model grid points that
!         overlap the image. The output is on the model grid. 'array_2d'
!         is the image terrain. 'ri_img' and 'rj_img' are dimensioned 
!         on the model grid and represent ri/rj image values.          
          call bilinear_laps_2d(ri_img,rj_img,iwidth,iheight,ni,nj 
     1                         ,array_2d,result_tile)

          write(6,*)' range of result_tile ',minval(result_tile)
     1                                      ,maxval(result_tile)

          where(result_tile(:,:) .ne. r_missing_data)
            result(:,:) = result_tile(:,:)
          endwhere

          write(6,*)' range of result',minval(result),maxval(result)

        else ! pixel average interpolation

          write(6,*)' Calls to latlon_to_rlapsgrid'
           
          if(.true.)then

!           Loop through image array
            do i = 1,iwidth
            do j = 1,iheight

!             Fill lat/lon of each image pixel
              j1 = iheight+1-j  
              rlon_img(i,j1) = rlon_start + float(i) * pix_latlon
              rlat_img(i,j1) = rlat_start + float(j) * pix_latlon

!             Fill LAPS ri/rj of each image pixel
              call latlon_to_rlapsgrid(rlat_img(i,j1),rlon_img(i,j1)
     1                                ,rlat_laps,rlon_laps,ni,nj
     1                                ,ri_mdl(i,j1),rj_mdl(i,j1)
     1                                ,istatus)

              if(itile .eq. 1)then
                if(j .eq. iheight .and. i .eq. 30*(i/30))then
                  write(6,*)i,j,j1,ri_mdl(i,j1),rj_mdl(i,j1)
                endif
              endif

            enddo ! j
            enddo ! i

          endif 

          write(6,*)' Pixel Ave Interpolate ',ni,nj,grid_spacing_m
     1                                       ,reslights_m
!         Check dimensions and purpose of each array
          call pixelave_interp_multi(ri_mdl,rj_mdl,iwidth,iheight 
     1         ,ni,nj
     1         ,itile,ntiles,pixsum,npix
     1         ,0.,array_2d,result)
          write(6,*)' interp terrain range is ',minval(result)
     1                                         ,maxval(result)
          i = ni/2
          do j = 1,nj
             write(6,*)i,j,rlat_laps(i,j),result(i,j)
          enddo ! j
        endif ! interpolation method

       enddo ! itiles

cdoc   Interpolate 2-d array to find the field values at fractional grid
cdoc   points.

      endif ! binary file exists

      topo = result

!     Normal end
      istatus = 1
      write(6,*)' normal return from get_topo_1s'
      goto 9999

!     Error condition
999   istatus = 0
      write(6,*)' error in get_topo_1s (check missing file)'

9999  if(allocated(rlat_img))deallocate(rlat_img)
      if(allocated(rlon_img))deallocate(rlon_img)
      if(allocated(img))deallocate(img)
      if(allocated(array_2d))deallocate(array_2d)
      if(allocated(ri_mdl))deallocate(ri_mdl)
      if(allocated(rj_mdl))deallocate(rj_mdl)

      return
      end

      subroutine pixelave_interp_multi(ri_mdl,rj_mdl,ni1,nj1,ni2,nj2
     1     ,itile,ntiles,pixsum,npix
     1     ,r_missing_data,array_2d,result)

      real ri_mdl(ni1,nj1)    ! image grid
      real rj_mdl(ni1,nj1)    ! image grid
      real pixsum(ni2,nj2)    ! model grid
      integer npix(ni2,nj2)   ! model grid
      real array_2d(ni1,nj1)  ! image grid
      real result(ni2,nj2)    ! model grid

      write(6,*)' ni2/nj2 (mdl grid) = ',ni2,nj2

      if(itile .eq. 1)then
        pixsum = 0.
        npix = 0
      endif

      rimin = 0.5
      rimax = float(ni2) + 0.5
      rjmin = 0.5
      rjmax = float(nj2) + 0.5

      jtest = nj1

      do i = 1,ni1
      do j = 1,nj1

        ri = ri_mdl(i,j)
        rj = rj_mdl(i,j)

        if(i .eq. 30*(i/30) .and. j .eq. jtest)then
          write(6,*)i,j,ri_mdl(i,j),rj_mdl(i,j),rimin,rimax,rjmin,rjmax
        endif

        if(ri .gt. rimin .and. ri .lt. rimax .and.
     1     rj .gt. rjmin .and. rj .lt. rjmax      )then

          iout = nint(ri)
          jout = nint(rj)
       
          npix(iout,jout) = npix(iout,jout) + 1
          pixsum(iout,jout) = pixsum(iout,jout) + array_2d(i,j)

          if(itile .eq. 1)then
            if(i .eq. 30*(i/30) .and. j .eq. jtest)then
              write(6,*)i,j,ri,rj,npix(iout,jout),pixsum(iout,jout)
            endif
          endif
 
        endif
      enddo ! j
      enddo ! i

      if(itile .eq. ntiles)then
        do iout = 1,ni2
        do jout = 1,nj2
          if(npix(iout,jout) .gt. 0)then
            result(iout,jout) = pixsum(iout,jout) 
     1                        / float(npix(iout,jout))
          else
            result(iout,jout) = r_missing_data
            write(6,*)' WARNING: npix = 0 at ',iout,jout
          endif
          if(i .eq. 30*(i/30) .and. jout .eq. jtest)then
            write(6,*)'iout/jout/result',iout,jout,pixsum(iout,jout)
     1                            ,npix(iout,jout),result(iout,jout)
          endif
        enddo ! jout
        enddo ! iout      
      endif

      return
      end 
