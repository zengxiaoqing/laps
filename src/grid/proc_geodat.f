!     subroutine proc_geodat(nx,ny,lat,lon,)

!1.     given lat/lon points of analysis domain
!2.     determine the domain subset corner points (lat0 and lon0
!       lat0=lat(1,1)
!         where real*4 lat0,lon0,dlat,dlon     !SW corner lat, lon, lat, lon spacing
!3.     adjust for staggering is done via lat/lon 
!
!4.     use latlon_2_llij(nx*ny,lat,lon,gri,grj)
!         with common /llgrid/ni,nj,nk,lat0,lon0,dlat,dlon
!         set appropriately for the raw (30s or whatever) data
!         ni/nj = num raw data gridpoints covering the domain
!5.     use read_usgs_geog (formerly read_usgs_veg.F)
!         to extract that part of the full-sized raw data
!         relevant to the domain.  This routine will need
!         to be written (starting with read_usgs_veg.F)
!6.     need function that determines the starting/ending i/j
!         values for the raw data surrounding the domain lat/lon point.
!
!         grid_ratio=(resolution input)/(resolution output)

!         function istartg(r_grid_ratio,gri)
!            is = nint(gri-(1./r_grid_ratio) * 0.5)
!            return
!            end

!         function iendg(r_grid_ratio,gri)
!            ie = int(gri+(1./r_grid_ratio) * 0.5)
!            return
!            end

!         function jstartg(r_grid_ratio,grj)
!            js = nint(grj-(1./r_grid_ratio) * 0.5)
!            return
!            end

!         function jendg(r_grid_ratio,grj)
!            je = int(grj+(1./r_grid_ratio) * 0.5)
!            return
!            end

!7.     now do this:
!
!       do j=1,ny
!          do i=1,nx
!             
!             is=istartg(r_grid_ratio,gri(i,j)) 
!             ie=iendg(r_grid_ratio,(gri(i,j))
!             js=jstartg(r_grid_ratio,gri(i,j))
!             je=jendg(r_grid-ratio,gri(i,j))
!             do jd=js,je
!                do id=is,ie
!                   dindex=rawdata(id,jd)
!                   dcnt(i,j,dindex)=dcnt(i,j,dindex)+1
!                   totnum(i,j)=totnum(i,j)+1 
!                enddo
!             enddo

!          enddo
!       enddo

c-----------------------------------------------


      subroutine proc_geodat(nx_dom,ny_dom,ncat
     1,path_to_tile_data,dom_lats,dom_lons,lmask_out
     1,geodat,istatus)
!    1,cat_pct)

      use horiz_interp

      implicit none 
      integer maxtiles 
      integer nx_dom 
      integer ny_dom
      integer ncat
      integer ntn
      integer itilesize_d
      integer lentd,lenp,lenf
      integer iblksizo
      integer isbego,iwbego
      integer no
      integer icurNS,icurEW
      integer itile
      integer itilex
      integer itiley
      integer,allocatable:: itile_lat(:)
      integer,allocatable:: itile_lon(:)
      integer itilelon
      integer it,jt
      integer is,js
      integer ie,je
      integer itx,ity
      integer i,j,k,ii,jj
      integer ix,iy
      integer icnta
      integer icat

      integer istatus

      integer nx_tile,ny_tile
      integer dom_i,dom_j
      integer point_value
      integer ilon,ilat
      integer method

      real    dom_lats(nx_dom,ny_dom)
      real    dom_lons(nx_dom,ny_dom)
      real    grx(nx_dom,ny_dom)
      real    gry(nx_dom,ny_dom)

!Define Array to keep count of number of raw tile data points found for each
!model grid point

!     real    cat_pct(nx_dom,ny_dom,ncat)

!     integer, allocatable :: num_raw_points_total (:,:)
!     integer, allocatable :: num_raw_points_cat(:,:,:)

      real,    allocatable :: raw_data(:,:,:)
      real,    allocatable :: data_proc(:,:,:)
      real,    allocatable :: lmask_src(:,:)

      real    geodat(nx_dom,ny_dom,ncat)
      real    lmask_out(nx_dom,ny_dom)
      real    amean(ncat)
      real    asum(ncat)

      real    dlat_tile
      real    dlon_tile
      real    tile_n_lat
      real    tile_s_lat
      real    tile_w_lon
      real    tile_e_lon
      real    deltallo
      real    min_lat
      real    max_lat
      real    min_lon
      real    max_lon
      real    raw_lon
      real    raw_lat
      real    start_lon
      real    start_lat
      real    rwoff,rsoff
      real    rlondif
      real    offset
      real    ri,rj
      real    min_val,max_val,def_val,val_mask
      real    r_missing_data
 
      logical  dem_data
      logical  lgotW,lgotE
      logical  make_srcmask
      logical  lexist,lopen

      character*(*) path_to_tile_data
      character*255 title
      character*255 cfname
      character(len=8),allocatable:: ctile_name_list(:)
      character*8   ctilename
      character*1   curNS,curEW
      character*1   ctiletype

      character*1   cgrddef
      integer nxst,nyst,nz
      real    lat0,lon0,dlat,dlon
      common /llgrid/nxst,nyst,nz,lat0,lon0,dlat,dlon,cgrddef

c     interface
c       subroutine read_dem(unit_no,unit_name,nn1,nn2,i1,i2
c    &,data)
c         character*(*) unit_name
c         integer  unit_no,nn1,nn2,nn3,nn4,nofr
c         integer  i1, i2
c         real     data(nn1,nn2,nn3,nn4)
c       end subroutine
c     end interface

c from Brent Shaw pseudocode

!LSM Field Processing for "Many-to-one" aggregation


!Define Model Grid - nx,ny, projection, etc.

! nx_dom = number of East-West points in anal/model domain   I
! ny_dom = number of North-south points      "               I
! dom_lats(nx_dom,ny_dom)  ! Latitude array                  I
! dom_lons(nx_dom,ny_dom)  ! Longitude array                 I

      call s_len(path_to_tile_data,lenp)
      ctiletype=path_to_tile_data(lenp:lenp)
      if(ctiletype.eq.'V'.or.ctiletype.eq.'O'.or.
     1   ctiletype.eq.'U'.or.ctiletype.eq.'T'.or.
     1   ctiletype.eq.'A'.or.ctiletype.eq.'G')then
         print*,'tile type to process = ',ctiletype
      else
         print*,'Unknown tile type in proc_geodat_tiles'
         print*,'path_to_tile_data: ',path_to_tile_data(1:lenp)
         stop
      endif

      TITLE=path_to_tile_data(1:lenp)//'HEADER'
      lentd=INDEX(TITLE,' ')-1
      CALL JCLGET(29,TITLE(1:lentd),'FORMATTED',1,istatus)
      if(istatus .ne. 1)then
         write(6,*)'Warning: proc_geodat_tiles opening HEADER: check'
     1            ,'geog paths and HEADER file'
         return
      endif

      READ(29,*)IBLKSIZO,NO,ISBEGO,IWBEGO,RWOFF,RSOFF
      print *,'title=',title(1:lentd)
      print *,'RWOFF,RSOFF = ',RWOFF,RSOFF
      print *,'isbego,iwbego=',isbego,iwbego
      print *,'iblksizo,no=',iblksizo,no
      CLOSE(29)

      maxtiles=(360/iblksizo)+(180/iblksizo)
      print*,'Maxtiles = ',maxtiles

      if(NO .GT. 0)then
         DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO)
      else
         print *,'NO <= 0 ... error '
         stop
      endif

      itilesize_d = iblksizo

! Define number of points in a raw tile and what the lat/lon increment
! between the points in a tile are

      nx_tile=no
      ny_tile=no
      dlat_tile=deltallo
      dlon_tile=dlat_tile

! ncat ! Number of categories (1 for topo, 24 for veg, etc.)
!Define Array to keep count of number of raw tile data points found for each
!model grid point

! Find min/max latitude and longitude so we can compute which tiles to 
! read

      min_lat = MINVAL(dom_lats)
      max_lat = MAXVAL(dom_lats)
      min_lon = MINVAL(dom_lons)
      max_lon = MAXVAL(dom_lons)

! apply offsets here since we need to concern ourselves with
! which data points in the tiles map to the domain.
      min_lat =  max(-89.9999,min(89.9999,min_lat + rsoff))
      max_lat =  max(-89.9999,min(89.9999,max_lat)+ rsoff)
      min_lon =  max(-179.9999,min(179.9999,min_lon + rwoff))
      max_lon =  max(-179.9999,min(179.9999,max_lon + rwoff))

!  Compute a list of tiles needed to fulfill lat/lon range just computed

      print*,'allocate ctile_name_list lat/lon arrays'
      print*,'with maxtiles = ',maxtiles

      allocate (ctile_name_list(maxtiles))
      allocate (itile_lat(maxtiles))
      allocate (itile_lon(maxtiles))

      call get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,iblksizo,ctiletype,itilesize_d
     1,ntn,ctile_name_list,istatus)

      if(istatus.ne.1)then
         print*,'Error:  returned from get_tile_list'
         return
      endif

      ALLOCATE(raw_data(nx_tile,ny_tile,ncat))

!     ALLOCATE(num_raw_points_total(nx_dom,ny_dom))
!     ALLOCATE(num_raw_points_cat(nx_dom,ny_dom,ncat))
!     num_raw_points_total(:,:)=0
!     num_raw_points_cat(:,:,:)=0

c preliminary ... change sign of lat/lon as necessary and
c find SW anchor point.

      min_lon=360
      max_lon=-360
      min_lat=90
      max_lat=0
      lgotW=.false.
      lgotE=.false.

      do itile=1,ntn
         read(ctile_name_list(itile)(1:2),'(i2.2)')itile_lat(itile)
         read(ctile_name_list(itile)(4:6),'(i3.3)')itile_lon(itile)

         if(ctile_name_list(itile)(3:3).eq.'S')then
            itile_lat(itile)=-1.0*itile_lat(itile)
         endif
         if(itile_lat(itile).lt.min_lat)min_lat=itile_lat(itile)
         if(itile_lat(itile).gt.max_lat)max_lat=itile_lat(itile)

         if(ctile_name_list(itile)(7:7).eq.'W')then
c           lgotW=.true.
            itile_lon(itile)=-1.0*itile_lon(itile)
         endif
         if(itile_lon(itile).lt.min_lon)min_lon=itile_lon(itile)
         if(itile_lon(itile).gt.max_lon)max_lon=itile_lon(itile)

c        if(ctile_name_list(itile)(7:7).eq.'W')then
c           itile_lon(itile)=360-itile_lon(itile)
c        endif
c        if(ctile_name_list(itile)(7:7).eq.'E')lgotE=.true.
      enddo

      tile_s_lat = float(minval(itile_lat(1:ntn)))+rsoff
      tile_n_lat = float(maxval(itile_lat(1:ntn)))+rsoff
      tile_w_lon = float(minval(itile_lon(1:ntn)))+rwoff
      tile_e_lon = float(maxval(itile_lon(1:ntn)))+rwoff

      print*,'S/N tile lat points = ',tile_s_lat,tile_n_lat
      print*,'W/E tile lon points = ',tile_w_lon,tile_e_lon

c determine the x/y size dimensions and allocate/initialize the super tile

      rlondif=tile_e_lon-tile_w_lon
      if(rlondif.lt.0.0)rlondif=360.+rlondif
      itx=nx_tile*(nint(rlondif)+itilesize_d)/float(itilesize_d)
      ity=ny_tile*nint((abs(tile_n_lat-tile_s_lat)+itilesize_d)
     ./float(itilesize_d))
      print*,'allocate data_proc: nx/ny/ncat ',itx,ity,ncat
      ALLOCATE(data_proc(itx,ity,ncat))
      call get_r_missing_data(r_missing_data,istatus)
      data_proc = r_missing_data

c for the supertile lat-lon "grid" definition:
c the LL setup must consider the offsets since this
c is relevant to the actual data points within the tile.
      nxst=itx
      nyst=ity
      dlat=dlat_tile
      dlon=dlat
      lat0=min_lat+rsoff 
c     if(lgotE .and. lgotW)min_lon=360+min_lon+rwoff
      lon0=min_lon+rwoff 
      cgrddef='S'

      print*,'generating supertile for domain'
      print*,'number of small tiles needed= ',ntn

      DO itile = 1, ntn  !number of tiles needed
 
       ctilename=ctile_name_list(itile)
       cfname = path_to_tile_data(1:lenp)//ctilename
       call s_len(cfname,lenf)
       print*,'processing tile: ',cfname(1:lenf)

       inquire(file=cfname,exist=lexist,opened=lopen)
       if(.not.lexist)then
          print*,'Error: file does not exist: ',cfname(1:lenp)
          goto 1000
       endif
       if(lopen)then
          print*,'Error: file already open: ',cfname(1:lenp)
          goto 1000
       endif


! Open the tile and read the points
       if( (ctiletype.eq.'U').and.(no.eq.1200) )then
           CALL READ_DEM(29,cfname,no,no,2,2,raw_data) ! world topo_30s
           dem_data=.true.
       elseif( ctiletype.eq.'O' )then      ! soiltype top and bot layer
           CALL READ_DEM(29,cfname,no,no,1,4,raw_data)
           dem_data=.true.
       elseif( ctiletype.eq.'V' )then      ! world USGS 30s landuse
           CALL READ_DEM(29,cfname,no,no,1,4,raw_data)
           dem_data=.true.
       elseif( (ctiletype.eq.'G').or.(ctiletype.eq.'A') )then      ! greenfrac/albedo
           CALL READ_DEM_G(29,cfname,no,no,1,ncat,itile,1,4,raw_data)
           dem_data=.true.
       elseif( ctiletype.eq.'T' )then      ! soiltemp
           CALL READ_DEM(29,cfname,no,no,2,2,raw_data)
           dem_data=.true.
       else                                ! other  like albedo
           CALL JCLGET(29,cfname,'FORMATTED',0,istatus)
!          CALL VFIREC(29,raw_data,NO*NO,'LIN')
           if ((ctiletype.eq.'U').and.(no.eq.121)) then
                dem_data=.false.           ! topo_30s
           endif
       endif
c
c make "super tile" from all smaller tiles
c
       read(ctilename(1:2),'(i2.2)')icurNS
       read(ctilename(3:3),'(a1)')curNS
       read(ctilename(4:6),'(i3.3)')icurEW
       read(ctilename(7:7),'(a1)')curEW

c      print*,'compute itiley/itilex'
       if(curNS .eq. 'S') icurNS=-1.0*icurNS
       itiley=nint(1.0+(float(icurNS)-tile_s_lat)/float(itilesize_d))
       if(curEW .eq. 'W') icurEW=-1.0*icurEW
       itilex=nint(1.0+(float(icurEW)-tile_w_lon)/float(itilesize_d))

       itx=1+(itilex-1)*nx_tile
       ity=1+(itiley-1)*ny_tile

c      print*,'itilex/itiley ',itilex,itiley
c      print*,'itx/ity ',itx,ity

       jj=0
       do iy=ity,ny_tile*itiley
          jj=jj+1
          ii=0
          do ix=itx,nx_tile*itilex
             ii=ii+1
             data_proc(ix,iy,:)=raw_data(ii,jj,:)
          enddo
       enddo

      ENDDO

      deallocate (raw_data,itile_lat,itile_lon,ctile_name_list)

      print*,'initializing hinterp grx/gry arrays'

      call init_hinterp(nxst,nyst,nx_dom,ny_dom,'LL',
     .dom_lats,dom_lons,grx,gry,1,'     ')

      print*,'grid rx/ry computed'
      print*,'SW: grx/gry    1,1: ',grx(1,1),gry(1,1)
      print*,'SE: grx/gry   nx,1: ',grx(nx_dom,1),gry(nx_dom,1)
      print*,'NW: grx/gry   1,ny: ',grx(1,ny_dom),gry(1,ny_dom)
      print*,'NE: grx/gry  nx,ny: ',grx(nx_dom,ny_dom)
     .,gry(nx_dom,ny_dom)

c compute mean value to use as def_value

      is=nint(MINVAL(grx))
      ie=nint(MAXVAL(grx))
      js=nint(MINVAL(gry))
      je=nint(MAXVAL(gry))
      do k=1,ncat
         asum(k)=0.0
         icnta=0
         do j=js,je
         do i=is,ie
            if(data_proc(i,j,k).ne.0.0)then
              icnta=icnta+1
              asum(k)=asum(k)+data_proc(i,j,k)
            endif
         enddo
         enddo
         amean(k)=asum(k)/icnta
      enddo

c     print*,'is/ie,js/je ',is,ie,js,je
c     if(ctiletype.eq.'A')then
c        write(10)data_proc(is:ie,js:je,1)
c     endif

      allocate (lmask_src(nxst,nyst))
      make_srcmask=.true.
      method  = 2
      val_mask = 1

      if(ctiletype.eq.'T')then

         min_val = 239.73
         max_val = 305.09

      elseif(ctiletype.eq.'A')then

         min_val = 2.0
         max_val = 100.

      elseif(ctiletype.eq.'G')then

         min_val = 1.0
         max_val = 100.

      endif

      do ii=1,ncat

         def_val = amean(ii)

         if(.false.)then
            call hinterp_field(nxst,nyst,nx_dom,ny_dom,1
     .,grx,gry,data_proc(1,1,ii),geodat(1,1,ii),1)
            geodat(:,:,ii)=geodat(:,:,ii)/100.

         else

            call interpolate_masked_val(nxst, nyst
     ., lmask_src, data_proc(1,1,ii), nx_dom, ny_dom, lmask_out
     ., geodat(1,1,ii), grx, gry, make_srcmask, min_val, max_val
     ., def_val, val_mask, method)
            if(ctiletype.eq.'A')then
               geodat(:,:,ii)=geodat(:,:,ii)/100.
            endif

         endif

      enddo

      deallocate (data_proc,lmask_src)

! Loop over the raw data points in the tile

!       DO jt = 1, ny_tile
!       DO it = 1, nx_tile

! Compute lat/lon of this raw data point using lat/lon from
! file name, the offset of the data, and the increment

!        raw_lat = float(ilat) + RSOFF + (jt-1)*dlat_tile
!        raw_lon = float(ilon) + RWOFF + (it-1)*dlon_tile

! Use lat/lon to ij to find wheter this tile point is in our domain

!        call latlon_to_rlapsgrid(raw_lat,raw_lon,dom_lats
!    &                ,dom_lons,nx_dom,ny_dom,ri,rj,istatus)
!        dom_i = NINT(ri)
!        dom_j = NINT(rj)
!        IF( (dom_i .GE. 1).AND.(dom_i .LE. nx_dom).AND.
!    &   (dom_j .GE. 1).AND.(dom_j .LE. ny_dom) ) THEN

!         if(ctiletype.eq.'V')then

!            point_value=int(raw_data(it,jt,1))
!            num_raw_points_cat(dom_i,dom_j,point_value) =
!    &       num_raw_points_cat(dom_i,dom_j,point_value) + 1

!            num_raw_points_total(dom_i, dom_j) =
!    &       num_raw_points_total(dom_i, dom_j)+ 1

!         elseif(ctiletype.eq.'T')then

!            num_raw_points_cat(dom_i,dom_j,1) =
!    &       num_raw_points_cat(dom_i,dom_j,1)+1

!            num_raw_points_total(dom_i, dom_j) =
!    &       num_raw_points_total(dom_i, dom_j)+raw_data(it,jt,1)

!         endif
         
!         ! The point is in our domain, so increment the counters.  The 
          ! value of the raw data point is equal to the category ID for
          ! things like veg-type and soil-type.  For topography, we would
          ! have to have additional arrays to keep the sum of the terrain
          ! values to compute the average value at the end. This example
          ! code only does category data

!        ENDIF
!       ENDDO
!       ENDDO 


! Now you can compute percentage of each category like so:
!     if(ctiletype.eq.'V')then
!        DO icat = 1, ncat
!           cat_pct(:,:,icat)=num_raw_points_cat(:,:,icat) / 
!    &                         num_raw_points_total(:,:)
!        ENDDO
!     elseif(ctiletype.eq.'T')then
!           cat_pct(:,:,1)= num_raw_points_total(:,:)/
!    &                       num_raw_points_cat(:,:,1)
!     endif

      return

1000  print*,'returning to main. no data processed'
      return
      end
