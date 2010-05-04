      subroutine proc_geodat(nx_dom,ny_dom,ncat
     1,path_to_tile_data,dom_lats_in,dom_lons_in,lmask_out
     1,geodat,istatus)
!    1,cat_pct)
!
! J. Smart (NOAA/FSL) : 2002            original version
! T. Hume/T. Simmers Meteorological Service of New Zealand
!                     : 2002
!                    Corrected problems with crossing the dateline
!                    and constrained albedo to one tile.
! J. Smart (NOAA/FSL) : 2003
!                    Further refined processing of tiles to accomodate
!                    dateline, Greenwich mean, and Poles. Added smoothing
!                    of temp field (subroutine one_two_one).

      use horiz_interp

      implicit none 
      integer maxtiles 
      integer nx_dom 
      integer ny_dom
      integer ncat
      integer ntn,nt
      integer bgmodel
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
      integer tile_n_lat
      integer tile_s_lat
      integer tile_w_lon
      integer tile_e_lon
      integer istatus
      integer istat

      integer nx_tile,ny_tile
      integer dom_i,dom_j
      integer point_value
      integer ilon,ilat
      integer method

      real    dom_lats_in(nx_dom,ny_dom)
      real    dom_lons_in(nx_dom,ny_dom)
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
      real,    allocatable,   dimension(:,:)::dom_lats,dom_lons

      real    geodat(nx_dom,ny_dom,ncat)
      real    lmask_out(nx_dom,ny_dom)
      real    amean(ncat)
      real    asum(ncat)

      real    dlat_tile
      real    dlon_tile
      real    sw(2),ne(2)
      real    deltallo
      real    min_lat
      real    max_lat
      real    min_lon
      real    max_lon
      real    minlat
      real    maxlat
      real    minlon
      real    maxlon
      real    rwoff,rsoff
      real    rlondif
      real    offset
      real    ri,rj
      real    rmult
      real    min_val,max_val,def_val,val_mask
      real    r_missing_data
 
      logical  dem_data
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

c original create from Brent Shaw pseudocode

      print*,'Start proc_geodat'

      if(.not. allocated(dom_lats))then
          allocate(dom_lats(nx_dom,ny_dom),dom_lons(nx_dom,ny_dom))
      endif
      dom_lats=dom_lats_in
      dom_lons=dom_lons_in

! TH: 8 Aug 2002 Now we may need to adjust the longitude values in 
! dom_lons so that they are always monotonically increasing, even if 
! we have the bad fortune to cross the date line.
      DO jj=1,ny_dom-1
         DO ii=1,nx_dom
            IF ((dom_lons(ii,jj+1) - dom_lons(ii,jj)) < -180.) THEN
               dom_lons(ii,jj+1) = dom_lons(ii,jj+1) + 360.
            ELSE IF ((dom_lons(ii,jj+1) - dom_lons(ii,jj)) > 180.) THEN
               dom_lons(ii,jj+1) = dom_lons(ii,jj+1) - 360.
            END IF
         END DO
      END DO
! TH: 8 Aug 2002 We also have to cope with spastic grids that have
! "horizontal" date lines (where we might cross the date line by
! moving up and down the grid rather than left and right).
      DO jj=1,ny_dom
         DO ii=1,nx_dom-1
            IF ((dom_lons(ii+1,jj) - dom_lons(ii,jj)) < -180.) THEN
               dom_lons(ii+1,jj) = dom_lons(ii+1,jj) + 360.
            ELSE IF ((dom_lons(ii+1,jj) - dom_lons(ii,jj)) > 180.) THEN
               dom_lons(ii+1,jj) = dom_lons(ii+1,jj) - 360.
            END IF
         END DO
      END DO
! TH: 8 Aug 2002 Now we want to make sure our longitudes don't fall
! outside the range (-360,360).
      IF (MAXVAL(dom_lons(:,:)) > 360.) dom_lons(:,:) = dom_lons(:,:)
     1   - 360
      max_lat = MAXVAL(dom_lats)
      IF (MINVAL(dom_lons(:,:)) < -360.) dom_lons(:,:) = dom_lons(:,:)
     1   + 360

! nx_dom = number of East-West points in anal/model domain   I
! ny_dom = number of North-south points      "               I
! dom_lats(nx_dom,ny_dom)  ! Latitude array                  I
! dom_lons(nx_dom,ny_dom)  ! Longitude array                 I

      call s_len(path_to_tile_data,lenp)
      ctiletype=path_to_tile_data(lenp:lenp)
      if(ctiletype.eq.'T'.or.
     1   ctiletype.eq.'A'.or.
     1   ctiletype.eq.'G'.or.
     1   ctiletype.eq.'U'.or.
     1   ctiletype.eq.'M' )then
         print*,'tile type to process = ',ctiletype
      else
         print*,'Unknown tile type in proc_geodat_tiles'
         print*,'path_to_tile_data: ',path_to_tile_data(1:lenp)
         stop
      endif

      bgmodel=6

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

      if(ctiletype.eq.'T'.and.iblksizo.eq.10)then
       print*
       print*,'Error: ****************************************'
       print*,'Error: ********** TERMINATING *****************'
       print*,'Error: Old Soil Temp database!!!'
       print*,'Error: You need to update these'
       print*,'Error: data: ',TRIM(path_to_tile_data)
       print*,'Error: by downloading new soil temp data'
       print*,'Error: ********** TERMINATING *****************'
       print*,'Error: ****************************************'
       print*
      endif

      maxtiles=(360/iblksizo)*(180/iblksizo)
      print*,'Maxtiles = ',maxtiles

      if(NO .GT. 0)then
         DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO)
      else
         print *,'ERROR:  NO <= 0 ... '
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

! Find min/max atitude and longitude so we can compute which tiles to 
! read

      minlat = MINVAL(dom_lats)
      maxlat = MAXVAL(dom_lats)
      minlon = MINVAL(dom_lons)
      maxlon = MAXVAL(dom_lons)

! no offsets needed here since these are the tiles, not the points
! in the tiles. However, since domains may need data from tiles with
! offset considered, add it in.
      min_lat = max(-89.9999, min(89.9999, minlat - abs(rsoff)))
      max_lat = max(-89.9999, min(89.9999, maxlat + abs(rsoff)))
      min_lon = max(-359.9999,min(359.9999,minlon - abs(rwoff)))
      max_lon = max(-359.9999,min(359.9999,maxlon + abs(rwoff)))

c     if(ctiletype .eq. 'G'.or.ctiletype.eq.'A')then
c        if(min_lon.lt.-180)min_lon=360.+min_lon
c     endif

      print*,'max/min lats: ',max_lat,min_lat
      print*,'max/min lons: ',max_lon,min_lon

      deallocate(dom_lats,dom_lons)

!  Compute a list of tiles needed to fulfill lat/lon range just computed

      print*,'allocate ctile_name_list lat/lon arrays'
      print*,'with maxtiles = ',maxtiles

      allocate (ctile_name_list(maxtiles))

      call get_tile_list(min_lat,max_lat,min_lon,max_lon
     1,maxtiles,isbego,iwbego,itilesize_d,ctiletype
     1,ntn,ctile_name_list,tile_w_lon,tile_e_lon
     1,tile_s_lat,tile_n_lat,istatus)

      if(istatus.ne.1)then
         print*,'ERROR:  returned from get_tile_list'
         return
      endif

      allocate (raw_data(nx_tile,ny_tile,ncat))
      allocate (itile_lat(ntn))
      allocate (itile_lon(ntn))

!     ALLOCATE(num_raw_points_total(nx_dom,ny_dom))
!     ALLOCATE(num_raw_points_cat(nx_dom,ny_dom,ncat))
!     num_raw_points_total(:,:)=0
!     num_raw_points_cat(:,:,:)=0

      do itile=1,ntn
         read(ctile_name_list(itile)(1:2),'(i2.2)')itile_lat(itile)
         read(ctile_name_list(itile)(4:6),'(i3.3)')itile_lon(itile)
         if(ctile_name_list(itile)(7:7).eq.'W')then
            if(itile_lon(itile).lt.180)then
               itile_lon(itile)=360-itile_lon(itile)
            endif
         endif
         if(ctile_name_list(itile)(3:3).eq.'S')then
            itile_lat(itile)=-1.0*itile_lat(itile)
         endif
      enddo

      print*,'S/N tile lat points = ',tile_s_lat,tile_n_lat
      print*,'W/E tile lon points = ',tile_w_lon,tile_e_lon

c determine the x/y size dimensions and allocate/initialize the super tile

      rlondif=abs(tile_e_lon-tile_w_lon)
      itx=nx_tile*nint((rlondif+itilesize_d)/float(itilesize_d))
      ity=ny_tile*nint((abs(tile_n_lat-tile_s_lat)+itilesize_d)
     ./float(itilesize_d))
      if(ctiletype.eq.'T' .or. ctiletype.eq.'M' .and.
     &itx.gt.360)itx=360
      if(ctiletype.eq.'G' .or. ctiletype.eq.'A' .and. 
     &itx.gt.2500)itx=2500
      print*,'allocate data_proc: nx/ny/ncat ',itx,ity,ncat
      allocate(data_proc(itx,ity,ncat))
      call get_r_missing_data(r_missing_data,istatus)
      data_proc = r_missing_data

c for the supertile lat-lon "grid" definition:
c the LL setup must consider the offsets since this
c is relevant to the actual data points within the tile.
      nxst=itx
      nyst=ity
      dlat=dlat_tile
      dlon=dlat
      lat0=tile_s_lat+rsoff
      if(tile_w_lon.gt.180)then
         lon0=tile_w_lon-360+rwoff 
      elseif(tile_w_lon.lt.-180)then
         lon0=360+tile_w_lon+rwoff
      else
         lon0=tile_w_lon+rwoff
      endif
      sw(1)=lat0
      sw(2)=lon0
      ne(1)=tile_n_lat+itilesize_d+rsoff
      ne(2)=tile_e_lon+itilesize_d+rwoff
      cgrddef='S'

      print*,'generate supertile for domain'
      print*,'number of small tiles needed= ',ntn
      print*

      DO itile = 1, ntn  !number of tiles needed
 
       ctilename=ctile_name_list(itile)
       cfname = path_to_tile_data(1:lenp)//ctilename(1:3)
       read(ctilename(4:6),'(i3.3)')icurEW
       read(ctilename(7:7),'(a1)')curEW
       IF (icurEW > 180)THEN    ! .and. icurEW < 360) THEN
          icurEW = 360 - icurEW
       ENDIF
       IF (icurEW == 180.and.curEW.ne.'W') THEN
          curEW = 'W' 
       END IF 

       write(cfname(lenp+4:lenp+6),'(i3.3)')icurEW
       write(cfname(lenp+7:lenp+7),'(a1)')curEW

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
           CALL READ_DEM(29,cfname,no,no,2,2,raw_data,istat) ! world topo_30s
           dem_data=.true.
       elseif( ctiletype.eq.'O' )then      ! soiltype top and bot layer
           CALL READ_DEM(29,cfname,no,no,1,4,raw_data,istat)
           dem_data=.true.
       elseif( ctiletype.eq.'V' )then      ! world USGS 30s landuse
           CALL READ_DEM(29,cfname,no,no,1,4,raw_data,istat)
           dem_data=.true.
       elseif( (ctiletype.eq.'G')
     1     .or.(ctiletype.eq.'A')
     1     .or.(ctiletype.eq.'M') )then      ! greenfrac/albedo/maxsnowalb
           CALL READ_DEM_G(29,cfname,no,no,1,ncat,1,1,4,raw_data,istat)
           dem_data=.true.
       elseif( ctiletype.eq.'T')then ! .or. ctiletype.eq.'M')then      ! soiltemp
           CALL READ_DEM(29,cfname,no,no,2,2,raw_data,istat)
           dem_data=.true.
       else                                ! other  like albedo
           CALL JCLGET(29,cfname,'FORMATTED',0,istatus)
           CALL VFIREC(29,raw_data,no*no,'LIN')
           if ((ctiletype.eq.'U').and.(no.eq.121)) then
                dem_data=.false.           ! topo_30s
           endif
       endif

       if(istat.ne.0)then
          print*,'Error returned: proc_geodat: READ_DEM'
          istatus=0
          return
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
c
       rlondif=float(itile_lon(itile))-tile_w_lon
       if(rlondif.lt.0)then
          rlondif=rlondif+360.
       elseif(rlondif.ge.360)then
          rlondif=rlondif-360
       endif

       itilex=nint(1.0+(rlondif/float(itilesize_d)))

       itx=1+(itilex-1)*nx_tile
       ity=1+(itiley-1)*ny_tile

       print*,'itilex/itiley ',itilex,itiley
       print*,'itx/ity ',itx,ity

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

      print*,'nxst= ',nxst
      print*,'nyst= ',nyst
      print*,'dlat= ',dlat
      print*,'dlon= ',dlon
      print*,'lat0= ',lat0
      print*,'lon0= ',lon0

      call init_hinterp(nxst,nyst,nx_dom,ny_dom,'LL',
     .dom_lats_in,dom_lons_in,grx,gry,bgmodel,'     ',.true.)

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

      print*,'is/ie/js/je ',is,ie,js,je

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
         if(icnta.gt.0)then
            amean(k)=asum(k)/icnta
         else
            amean(k)=r_missing_data
         endif
      enddo

      allocate (lmask_src(nxst,nyst))
      make_srcmask=.true.
      val_mask = 1
      nt=1

      if(ctiletype.eq.'T')then

         min_val = 239.73
         max_val = 305.09
         method  = 1
         geodat(:,:,:)=r_missing_data

      elseif(ctiletype.eq.'A')then

         min_val = 2.0
         max_val = 100.
         method  = 2
         where(amean .eq. r_missing_data)amean=18.

      elseif(ctiletype.eq.'G')then

         min_val = 1.0
         max_val = 100.
         method  = 2
         where(amean .eq. r_missing_data)amean=65.

      elseif(ctiletype.eq.'M')then

         min_val = 1.0
         max_val = 100.
         method  = 1
         geodat(:,:,:)=65.

      endif

      do ii=1,ncat

         def_val = amean(ii)

         call interpolate_masked_val(nxst, nyst
     ., lmask_src, data_proc(1,1,ii), nx_dom, ny_dom, lmask_out
     ., geodat(1,1,ii), grx, gry, make_srcmask, min_val, max_val
     ., def_val, val_mask, method)

         if(ctiletype.eq.'A'.or.ctiletype.eq.'M')then
            geodat(:,:,ii)=geodat(:,:,ii)/100.
         endif

         if(ctiletype.eq.'T')then  !.or.ctiletype.eq.'M')then
            print*,'Filering Deep Soil Temp with 1-2-1'
            call one_two_one(nx_dom,ny_dom,nt,geodat(1,1,ii))
         endif

      enddo

      deallocate (data_proc,lmask_src)

! --------------------------------------------------------
! ---- we are not yet processing categorical data (like
! ---- landuse or soil type, so this section is commented
! ---- for now.
! --------------------------------------------------------
! Loop over the raw data points in the tile

!       DO jt = 1, ny_tile
!       DO it = 1, nx_tile

! Compute lat/lon of this raw data point using lat/lon from
! file name, the offset of the data, and the increment

!        raw_lat = float(ilat) + RSOFF + (jt-1)*dlat_tile
!        raw_lon = float(ilon) + RWOFF + (it-1)*dlon_tile

! Use lat/lon to ij to find wheter this tile point is in our domain

!        call latlon_to_rlapsgrid(raw_lat,raw_lon,dom_lats_in
!    &                ,dom_lons_in,nx_dom,ny_dom,ri,rj,istatus)
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
c
c-----------------------------------------------------------
c
      subroutine one_two_one(nx,ny,nt,data)

      implicit none
      integer i,j,ii,jj,icnt,n
      integer nx,ny,nt
      integer istatus
      real  r_missing_data
      real  fact,sum
      real  data(nx,ny)
      real, allocatable :: temp(:,:)

      allocate(temp(nx,ny))

      call get_r_missing_data(r_missing_data,istatus)

      temp=data

      do n=1,nt

         do j=2,ny-1
         do i=2,nx-1

            sum=0.0
            icnt=0
            do jj=j-1,j+1
            do ii=i-1,i+1

               fact=1
               if(jj==j.and.ii==i)fact=2
               if(data(ii,jj).lt.r_missing_data)then
                  icnt=icnt+fact
                  sum=sum+data(ii,jj)*fact
               endif

            enddo
            enddo

            if(icnt.gt.0)then
               temp(i,j)=sum/float(icnt)
            else
               temp(i,j)=data(i,j)
            endif

         enddo
         enddo

         data=temp

      enddo

      deallocate (temp)

      return
      end
