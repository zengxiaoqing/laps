
      subroutine compute_latlon(nnxp,nnyp,n_staggers
     + ,deltax,xtn_ret,ytn_ret,lats,lons,istatus)


      implicit none

      integer nnxp,nnyp
      integer n_staggers
      integer i,j,k

      real    deltax,deltay
      real    mdlat,mdlon
      real    stagger_ew,stagger_ns
      real    erad

c these used internally
      real    xtn(nnxp)
      real    ytn(nnyp)

c these arrays returned for geodat subroutine
      real    xtn_ret(nnxp)
      real    ytn_ret(nnyp)

c these are the anchor points (SW corner) of the grid
      real    xmn1,ymn1

c     real    lat(nnxp,nnyp)
c     real    lon(nnxp,nnyp)

c A-c staggers contained within these arrays.
      real    lats(nnxp,nnyp,n_staggers)
      real    lons(nnxp,nnyp,n_staggers)

      integer istatus

      print*,'calculate lat/lon at stagger grid points.'

      deltay=deltax

      call get_grid_center(mdlat,mdlon,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error returned: get_grid_center'
         return
      endif

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error calling get_earth_radius'
         return
      endif

      do k=1,n_staggers

         if(k.eq.1)then

c           Get X/Y for lower left corner
            CALL POLAR_GP(mdlat,mdlon,XMN1,YMN1,DELTAX,DELTAY,
     +  NNXP,NNYP)
            stagger_ns=0
            stagger_ew=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)
c
c Return to main the non-staggered x/y.
            do i=1,nnxp
               xtn_ret(i)=xtn(i)
            enddo
            do j=1,nnyp
               ytn_ret(j)=ytn(j)
            enddo

         elseif(k.eq.2)then    !this is the B stagger (E-W stagger)

            stagger_ew=0.5*deltax
            stagger_ns=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

         elseif(k.eq.3)then !this is the C stagger (N-S stagger)

            call get_grid_center(mdlat,mdlon,istatus)
            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn,ytn)
            stagger_ew=0
            stagger_ns=0.5*deltay 
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn) 

         endif
C
C*****************************************************************
C*  Convert it to lat/lon using the library routines.            *

         Do J = 1,nnyp
         Do I = 1,nnxp

            call xy_to_latlon(xtn(i),ytn(j),erad ! ,90.,std_lon
     1                                     ,lats(I,J,k),lons(I,J,k))

c             print *,'i,j,xtn,ytn,pla,lplo=',i,j,xtn,ytn,pla,plo
         enddo
         enddo

      enddo

      istatus = 1

      return
      end

c ----------------------------------------------

      subroutine get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

      implicit none

      integer  nnxp,nnyp
      integer  i,j

      real    deltax,deltay
      real    stagger_ew,stagger_ns
      real    xmn1,ymn1
      real    xmn(nnxp)
      real    ymn(nnyp)
      real    xtn(nnxp)
      real    ytn(nnyp)

      xmn(1)=xmn1+stagger_ew
      ymn(1)=ymn1+stagger_ns

      DO 600 I=2,NNXP
         XMN(I)=XMN(I-1)+DELTAX
 600  CONTINUE
      XMN(NNXP)=2*XMN(NNXP-1)-XMN(NNXP-2)
 
      DO 610 J=2,NNYP
         YMN(J)=YMN(J-1)+DELTAY
 610  CONTINUE
      YMN(NNYP)=2*YMN(NNYP-1)-YMN(NNYP-2)
 
      DO 650 I=2,NNXP
         XTN(I)=.5*(XMN(I)+XMN(I-1))
 650  CONTINUE
      XTN(1)=1.5*XMN(1)-.5*XMN(2)

      DO 660 J=2,NNYP
         YTN(J)=.5*(YMN(J)+YMN(J-1))
 660  CONTINUE
      YTN(1)=1.5*YMN(1)-.5*YMN(2)

      return
      end
c
c -------------------------------------------------------------
c
      subroutine get_map_factor_grid(nx,ny,n_staggers
     +,rlats,rlons ,rmap_factors,istatus)

      implicit none

      integer nx,ny
      integer n_staggers
      integer i,j,k
      integer istatus

      real    rlats(nx,ny,n_staggers)
      real    rlons(nx,ny,n_staggers)
      real    rmap_factors(nx,ny,n_staggers)
      real    sigma

      do k=1,n_staggers
      do j=1,ny
      do i=1,nx
         call get_sigma(rlats(i,j,k),rlons(i,j,k)
     +,sigma,istatus)
         if(istatus.ne.1)goto 999
         rmap_factors(i,j,k)=sigma
      enddo
      enddo
      enddo

      return

999   print*,'Error returned: get_sigma'

      return
      end
c
c ------------------------------------------------------------
c
      subroutine get_coriolis_components(nx,ny,lat,coriolis_parms)
c

      implicit none

      include 'trigd.inc'

      integer nx,ny
      integer i,j

      real    lat(nx,ny)
      real    coriolis_parms(nx,ny,2)

      real    omega_ear
      data    omega_ear/7.292e-5/

      do j=1,ny
      do i=1,nx
         coriolis_parms(i,j,1)=2*omega_ear*sind(lat(i,j))
         coriolis_parms(i,j,2)=2*omega_ear*cosd(lat(i,j))
      enddo
      enddo

      return
      end

      subroutine get_projrot_grid(nx,ny,lat,lon,projrot_grid
     +,istatus)

      implicit none

      integer nx,ny
      integer istatus
      integer i,j
      real    r_missing_data
      real    projrot_deg
      real    projrot_latlon
      real    projrot_grid(nx,ny,2)
      real    lat(nx,ny),lon(nx,ny)

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: get_r_missing_data'
         return
      endif

      do j=1,ny
      do i=1,nx
         projrot_grid(i,j,1)=r_missing_data
         projrot_grid(i,j,2)=r_missing_data
         projrot_deg=projrot_latlon(lat(i,j),lon(i,j),istatus)
         if(istatus.eq.1)then
            projrot_grid(i,j,1)=sind(projrot_deg)
            projrot_grid(i,j,2)=cosd(projrot_deg)
         endif
      enddo
      enddo

      return
      end
c
c--------------------------------------------------------
c
      subroutine get_gridgen_var(nf,ngrids,var,comment)

      implicit none

      integer        nf,ngrids
      character*(*)  var(nf)
      character*(*)  comment(nf)

      if(ngrids.eq.6)then

         var(1)    = 'LAT'
         var(2)    = 'LON'
         var(3)    = 'AVG'
         var(4)    = 'LDF'
         var(5)    = 'USE'
         var(6)    = 'ZIN'

         comment(1) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(2) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(3) = 'Average terrain elevation (m) \0'
         comment(4) = '\0'
         comment(5) = '\0'

      elseif(ngrids.eq.17)then

         var(1)    = 'LAT'  ! non-staggered (A-grid) lats
         var(2)    = 'LON'  ! non-staggered (A-grid) lons
         var(3)    = 'LAB'  ! b-stagger (.5*deltax (e-w)) lats
         var(4)    = 'LOB'  ! b-stagger (.5*deltax (e-w)) lons
         var(5)    = 'LAC'  ! c-stagger (.5*deltay (n-s)) lats
         var(6)    = 'LOC'  ! c-stagger (.5*deltay (n-s)) lons
         var(7)    = 'AVG'  ! Topo (m) on A-grid
         var(8)    = 'LDF'  ! Land Fraction
         var(9)    = 'USE'  !
         var(10)   = 'SPR'  ! Sin(projection rotation) from true
         var(11)   = 'CPR'  ! Cos(projection rotation) from true
         var(12)   = 'MFA'  ! Map factor A grid
         var(13)   = 'MFB'  ! Map factor b-stagger grid
         var(14)   = 'MFC'  ! Map factor c-stagger grid
         var(15)   = 'CPH'  ! Horizontal component of coriolis parameter
         var(16)   = 'CPV'  ! Vertical component of coriolis parameter
         var(17)   = 'ZIN'

         comment(1) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(2) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(3) = 'B-stagger grid latitudes \0'
         comment(4) = 'B-stagger grid longitudes for WRF_SI \0'
         comment(5) = 'C-stagger grid latitudes for WRF_SI \0'
         comment(6) = 'C-stagger grid longitudes for WRF_SI \0'
         comment(7) = 'Average terrain elevation (m) \0'
         comment(8) = 'Land Fraction A-grid \0'
         comment(9) = 'Land Use Categories \0'
         comment(10)= 'Sin of projection rotation (rad) \0'
         comment(11)= 'Cosine of projection rotation (rad) \0'
         comment(12)= 'Map Factor A-grid \0'
         comment(13)= 'Map Factor b-stagger grid \0'
         comment(14)= 'Map Factor c-stagger grid \0'
         comment(15)= 'Horizontal component coriolis parameter \0'
         comment(16)= 'Vertical component coriolis parameter \0'
         comment(17)= '\0'

      endif
      return
      end
