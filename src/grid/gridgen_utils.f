
      subroutine compute_latlon(nnxp,nnyp,n_staggers
     + ,deltax,xtn_ret,ytn_ret,lats,lons,istatus)


      implicit none

      integer nnxp,nnyp
      integer n_staggers
      integer i,j,k,nc

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

      character c_dataroot*255
      character c10_grid_fname*10

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

c     call find_domain_name(c_dataroot,c10_grid_fname
c    1,istatus)
c     if(istatus .ne. 1)then
c        write(6,*)' Error calling find_domain_name'
c        return
c     endif
c     call s_len(c10_grid_fname,nc)


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
c Return to main the non-staggered x/y for terrain analysis.
            
c           if(c10_grid_fname(1:nc).eq.'nest7grid')then

               do i=1,nnxp
                  xtn_ret(i)=xtn(i)
               enddo
               do j=1,nnyp
                  ytn_ret(j)=ytn(j)
               enddo
c           endif

         elseif(k.eq.2)then    !this is A-stagger (.5 E-W stagger) 

            stagger_ew=0.5*deltax
            stagger_ns=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

         elseif(k.eq.3)then !this is the B-stagger (.5 N-S stagger)

            call get_grid_center(mdlat,mdlon,istatus)
            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn,ytn)
            stagger_ew=0
            stagger_ns=0.5*deltay 
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn) 

         elseif(k.eq.4)then !this is the C-stagger (.5 both N-S and E-W)

            call get_grid_center(mdlat,mdlon,istatus)
            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn,ytn)
            stagger_ew=0.5*deltax
            stagger_ns=0.5*deltay
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)


c           if(c10_grid_fname(1:nc).eq.'wrfsi')then
c              do i=1,nnxp
c                 xtn_ret(i)=xtn(i)
c              enddo
c              do j=1,nnyp
c                 ytn_ret(j)=ytn(j)
c              enddo
c           endif

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
      include 'trigd.inc'
      implicit none

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
c
c ----------------------------------------------------------
c
      subroutine get_projrot_grid(nx,ny,lat,lon,projrot_grid
     +,istatus)

      include 'trigd.inc'
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
      subroutine get_static_albedo(nx,ny,lat,lon,landfrac
     +,static_albedo,istatus)

      implicit none

      integer nx,ny
      real    lat(nx,ny),lon(nx,ny)
      real    landfrac(nx,ny)
      real    static_albedo(nx,ny)
      real    water_albedo_cmn
      real    r_missing_data
      integer istatus
      integer i,j
      integer nwater

      data water_albedo_cmn/0.04/

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: get_r_missing_data'
         return
      endif

      nwater = 0
      do j=1,ny
      do i=1,nx
         if(landfrac(i,j).le.0.01)then
            nwater = nwater+1
            static_albedo(i,j)=water_albedo_cmn
         else
            static_albedo(i,j)=r_missing_data
         endif
      enddo
      enddo
      if(nwater .gt. 0)then
         print*,'static_albedo ',water_albedo_cmn,' used ',
     +'at ',nwater,' grid points '
      else
         print*,'No water grid points for water_albedo',
     +' in this domain'
      endif

      return
      end
c
c --------------------------------------------------------
c
      subroutine bilinear_interp(i,j,imax,jmax,array_2d,result)
c
c this used only for getting topography on the c-staggered grid

      implicit none
      integer i,j,imax,jmax
      real*4 result
      real*4 array_2d(imax,jmax)
      real*4 Z1,Z2,Z3,Z4
      real*4 fraci,fracj

      fraci = 0.5
      fracj = 0.5

      Z1=array_2d(i  , j  )
      Z2=array_2d(i-1, j  )
      Z3=array_2d(i-1, j-1)
      Z4=array_2d(i  , j-1)

      result= Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1     - (Z2+Z4-Z3-Z1)*fraci*fracj

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

      if(ngrids.eq.10)then

         var(1)    = 'LAT'
         var(2)    = 'LON'
         var(3)    = 'AVG'
         var(4)    = 'LDF'
         var(5)    = 'USE'
         var(6)    = 'ALB'
         var(7)    = 'STD'
         var(8)    = 'SLN'
         var(9)    = 'SLT'
         var(10)   = 'ZIN'

         comment(1) = 'Lat: From MODEL by J. Snook/ S. Albers 1-95\0'
         comment(2) = 'Lon: From MODEL by J. Snook/ S. Albers 1-95\0'
         comment(3) = 'Average terrain elevation (m) \0'
         comment(4) = 'Land Fraction \0'
         comment(5) = 'Soil Type (Missing) \0'
         comment(6) = 'Clear Sky Albedo - fixed at .04 over water\0'
         comment(7) = 'Standard Deviation of Elevation data (m)\0'
         comment(8) = 'Mean longitudinal terrain slope (m/m)\0'
         comment(9) = 'Mean latitudinal terrain slope (m/m)\0'
         comment(10)='\0'


      elseif(ngrids.eq.25)then

         var(1)    = 'LAT'  ! non-staggered (Analysis-grid) lats
         var(2)    = 'LON'  ! non-staggered (Analysis-grid) lons
         var(3)    = 'LAA'  ! a-stagger (.5*deltax (e-w)) lats
         var(4)    = 'LOA'  ! a-stagger (.5*deltax (e-w)) lons
         var(5)    = 'LAB'  ! b-stagger (.5*deltay (n-s)) lats
         var(6)    = 'LOB'  ! b-stagger (.5*deltay (n-s)) lons
         var(7)    = 'LAC'  ! c-stagger (.5*deltax and .5*deltay) lats
         var(8)    = 'LOC'  ! c-stagger (.5*deltay and .5*deltay) lons

         var(9)    = 'AVG'  ! Topo (m) on Analysis-grid

         var(10)   = 'LDF'  ! Land Fraction
         var(11)   = 'USE'  !
         var(12)   = 'SPR'  ! Sin(projection rotation) from true
         var(13)   = 'CPR'  ! Cos(projection rotation) from true
         var(14)   = 'MFL'  ! Map factor Analysis grid
         var(15)   = 'MFA'  ! Map factor a-stagger grid
         var(16)   = 'MFB'  ! Map factor b-stagger grid
         var(17)   = 'MFC'  ! Map factor c-stagger grid
         var(18)   = 'CPH'  ! Horizontal component of coriolis parameter
         var(19)   = 'CPV'  ! Vertical component of coriolis parameter
         var(20)   = 'ALB'  ! Static (climatological) albedo
         var(21)   = 'STD'  ! Standard Deviation of Elevation Data (m)
         var(22)   = 'SLN'  ! Terrain Slope; Longitudinal Component (m/m)
         var(23)   = 'SLT'  ! Terrain Slope; Latitudinal Component (m/m)
         var(24)   = 'AVC'  ! Topo (m) on c-stagger grid
         var(25)   = 'ZIN'

         comment(1) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(2) = 'Made from MODEL by J. Snook/ S. Albers 1-95\0'
         comment(3) = 'a-stagger grid latitudes \0'
         comment(4) = 'a-stagger grid longitudes for WRF_SI \0'
         comment(5) = 'b-stagger grid latitudes for WRF_SI \0'
         comment(6) = 'b-stagger grid longitudes for WRF_SI \0'
         comment(7) = 'c-stagger grid latitudes for WRF_SI \0'
         comment(8) = 'c-stagger grid longitudes for WRF_SI \0'

         comment(9) = 'Average terrain elevation (m) \0'

         comment(10)= 'Land Fraction A-grid \0'
         comment(11)= 'Land Use Categories \0'
         comment(12)= 'Sin of projection rotation (rad) \0'
         comment(13)= 'Cosine of projection rotation (rad) \0'
         comment(14)= 'Map Factor Analysis grid \0'
         comment(15)= 'Map Factor a-stagger grid \0'
         comment(16)= 'Map Factor b-stagger grid \0'
         comment(17)= 'Map Factor c-stagger grid \0'
         comment(18)= 'Horizontal component coriolis parameter \0'
         comment(19)= 'Vertical component coriolis parameter \0'
         comment(20)= 'Static Albedo (%) valid only over water atm \0'
         comment(21) = 'Standard Deviation of Elevation data (m)\0'
         comment(22) = 'Mean longitudinal terrain slope (m/m)\0'
         comment(23) = 'Mean latitudinal terrain slope (m/m)\0'
         comment(24) = 'Average terrain elevation (c-stagger) (m)\0'
         comment(25)= '\0'

      endif
      return
      end
