      subroutine uvgrid_to_uvtrue_a(u,v,lon,std_lon,nx,ny,nz)
c
c *** Convert grid north winds to true north winds.
c
      
c
      integer nx,ny,nz,i,j
c
      real*4    u(nx,ny,nz),
     .          v(nx,ny,nz),
     .          lon(nx,ny),
     .          std_lon,
     .          angle(nx,ny)
c_______________________________________________________________________________
c
      do j=1,ny
      do i=1,nx
         angle(i,j)=lon(i,j)-std_lon
      enddo
      enddo
c
      call rotate_vec_a(u,v,angle,nx,ny,nz)
c
      return
      end
c
c===============================================================================
c
      subroutine uvtrue_to_uvgrid_a(u,v,lon,std_lon,nx,ny,nz)
c
c *** Convert true north winds to grid north winds.
c
c
      
c
      integer nx,ny,nz,i,j
c
      real*4    u(nx,ny,nz),
     .          v(nx,ny,nz),
     .          lon(nx,ny),
     .          std_lon,
     .          angle(nx,ny)
c_______________________________________________________________________________
c
      do j=1,ny
      do i=1,nx
         angle(i,j)=std_lon-lon(i,j)
      enddo
      enddo
c
      call rotate_vec_a(u,v,angle,nx,ny,nz)
c
      return
      end
c
c===============================================================================
c
      subroutine rotate_vec_a(u,v,angle,nx,ny,nz)
c
      include 'trigd.inc'
      
c
      integer nx,ny,nz,i,j,k
c
      real*4    u(nx,ny,nz),
     .          v(nx,ny,nz),
     .          angle(nx,ny),
     .          utmp
c_______________________________________________________________________________
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if (u(i,j,k) .gt. -400. .and. u(i,j,k) .lt. 400. .and.
     .       v(i,j,k) .gt. -400. .and. v(i,j,k) .lt. 400.) then
            utmp    =+u(i,j,k)*cosd(angle(i,j))
     .               +v(i,j,k)*sind(angle(i,j))
            v(i,j,k)=-u(i,j,k)*sind(angle(i,j))
     .               +v(i,j,k)*cosd(angle(i,j))
            u(i,j,k)=utmp
         endif
      enddo
      enddo
      enddo
c
      return
      end
c
c==================================================================================
c
      subroutine rotate_lga_winds(ldir,nx,ny,nz,lon
     1,uw3d,vw3d,uw2d,vw2d)


      implicit none
      integer nx,ny,nz
      integer i,j,k
      logical ldir
      real*4  u_rot,v_rot
      real*4  lon(nx,ny)
      real*4  uw3d(nx,ny,nz)
      real*4  vw3d(nx,ny,nz)
      real*4  uw2d(nx,ny)
      real*4  vw2d(nx,ny)

      if(ldir)then   !from grid north to true north
c 3d
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
            call uvgrid_to_uvtrue(
     1            uw3d(i,j,k),vw3d(i,j,k)
     1           ,u_rot ,v_rot
     1           ,lon(i,j)       )
            uw3d(i,j,k) = u_rot
            vw3d(i,j,k) = v_rot
         enddo
         enddo
         enddo
c 2d
         do j = 1, ny
         do i = 1, nx
            call uvgrid_to_uvtrue(
     1            uw2d(i,j),vw2d(i,j)
     1           ,u_rot   ,v_rot
     1           ,lon(i,j)       )
            uw2d(i,j) = u_rot
            vw2d(i,j) = v_rot
         enddo
         enddo

      else
c 3d
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
            call uvtrue_to_uvgrid(
     1            uw3d(i,j,k),vw3d(i,j,k)
     1           ,u_rot   ,v_rot
     1           ,lon(i,j)       )
            uw3d(i,j,k) = u_rot
            vw3d(i,j,k) = v_rot
         enddo
         enddo
         enddo
c 2d
         do j = 1, ny
         do i = 1, nx
            call uvtrue_to_uvgrid(
     1            uw2d(i,j),vw2d(i,j)
     1           ,u_rot   ,v_rot
     1           ,lon(i,j)       )
            uw2d(i,j) = u_rot
            vw2d(i,j) = v_rot
         enddo
         enddo

      endif
      return
      end
c
c===============================================================================
c
      subroutine rotate_background_uv(nx,ny,nz,lon,gproj,slon0
     +               ,slat1,slat2,uw,vw,uw_sfc,vw_sfc,istatus)
c
c
c

      implicit none
      integer nx,ny,nz
      integer istatus

      real    slon0,slat1,slat2
      real    std_lon,std_lat1,std_lat2

      real    uw(nx,ny,nz)
      real    vw(nx,ny,nz)
      real    uw_sfc(nx,ny)
      real    vw_sfc(nx,ny)
      real    lon(nx,ny)

      character  gproj*2
      character  c6_maproj*6

      call get_c6_maproj(c6_maproj,istatus)
      call get_standard_longitude(std_lon,istatus)
      call get_standard_latitudes(std_lat1,std_lat2,istatus)

      print*,'Rotate u/v components'
      print*

      if(c6_maproj.eq.'merctr')then
         if(gproj.eq.'MC'.or.gproj.eq.'LL')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'no rotation required for background'
         else
c rotate from grid north to LAPS true north (this subroutine uses
c std_lon by virtue of library routines.
            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'Rotate from grid-north to true-north only'

            call rotate_lga_winds(.true.,nx,ny,nz,lon
     +                   ,uw,vw,uw_sfc,vw_sfc)
         endif

      elseif(c6_maproj.eq.'plrstr')then

         if(gproj.eq.'MC'.or.gproj.eq.'LL'.or.gproj.eq.'LE')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'Rotate from true-north to grid-north only'

            call rotate_lga_winds(.false.,nx,ny,nz,lon
     +                   ,uw,vw,uw_sfc,vw_sfc)

         elseif( (gproj.eq.'PS').and.(slon0.eq.std_lon) )then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'no rotation required for background'

         else

c rotate from native grid north to true north
            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'Rotate from bkgd grid-north to LAPS true-north'
            print*,'Rotate from LAPS true-north to LAPS grid-north'

            call uvgrid_to_uvtrue_a(uw,vw,lon,slon0
     +,nx,ny,nz)
            call uvgrid_to_uvtrue_a(uw_sfc,vw_sfc,lon,slon0
     +,nx,ny,1)

c rotate from true north to LAPS grid north
            call rotate_lga_winds(.false.,nx,ny,nz,lon
     +                   ,uw,vw,uw_sfc,vw_sfc)
         endif

      else                         !if(c6_maproj.eq.'lambrt')then

         if(gproj.eq.'MC' .or.
     +      gproj.eq.'LL' .or.
     +      gproj.eq.'LE')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'Rotate from true-north to grid-north only'

            call rotate_lga_winds(.false.,nx,ny,nz,lon
     +                   ,uw,vw,uw_sfc,vw_sfc)

         elseif(gproj.eq.'LC' .and.
     +         (std_lon.eq.slon0).and.
     +         (std_lat1.eq.slat1).and.
     +         (std_lat2.eq.slat2) )then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'no rotation required for background'

         else   !background model is polar stereo

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'Rotate from bkgd grid-north to LAPS true-north'
            print*,'Rotate from LAPS true-north to LAPS grid-north'

c rotate from native grid-north to true north
            call uvgrid_to_uvtrue_a(uw,vw,lon,slon0
     +,nx,ny,nz)
            call uvgrid_to_uvtrue_a(uw_sfc,vw_sfc,lon,slon0
     +,nx,ny,1)

c rotate from true north to LAPS grid-north
            call rotate_lga_winds(.false.,nx,ny,nz,lon
     +                   ,uw,vw,uw_sfc,vw_sfc)
         endif

      endif

      return
      end
c
c =======================================================================
c
      subroutine print_rotproj(gproj,c6_maproj,lon0,lat1,lat2
     +,std_lon,std_lat1,std_lat2)

      character*2 gproj
      character*6 c6_maproj
      real lon0,lat1,lat2
      real std_lon,std_lat1,std_lat2

      print*,'BACKGROUND: /gproj/lon0/lat1/lat2: '
      print*,'            ',gproj,' ',lon0,lat1,lat2

      print*,'ANALYSIS: /c6_maproj/std_lon/std_lat1/std_lat2: '
      print*,'          ',c6_maproj,' ',std_lon,std_lat1,std_lat2

      return
      end
c
c===============================================================================
c
      subroutine thvpc2tq(thv,pc,p,t,q)
c
c *** Subprogram:  thvpc2tq - Calculates temperature (K) and specific 
c                             humidity (kg/kg) given virtual potential 
c                             temperature (K), pressure (mb), and
c                             condensation pressure (mb).
c
c *** Program history log:
c        93-12-20  S. Benjamin - Original version 
c        96-09-17  J. Snook    - es calculated in a table
c
c *** Usage:  call thvpc2tq(thv,pc,p,t,q)
c
c *** Input argument list:
c        thv    - real  virtual potential temperature (K)
c        pc     - real  condensation pressure (mb)
c        p      - real  pressure (mb)
c
c *** Output argument list:
c        t      - real  temperature (K)
c        q      - real  specific humidity (kg/kg)
c
c *** Subprograms called:
c        tv2tq  - calculate temp and spec. hum. from virtual
c                    temp and relative humidity
c        es   - calculate saturation vapor pressure (from a table)
c_______________________________________________________________________________
c
      
c
      real tv,rh,p,t,q,thv,kappa,templcl,x,x1,pc
      integer it
      data kappa/0.285714/
c
      include 'bgdata.inc'
c_______________________________________________________________________________
c
      tv=thv*(p*0.001)**kappa
      templcl=thv*(pc*0.001)**kappa
      it=tv*100
      it=min(45000,max(15000,it))
      x =es(it)
      it=templcl*100
      it=min(45000,max(15000,it))
      x1=es(it)
      rh=x1/x * (p-x) / (pc-x1)
      call tv2tq(tv,rh,p,t,q)
      return
      end
c
c===============================================================================
c
      subroutine tv2tq(tv,rh,p,t,q)
c
c *** Subprogram:  tv2tq - Calculates temperature (K) and specific 
c                          humidity (kg/kg) given virtual temperature (K),
c                          pressure (mb), and relative humidity.
c
c *** Program history log:
c        93-01-12  S. Benjamin - Original version
c
c *** Usage:  call tv2tq(tv,rh,p,t,q)
c
c *** Input argument list:
c        tv     - real  virtual temperature (K)
c        rh     - real  relative humidity (range 0.0-1.0)
c        p      - real  pressure (mb)
c
c *** Output argument list:
c        t      - real  temperature (K) 
c        q      - real  specific humidity (kg/kg)
c
c *** Reamrks:
c        It uses an iterative newton-raphson technique.  Four iterations are
c        generally adequate to provide convergence to 5 decimal places.
c        the wobus function for saturation vapor pressure over liquid water
c        is used.
c_______________________________________________________________________________
c
      
c
      real tv,rh,p,t,q,t1,estv1,etv,t2,estv2,dt,dum
c
      integer j,it
c
      include 'bgdata.inc'
      
c_______________________________________________________________________________
c
      t1 = tv
c
c *** estv = saturation vapor pressure (mb) for tv.
c
      it=t1*100
      it=min(45000,max(15000,it))
      estv1=es(it)
c
      do j=1,3
c
c ****** etv = vapor pressure (mb) for tv and rh*
c
         etv=estv1*rh
c
c ****** q = mixing ratio for tv (kg/kg).
c
         q=0.62197*etv/(p-etv)
         t2=tv/(1.+0.608*q)
         if (abs(t2-t1) .lt. 0.001) goto 77
c
c ****** estv2 = saturation vapor pressure (mb) for estimated t (=t2).
c
         it=t2*100
         it=min(45000,max(15000,it))
         estv2=es(it)
c
c ****** etv = vapor pressure (mb) for estimated t and rh.
c
         etv=estv2*rh
c
c ****** q = mixing ratio for estimated t and rh (kg/kg).
c
         q=0.62197*etv/(p-etv)
         t=t2*(1.+0.608*q)
c
c ****** Recalc. tv.
c
         dt=tv-t
         dum=(estv2-estv1)/(t2-t1)
         etv=estv2+dum*dt*rh
c
c ****** Reset t1 and estv1 before next iteration.
c
         t1=t2
         estv1=estv2
c
      enddo
77    continue
      t=t2
c
      return
      end
