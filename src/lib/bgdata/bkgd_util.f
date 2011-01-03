      subroutine rotate_lga_winds(ldir,bgmodel,cmodel,fullname
     1,gproj,nx,ny,nz,lon,uw3d,vw3d,uw2d,vw2d,lgb_only)

      use mem_namelist
      implicit none

      character fullname*200
      character cmodel*132
      character gproj*2
      integer nx,ny,nz
      integer i,j,k
      integer bgmodel
      integer istatus,istatus_rot,ishow_timer
      logical ldir,lgb_only
      real    u_true,v_true
      real    u_grid,v_grid
      real    u_true_2d(nx,ny),v_true_2d(nx,ny)
      real    u_grid_2d(nx,ny),v_grid_2d(nx,ny)
      real    lon(nx,ny)
      real    uw3d(nx,ny,nz)
      real    vw3d(nx,ny,nz)
      real    uw2d(nx,ny)
      real    vw2d(nx,ny)
      real    angle(nx,ny,2)  !both grid to true N and true to grid N.
      real    latitude(nx,ny),projrot_latlon(nx,ny)

c
c reset or restore common projection parameters
c
      call reset_lapsparms_common(ldir,bgmodel,cmodel,fullname
     1,gproj,istatus)
c
c build look-up-tables for rotation angles in 2D grid and
c apply this to grid winds
c
!     istatus_rot=ishow_timer()

      latitude = -999. ! Since lat is not yet passed in
      call projrot_latlon_2d(latitude,lon,nx,ny,projrot_latlon,istatus)

      do j = 1, ny
      do i = 1, nx
         angle(i,j,1)= -projrot_latlon(i,j)                        !grid 2 true
         angle(i,j,2)= -angle(i,j,1)                               !true 2 grid
      enddo
      enddo

      write(6,*)'rotate_lga_winds - built lookup table'
!     istatus_rot=ishow_timer()

      if(ldir)then   !from grid to true North
c 3d
        if(.not. lgb_only)then
          do k = 1, nz
          do j = 1, ny
          do i = 1, nx

            call     rotate_vec(uw3d(i,j,k),
     1                          vw3d(i,j,k),
     1                          u_true,
     1                          v_true,
     1                          angle(i,j,1))
            uw3d(i,j,k) = u_true
            vw3d(i,j,k) = v_true

          enddo
          enddo
          enddo
        endif
c 2d
        call   rotate_vec_2d(uw2d,
     1                       vw2d,
     1                       u_true_2d,
     1                       v_true_2d,
     1                       angle(1,1,1),nx,ny)

        uw2d = u_true_2d
        vw2d = v_true_2d

      else !rotate from true to grid N.
c 3d
        if(.not. lgb_only)then
          do k = 1, nz
          do j = 1, ny
          do i = 1, nx

            call     rotate_vec(uw3d(i,j,k),
     1                          vw3d(i,j,k),
     1                          u_grid,
     1                          v_grid,
     1                          angle(i,j,2))
            uw3d(i,j,k) = u_grid
            vw3d(i,j,k) = v_grid

          enddo
          enddo
          enddo
        endif
c 2d
        call rotate_vec_2d(uw2d,
     1                     vw2d,
     1                     u_grid_2d,
     1                     v_grid_2d,
     1                     angle(1,1,2),nx,ny)

        uw2d = u_grid_2d
        vw2d = v_grid_2d

      endif

      write(6,*)'end of rotate_lga_winds, lgb_only = ',lgb_only
!     istatus_rot=ishow_timer()

      return
      end
c
c=======================================================================
c
      subroutine rotate_background_uv(nx,ny,nz,lon,bgmodel,cmodel
     +,fullname,gproj,slon0,slat1,slat2,uw,vw,uw_sfc,vw_sfc,lgb_only
     +,istatus)
c
c
c

      use mem_namelist
      implicit none

      integer nx,ny,nz
      integer bgmodel
      integer istatus

      real    slon0,slat1,slat2
      real    std_lon,std_lat1,std_lat2

      real    uw(nx,ny,nz)
      real    vw(nx,ny,nz)
      real    uw_sfc(nx,ny)
      real    vw_sfc(nx,ny)
      real    lon(nx,ny)

      character  gproj*2
      character  fullname*200
      character  cmodel*132

c     character  c6_maproj*6

c     call get_c6_maproj(c6_maproj,istatus)
c     call get_standard_longitude(std_lon,istatus)
c     call get_standard_latitudes(std_lat1,std_lat2,istatus)

      logical lgb_only

      print*,'Rotate u/v components'
      print*

c For the case when rotating the bkgd grid winds, we use subroutine
c reset_lapsparms_common to reset the appropriate projection parameters
c within lapsparms.cmn.

      std_lon =standard_longitude
      std_lat1=standard_latitude
      std_lat2=standard_latitude2

c ----------------------------------------------------------------
      if(c6_maproj.eq.'merctr')then
c ----------------------------------------------------------------

         if(gproj.eq.'MC'.or.gproj.eq.'LL')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
            print*,'no rotation required for background'

         else

c rotate from grid north to LAPS true north (this subroutine uses
c std_lon by virtue of library routines.

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'Rotate grid-north (bkgd) to true-north (LAPS)'

            call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         endif
c ----------------------------------------------------------------
      elseif(c6_maproj.eq.'plrstr')then
c ----------------------------------------------------------------

         if(gproj.eq.'MC'.or.gproj.eq.'LL'.or.gproj.eq.'LE')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'Rotate true-north (bkgd) to grid-north (LAPS)'

            call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         elseif(gproj.eq.'PS')then

             if(slon0.eq.std_lon)then

c also, this only good if both domains share a common pole point.
c currently we are only considering polar stereo, not local stereo.

                call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)
                print*,'no rotation required for background'

             else

                call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

                print*,'Rotate grid-north (bkgd) to true-north'
                call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

                print*,'Rotate true-north to grid-north (LAPS)'
                call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

             endif

         else

c background is lambert

             call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

             print*,'Rotate grid-north (bkgd) to true-north'

             call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

             print*,'Rotate true-north to grid-north (LAPS)'

             call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         endif

c ----------------------------------------------------------------
      else     !LAPS is Lambert --- if(c6_maproj.eq.'lambrt')then
c ----------------------------------------------------------------

         if(gproj.eq.'MC' .or.
     +      gproj.eq.'LL' .or.
     +      gproj.eq.'LE')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'Rotate true-north (bkgd) to grid-north (LAPS)'   ! because MC/LL/LE grids are true north

            call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

         elseif(gproj.eq.'PS'.or.gproj.eq.'LC')then

            call print_rotproj(gproj,c6_maproj,slon0,slat1,slat2
     +,std_lon,std_lat1,std_lat2)

            print*,'Rotate grid-north (bkgd) to true-north'

            call rotate_lga_winds(.true.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

             print*,'Rotate true-north to grid-north (LAPS)'

             call rotate_lga_winds(.false.,bgmodel,cmodel,fullname
     +,gproj,nx,ny,nz,lon,uw,vw,uw_sfc,vw_sfc,lgb_only)

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
c
c ---------------------------------
c
      subroutine reset_lapsparms_common(ldir,bgmodel,cmodel
     1,fullname,gproj,istatus)
c
c acquires background projection parameters and uses them to
c temporarily replace lapsparms parameters. It is then possible
c to use the library projection routines to calculate wind rotation
c angles for the background.
c
      use mem_namelist
      implicit none

      character*200 fullname
      character*132 cmodel
      character*2   gproj
      character*1   cgrddef
      integer       istatus
      integer       bgmodel
      integer       nxbg,nybg
      integer       nzbg_ht
      integer       nzbg_tp
      integer       nzbg_sh
      integer       nzbg_uv
      integer       nzbg_ww
      logical       ldir
      real          dlat,dlon
      real          sw(2),ne(2)
      real          dxbg,dybg
      real          La1,Lo1,La2,Lo2
      real          Xcen,Ycen
      real          Xsw,Ysw,Xne,Yne
      real          erad
      real          rlatcen,rloncen

      character*6   c6_maproj_save
      save          c6_maproj_save
      real          centrallat_save
      save          centrallat_save
      real          centrallon_save
      save          centrallon_save
      real          Lat0_save
      save          Lat0_save
      real          Lat1_save
      save          Lat1_save
      real          Lon0_save
      save          Lon0_save

      integer       isave
      data          isave/0/
      save          isave



      if(isave.eq.0)then
         if(ldir)then
            c6_maproj_save=c6_maproj
            centrallat_save=grid_cen_lat
            centrallon_save=grid_cen_lon
            Lat0_save=standard_latitude
            Lat1_save=standard_latitude2
            Lon0_save=standard_longitude
            isave=1
         else
            print*,'Return to rotate_background_uv: isave=0; ldir=false'
            return 
         endif
      elseif(isave.eq.1)then        
c Restore original nest7grid.parms settings
         c6_maproj=c6_maproj_save
         grid_cen_lat=centrallat_save
         grid_cen_lon=centrallon_save
         standard_latitude=Lat0_save
         standard_latitude2=Lat1_save
         standard_longitude=Lon0_save
         isave=0
         return
      endif

      print*
      print*,'reset_lapsparms_common: isave = ',isave
      print*
      print*,'bgmodel: ',bgmodel
      print*,'cmodel: ', TRIM(cmodel)
      print*,'fullname:',TRIM(fullname)
      print*

      call get_bkgd_mdl_info(bgmodel,cmodel,fullname
     &,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,grid_cen_lat,grid_cen_lon
     &,dxbg,dybg,standard_latitude,standard_latitude2
     &,standard_longitude,sw,ne,cgrddef,istatus)

c Temporarily set c6_maproj.
      if(gproj.eq.'LC')c6_maproj='lambrt'
      if(gproj.eq.'PS')c6_maproj='plrstr'
      if(gproj.eq.'MC')c6_maproj='merctr' 
cc
c if RUC_NATIVE
      if(bgmodel.eq.5 .or. bgmodel.eq.3 .or. bgmodel.eq.13)then
         if(TRIM(cmodel).eq.'CWB_20FA_LAMBERT_NF'.or.
     +      TRIM(cmodel).eq.'CWB_20FA_LAMBERT_RE'.or.
     +      TRIM(cmodel).eq.'RUC40_NATIVE'.or.
     +      TRIM(cmodel).eq.'RUC')then  
            La1=sw(1)
            Lo1=sw(2)
            La2=ne(1)
            Lo2=ne(2)
            print*,'*** SW lat/lon = ',La1,Lo1,'***'
            print*,'*** NE lat/lon = ',La2,Lo2,'***'
            call get_earth_radius(erad,istatus)
c Get X/Y for grid center
            call latlon_to_xy(La1,Lo1,ERAD,Xsw,Ysw) 
            call latlon_to_xy(La2,Lo2,ERAD,Xne,Yne) 
            Xcen=(Xne+Xsw)/2.
            Ycen=(Yne+Ysw)/2.
            call xy_to_latlon(Xcen,Ycen,erad,rlatcen,rloncen)
            print*,'*** Center lat/lon= ',rlatcen,rloncen,' ***'
         endif
      endif

      return
      end
