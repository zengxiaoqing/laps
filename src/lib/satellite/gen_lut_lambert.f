cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
       subroutine gen_lut_lambert(isat,jtype,kchl,
     &nx_l,ny_l,lat,lon,ri_laps,rj_laps,jstatus)
c
c
      implicit none

      integer   nx_l,ny_l

      real    lat(nx_l,ny_l)
      real    lon(nx_l,ny_l)   ! laps lat/lon data  -  input

      real    xlat(nx_l+2,ny_l+2)
      real    xlon(nx_l+2,ny_l+2) !expanded domain lats/lons
 
      real    dx,dy  !both in km from /public and /WFO SBN
      real    du,dv
      real    u_orig,v_orig,sioffset,sjoffset
      real    reflat,reflon,refi,refj,refu,refv
      real    pi
      real    u1,v1,u2,v2,rlats,rlons

      real    lapterm
      real    lovterm
      real    latterm
      real    lonterm
      real    r_missing_data
      real    rls,rle,res,ree
      real    rlatin,rlap,rlov
      real    centerlat,centerlon
      real    rla1,rlo1
      real    rla1nxny,rlo1nxny
      real    rla100,rlo100
      real    rlatdxdy,rlondxdy

      real    ri(nx_l+2,ny_l+2)
      real    rj(nx_l+2,ny_l+2)
      real    rel_ri(nx_l+2,ny_l+2)
      real    rel_rj(nx_l+2,ny_l+2)
      real    ri_laps(nx_l,ny_l)
      real    rj_laps(nx_l,ny_l)
      real    u,v,u0,v0,uscale,vscale
      real    usmin,usmax,vsmin,vsmax
      real    usat,vsat

      integer center_id
      integer process_id,wmo_sat_id
      double precision reftime, valtime
      character*132 earth_shape, grid_name, grid_type,
     +     origin_name, process_name, wavelength, x_dim, y_dim

      integer isat,jtype,kchl
c     integer i,j,n,nc
      integer i,j,nc,nt
      integer ii,jj
      integer indx
      integer n1,np
      integer i1,j1
      integer nx,ny
c     integer lend
      integer istatus
      integer jstatus
      integer istatus_wp
      integer nx3,ny3
      integer nx3mx,ny3mx
      integer linestart,lineend
      integer elemstart,elemend
      integer i4time_latest,ifstat
      integer nijout
      integer ierr

      logical lpoint

      character*200 table_path
      character*255 path
      character*255 fname
      character*200 cname
      character*9   cftime9
c     character*6   csatid  !satellite data identifier {goes08, goes09, meteos, etc}
      character*3   cdtype !satellite data type       {'cdf', 'gvr', 'gwc', 'asc', etc}
      character*3   ct     !                         {'vis', 'ir', or 'wv'}
c
c ===========================================================
c the following include file contains static navigation parameters for
c all possible satellites.  It contains fortran logic that uses csatid,
c and csattyp to set the appropriate navigation parameters.
c
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      logical     lwrite
      data lwrite /.false./

c     logical     lfirst(maxtype,maxsat)              !4 types x 3 sats (5-12-98) 6 sats (1-31-03)
c     data lfirst /.false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false./
c     save lfirst

      integer nxlc,nylc,nzlc
      real    lat1,lat2,lon0,sw(2),ne(2)
      common /lcgrid/nxlc,nylc,nzlc,lat1,lat2,lon0,sw,ne


      jstatus = -1
      cdtype=c_sat_types(jtype,isat)

      ct=c_channel_types(kchl,jtype,isat)
      nc=index(ct,' ')-1
      if(nc.le.0)nc=3
      call lvd_file_specifier(ct,indx,istatus)
c
c retrieve the latest nav info from a file if it exists
c
      if(cdtype.eq.'wfo')then

         call get_wfo_nav_parms(path_to_raw_sat(kchl,jtype,isat),
     &ct,centerlat,centerlon,rla100,rlo100,rla1nxny,rlo1nxny,
     &rlatdxdy,rlondxdy,dx,dy,nx3mx,ny3mx,istatus_wp)

         if(istatus_wp.eq.1)then

            print*,'WARNING:, using namelist values for wfo mapping'
            print*,'because new wfo nav parameters were not obtained'
            rla100 = r_la1(jtype,isat)
            rlo100 = r_lo1(jtype,isat)
            centerlat = r_lap(jtype,isat)
            centerlon = r_lov(jtype,isat)

            if(indx.eq.1)then

              dx = r_resolution_x_vis(jtype,isat)
              dy = r_resolution_y_vis(jtype,isat)
              nx3= n_pixels_vis(jtype,isat)
              ny3= n_lines_vis(jtype,isat)

            elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

c             if(.not.lfirst(jtype,isat))then
c --- experimental: go ahead and compute these for all satellite types
                 dx = r_resolution_x_ir(jtype,isat)
                 dy = r_resolution_y_ir(jtype,isat)
                 nx3 = n_pixels_ir(jtype,isat)
                 ny3 = n_lines_ir(jtype,isat)

c                lfirst(jtype,isat)=.true.
c             else
c                jstatus = 0
c                goto 1000
c             endif

            elseif(indx.eq.3)then

              dx = r_resolution_x_wv(jtype,isat)
              dy = r_resolution_y_wv(jtype,isat)
              nx3= n_pixels_wv(jtype,isat)
              ny3= n_lines_wv(jtype,isat)

            endif

         elseif(istatus_wp.eq.0)then

c --- experimental: set rlap/rlov for all satellite types.
c     Note: types 2, 4, and 5 should all be the same lambert projections.
c           if(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then
c           if(.not.lfirst(jtype,isat))then
c               lfirst(jtype,isat)=.true.
c           else
c               jstatus = 0
c               goto 1000
c           endif
c           endif

            rlap=centerlat ! assume standard lat = central lat
            rlov=centerlon ! assume standard lon = central lon

         else
            print*,'Error returned from get_wfo_nav_parms'
            goto 1000
         endif

c goes noaaport
      elseif(cdtype.eq.'gnp')then

         if(indx.eq.1)then    
            dx = r_resolution_x_vis(jtype,isat) / 1000. ! km
            dy = r_resolution_y_vis(jtype,isat) / 1000. ! km
            nx3 = n_pixels_vis(jtype,isat)
            ny3 = n_lines_vis(jtype,isat)

         elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then
            dx = r_resolution_x_ir(jtype,isat) / 1000. ! km
            dy = r_resolution_y_ir(jtype,isat) / 1000. ! km
            nx3 = n_pixels_ir(jtype,isat)
            ny3 = n_lines_ir(jtype,isat)

         elseif(indx.eq.3)then

            dx = r_resolution_x_wv(jtype,isat) / 1000. ! km
            dy = r_resolution_y_wv(jtype,isat) / 1000. ! km
            nx3 = n_pixels_wv(jtype,isat)
            ny3 = n_lines_wv(jtype,isat)

         endif

         nx3mx=nx3
         ny3mx=ny3

         rlap = r_lap(jtype,isat)
         rlov = r_lov(jtype,isat)
!        call get_attribute_gnp()

c not wfo data type
      elseif(cdtype.eq.'cdf')then

         if(indx.eq.1)then

            dx = r_resolution_x_vis(jtype,isat)
            dy = r_resolution_y_vis(jtype,isat)
            nx3 = n_pixels_vis(jtype,isat)
            ny3 = n_lines_vis(jtype,isat)

         elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

c See experimental comment above.
c           if(.not.lfirst(jtype,isat))then

               dx = r_resolution_x_ir(jtype,isat)
               dy = r_resolution_y_ir(jtype,isat)
               nx3 = n_pixels_ir(jtype,isat)
               ny3 = n_lines_ir(jtype,isat)

c              lfirst(jtype,isat)=.true.
c           else
c              jstatus = 0
c              goto 1000
c           endif

         elseif(indx.eq.3)then

            dx = r_resolution_x_wv(jtype,isat)
            dy = r_resolution_y_wv(jtype,isat)
            nx3 = n_pixels_wv(jtype,isat)
            ny3 = n_lines_wv(jtype,isat)

         endif

         call s_len(path_to_raw_sat(kchl,jtype,isat),np)
         fname=path_to_raw_sat(kchl,jtype,isat)(1:np)
         fname=fname(1:np)//'*_'//ct
         call get_latest_file_time(fname,i4time_latest)
         call make_fnam_lp(i4time_latest,cftime9,ifstat)
         fname=path_to_raw_sat(kchl,jtype,isat)(1:np)
         fname=fname(1:np)//cftime9//'_'//ct
         call  rdcdfhead(fname, nx3mx, ny3mx, center_id, 
     +     process_id, wmo_sat_id, Dx, Dy, rla100, rLatin, rlo100, 
     +     rlov, reftime, valtime, earth_shape, grid_name, grid_type, 
     +     origin_name, process_name, wavelength, x_dim, y_dim)

         rlap=rlatin

      else  !taiwan gms satellite data

         nxlc = 512 
         nylc = 512
         nx3mx=nxlc
         ny3mx=nylc
         lat1=10.
         lat2=40.
         lon0=+120.
         rla100=lat1
         rlo100=lon0
         sw(1)=5.119949
         sw(2)=+104.1190
         ne(1)=44.09481
         ne(2)=+153.3814
         dx = 8.77
         dy = dx

      endif

      write(6,*)'Satellite Nav Parameters (gen_lut_lambert) ',ct
      write(6,*)'dx (km)',dx
      write(6,*)'dy (km)',dy
      write(6,*)'nx3    ',nx3
      write(6,*)'ny3    ',ny3
      write(6,*)'nx3mx  ',nx3mx
      write(6,*)'ny3mx  ',ny3mx
      write(6,*)'rla100 ',rla100   ! Lower left latitude
      write(6,*)'rlo100 ',rlo100   ! Lower left longitude
      write(6,*)'rlatin ',rlatin
      write(6,*)'rlov   ',rlov     ! Standard longitude
      write(6,*)'rlap   ',rlap     ! Standard latitude
c
c expand domain lats/lons: extra row and column "laps domain" used
c to build the domain relative look up table
c
      nx=nx_l+2
      ny=ny_l+2
      pi = acos(-1.0)
 
      call expand_domain(nx_l,ny_l,lat,lon,nx,ny,xlat,xlon,
     &istatus)
c
c new code to compute lambert projected satellite data for
c laps domain as specified in laps lat/lon arrays. Follows
c method in lib/latlon_to_rlaps.f by Albers.
c
      if(cdtype.eq.'wfo')then

         call latlon_to_uv_lc(rla100,rlo100,rlap,rlap
     &,rlov,usmin,vsmin)
         call latlon_to_uv_lc(rla1nxny,rlo1nxny,rlap
     &,rlap,rlov,usmax,vsmax)
         uscale = (usmax - usmin)/(float(nx3mx)-1.)
         vscale = (vsmax - vsmin)/(float(ny3mx)-1.)
         u0 = usmin - uscale
         v0 = vsmin - vscale

         do j=1,ny
         do i=1,nx

            call latlon_to_uv_lc(xlat(i,j),xlon(i,j),rlap,rlap
     &,rlov,usat,vsat)

            ri(i,j)=(usat - u0)/uscale
            rj(i,j)=(vsat - v0)/vscale

         enddo
         enddo

      elseif(cdtype.eq.'cdf')then

c use original lambert software for fsl-conus in FSL's /public
         write(6,*)' calling getdudv_lam'

         call getdudv_lam(rlov,rlap,dx,dy,
     &rla100,rlo100,du,dv,u_orig,v_orig)
         lapterm=(rlap*pi)/180.
         lovterm=(rlov*pi)/180.
         do j = 1, ny
         do i = 1, nx
            latterm=(xlat(i,j)*pi)/180.
            lonterm=(xlon(i,j)*pi)/180.
            call getuv_lam(lapterm,lovterm,latterm,lonterm,u,v)
            call uv_ij(ny3mx,u_orig,v_orig,du,dv,u,v
     +,ri(i,j),rj(i,j))
         enddo
         enddo

      elseif(cdtype.eq.'gnp')then

c use original lambert software for fsl-conus in FSL's /public
         if(indx.eq.1)then                               ! vis
            reflat = 42.09345
            reflon = -106.13704
            refi = 512.5 + 1024. ! center of center tile on upper row
            refj = 512.5 + 1024. 
         elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then ! ir
!           reflat = 45.693
!           reflon = -115.980
            reflat = 45.8778
            reflon = -114.5280
            refi = 512.5          ! center of left tile
            refj = 512.5
         endif

         write(6,*)' calling getdudv_lam'

         lapterm=(rlap*pi)/180.
         lovterm=(rlov*pi)/180.

         if(.false.)then ! use center point
!          Initial calculation of 'u_orig' and 'v_orig' for center
           call getdudv_lam(rlov,rlap,dx,dy,
     &                      centerlat,centerlon,du,dv,u_orig,v_orig)

           write(6,*)' u_orig/v_orig at center ',u_orig,v_orig

!          Recalculate 'u_orig' and 'v_orig' for LL sat corner
           sioffset = float(nx3mx-1) / 2.
           sjoffset = float(ny3mx-1) / 2.
           u_orig = u_orig - du * sioffset 
           v_orig = v_orig - dv * sjoffset

         else ! use reference point
           call getdudv_lam(rlov,rlap,dx,dy,
     &                      reflat,reflon,du,dv,refu,refv)

           write(6,*)' refu/refv at reference point ',refu,refv
           
!          Recalculate 'u_orig' and 'v_orig' for LL sat corner
           sioffset = refi - 1.
           sjoffset = refj - 1.
           u_orig = refu - du * sioffset 
           v_orig = refv - dv * sjoffset

         endif

         write(6,*)' nx/ny (model grid) = ',nx,ny
         write(6,*)' u_orig/v_orig at corner ',u_orig,v_orig
         u1 = u_orig
         u2 = u_orig + du * float(nx3mx-1)
         v1 = v_orig
         v2 = v_orig + dv * float(ny3mx-1)
         call uv_to_latlon_lc(u1,v1,rlap,rlap,rlov,rlats,rlons)
         write(6,*)' lat/lon at ll sat corner ',rlats,rlons

         call uv_to_latlon_lc(u2,v1,rlap,rlap,rlov,rlats,rlons)
         write(6,*)' lat/lon at lr sat corner ',rlats,rlons

         call uv_to_latlon_lc(u1,v2,rlap,rlap,rlov,rlats,rlons)
         write(6,*)' lat/lon at ul sat corner ',rlats,rlons

         call uv_to_latlon_lc(u2,v2,rlap,rlap,rlov,rlats,rlons)
         write(6,*)' lat/lon at ur sat corner ',rlats,rlons

!        For each model grid point obtain the satellite pixel         
         do j = 1, ny
         do i = 1, nx
            latterm=(xlat(i,j)*pi)/180.
            lonterm=(xlon(i,j)*pi)/180.
            call getuv_lam(lapterm,lovterm,latterm,lonterm,u,v)
            call uv_ij(ny3mx,u_orig,v_orig,du,dv,u,v
     +                ,ri(i,j),rj(i,j))
         enddo
         enddo

      elseif(cdtype.eq.'twn')then           !only other type is taiwan gms

         call  latlon_2_lcij(nx*ny,xlat,xlon,ri,rj)

      endif ! cdtype

      write(6,*)'Sat ri/rj for corners of expanded model domain'
      write(6,*)'ri1/rj1 (SW) ',ri(1,1),rj(1,1)
      write(6,*)'ri2/rj2 (SE) ',ri(nx,1),rj(nx,1)
      write(6,*)'ri3/rj3 (NW) ',ri(1,ny),rj(1,ny)
      write(6,*)'ri4/rj4 (NE) ',ri(nx,ny),rj(nx,ny)
c
c get new i/j start/end values for this domain
c
      call get_sat_boundary(nx,ny,nx,ny,0,ny3mx,nx3mx,
     &ri,rj,linestart,lineend,elemstart,elemend,
     &rls,rle,res,ree,istatus)

      if(istatus.ne.1)then
         write(6,*)'WARNING: Laps domain outside sat data cover!'
      endif
c
c compute ri, rj relative look up table for the block of data surrounding
c the laps domain.
c
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error getting r_missing_data'
         return
      endif

      do j = 1,ny
      do i = 1,nx
       if(ri(i,j).ne.r_missing_data.and.rj(i,j).ne.r_missing_data)then
          rel_ri(i,j) = ri(i,j) - res
          rel_rj(i,j) = rj(i,j) - rls 
          if(lpoint)then
             i1=i
             j1=j
             lpoint=.false.
          endif
       else
          rel_ri(i,j) = r_missing_data
          rel_rj(i,j) = r_missing_data
          nijout=nijout+1
       endif

      enddo
      enddo

      if(nijout.gt.0)then
         print*,'Found ',nijout,' points outside domain'
      endif
c
c put the expanded domain ri/rj's into the original laps domain
c
      jj = 0
      do j = 2,ny-1
      jj = jj+1
      ii = 0
      do i = 2,nx-1
         ii = ii+1
         ri_laps(ii,jj) = rel_ri(i,j)
         rj_laps(ii,jj) = rel_rj(i,j)
      enddo
      enddo

      if(lwrite)then
        do i = 1,nx_l,10
        do j = 1,ny_l,10
          print*,'i,j,ri,rj: ',i,j,ri_laps(i,j),rj_laps(i,j)
        enddo
        enddo
      endif

      call get_directory('static',path,n1)
      path=path(1:n1)//'/lvd/'
      n1=index(path,' ')-1

      nt=3
      if(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then
         ct='ir'
         nt=2
      endif

c     cname=path(1:n1)//c_sat_id(isat)//'-llij-'//ct(1:nt)
c     n1=index(cname,' ')-1
c     table_path = cname(1:n1)//'-'//cdtype//'.lut'
c     n1=index(table_path,' ')
c     write(6,*)'Write lat/lon to i/j look up table'
c     write(6,*)table_path(1:n1)
c     call write_table (table_path,nx_l,ny_l,xlat,xlon,
c    &ri_laps,rj_laps,istatus)
c     if(istatus .ne. 1)then
c        write(6,*)'Error writing look-up table'
c        goto 1000
c     endif

      if(elemstart.le.0)elemstart=1
      if(elemend.gt.nx3mx)elemend=nx3mx
      if(linestart.le.0)linestart=1
      if(lineend.gt.ny3mx)lineend=ny3mx


      r_lap(jtype,isat) = rlap
      r_lov(jtype,isat) = rlov
      r_latin(jtype,isat) = rlatin

      if(indx.eq.1)then

         i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_la1(jtype,isat) = rla100
         r_lo1(jtype,isat) = rlo100
         r_resolution_x_vis(jtype,isat) = dx*1000.
         r_resolution_y_vis(jtype,isat) = dy*1000.
         n_pixels_vis(jtype,isat) = nx3mx
         n_lines_vis(jtype,isat)  = ny3mx

      elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

         i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_la1(jtype,isat) = rla100
         r_lo1(jtype,isat) = rlo100
         r_resolution_x_ir(jtype,isat) = dx*1000.
         r_resolution_y_ir(jtype,isat) = dy*1000.
         n_pixels_ir(jtype,isat) = nx3mx
         n_lines_ir(jtype,isat)  = ny3mx

      elseif(indx.eq.3)then

         i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_la1(jtype,isat) = rla100
         r_lo1(jtype,isat) = rlo100
         r_resolution_x_wv(jtype,isat) = dx*1000.
         r_resolution_y_wv(jtype,isat) = dy*1000.
         n_pixels_wv(jtype,isat) = nx3mx
         n_lines_wv(jtype,isat)  = ny3mx

      endif

      write(6,1)ct,indx,elemstart,elemend,linestart,lineend
1     format(' gen_lut_lambert indx i/j start/end ',a,5i6)
         
      jstatus = 1

1000  return
      end
