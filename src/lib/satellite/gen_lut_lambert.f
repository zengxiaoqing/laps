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
     &nx_l,ny_l,xlat,xlon,jstatus)
c
c
      implicit none

      integer   nx_l,ny_l

      real*4    xlat(nx_l,ny_l)
      real*4    xlon(nx_l,ny_l)   ! laps lat/lon data  -  input

      real*4    lat(nx_l+2,ny_l+2)
      real*4    lon(nx_l+2,ny_l+2) !expanded domain lats/lons
 
      real*4    dx,dy
      real*4    du,dv
      real*4    u_orig,v_orig
      real*4    pi
      real*4    u1,v1
      real*4    resx,resy

      real*4    ri1,rj1
      real*4    ri2,rj2
      real*4    ri3,rj3
      real*4    ri4,rj4

      real*4    lapterm
      real*4    lovterm
      real*4    latterm
      real*4    lonterm
      real*4    dxterm
      real*4    dyterm
      real*4    wdw_lat
      real*4    wdw_lon
      real*4    rls,rle,res,ree
      real*4    rlatin,rlap,rlov
      real*4    rla1,rlo1

      real*4    ri(nx_l+2,ny_l+2)
      real*4    rj(nx_l+2,ny_l+2)
      real*4    rel_ri(nx_l+2,ny_l+2)
      real*4    rel_rj(nx_l+2,ny_l+2)
      real*4    ri_laps(nx_l,ny_l)
      real*4    rj_laps(nx_l,ny_l)
      real*4    u
      real*4    v

      integer isat,jtype,kchl
      integer i,j,n,nc
      integer ii,jj
      integer indx
      integer n1
      integer nx,ny
      integer lend
      integer istatus
      integer jstatus
      integer istatus_wp
      integer nx3,ny3
      integer nx3mx,ny3mx
      integer linestart,lineend
      integer elemstart,elemend
      integer idum

      character*200 table_path
      character*255 path
      character*200 cname
      character*6   csatid  !satellite data identifier {goes08, goes09, meteos, etc}
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

      logical     lfirst(maxtype,maxsat)              !4 types x 3 sats (5-12-98)
      data lfirst /.false.,.false.,.false.,.false.,
     &             .false.,.false.,.false.,.false.,
     &             .false.,.false.,.false.,.false./
      save

      jstatus = -1
      cdtype=c_sat_types(jtype,isat)

      ct=c_channel_types(kchl,jtype,isat)
      nc=index(ct,' ')-1
      if(nc.le.0)nc=3
      call lvd_file_specifier(ct,indx,istatus)
c
c retrieve the latest nav info from a file if it exists
c Note: for wfo data we are using the NW corner and not the SW
c as our reference point.
c
      if(cdtype.eq.'wfo')then

         call get_wfo_nav_parms(path_to_raw_sat(kchl,jtype,isat),
     &ct,rla1,rlo1,dx,dy,nx3,ny3,istatus_wp)

         if(istatus_wp.le.0)then

            print*,'Warning, using namelist values for wfo mapping'
            rla1 = r_la1(jtype,isat)
            rlo1 = r_lo1(jtype,isat)

            goto(52,53,54,53,53)indx
52            dx = r_resolution_x_vis(jtype,isat)
              dy = r_resolution_y_vis(jtype,isat)
              nx3= n_pixels_vis(jtype,isat)
              ny3= n_lines_vis(jtype,isat)
              goto 55

53            if(.not.lfirst(jtype,isat))then
                 ct='ir'
                 dx = r_resolution_x_ir(jtype,isat)
                 dy = r_resolution_y_ir(jtype,isat)
                 nx3 = n_pixels_ir(jtype,isat)
                 ny3 = n_lines_ir(jtype,isat)
                 lfirst(jtype,isat)=.true.
              else
                 jstatus = 0
                 goto 1000
              endif
              goto 55

54            dx = r_resolution_x_wv(jtype,isat)
              dy = r_resolution_y_wv(jtype,isat)
              nx3= n_pixels_wv(jtype,isat)
              ny3= n_lines_wv(jtype,isat)

55          continue

         elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

            if(.not.lfirst(jtype,isat))then
                ct='ir'
                lfirst(jtype,isat)=.true.
            else
                jstatus = 0
                goto 1000
            endif

         endif

c not wfo data type
      else

         goto(61,62,63,62,62)indx

61          dx = r_resolution_x_vis(jtype,isat)
            dy = r_resolution_y_vis(jtype,isat)
            nx3 = n_pixels_vis(jtype,isat)
            ny3 = n_lines_vis(jtype,isat)
            goto 65

62          if(.not.lfirst(jtype,isat))then
               ct='ir'
               dx = r_resolution_x_ir(jtype,isat)
               dy = r_resolution_y_ir(jtype,isat)
               nx3 = n_pixels_ir(jtype,isat)
               ny3 = n_lines_ir(jtype,isat)
               lfirst(jtype,isat)=.true.
            else
               jstatus = 0
               goto 1000
            endif
            goto 65

63          dx = r_resolution_x_wv(jtype,isat)
            dy = r_resolution_y_wv(jtype,isat)
            nx3 = n_pixels_wv(jtype,isat)
            ny3 = n_lines_wv(jtype,isat)

65       continue

         rla1=r_la1(jtype,isat)
         rlo1=r_lo1(jtype,isat)
 
      endif

      rlatin = r_latin(jtype,isat)
      rlov = r_lov(jtype,isat)
      rlap = r_lap(jtype,isat)
 
      if(cdtype.eq.'wfo')then
         nx3mx=nx3
         ny3mx=ny3
         nx3=1
         ny3=1
      else
         nx3mx=nx3
         ny3mx=ny3
      endif

      write(6,*)'Satellite Nav Parameters'
      write(6,*)'dx     ',dx
      write(6,*)'dy     ',dy
      write(6,*)'nx3    ',nx3
      write(6,*)'ny3    ',ny3
      write(6,*)'nx3mx  ',nx3mx
      write(6,*)'ny3mx  ',ny3mx
      write(6,*)'la1    ',rla1
      write(6,*)'lo1    ',rlo1
      write(6,*)'latin  ',rlatin
      write(6,*)'Lov    ',rlov
      write(6,*)'lap    ',rlap
c
c expand domain lats/lons: extra row and column "laps domain" used
c to build the domain relative look up table
c
      nx=nx_l+2
      ny=ny_l+2
      pi = acos(-1.0)
 
      call expand_domain(nx_l,ny_l,xlat,xlon,nx,ny,lat,lon,
     &istatus)
c
c build the absolute lut for the slightly larger domain
c
      dxterm = dx/1000.
      dyterm = dy/1000.
      call getdudv_lam(rlov,rlap,dxterm,dyterm,
     &rla1,rlo1,du,dv,u_orig,v_orig)
c
c compute corner point ri/rj pairs used to determine # of i and j
c
      write(6,*)'Compute corner points'

      lapterm=(rlap*pi)/180.
      lovterm=(rlov*pi)/180.
      latterm=(lat(1,1)*pi)/180.
      lonterm=(lon(1,1)*pi)/180.

      call getuv_lam(lapterm,lovterm,latterm,lonterm,u1,v1)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u1,v1,ri1,rj1)

      latterm=(lat(nx,1)*pi)/180.
      lonterm=(lon(nx,1)*pi)/180.

      call getuv_lam(lapterm,lovterm,latterm,lonterm,u1,v1)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u1,v1,ri2,rj2)

      latterm=(lat(1,ny)*pi)/180.
      lonterm=(lon(1,ny)*pi)/180.

      call getuv_lam(lapterm,lovterm,latterm,lonterm,u1,v1)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u1,v1,ri3,rj3)

      latterm=(lat(nx,ny)*pi)/180.
      lonterm=(lon(nx,ny)*pi)/180.

      call getuv_lam(lapterm,lovterm,latterm,lonterm,u1,v1)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u1,v1,ri4,rj4)
c
      write(6,*)'Sat ri/rj corners for domain'
      write(6,*)'ri1/rj1 (SW) ',ri1,rj1
      write(6,*)'ri2/rj2 (SE) ',ri2,rj2
      write(6,*)'ri3/rj3 (NW) ',ri3,rj3
      write(6,*)'ri4/rj4 (NE) ',ri4,rj4
      write(6,*)
c
c compute ri, rj look up tables for satellite data file (for each point in
c LAPS domain). This is absolute lut (relative to the full satellite data file).
c
      do j = 1, ny
      do i = 1, nx
         latterm=(lat(i,j)*pi)/180.
         lonterm=(lon(i,j)*pi)/180.
         call getuv_lam(lapterm,lovterm,latterm,lonterm,u,v)
         call uv_ij(ny3,u_orig,v_orig,du,dv,u,v,ri(i,j),rj(i,j))
      enddo
      enddo
c
c get new i/j start/end values for this domain
c
      call get_sat_boundary(nx,ny,ny3mx,nx3mx,ri,rj,
     &linestart,lineend,elemstart,elemend,
     &rls,rle,res,ree,istatus)
      if(istatus.ne.1)then
         write(6,*)'Laps boundary outside satellite data!'
         write(6,*)'Terminating this lut generation!'
         goto 1000
      endif
c
c compute ri, rj relative look up table for the block of data surrounding
c the laps domain.
c
      do j = 1,ny
      do i = 1,nx
         rel_ri(i,j) = ri(i,j) - res + 1.
         rel_rj(i,j) = rj(i,j) - rls + 1.
      enddo
      enddo
c
c put the rel ri/rj lut into the laps domain size
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

      do i = 1,nx_l,10
      do j = 1,ny_l,10
         write(6,*)'i,j,ri,rj: ',i,j,ri_laps(i,j),rj_laps(i,j)
      enddo
      enddo

      call get_directory('static',path,n1)
      path=path(1:n1)//'/lvd/'
      n1=index(path,' ')-1

      cname=path(1:n1)//c_sat_id(isat)//'-llij-'//ct(1:nc)
      n1=index(cname,' ')-1
      table_path = cname(1:n1)//'-'//cdtype//'.lut'

      n1=index(table_path,' ')
      write(6,*)'Write lat/lon to i/j look up table'
      write(6,*)table_path(1:n1)

      call write_table (table_path,nx_l,ny_l,xlat,xlon,
     &ri_laps,rj_laps,istatus)
      if(istatus .ne. 1)then
         write(6,*)'Error writing look-up table'
         goto 1000
      endif

      r_lap(jtype,isat) = rlap
      r_lov(jtype,isat) = rlov
      r_latin(jtype,isat) = rlatin

      goto(71,72,73,72,72)indx

71       i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_la1(jtype,isat) = rla1
         r_lo1(jtype,isat) = rlo1
         r_resolution_x_vis(jtype,isat) = dx
         r_resolution_y_vis(jtype,isat) = dy
         n_pixels_vis(jtype,isat) = nx3mx
         n_lines_vis(jtype,isat)  = ny3mx

         goto 75

72       i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_la1(jtype,isat) = rla1
         r_lo1(jtype,isat) = rlo1
         r_resolution_x_ir(jtype,isat) = dx
         r_resolution_y_ir(jtype,isat) = dy
         n_pixels_ir(jtype,isat) = nx3mx
         n_lines_ir(jtype,isat)  = ny3mx

         goto 75

73       i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_la1(jtype,isat) = rla1
         r_lo1(jtype,isat) = rlo1
         r_resolution_x_wv(jtype,isat) = dx
         r_resolution_y_wv(jtype,isat) = dy
         n_pixels_wv(jtype,isat) = nx3mx
         n_lines_wv(jtype,isat)  = ny3mx

75    continue
      jstatus = 1

1000  return
      end
