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
       subroutine gen_lut_ps(isat,jtype,kchl,
     &nx_l,ny_l,lat,lon,ri_laps,rj_laps,istatus)
c
c
c     Returns satellite pixel locations for each point on model (LAPS) domain
c
      implicit none

      integer   nx_l,ny_l

      real    lat(nx_l,ny_l)
      real    lon(nx_l,ny_l)   ! laps lat/lon data  -  input

      real    xlat(nx_l+2,ny_l+2)
      real    xlon(nx_l+2,ny_l+2) !expanded domain lats/lons
 
      real    dx,dy  !both in km from /public and /WFO SBN
      real    du,dv
      real    u_orig,v_orig
      real    pi
      real    u1,v1

      real    lapterm
      real    lovterm
      real    latterm
      real    lonterm
      real    r_missing_data
      real    rls,rle,res,ree
      real    polat,slat,slon
      real    rLa1,rLo1
      real    rLa2,rLo2
      real    level

      real    ri(nx_l+2,ny_l+2)
      real    rj(nx_l+2,ny_l+2)
      real    rel_ri(nx_l+2,ny_l+2)
      real    rel_rj(nx_l+2,ny_l+2)
      real    ri_laps(nx_l,ny_l)       ! Output
      real    rj_laps(nx_l,ny_l)       ! Output
      real    u,v,u0,v0,uscale,vscale
      real    usmin,usmax,vsmin,vsmax
      real    usat,vsat

      integer center_id
      integer process_id,wmo_sat_id
      double precision reftime, valtime
      character*132 earth_shape, grid_name, grid_type,
     +     origin_name, process_name, wavelength,
     +     channel_comment_, x_dim, y_dim
      character*132 asctime

      integer isat,jtype,kchl
c     integer i,j,n,nc
      integer i,j,nc,nt
      integer ii,jj
      integer indx
      integer n1,np
      integer i1,j1
      integer nx,ny,nz
      integer istatus
      integer iostatus
      integer istatus_wp
      integer nf_status
      integer nxs,nys
      integer linestart,lineend
      integer elemstart,elemend
      integer i4time_latest,ifstat
      integer nijout
      integer ierr
      integer imax,jmax,kmax
      integer kdim
      integer nav
      integer channel_fcinv
      integer record
      integer nf_fid

      logical lpoint

      character*255 path
      character*255 fname
      character*200 cname

      character*13  cftime13,cvt_i4time_wfo_fname13
c     character*6   csatid  !satellite data identifier {goes08, goes09, meteos, etc}
      character*3   cdtype !satellite data type       {'cdf', 'gvr', 'gwc', 'asc', etc}
      character*3   ct     !                         {'vis', 'ir', or 'wv'}

      real    lat0,lon0,rota
      real    sw(2),ne(2)
      common /psgrid/nxs,nys,nz,lat0,lon0,rota,sw,ne
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

      istatus = -1
      cdtype=c_sat_types(jtype,isat)

      ct=c_channel_types(kchl,jtype,isat)
      nc=index(ct,' ')-1
      if(nc.le.0)nc=3
      call lvd_file_specifier(ct,indx,istatus)
c
c retrieve the latest nav info from a file if it exists
c
      call s_len(path_to_raw_sat(kchl,jtype,isat),np)
      fname=path_to_raw_sat(kchl,jtype,isat)(1:np)
      fname=fname(1:np)//'/*'  !!!.'//ct(1:nc)

c     fname=fname(1:np)//'/20070523_1036'  !!!.'//ct(1:nc)

c test/develop this later
      call get_latest_file_time(fname,i4time_latest)
      cftime13=cvt_i4time_wfo_fname13(i4time_latest) 
      fname=path_to_raw_sat(kchl,jtype,isat)(1:np)
      fname=fname(1:np)//cftime13

      call s_len(fname,np)
C
C  Open netcdf File for reading
C
      call nfopenfile(fname,nf_fid,iostatus)
      if(iostatus.ne.0)then
         print*,'!!! Error opening file: ',TRIM(fname)
         return
      endif

      call read_FMI_sat_cdfheader(nf_fid,nav, record,nx,ny,nz,
     + Nxs,Nys,channel_fcinv, imax, jmax, kdim, kmax, Dx, Dy,
     + rLa1, rLo1, rLo2, rLa2, sLat, sLon, poLat, level, asctime,
     + channel_comment_, earth_shape, grid_name, grid_type,
     + origin_name, process_name, x_dim, y_dim, reftime, valtime,
     + istatus)

      if(istatus.ne.0)then
         print*,'!!!  Error in rd_FMI_sat_cdfheader',istatus
         print*,'fname: ',fname
         print*,'!!!  Aborting'
         print*
         return
      else
         print*,'FMI cdf Header successfully read'
         print*,'---------------------------------'
         print*,'Grid Name:    ',TRIM(grid_name)
         print*,'Grid Type:    ',TRIM(grid_type)
         print*,'Origin Name:  ',TRIM(origin_name)
         print*,'Process Name: ',TRIM(process_name)
         print*
      endif

      call nfclosefile(nf_fid,iostatus)
      if(iostatus.ne.0)then
         print*,'!!! Error closing file: ',TRIM(fname)
         return
      endif

      print*,'PS Satellite Nav Parameters'
      print*,'---------------------------------'
      print*,'dx     ',dx
      print*,'dy     ',dy
      print*,'nx     ',nxs
      print*,'ny     ',nys
      print*,'rLa1   ',rLa1
      print*,'rLo1   ',rLo1
      print*,'rLa2   ',rLa2
      print*,'rLo2   ',rLo2
      print*,'slat   ',slat
      print*,'sLon   ',slon
      print*,'polat  ',poLat
      print*

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
c code to compute polar stereo  projected satellite data for
c laps domain as specified in laps lat/lon arrays. Follows
c method in lib/latlon_to_rlaps.f by Albers.
c
      call latlon_to_uv_ps(rla1,rlo1,sLat,poLat
     &,sLon,usmin,vsmin)
      call latlon_to_uv_ps(rla2,rlo2,sLat
     &,poLat,sLon,usmax,vsmax)
      uscale = (usmax - usmin)/(float(nxs)-1.)
      vscale = (vsmax - vsmin)/(float(nys)-1.)
      u0 = usmin - uscale
      v0 = vsmin - vscale

      do j=1,ny
      do i=1,nx
         call latlon_to_uv_ps(xlat(i,j),xlon(i,j),sLat,poLat
     &,sLon,usat,vsat)
         ri(i,j)=(usat - u0)/uscale
         rj(i,j)=(vsat - v0)/vscale
      enddo
      enddo

      write(6,*)'Sat ri/rj corners for expanded domain'
      print*,'LAPS Map proj software results'
      write(6,*)'ri1/rj1 (SW) ',ri(1,1),rj(1,1)
      write(6,*)'ri2/rj2 (SE) ',ri(nx,1),rj(nx,1)
      write(6,*)'ri3/rj3 (NW) ',ri(1,ny),rj(1,ny)
      write(6,*)'ri4/rj4 (NE) ',ri(nx,ny),rj(nx,ny)
      write(6,*)

c     sw(1)=rLa1
c     sw(2)=rLo1
c     ne(1)=rLa2
c     ne(2)=rLo2
c     nz=1
c     lat0=sLat
c     lon0=sLon

c     call latlon_2_psij(nx*ny,xlat,xlon,ri,rj)

c     write(6,*)'Sat ri/rj corners for expanded domain'
c     print*,'LAPS gridconv software results'
c     write(6,*)'ri1/rj1 (SW) ',ri(1,1),rj(1,1)
c     write(6,*)'ri2/rj2 (SE) ',ri(nx,1),rj(nx,1)
c     write(6,*)'ri3/rj3 (NW) ',ri(1,ny),rj(1,ny)
c     write(6,*)'ri4/rj4 (NE) ',ri(nx,ny),rj(nx,ny)
c     write(6,*)
c
c get new i/j start/end values for this domain
c
      call get_sat_boundary(nx,ny,nx,ny,0,nxs,nys,
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
          rel_ri(i,j) = ri(i,j) - res + 1.
          rel_rj(i,j) = rj(i,j) - rls + 1.
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

      if(elemstart.le.0)elemstart=1
      if(elemend.gt.nxs)elemend=nxs
      if(linestart.le.0)linestart=1
      if(lineend.gt.nys)lineend=nys


      r_lap(jtype,isat) = sLat
      r_lov(jtype,isat) = sLon
      r_latin(jtype,isat) = poLat

      if(indx.eq.1)then

         i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_la1(jtype,isat) = rLa1
         r_lo1(jtype,isat) = rLo1
         r_resolution_x_vis(jtype,isat) = dx*1000.
         r_resolution_y_vis(jtype,isat) = dy*1000.
         n_pixels_vis(jtype,isat) = nxs
         n_lines_vis(jtype,isat)  = nys

      elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

         i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_la1(jtype,isat) = rLa1
         r_lo1(jtype,isat) = rLo1
         r_resolution_x_ir(jtype,isat) = dx*1000.
         r_resolution_y_ir(jtype,isat) = dy*1000.
         n_pixels_ir(jtype,isat) = nxs
         n_lines_ir(jtype,isat)  = nys

      elseif(indx.eq.3)then

         i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_la1(jtype,isat) = rLa1
         r_lo1(jtype,isat) = rLo1
         r_resolution_x_wv(jtype,isat) = dx*1000.
         r_resolution_y_wv(jtype,isat) = dy*1000.
         n_pixels_wv(jtype,isat) = nxs
         n_lines_wv(jtype,isat)  = nys

      endif

      istatus = 1

1000  return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_psFMI_dims(nf_fid, x, y, z, record, nav,
     &istatus)

      implicit none

      include 'netcdf.inc'

      integer nav, record, x, y, z,nf_fid, nf_vid, nf_status
      integer istatus

      istatus = 1
C
C  Fill all dimension values
C
C
C Get size of nav
C
      nf_status=NF_INQ_DIMID(nf_fid,'nav',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nav'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nav)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nav'
      endif
C
C Get size of record
C
      nf_status=NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,record)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
C
C Get size of x
C
      nf_status=NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,x)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
C
C Get size of y
C
      nf_status=NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,y)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
C
C Get size of z
C
      nf_status=NF_INQ_DIMID(nf_fid,'z',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,z)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
C
C
      subroutine read_FMI_sat_data(file_name, nav, record, x, y, z,
     +     i4time_sys, ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, istatus)


      implicit none

      include 'netcdf.inc'
      integer nav, record, x, y, z,nf_fid, nf_vid, nf_status
      integer Nxs, Nys, channel_fcinv, imax, jmax,
     +     kdim, kmax
      integer lun_out
      integer istatus
      integer i4time_latest,i4time_earliest
      real Dx, Dy, La1, Latin1, Latin2,
     +     Lo1, LoV, La2, Lo2,
     +     level
      real north,south,east,west
      real r_missing_data
      double precision reftime, valtime
      character*132 channel_comment_
      character*132 earth_shape
      character*132 x_dim
      character*132 origin_name
      character*132 grid_name
      character*132 grid_type
      character*132 process_name
      character*132 y_dim
      character*132 asctime
      character*(*) file_name

!     Declarations for 'write_lvd' call
!     integer iwmostanum(recNum)
      integer nx_l,ny_l
      integer ilaps_cycle_time
      integer i4time_sys
      logical l_closest_time, l_closest_time_i, l_in_domain
      real   lat_a(NX_L,NY_L)
      real   lon_a(NX_L,NY_L)
      real   topo_a(NX_L,NY_L)
      real   psi(nx_l,ny_l),psj(nx_l,ny_l)
      

      integer nx,ny,nz
      real    lat0,lon0,rota
      real    sw(2),ne(2)
      common /psgrid/nx,ny,nz,lat0,lon0,rota,sw,ne

      print*,'NX_L/NY_L',NX_L,NY_L
      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif
      call get_domain_perimeter(NX_L,NY_L,'nest7grid',lat_a,lon_a,
     1            topo_a,1.0,north,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error in get_domain_perimeter'
          return
      endif

      call read_FMI_sat_cdfheader(nf_fid,nav,record,x,y,z,
     +     Nxs,Nys, 
     +     channel_fcinv, imax, jmax, kdim, kmax, Dx, Dy, 
     + La1, Lo1, Lo2, La2, Latin1, LoV, Latin2, level, asctime,
     +     channel_comment_, earth_shape, grid_name, grid_type, 
     +     origin_name, process_name, x_dim, y_dim, reftime,
     +     valtime, istatus)
C
C The netcdf variables are filled - your lvd write call may go here
C
!     Initial loop through obs to get times and stanums
c     do iob = 1,recNum
c         iwmostanum(iob) = 0
c         if(abs(observationTime(iob)) .le. 1e10)then
c             i4time_ob = idint(observationTime(iob))+315619200
c             call make_fnam_lp(i4time_ob,a9time_ob_r(iob),istatus)
c         endif

c     enddo ! iob

c     c8_obstype = 

c     do iob = 1,recNum
c         height_m = r_missing_data
c         l_closest_time = .true.

c     enddo ! iob

c Ok, here's where we set up the ps nav parameters for common block
c used in latlon_2_psij

      print*,'Results of Read'
      print*,'---------------'
      print*,'x/y ',x,y
      print*,'La1/Latin1/Latin2/Lo1/LoV ',La1,Latin1,Latin2,Lo1,LoV

      sw(1)=La1
      sw(2)=Lo1
      ne(1)=61.162
      ne(2)=28.557
      nx=Nxs
      ny=Nys
      nz=1
      lat0=Latin1
      lon0=LoV

      call latlon_2_psij(nx_l*ny_l,lat_a,lon_a,psi,psj)

      print*,'psi/psj(1,1) = ',psi(1,1),psj(1,1)
      print*,'psi/psj(nx_l,1) = ',psi(nx_l,1),psj(nx_l,1)
      print*,'psi/psj(1,ny_l) = ',psi(1,ny_l),psj(1,ny_l)
      print*,'psi/psj(nx_l,ny_l) = ',psi(nx_l,ny_l),psj(nx_l,ny_l)


      print*,'Return to main'
      return
      end
C
C  Subroutine to read the file "LAPS lvd file - satellite data" 
C
      subroutine read_FMI_sat_cdfheader(nf_fid, nav, record,
     +       x, y, z, Nx, 
     +     Ny, channel_fcinv, imax, jmax, kdim, kmax, Dx, Dy, 
     + La1, Lo1, Lo2, La2, Latin1, LoV, Latin2, level, asctime,
     +     channel_comment_, earth_shape, grid_name, grid_type, 
     +     origin_name, process_name, x_dim, y_dim, reftime, valtime,
     +     istatus)
C
      implicit none

      include 'netcdf.inc'
      integer nav, record, x, y, z,nf_fid, nf_vid, nf_status
      integer Nx, Ny, channel_fcinv, imax, jmax,
     +     kdim, kmax
      integer istatus
      real Dx, Dy, La1, Latin1, Latin2,
     +     Lo1, LoV,La2,Lo2,
     +     level
c not reading the data here
c     real channel( x,  y,  z, record)

      double precision reftime, valtime
      character*132 channel_comment_
      character*132 earth_shape
      character*132 asctime
      character*132 x_dim
      character*132 origin_name
      character*132 grid_name
      character*132 grid_type
      character*132 y_dim
      character*132 process_name
      character*255 file_name

      istatus = 0

C
C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     Dx            "x grid increment"
C
      nf_status=NF_INQ_VARID(nf_fid,'Dx',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Dx'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,Dx)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Dx'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     Dy            "y grid increment"
C
      nf_status=NF_INQ_VARID(nf_fid,'Dy',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Dy'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,Dy)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Dy'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     La1           "first latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'La1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for La1'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,La1)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for La1'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     La2           "last latitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'La2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for La2'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,La2)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for La2'
       endif
      endif

C
C     Variable        NETCDF Long Name
C     Latin1        "orientation of grid"
C
      nf_status=NF_INQ_VARID(nf_fid,'Latin1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Latin1'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,Latin1)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Latin1'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     Latin2        "orientation of grid"
C
      nf_status=NF_INQ_VARID(nf_fid,'Latin2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Latin2'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,Latin2)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Latin2'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     Lo1           "first longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'Lo1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Lo1'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,Lo1)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Lo1'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     Lo2           "last longitude"
C
      nf_status=NF_INQ_VARID(nf_fid,'Lo2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Lo2'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,Lo2)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Lo2'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     LoV           "orientation of grid"
C
      nf_status=NF_INQ_VARID(nf_fid,'LoV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for LoV'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,LoV)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for LoV'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     channel       "Channel 1 (0.58 - 0.68) Daytime cloud and surface mapping"
C
c     nf_status=NF_INQ_VARID(nf_fid,'channel',nf_vid)
c     if(nf_status.ne.NF_NOERR) then
c      print *, NF_STRERROR(nf_status),' for channel'
c     else
c      nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,channel)
c      if(nf_status.ne.NF_NOERR) then
c       print *, NF_STRERROR(nf_status),' for channel'
c      endif
c     endif
C
C     Variable        NETCDF Long Name
C     level         "level of data"
C
      nf_status=NF_INQ_VARID(nf_fid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for level'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,level)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for level'
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     Nx            "number of x points"
C
      nf_status=NF_INQ_VARID(nf_fid,'Nx',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Nx'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,Nx)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Nx'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     Ny            "number of y points"
C
      nf_status=NF_INQ_VARID(nf_fid,'Ny',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for Ny'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,Ny)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for Ny'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     channel_fcinv 
C
      nf_status=NF_INQ_VARID(nf_fid,'channel_fcinv',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for channel_fcinv'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,channel_fcinv)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for channel_fcinv'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     imax          
C
      nf_status=NF_INQ_VARID(nf_fid,'imax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for imax'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,imax)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for imax'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     jmax          
C
      nf_status=NF_INQ_VARID(nf_fid,'jmax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for jmax'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,jmax)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for jmax'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     kdim          
C
      nf_status=NF_INQ_VARID(nf_fid,'kdim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for kdim'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,kdim)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for kdim'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     kmax          
C
      nf_status=NF_INQ_VARID(nf_fid,'kmax',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for kmax'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,kmax)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for kmax'
       endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     reftime       "reference time"
C
      nf_status=NF_INQ_VARID(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for reftime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reftime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for reftime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     valtime       "valid time"
C
      nf_status=NF_INQ_VARID(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for valtime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,valtime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for valtime'
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     asctime       
C
      nf_status=NF_INQ_VARID(nf_fid,'asctime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for asctime'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,asctime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for asctime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     channel_comment_
C
      nf_status=NF_INQ_VARID(nf_fid,'channel_comment_',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for channel_comment_'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,channel_comment_)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for channel_comment_'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     earth_shape   
C
      nf_status=NF_INQ_VARID(nf_fid,'earth_shape',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for earth_shape'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,earth_shape)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for earth_shape'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     grid_name     
C
      nf_status=NF_INQ_VARID(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for grid_name'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_name)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for grid_name'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     grid_type     "GRIB-1 grid type"
C
      nf_status=NF_INQ_VARID(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for grid_type'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_type)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for grid_type'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     origin_name   
C
      nf_status=NF_INQ_VARID(nf_fid,'origin_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for origin_name'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,origin_name)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for origin_name'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     process_name  
C
      nf_status=NF_INQ_VARID(nf_fid,'process_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for process_name'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,process_name)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for process_name'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     x_dim         "longitude dimension"
C
      nf_status=NF_INQ_VARID(nf_fid,'x_dim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for x_dim'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,x_dim)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for x_dim'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     y_dim         "latitude dimension"
C
      nf_status=NF_INQ_VARID(nf_fid,'y_dim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for y_dim'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,y_dim)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for y_dim'
       endif
      endif

      return
      end
c ----------------------------------------------------------
      subroutine nfopenfile(filename, nf_fid, istatus)

      implicit none

      include 'netcdf.inc'

      character*(*) filename
      integer nf_fid
      integer nf_status
      integer istatus

      istatus = 0
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),TRIM(filename)
        istatus=1
        return
      endif
      return
      end
c ----------------------------------------------------------
      subroutine nfclosefile(nf_fid, istatus)

      implicit none

      include 'netcdf.inc'

      integer nf_fid
      integer nf_status
      integer istatus

      istatus = 0
C
C  Open netcdf File for reading
C
      nf_status=NF_CLOSE(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print*,'nf_close'
        istatus=1
        return
      endif

      return
      end
