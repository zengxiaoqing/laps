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
       subroutine gen_llij_lut_v01(cid,radlat,radlon,
     &nx,ny,lat,lon,la1,lo1,nx3,ny3,dx,dy,ri,rj,jstatus)
c
c
      implicit none

      integer*4 nx,ny
      integer*4 nxv01,nyv01

      real*4    lat(nx,ny)
      real*4    lon(nx,ny)   ! laps lat/lon data  -  input
 
      real*4    la1, lo1, lap
      real*4    dx,dy
      real*4    du,dv
      real*4    u_orig,v_orig
      real*4    latin,lov
      real*4    pi

      real*4    ri1,rj1
      real*4    ri2,rj2
      real*4    ri3,rj3
      real*4    ri4,rj4
      real*4    radri,radrj

      real*4    lapterm
      real*4    lovterm
      real*4    latterm
      real*4    lonterm
      real*4    dxterm
      real*4    dyterm
      real*4    grid_spacing_deg
      real*4    radlat
      real*4    radlon
      real*4    rls,rle,res,ree

      real*4    ri(nx,ny)
      real*4    rj(nx,ny)
      real*4    u
      real*4    v

      integer*4 i,j,n
      integer*4 ii,jj
      integer*4 n1,n2
      integer*4 lend
      integer*4 istatus
      integer*4 jstatus
      integer*4 nx3,ny3

      logical   lwrite
      data      lwrite/.false./

      character*200 table_path
      character*200 file
      character*255 path
      character*200 cname
      character*150 cdir
      character*4   cid     !Radar id {KFTG, KGLD, KCYS, etc}
c
c =================================
c
      jstatus = -1
      pi = acos(-1.0)
c
c     call get_directory('static',cdir,lend)
c     path=cdir(1:lend)//'lvd/'
c     n1=index(path,' ')-1
c     cname=cid//'_'//cdtyp//'_'//cchtyp
c     n2=index(cname,' ')-1
c     file=path(1:n1)//cname(1:n2)//'.parms'
c     n=index(file,' ')

c     call get_satnav_parms(cid,cdtyp,cchtyp,'install   ',
c    &dx,dy,nx3,ny3,la1,lo1,latin,lov,lap,
c    &elemstart,elemend,linestart,lineend,
c    &istatus)
c     if(istatus.ne.1)then
c        write(6,*)'Error getting sat parms - gen_cdf_lut'
c        goto 1000
c     endif
c
c fsl-conus fixed parameters
c ----------------------
      latin=25.00000
      Lov=-95.00000
      lap=25.00000

      write(6,*)'Parameters for ',cid
      write(6,*)'Radar lat/lon  ',radlat,radlon
      write(6,*)'dx     ',dx
      write(6,*)'dy     ',dy
      write(6,*)'nx3    ',nx3
      write(6,*)'ny3    ',ny3
      write(6,*)'la1    ',la1
      write(6,*)'lo1    ',lo1
      write(6,*)'latin  ',latin
      write(6,*)'Lov    ',Lov
      write(6,*)'lap    ',lap
c
      if(la1.ne.0.0 .and. lo1.ne.0.0 .and.
     &    dx.ne.0.0 .and.  dy.ne.0.0 .and.
     &    radlat.ne.-999.99          .and.
     &    radlon.ne.-999.99)then

      nyv01=ny3
      dxterm = dx/1000.
      dyterm = dy/1000.
      call getdudv_lam(lov,lap,dxterm,dyterm,
     &    la1,lo1,du,dv,u_orig,v_orig)
c
c compute corner point ri/rj pairs used to determine # of i and j
c
      write(6,*)'Compute corner points'

      lapterm=(lap*pi)/180.
      lovterm=(lov*pi)/180.
      latterm=(lat(1,1)*pi)/180.
      lonterm=(lon(1,1)*pi)/180.

      call getuv_lam (lapterm,lovterm,latterm,lonterm,u,v)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u,v,ri1,rj1)
      rj1=nyv01-rj1

      latterm=(lat(nx,1)*pi)/180.
      lonterm=(lon(nx,1)*pi)/180.

      call getuv_lam (lapterm,lovterm,latterm,lonterm,u,v)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u,v,ri2,rj2)
      rj2=nyv01-rj2

      latterm=(lat(1,ny)*pi)/180.
      lonterm=(lon(1,ny)*pi)/180.

      call getuv_lam (lapterm,lovterm,latterm,lonterm,u,v)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u,v,ri3,rj3)
      rj3=nyv01-rj3

      latterm=(lat(nx,ny)*pi)/180.
      lonterm=(lon(nx,ny)*pi)/180.

      call getuv_lam (lapterm,lovterm,latterm,lonterm,u,v)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u,v,ri4,rj4)
      rj4=nyv01-rj4

      latterm=(radlat*pi)/180.
      lonterm=(radlon*pi)/180.

      call getuv_lam (lapterm,lovterm,latterm,lonterm,u,v)
      call uv_ij (ny3,u_orig,v_orig,du,dv,u,v,radri,radrj)
c     radrj=nyv01-radrj
c
      write(6,*)'Radar ri/rj corners for domain'
      write(6,*)'ri1/rj1 (SW) ',ri1,rj1
      write(6,*)'ri2/rj2 (SE) ',ri2,rj2
      write(6,*)'ri3/rj3 (NW) ',ri3,rj3
      write(6,*)'ri4/rj4 (NE) ',ri4,rj4
      write(6,*)'Radar ri/rj  ',radri,radrj
c
c first get uv in lambert grid
c
      do j=1,ny
      do i=1,nx
         latterm=(lat(i,j)*pi)/180.
         lonterm=(lon(i,j)*pi)/180.
         call getuv_lam(lapterm,lovterm,latterm,lonterm,u,v)
         call uv_ij(ny3,u_orig,v_orig,du,dv,u,v,ri(i,j),rj(i,j))
         rj(i,j)=nyv01-rj(i,j)
      enddo
      enddo
c
      if(lwrite)then
      do i=1,nx,10
      do j=1,ny,10
         write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)
      enddo
      enddo
      endif

      else
        print*,'Not able to compute lut given those parms'
      endif

c     if(.false)then
c     cname=path(1:n1)//cid//'-llij-'//cchtyp
c     n2=index(cname,' ')-1
c     table_path = cname(1:n2)//'-'//cdtyp//'.lut'
c     n2=index(table_path,' ')
c     write(6,*)'Write lat/lon to i/j look up table'
c     write(6,*)table_path(1:n2)

c     call write_table (table_path,nx,ny,lat,lon,ri,rj,istatus)
c     if(istatus .ne. 1)then
c        write(6,*)'Error writing look-up table'
c        goto 1000
c     endif

c     endif
 
      jstatus = 1

1000  return
      end
