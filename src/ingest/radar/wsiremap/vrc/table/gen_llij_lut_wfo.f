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
       subroutine gen_vrc_wfo_cdf_lut(ich,imax,jmax,lat,lon,istatus)
c
      implicit none

ccc      include 'lapsparms.for'

      integer*4 imax,jmax

      real*4    lat(imax,jmax)
      real*4    lon(imax,jmax)
      real*4    la1, lo1, lap
      real*4    dx,dy
      real*4    du,dv
      real*4    u_orig,v_orig
      real*4    latin,lov
      real*4    pi
      real*4    u1,v1

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

      real*4    ri(imax,jmax)
      real*4    rj(imax,jmax)
      real*4    u(imax,jmax)
      real*4    v(imax,jmax)

      integer*4 i,j,n,nn
      integer*4 n1,n2
      integer*4 ich
      integer*4 isat
      integer*4 istatus
      integer*4 nx3,ny3

      character*200 table_path
      character*200 file
      character*255 path
      character*200 cname
      character*3   ctype
c
c =================================
      pi = acos(-1.0)
c
c parameters needed to generate the liberal goes8 lookup table
c
c     nn=index(c_channel_type(ich,isat),' ')-1
c     if(nn.lt.0)nn=3
c     ctype=c_channel_type(ich,isat)
      call get_directory('static',path,n1)
      path=path(1:n1)//'vrc/ '
      n1=index(path,' ')-1
      cname='wsi_cdf_lut_wfo'
      n2=index(cname,' ')-1
      file=path(1:n1)//cname(1:n2)//'.parms'
      n=index(file,' ')-1

      open(22,file=file(1:n),
     &     form='formatted',status='old',err=901)
      read(22,*)
      read(22,*)
      read(22,*)
      read(22,*)
      read(22,*)
      read(22,50)dx
      read(22,50)dy
      read(22,51)nx3
      read(22,51)ny3
      read(22,50)la1
      read(22,50)lo1
      read(22,50)latin
      read(22,50)Lov
      read(22,50)lap
50    format(f10.5)
51    format(i4)
      close(22)

      write(6,*)'Parameters from ',file(1:n)
      write(6,*)'dx    ',dx
      write(6,*)'dy    ',dy
      write(6,*)'nx3   ',nx3
      write(6,*)'ny3   ',ny3
      write(6,*)'la1   ',la1
      write(6,*)'lo1   ',lo1
      write(6,*)'latin ',latin
      write(6,*)'Lov   ',Lov
      write(6,*)'lap   ',lap
  
c compute lat lons  fsl conus grid (or read fm disk)
c determine ri/rj pair for the four laps domain corners.
c compute i,j  in lambert grid for desired points.
c get delta u and delta v in the lambert grid

         dxterm = dx/1000.
         dyterm = dy/1000.
         call getdudv_lam(lov,lap,dxterm,dyterm,
     &       la1,lo1,du,dv,u_orig,v_orig)
c
c compute corner point ri/rj pairs used to determine # of i and j
c
         write(6,*)'Compute corner points'

         lapterm=(lap*pi)/180.
         lovterm=(lov*pi)/180.
         latterm=(lat(1,1)*pi)/180.
         lonterm=(lon(1,1)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u1,v1)
            call uv_ij (ny3,u_orig,v_orig,du,dv,
     &       u1,v1,ri1,rj1)

         latterm=(lat(imax,1)*pi)/180.
         lonterm=(lon(imax,1)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u1,v1)
            call uv_ij (ny3,u_orig,v_orig,du,dv,
     &       u1,v1,ri2,rj2)

         latterm=(lat(1,jmax)*pi)/180.
         lonterm=(lon(1,jmax)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u1,v1)
            call uv_ij (ny3,u_orig,v_orig,du,dv,
     &       u1,v1,ri3,rj3)

         latterm=(lat(imax,jmax)*pi)/180.
         lonterm=(lon(imax,jmax)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u1,v1)
            call uv_ij (ny3,u_orig,v_orig,du,dv,
     &       u1,v1,ri4,rj4)
c
         write(6,*)'Sat ri/rj corners for domain: '
ccc     &laps_domain_file
         write(6,*)'ri1/rj1 (SW) ',ri1,rj1
         write(6,*)'ri2/rj2 (SE) ',ri2,rj2
         write(6,*)'ri3/rj3 (NW) ',ri3,rj3
         write(6,*)'ri4/rj4 (NE) ',ri4,rj4
         write(6,*)
c
c first get uv in lambert grid
c
         do j = 1, jmax
         do i = 1, imax

            latterm=(lat(i,j)*pi)/180.
            lonterm=(lon(i,j)*pi)/180.

            call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u(i,j),v(i,j) )

         enddo
         enddo
c
c compute ri, rj look uptables for ispan file (for each point in 
c LAPS domain)
c
         do i = 1, imax
         do j = 1, jmax

            call uv_ij (ny3,u_orig,v_orig,du,dv,
     &       u(i,j),v(i,j),ri(i,j),rj(i,j))

         enddo
         enddo

         do i = 1,imax,10
         do j = 1,jmax,10

            write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)

         enddo
         enddo
         call get_directory('static',table_path,n1)
         table_path = table_path(1:n1)//'vrc/wsi_llij_lut_wfo.lut'
         write(6,*)'Write lat/lon to i/j look up table'
         n1=index(table_path,' ')
         write(6,*)table_path(1:n1)

         call write_table (table_path,imax,jmax,lat,lon,ri,rj,istatus)
         if(istatus .ne. 1)then
            write(6,*)'Error writing look-up table'
            goto 900
         endif
c
         goto 900

901      write(6,*)'Error openning file ',file(1:n)

900      return
         end
