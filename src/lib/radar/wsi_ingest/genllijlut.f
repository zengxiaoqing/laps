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
      subroutine gen_llij_lut_wsi(irad,imax,jmax,lat,lon
     +    ,istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Routine reads netCDF nowrad (WSI) data using subroutine read_wsi_cdf.
c Nowrad data is then remapped to LAPS domain given lat/lon of domain.
c Routine automatically moves the boundary for any domain with the nowrad
c confines.
c
       implicit none

       integer*4 imax,jmax
       real*4 lat(imax,jmax)
       real*4 lon(imax,jmax)
       real*4 ri(imax,jmax)
       real*4 rj(imax,jmax)

       real*4 dx,dy
       real*4 rla1,rlo1
       real*4 rla2,rlo2
       real*4 rlatin,rlap
       real*4 rlat,rlon
       real*4 rlatc,rlonc
       real*4 dlat,dlon
       real*4 nw(2),se(2)
       real*4 r_missing_data

       integer i,j
       integer n,n1
       integer nx,ny,nz
       integer irad
       integer istatus
       integer nlines
       integer nelems
       integer i_out
       integer j_out

       character table_path*255
       character path*100
       character file*255

       common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c ***************************************************************************
c
      istatus=1
      call get_wsi_parms_vrc(irad,nlines,nelems,
     +dx,dy,rla1,rlo1,rla2,rlo2,rlat,rlon,rlatin,
     +istatus)

      write(6,*)'Parameters from vrc_nl namelist'
      write(6,*)'dx     ',dx
      write(6,*)'dy     ',dy
      write(6,*)'nelems ',nelems
      write(6,*)'nlines ',nlines
      write(6,*)'rla1   ',rla1
      write(6,*)'rlo1   ',rlo1
      write(6,*)'rla2   ',rla2
      write(6,*)'rlo2   ',rlo2
      write(6,*)'rlon   ',rlon
      write(6,*)'rlat   ',rlat
      write(6,*)'rlatin ',rlatin
c
      dlon = 0.019119
      dlat = 0.017966
      rlatc = rlat
      rlonc = rlon
      nx = nelems
      ny = nlines
      nw(1)=rla1
      nw(2)=rlo1
      se(1)=rla2
      se(2)=rlo2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Build ri/rj look up for laps domain
c
      call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)



      write(6,*)'ri/rj corners for domain: '
      write(6,*)'--------------------------'
      write(6,*)'(SW) ',ri(1,1),rj(1,1)
      write(6,*)'(SE) ',ri(imax,1),rj(imax,1)
      write(6,*)'(NW) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'(NE) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)
      write(6,*)'Check for domain coverage'
c
c check for laps domain outside data domain
c
       call get_r_missing_data(r_missing_data,istatus)
       i_out=0
       j_out=0
       do j=1,jmax
         do i=1,imax
           if(ri(i,j).gt.nx)then
              i_out=i_out+1
              ri(i,j)=r_missing_data
           endif
           if(ri(i,j).le.0.)then
              i_out=i_out+1
              ri(i,j)=r_missing_data
           endif
           if(rj(i,j).gt.ny)then
              j_out=j_out+1
              rj(i,j)=r_missing_data
           endif
           if(rj(i,j).le.0.0)then
              j_out=j_out+1
              rj(i,j)=r_missing_data
           endif
         enddo
       enddo

       if(i_out.gt.0)then
          print*,'Warning! LAPS domain outside of data domain'
          print*,'i out =',i_out
       endif
       if(j_out.gt.0)then
          print*,'Warning! LAPS domain outside of data domain'
          print*,'j out =',j_out
       endif
       if(i_out.gt.0.or.j_out.gt.0)then
          print*,'Check yor laps domain (nest7grid.parms) settings'
          print*,'if you do not want partial coverage'
       endif
c
c output
c
      write(6,*)'ri/rj corners for domain after check: '
      write(*,*)'--------------------------------------'
      write(6,*)'(SW) ',ri(1,1),rj(1,1)
      write(6,*)'(SE) ',ri(imax,1),rj(imax,1)
      write(6,*)'(NW) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'(NE) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)

      if(.false.)then
        do i = 1,imax,10
        do j = 1,jmax,10
           write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)
        enddo
        enddo
      endif

      call get_directory('static',table_path,n1)
      table_path = table_path(1:n1)//'vrc/wsi_llij_lut_wsi.lut'
      write(6,*)'Write lat/lon to i/j look up table'
      n1=index(table_path,' ')
      write(6,*)table_path(1:n1)

      call write_table (table_path,imax,jmax,lat,lon,ri,rj,istatus)
      if(istatus .ne. 1)then
         write(6,*)'Error writing look-up table'
         goto 900
      endif

      goto 16

900   write(6,*)'Error writting table ',file(1:n)
      goto 16

901   write(6,*)'Error reading parm file ',file(1:n)

16    write(6,*)'Finished in get_llij_lut_wsi'
      return
      end
c
c *******************************************************************
c *******************************************************************
c
      subroutine gen_vrc_wfo_cdf_lut(irad,imax,jmax,lat,lon,istatus)
c
      implicit none

      integer*4 imax,jmax

      real*4    lat(imax,jmax)
      real*4    lon(imax,jmax)
      real*4    rla1,rlo1
      real*4    rla2,rlo2
      real*4    rlat,rlon
      real*4    dx,dy
      real*4    du,dv
      real*4    u_orig,v_orig
      real*4    rlatin
      real*4    pi
      real*4    u,v

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

      integer   i,j,n,nn
      integer   n1,n2
      integer   nlines,nelems
      integer   irad
      integer   istatus

      character*200 table_path
      character*200 file
      character*255 path
      character*200 cname
      character*3   ctype
c
c =================================
      pi = acos(-1.0)
c

      call get_wsi_parms_vrc(irad,nlines,nelems,
     +dx,dy,rla1,rlo1,rla2,rlo2,rlat,rlon,rlatin,
     +istatus)

      write(6,*)'Parameters from vrc_nl namelist'
      write(6,*)'dx     ',dx
      write(6,*)'dy     ',dy
      write(6,*)'nelems ',nelems
      write(6,*)'nlines ',nlines
      write(6,*)'rla1   ',rla1
      write(6,*)'rlo1   ',rlo1
      write(6,*)'rla2   ',rla2
      write(6,*)'rlo2   ',rlo2
      write(6,*)'rlon   ',rlon
      write(6,*)'rlat   ',rlat
      write(6,*)'rlatin ',rlatin

c compute lat lons  fsl conus grid (or read fm disk)
c determine ri/rj pair for the four laps domain corners.
c compute i,j  in lambert grid for desired points.
c get delta u and delta v in the lambert grid

         dxterm = dx/1000.
         dyterm = dy/1000.
         call getdudv_lam(rlon,rlat,dxterm,dyterm,
     &       rla1,rlo1,du,dv,u_orig,v_orig)
c
c compute corner point ri/rj pairs used to determine # of i and j
c
         write(6,*)'Compute corner points'

         lapterm=(rlat*pi)/180.
         lovterm=(rlon*pi)/180.
         latterm=(lat(1,1)*pi)/180.
         lonterm=(lon(1,1)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u,v)
            call uv_ij (nlines,u_orig,v_orig,du,dv,
     &       u,v,ri1,rj1)

         latterm=(lat(imax,1)*pi)/180.
         lonterm=(lon(imax,1)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u,v)
            call uv_ij (nlines,u_orig,v_orig,du,dv,
     &       u,v,ri2,rj2)

         latterm=(lat(1,jmax)*pi)/180.
         lonterm=(lon(1,jmax)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u,v)
            call uv_ij (nlines,u_orig,v_orig,du,dv,
     &       u,v,ri3,rj3)

         latterm=(lat(imax,jmax)*pi)/180.
         lonterm=(lon(imax,jmax)*pi)/180.

         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u,v)
            call uv_ij (nlines,u_orig,v_orig,du,dv,
     &       u,v,ri4,rj4)
c
         write(6,*)'WFO/WSI ri/rj corners for domain: '
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
     &       latterm,lonterm,u,v)
            call uv_ij (nlines,u_orig,v_orig,du,dv,
     &       u,v,ri(i,j),rj(i,j))

         enddo
         enddo
c
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
