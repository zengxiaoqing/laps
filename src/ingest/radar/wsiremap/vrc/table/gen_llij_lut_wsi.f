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
     +    ,c_raddat_type,istatus)
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
       character*3 c_raddat_type
       real*4 lat(imax,jmax)
       real*4 lon(imax,jmax)
       real*4 ri(imax,jmax)
       real*4 rj(imax,jmax)

       real*4 pi
       real*4 rdtodg
       real*4 dgtord
       real*4 dgtokm
       real*4 dx,dy
       real*4 du,dv
       real*4 rla1,rlo1
       real*4 rla2,rlo2
       real*4 rlatin,rlap
       real*4 rlat,rlon
       real*4 ri1,ri2,ri3,ri4
       real*4 rj1,rj2,rj3,rj4
       real*4 dxterm,dyterm
       real*4 lat0,lon0
       real*4 dlat,dlon

       real*4    fraclat
       real*4    fraclon
       real*4    rlat_diff_deg
       real*4    rlon_diff_deg
       real*4    iline,jline
       real*4    idiff,jdiff

       integer i,j
       integer k,l
       integer n,n1,n2
       integer nx,ny,nz
       integer irad
       integer istart,jstart
       integer iend,jend
       integer ishow_timer
       integer init_timer
       integer itstatus
       integer istatus
       integer nlines
       integer nelems

       logical   found_line
       logical   found_elem

       character path*100
       character cname*100
       character file*255

       integer lines
       integer nelements

       common /cegrid/nx,ny,nz,lat0,lon0,dlat,dlon
c
c ***************************************************************************
c
      pi=acos(-1.)
      rdtodg=180.0/pi
      dgtord=1./rdtodg
      dgtokm=111.1
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
      lat0 = rla1
      lon0 = rlo1
      nx = nelems
      ny = nlines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Build ri/rj look up for laps domain
c
       call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)
c
c output
c
      write(6,*)'WFO/WSI ri/rj corners for domain: '
      write(6,*)'ri1/rj1 (SW) ',ri(1,1),rj(1,1)
      write(6,*)'ri2/rj2 (SE) ',ri(imax,1),rj(imax,1)
      write(6,*)'ri3/rj3 (NW) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'ri4/rj4 (NE) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)

       do i = 1,imax,10
       do j = 1,jmax,10

          write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)

       enddo
       enddo

       cname='wsi_llij_lut_'//c_raddat_type
       n2=index(cname,' ')-1
       file = path(1:n1)//cname(1:n2)//'.lut'
       n=index(file,' ')
       write(6,*)'Write lat/lon to i/j look up table'
       write(6,*)file(1:n)

       call write_table (file,imax,jmax,lat,lon,ri,rj,istatus)
       if(istatus .ne. 1)then
          write(6,*)'Error writing look-up table'
          goto 900
       endif

       goto 16

900    write(6,*)'Error writting table ',file(1:n)
       goto 16

901    write(6,*)'Error reading parm file ',file(1:n)

16     write(6,*)'Finished in get_llij_lut_polar'
       return
       end
