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
       character*3 c_raddat_type
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

       integer i,j
       integer n,n1
       integer nx,ny,nz
       integer irad
       integer istatus
       integer nlines
       integer nelems

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
c
c output
c
      write(6,*)'WFO/WSI ri/rj corners for domain: '
      write(6,*)'(SW) ',ri(1,1),rj(1,1)
      write(6,*)'(SE) ',ri(imax,1),rj(imax,1)
      write(6,*)'(NW) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'(NE) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)

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

      goto 16

900   write(6,*)'Error writting table ',file(1:n)
      goto 16

901   write(6,*)'Error reading parm file ',file(1:n)

16    write(6,*)'Finished in get_llij_lut_wsi'
      return
      end
