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
       subroutine expand_domain(imax,jmax,xlat,xlon,kmax,lmax,lat,lon,
     &istatus)
c
      include 'trigd.inc'
      implicit none

      integer imax,jmax
      integer kmax,lmax

c these are input lat/lon
      real*4    xlat(imax,jmax)
      real*4    xlon(imax,jmax)

c these are output (expanded) lat/lon
      real*4    lat(kmax,lmax)
      real*4    lon(kmax,lmax)

      real*4    g_space_deg
      real*4    wdw_lat_ns
      real*4    wdw_lon_ns
      integer i,j
      integer ii,jj
      integer nxl,nyl
      integer istatus
c
c =================================
c
      istatus=-1
      nxl = kmax
      nyl = lmax
c
      g_space_deg = sqrt( 
     1       (  xlat(1,2) - xlat(1,1)                   )**2
     1     + ( (xlon(1,2) - xlon(1,1))*cosd(xlat(1,1))  )**2 
     1                         )
c.....       Define north-south window dimensions

      wdw_lat_ns =  g_space_deg
      wdw_lon_ns =  g_space_deg  / cosd(xlat(1,1))

c     g_space_deg = sqrt(
c    1       (  xlon(1,2) - xlon(1,1)                   )**2
c    1     + ( (xlat(1,2) - xlat(1,1))*cosd(xlat(1,1))  )**2
c    1                         )
c.....       Define east-west window dimensions

c     wdw_lat_ew =  g_space_deg
c     wdw_lon_ew =  g_space_deg  / cosd(xlat(1,1))

      do j = 1,jmax
      jj = j+1
      do i = 1,imax
         ii = i+1
         lat(ii,jj) = xlat(i,j)
         lon(ii,jj) = xlon(i,j)
      enddo
      enddo

c across the bottom and top
      do i = 2,nxl-1
         lon(i,1) = lon(i,2)
         lat(i,1) = lat(i,2)-wdw_lat_ns
         lon(i,nyl)= lon(i,nyl-1)
         lat(i,nyl)= lat(i,nyl-1)+wdw_lat_ns
         if(lat(i,nyl).gt.90.0)lat(i,nyl)=lat(i,nyl)-180.0
         if(lat(i,nyl).lt.-90.0)lat(i,nyl)=180.0+lat(i,nyl)
      enddo
c on the sides
      do j = 2,nyl-1
         lon(1,j) = lon(2,j)-wdw_lon_ns
         lat(1,j) = lat(2,j)
         lon(nxl,j)= lon(nxl-1,j)+wdw_lon_ns
         lat(nxl,j)= lat(nxl-1,j)
         if(lon(1,j).lt.-180.0)lon(1,j)=lon(1,j)+360.
         if(lon(nxl,j).gt.180.0)lon(nxl,j)=lon(nxl,j)-360.0
      enddo
c now the corners
      lat(1,1)=lat(1,2)-(lat(1,3)-lat(1,2))
      lon(1,1)=lon(2,1)-(lon(3,1)-lon(2,1))
      lat(1,nyl)=lat(2,nyl)-(lat(3,nyl)-lat(2,nyl))
      lon(1,nyl)=lon(2,nyl)-(lon(3,nyl)-lon(2,nyl))
      lat(nxl,1)=lat(nxl,2)-(lat(nxl,3)-lat(nxl,2)) 
      lon(nxl,1)=lon(nxl-1,1)-(lon(nxl-2,1)-lon(nxl-1,1))
      lat(nxl,nyl)=lat(nxl,nyl-1)-
     &(lat(nxl,nyl-2)-lat(nxl,nyl-1))
      lon(nxl,nyl)=lon(nxl-1,nyl)-
     &(lon(nxl-2,nyl)-lon(nxl-1,nyl))

      if(lat(1,1).gt.90.0)lat(1,1)=lat(1,1)-180.0
      if(lat(1,1).lt.-90.0)lat(1,1)=-180.0-lat(1,1)
      if(lat(nxl,1).gt.90.0)lat(nxl,1)=lat(nxl,1)-180.0
      if(lat(nxl,1).lt.-90.0)lat(nxl,1)=180.0-lat(nxl,1)
      if(lat(1,nyl).gt.90.0)lat(1,nyl)=lat(1,nyl)-180.0
      if(lat(1,nyl).lt.-90.0)lat(1,nyl)=-180.0-lat(1,nyl)
      if(lat(nxl,nyl).gt.90.0)lat(nxl,nyl)=lat(nxl,nyl)-180.0
      if(lat(nxl,nyl).lt.-90.0)lat(nxl,nyl)=-180.0-lat(nxl,nyl)

      if(lon(1,1).gt.180.0)lon(1,1)=lon(1,1)-360.0
      if(lon(1,1).lt.-180.0)lon(1,1)=lon(1,1)+360.0
      if(lon(nxl,1).gt.180.0)lon(nxl,1)=lon(nxl,1)-360.0
      if(lon(nxl,1).lt.-180.0)lon(nxl,1)=lon(nxl,1)+360.0
      if(lon(1,nyl).gt.180.0)lon(1,nyl)=lon(1,nyl)-360.0
      if(lon(1,nyl).lt.-180.0)lon(1,nyl)=lon(1,nyl)+360.0
      if(lon(nxl,nyl).gt.180.0)lon(nxl,nyl)=lon(nxl,nyl)-360.0
      if(lon(nxl,nyl).lt.-180.0)lon(nxl,nyl)=lon(nxl,nyl)+360.0
c  
      istatus = 1

1000  return
      end
