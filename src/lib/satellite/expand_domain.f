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
       subroutine expand_domain(imax,jmax,lat,lon,kmax,lmax,xlat,xlon,
     &istatus)
c
      include 'trigd.inc'
      implicit none

      integer imax,jmax
      integer kmax,lmax

c these are input lat/lon
      real    lat(imax,jmax)
      real    lon(imax,jmax)

c these are output (expanded) lat/lon
      real    xlat(kmax,lmax)
      real    xlon(kmax,lmax)

      real    g_space_mx,g_space_my,g_space_deg_x,g_space_deg_y
      real    wdw_lat_ns
      real    wdw_lon_ns
      real    wdw_lon_ns_test
      integer i,j
      integer ii,jj
      integer nxl,nyl
      integer icen,jcen
      integer istatus
c
c =================================
c
      istatus=-1
      nxl = kmax
      nyl = lmax
      icen = imax/2+1
      jcen = jmax/2+1

      call get_grid_spacing_actual_xy(lat(icen,jcen),lon(icen,jcen)
     &                               ,g_space_mx,g_space_my,istatus)
      g_space_deg_x = g_space_mx/111100.
      g_space_deg_y = g_space_my/111100.
c
c     g_space_deg = sqrt( 
c    1       (  lat(1,2) - lat(1,1)                   )**2
c    1     + ( (lon(1,2) - lon(1,1))*cosd(lat(1,1))  )**2 
c    1                         )
c.....       Define north-south window dimensions

      wdw_lat_ns =  g_space_deg_y
      wdw_lon_ns =  wdw_lat_ns           !g_space_deg  / cosd(xlat(1,1))

c     wdw_lon_ns_test = g_space_m/111100./cosd(lat(icen,jcen))

c     g_space_deg = sqrt(
c    1       (  lon(1,2) - lon(1,1)                   )**2
c    1     + ( (lat(1,2) - lat(1,1))*cosd(lat(1,1))  )**2
c    1                         )
c.....       Define east-west window dimensions

c     wdw_lat_ew =  g_space_deg
c     wdw_lon_ew =  g_space_deg  / cosd(lat(1,1))

      do j = 1,jmax
      jj = j+1
      do i = 1,imax
         ii = i+1
         xlat(ii,jj) = lat(i,j)
         xlon(ii,jj) = lon(i,j)
      enddo
      enddo

c across the bottom and top
      do i = 2,nxl-1
         xlon(i,1) = xlon(i,2)-(xlon(i,3)-xlon(i,2))
         xlat(i,1) = xlat(i,2)-wdw_lat_ns
         xlon(i,nyl)= xlon(i,nyl-1)-(xlon(i,nyl-2)-xlon(i,nyl-1))
         xlat(i,nyl)= xlat(i,nyl-1)+wdw_lat_ns
         if(xlat(i,nyl).gt.90.0)xlat(i,nyl)=xlat(i,nyl)-180.0
         if(xlat(i,nyl).lt.-90.0)xlat(i,nyl)=180.0+xlat(i,nyl)
      enddo
c on the sides
      do j = 2,nyl-1
         xlon(1,j) = xlon(2,j)-wdw_lon_ns
         xlat(1,j) = xlat(2,j)-(xlat(3,j)-xlat(2,j))
         xlon(nxl,j)= xlon(nxl-1,j)+wdw_lon_ns
         xlat(nxl,j)= xlat(nxl-1,j)-(xlat(nxl-2,j)-xlat(nxl-1,j))
         if(xlon(1,j).lt.-180.0)xlon(1,j)=xlon(1,j)+360.
         if(xlon(nxl,j).gt.180.0)xlon(nxl,j)=xlon(nxl,j)-360.0
      enddo
c now the corners
      xlat(1,1)=xlat(1,2)-wdw_lat_ns                      !(xlat(1,3)-xlat(1,2))
      xlon(1,1)=xlon(2,1)-wdw_lon_ns                      !(xlon(3,1)-xlon(2,1))
      xlat(1,nyl)=xlat(1,nyl-1)+wdw_lat_ns                !(xlat(3,nyl)-xlat(2,nyl))
      xlon(1,nyl)=xlon(2,nyl)-wdw_lon_ns                  !(xlon(3,nyl)-xlon(2,nyl))
      xlat(nxl,1)=xlat(nxl,2)-wdw_lat_ns                  !(xlat(nxl,3)-xlat(nxl,2)) 
      xlon(nxl,1)=xlon(nxl-1,1)+wdw_lon_ns                !-(xlon(nxl-2,1)-xlon(nxl-1,1))
      xlat(nxl,nyl)=xlat(nxl,nyl-1)+wdw_lat_ns
                                                          !&(xlat(nxl,nyl-2)-xlat(nxl,nyl-1))
      xlon(nxl,nyl)=xlon(nxl-1,nyl)+wdw_lon_ns
                                                          !&(xlon(nxl-2,nyl)-xlon(nxl-1,nyl))

      if(xlat(1,1).gt.90.0)xlat(1,1)=xlat(1,1)-180.0
      if(xlat(1,1).lt.-90.0)xlat(1,1)=-180.0-xlat(1,1)
      if(xlat(nxl,1).gt.90.0)xlat(nxl,1)=xlat(nxl,1)-180.0
      if(xlat(nxl,1).lt.-90.0)xlat(nxl,1)=180.0-xlat(nxl,1)
      if(xlat(1,nyl).gt.90.0)xlat(1,nyl)=xlat(1,nyl)-180.0
      if(xlat(1,nyl).lt.-90.0)xlat(1,nyl)=-180.0-xlat(1,nyl)
      if(xlat(nxl,nyl).gt.90.0)xlat(nxl,nyl)=xlat(nxl,nyl)-180.0
      if(xlat(nxl,nyl).lt.-90.0)xlat(nxl,nyl)=-180.0-xlat(nxl,nyl)

      if(xlon(1,1).gt.180.0)xlon(1,1)=xlon(1,1)-360.0
      if(xlon(1,1).lt.-180.0)xlon(1,1)=xlon(1,1)+360.0
      if(xlon(nxl,1).gt.180.0)xlon(nxl,1)=xlon(nxl,1)-360.0
      if(xlon(nxl,1).lt.-180.0)xlon(nxl,1)=xlon(nxl,1)+360.0
      if(xlon(1,nyl).gt.180.0)xlon(1,nyl)=xlon(1,nyl)-360.0
      if(xlon(1,nyl).lt.-180.0)xlon(1,nyl)=xlon(1,nyl)+360.0
      if(xlon(nxl,nyl).gt.180.0)xlon(nxl,nyl)=xlon(nxl,nyl)-360.0
      if(xlon(nxl,nyl).lt.-180.0)xlon(nxl,nyl)=xlon(nxl,nyl)+360.0
c  
      istatus = 1

1000  return
      end
