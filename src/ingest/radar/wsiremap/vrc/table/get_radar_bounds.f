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
        subroutine get_radar_bounds(lat,lon,nx,ny,nlines,nelements,
     &wsi_lat,wsi_lon,istart,jstart,iend,jend)
cc
c
c program builds static file for WSI radar processor. Requires lat/lon static
c file for laps. Output is i and j starting and ending points for processing the
c WSI radar file.  This program must be run before installing laps
c
c This program should run out of ../exe
c

      Implicit NONE

      include 'lapsparms.for'

      integer*4 lines,elements
      integer*4 nx,ny

      real*4 lat(nx,ny)
      real*4 lon(nx,ny)

      integer*4 istart,iend,jstart,jend
      integer*4 nlines,nelements

      real*4 dx,dy
      real*4 la1,lo1,Lov
      real*4 latin,lap

      real*4 rdtodg
      real*4 dlat_deg
      real*4 dlon_deg
      real*4 e_radius
      real*4 wsi_lat(nelements,nlines)
      real*4 wsi_lon(nelements,nlines)
      real*4 ri,rj
      real*4 pi
      integer*4 i,j
      integer*4 ii,jj
      integer*4 n,n1,n2
      integer*4 istatus

      character path*100
      character cname*100
      character file*255
      integer len 
      data e_radius/6367000./
c
c  BEGIN
c
      call get_directory('static',path,len)
      path = path(1:len)//'vrc/'

c      path='../static/vrc/'
      n1=index(path,' ')-1
      cname='wsi_cdf_lut_wsi'
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
      read(22,50)dx           !meters
      read(22,50)dy           !meters
      read(22,51)elements     !integer
      read(22,51)lines        !integer
      read(22,50)la1          !degrees
      read(22,50)lo1          !degrees
      read(22,50)latin        !not used
      read(22,50)Lov          !degrees
      read(22,50)lap          !not used
50    format(f10.5)
51    format(i4)

      close(22)

      pi=acos(-1.)
      rdtodg=180.0/pi

      dlat_deg=tan(dy/e_radius)*rdtodg
      dlon_deg=tan(dx/e_radius)*rdtodg

      istart=elements
      jstart=lines
      iend  =0
      jend  =0
c
c some caution here. For the wsi raw data (polar stereo) grid the first point
c is NW corner and thus wsi_lat is decremented.  If the
c first point is the SW corner, then wsi_lat computation must be incremented.
c
      do j =1,lines
         do i=1,elements

         wsi_lat(i,j)=la1-(dlat_deg*(j-1))
         wsi_lon(i,j)=lo1+(dlon_deg*(i-1))

         enddo
      enddo

      do j = 1,lines
         do i = 1,elements

            call latlon_to_rlapsgrid(wsi_lat(i,j),wsi_lon(i,j),
     &                                lat,lon,nx_l,ny_l,ri,rj,
     &                                istatus)
            ii = nint(ri)
            jj = nint(rj)

            if(ii .gt. 0 .and. ii .le. nx_l)then
               if(jj .gt. 0 .and. jj .le. ny_l)then
                  jstart = min(jstart,j)
                  jend   = max(jend,j)
                  istart = min(istart,i)
                  iend   = max(iend,i)
               end if
            end if

         end do       !all i within window
      end do       !all j within window.
c
c allow a few lines and elements near the border so that a window around
c each laps grid point will be available. Note that wsi = 2 km.
c
      if(istart .gt. 3)then
         istart = istart - 3
      elseif(istart.gt.0)then
         write(6,*)'Could not adjust istart.'
         write(6,*)'Domain too close to wsi western boundary'
      elseif(istart.le.0)then
         write(6,*)'***WARNING*** istart <= 0 !'
      endif
c
      if(iend .lt. elements-3)then
         iend = iend + 3
      elseif(iend .le. elements)then
         write(6,*)'Could not adjust iend.'
         write(6,*)'Domain too close to wsi eastern boundary'
      elseif(iend .gt. elements)then
         write(6,*)'***WARNING*** iend > # of wsi elements!'
      endif
c
      if(jstart .gt. 3)then
         jstart = jstart - 3
      elseif(jstart.gt.0)then
         write(6,*)'Could not adjust jstart.'
         write(6,*)'Domain too close to wsi northern boundary'
      elseif(jstart.le.0)then
         write(6,*)'***WARNING*** jstart <= 0 !'
      endif

      if(jend .lt. lines-3)then
         jend = jend + 3
      elseif(jend .le. lines)then
         write(6,*)'Could not adjust jend.'
         write(6,*)'Domain too close to wsi southern boundary'
      elseif(jend .gt. lines)then
         write(6,*)'***WARNING*** jend > # of wsi lines!'
      endif
c
c output for the user
c
      write(6,*)'istart/iend ',istart,iend
      write(6,*)'jstart/jend ',jstart,jend
c
      goto 999

901   write(6,*)'Error opening parameter file'
      goto 1000

995   write(6,*)'Error opening radar_bounds.inc file'
      write(6,*)'Not generating radar_bounds.inc file'
      goto 1000

999   write(6,*)'generated wsi radar start/end '
      return
1000  end
