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
       subroutine gen_llij_lut_polar(irad,imax,jmax,lat,lon,istatus)
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

       include 'lapsparms.for'
       include 'vrc.inc'

       integer*4 imax,jmax

       real*4 lat(imax,jmax)
       real*4 lon(imax,jmax)
       real*4 ri(nx_l,ny_l)
       real*4 rj(nx_l,ny_l)
       real*4 wsi_lat(nelems,nlines)
       real*4 wsi_lon(nelems,nlines)

       real*4 pi
       real*4 rdtodg
       real*4 dgtord
       real*4 dgtokm
       real*4 dx,dy
       real*4 la1,lo1,Lov
       real*4 latin,lap

       real*4    fraclat
       real*4    fraclon
       real*4    rlat_diff_deg
       real*4    rlon_diff_deg
       real*4    iline,jline
       real*4    idiff,jdiff

       integer*4 i,j
       integer*4 k,l
       integer*4 n,n1,n2
       integer*4 irad
       integer*4 istart,jstart
       integer*4 iend,jend
       integer*4 ishow_timer
       integer*4 init_timer
       integer*4 itstatus
       integer*4 istatus

       logical   found_line
       logical   found_elem

       character path*100
       character cname*100
       character file*255

       integer lines
       integer nelements
c
c ***************************************************************************
c
      pi=acos(-1.)
      rdtodg=180.0/pi
      dgtord=1./rdtodg
      dgtokm=111.1
      istatus=1
      call get_directory('static',path,n1)
      path=path(1:n1)//'vrc/ '
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
      read(22,51)nelements    !integer
      read(22,51)lines        !integer
      read(22,50)la1          !degrees
      read(22,50)lo1          !degrees
      read(22,50)latin        !not used
      read(22,50)Lov          !degrees
      read(22,50)lap          !not used
50    format(f10.5)
51    format(i4)

      close(22)

      write(6,*)'Parameters from ',file(1:n)
      write(6,*)'dx    ',dx
      write(6,*)'dy    ',dy
      write(6,*)'nelems',nelements
      write(6,*)'nlines',lines
      write(6,*)'la1   ',la1
      write(6,*)'lo1   ',lo1
      write(6,*)'Lov   ',Lov
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Build ri/rj look up for laps domain
c
c  First guess so that all 3661 x 1837 points do not have to be
c  searched.
c
      call  get_radar_bounds(lat,lon,nx_l,ny_l,nlines,nelems,
     &wsi_lat,wsi_lon,istart,jstart,iend,jend)

      write(6,*)'Got the wsi lat/lon get_radar_bounds'
      write(6,*)'This takes awhile, maybe an hour or more.'
      write(6,*)'Some output is coming, and timing stats'
      write(6,*)

      itstatus=init_timer()
      itstatus=ishow_timer()

      do j = ny_l,1,-1
         jline = float(j)/10.
         jdiff = jline - int(jline)

      do i = 1,nx_l
         iline = float(i)/10.
         idiff = iline - int(iline)

         if(idiff.eq.0.00 .and. jdiff.eq.0.00)then
            write(6,29)i,j
29          format(1x,'LAPS(i,j) ',2i6)
            itstatus=ishow_timer()
         end if

         l=jstart
         found_line=.false.
         found_elem=.false.

         do while(.not.found_line)
         if(l.lt.jend)then
            k=istart
            found_elem=.false.
            do while(.not.found_elem)
            if(k.le.iend)then

            if( (wsi_lat(k,l).ge.lat(i,j))        .and.
     &          (wsi_lat(k+1,l+1).le.lat(i,j)) )then

                if( (wsi_lon(k,l).le.lon(i,j))    .and.
     &              (wsi_lon(k+1,l).gt.lon(i,j)) )then

               rlat_diff_deg=abs(wsi_lat(k,l)-wsi_lat(k,l+1))
               rlon_diff_deg=abs(wsi_lon(k,l)-wsi_lon(k+1,l))
               if(rlat_diff_deg.eq.0.0000000)then
                  write(6,*)i,j,k,l,'rlat_diff=0.0'
               endif
               if(rlon_diff_deg.eq.0.0000000)then
                  write(6,*)i,j,k,l,'rlon_diff=0.0'
               endif

               fraclat=(wsi_lat(k,l)-lat(i,j))/rlat_diff_deg
               fraclon=(abs(wsi_lon(k,l))-abs(lon(i,j)))
     &/rlon_diff_deg

               ri(i,j) = float(k)+fraclon
               rj(i,j) = float(l)+fraclat
               found_line=.true.
               found_elem=.true.

                endif

            endif

            k=k+1

            elseif(k.gt.iend)then
               found_elem=.true.
            endif

            enddo

         elseif(l.gt.jend)then
            found_line=.true.
         endif
         l=l+1

         enddo

c            call get_rirj_pol(lat(i,j),
c    &                         lon(i,j),
c    &                         nlines,
c    &                         Lov,dx,dy,la1,lo1,
c    &                         rj(i,j),
c    &                         ri(i,j),
c    &                         istatus)

c
c

       end do       !all i within window
       end do       !all j within window.
c
c output
c
       do i = 1,nx_l,10
       do j = 1,ny_l,10

          write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)

       enddo
       enddo

       cname='wsi_llij_lut_'//c_raddat_types(irad)
       n2=index(cname,' ')-1
       file = path(1:n1)//cname(1:n2)//'.lut'
       n=index(file,' ')
       write(6,*)'Write lat/lon to i/j look up table'
       write(6,*)file(1:n)

       call write_table (file,nx_l,ny_l,lat,lon,ri,rj,istatus)
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
