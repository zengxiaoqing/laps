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
      subroutine set_missing_sat(csatid,csattype,chtype,
     &       image_in,nx,ny,smsng,r_missing_data,istatus)
c
c
c J. Smart    11-13-95   include technique to determine bad (outliers) data in
c                        addition to missing sat pixels
c
      implicit none
c
      integer nx,ny,i,j,n,nn,n1,n2
      integer ii,jj
      integer istatus
      integer nstatus
      integer mstatus
      integer istat_status
      integer imiss_status
      integer linestart,lineend
      integer elemstart,elemend
      integer ibnd,jbnd
      integer idum
      integer i2dum

      real*4    image_in(nx,ny)
      real*4    image_temp(nx,ny)
      real*4    data(125)
      real*4    ave,adev,sdev,var,skew,curt
      real*4    smsng
      real*4    r_missing_data
      real*4    rlow,rhigh
      character csattype*(*)
      character chtype*(*)
      character csatid*(*)
      character cid*2
      character path*100
      character cname*255
      character filename*255
c
c note that this quality control step is performed with satellite counts
c both ir and vis.
c
      ibnd=nx
      jbnd=ny
c
c the test is different depending on the value of smsng (satellite missing).
c the missing satellite data values are defined in data/static/satellite_lvd.nl
c
      rlow=0.0
      if(smsng.gt.0.0)then
         rhigh=smsng
      else
         rhigh=255.
      endif

      if(csatid.ne.'gmssat'.and.csatid.ne.'meteos')then
         if(csattype.eq.'gvr'.or.csattype.eq.'gwc')then
            rhigh = 1023.
            if(chtype.eq.'4u ')rlow=69.
            if(chtype.eq.'wv ')rlow=30.
            if(chtype.eq.'11u'.or.chtype.eq.'12u')rlow=16.
            if(chtype.eq.'vis')rlow=28.
         endif
      endif

      istat_status=0
      imiss_status=0

      do j=2,jbnd-1
      do i=2,ibnd-1

         if(image_in(i,j).eq.smsng.or.
     &        image_in(i,j).lt.rlow .or.
     &        image_in(i,j).gt.rhigh   )then

             n=0
             do jj=j-1,j+1
             do ii=i-1,i+1
                if(image_in(ii,jj).ne.smsng .and.
     &              image_in(ii,jj).gt.rlow  .and.
     &              image_in(ii,jj).lt.rhigh )then
                   n=n+1
                   data(n)=image_in(ii,jj)
                endif
             enddo
             enddo

             if(n.ge.2)then
                call moment(data,n,
     &              ave,adev,sdev,var,skew,curt,
     &              mstatus)

                if(abs(image_in(i,j)-ave).gt.3.0*sdev.or.
     &              image_in(i,j).eq.smsng.or.
     &              image_in(i,j).lt.rlow .or.
     &              image_in(i,j).gt.rhigh)then

                   image_temp(i,j) = ave
                   istat_status=istat_status-1 
                else
                   image_temp(i,j) = r_missing_data
                   imiss_status=imiss_status-1
                endif
             else
                image_temp(i,j)=r_missing_data
                imiss_status=imiss_status-1
             endif

         else

            image_temp(i,j) = image_in(i,j)

         endif

      end do
      end do
c
c take of boundaries
c
      do i=1,ibnd
         if(image_in(i,1).eq.smsng.or.
     &        image_in(i,1).lt.rlow .or.
     &        image_in(i,1).gt.rhigh)then

            image_temp(i,1)=r_missing_data
            imiss_status=imiss_status-1
         else
            image_temp(i,1)=image_in(i,1)
         endif

         if(image_in(i,jbnd).eq.smsng.or.
     &        image_in(i,jbnd).lt.rlow .or.
     &        image_in(i,jbnd).gt.rhigh)then

            image_temp(i,jbnd)=r_missing_data
            imiss_status=imiss_status-1
         else
            image_temp(i,jbnd)=image_in(i,jbnd)
         endif
      enddo
      do j=1,jbnd
         if(image_in(1,j).eq.smsng.or.
     &        image_in(1,j).lt.rlow .or.
     &        image_in(1,j).gt.rhigh)then

            image_temp(1,j)=r_missing_data
            imiss_status=imiss_status-1
         else
            image_temp(1,j)=image_in(1,j)
         endif

         if(image_in(ibnd,j).eq.smsng.or.
     &        image_in(ibnd,j).lt.rlow .or.
     &        image_in(ibnd,j).gt.rhigh)then

            image_temp(ibnd,j)=r_missing_data
            imiss_status=imiss_status-1
         else
            image_temp(ibnd,j)=image_in(ibnd,j)
         endif
      enddo

      do j=1,jbnd
      do i=1,ibnd
         image_in(i,j)=image_temp(i,j)
      enddo
      enddo
      write(6,*)'   # reset to r_missing: ',imiss_status
      write(6,*)'   # reset to average  : ',istat_status
      istatus=imiss_status+istat_status

1000  return
      end
