c
c
      subroutine writev(a,imax,jmax,m,cs,atime,onoff,offset)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      parameter(im=500,jm=500)
      parameter (badflag = -99.9)
      real a(m,m)
      integer ii(im,jm)   ,onoff,m 
      real smax,smin,s,sum,rms
      character*12 cs,atime*24
      integer i,j
c
      smin=1.e25
      smax=-1.e25
      sum=0
      rms=0
      cnt=0
      do j=1,jmax
      do i=1,imax
         ii(i,j)=0 
         if(a(i,j).ne.badflag) then
            sum=sum+a(i,j)-offset
            if((a(i,j)-offset).gt.smax) smax=a(i,j)-offset
            if((a(i,j)-offset).lt.smin) smin=a(i,j)-offset
            cnt=cnt+1.
         endif
      enddo !i
      enddo !j
      sum=sum/cnt              
      s=amax1(abs(smax),abs(smin)) 
      do i=1,25
         if((s/10.**(i-12)).lt.1000.) then
            iii=i-12
            go to 1
         endif
      enddo !i
    1 write(8,999) cs
 999  format(//1x,a12)
      write(8,1000) atime,iii, smax,smin
 1000 format(13x,a24,' x 10**',i3,' Max: ',e10.4,'  Min: ',e10.4)
      do j=1,jmax
      do i=1,imax
         if(a(i,j).ne.badflag) then
            rms=rms+((a(i,j)-offset)-sum)**2 
            ii(i,j)=(a(i,j)-offset)*10.**(-iii)+.5
         endif
      enddo !i
      enddo !j
      rms=sqrt(rms/cnt)                
      write(8,1004) rms,sum,offset
 1004 format(1x,'RMS  = ',e12.4,' MEAN = ',e12.4,' OFFSET = ',e12.4)
 1001 format(1x,'transposed vector')
      if (jmax.eq.1) then
         write (8,1001)
         if(ONOFF.eq.1)write (8,1002) (ii(i,1),i=1,imax)
      else
         continue
         if(onoff.eq.1) then
            do i=1,imax
               write(8,1005) i
 1005          format(/i3)
               write(8,1002) (j,ii(i,j),j=1,jmax)
            enddo !i
         endif
      endif
 1002 format(10i7)
c
      return
      end 
c
c
      subroutine writei(a,imax,m,cs,atime)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      integer a(m)
      character*72 cs,atime*24
      integer i
c
      write(8,1000) cs,atime
 1000 format(//1x,a72)
       write (8,1002) (a(i)  ,i=1,imax)
 1002 format(10i7)
c
      return
      end 
c
c
      subroutine writec(a,imax,m,cs,atime)
c
c*********************************************************************
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
      character*5   a(m)
      character*14 cs,atime*24
c
      write(8,1000) cs,atime
 1000 format(//1x,a14)
       write (8,1002) (a(i)(1:5)  ,i=1,imax)
 1002 format(10(2x,a5))
c
      return
      end 
c
c
       Subroutine writemon(ta,tda,ua,va,pmsla,alta,
     &   nvar,maxsta,m,monster,itm1) 
c
c*********************************************************************
c     Subroutine puts the latest sets of obs into monster for
c     fourier processing.
c
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       24 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c*********************************************************************
c
       real ta(m),tda(m),ua(m),va(m),pmsla(m),alta(m)
       real monster(m,m,nvar)
c
       do k=1,nvar
          do l=m-1,1,-1
             do i=1,maxsta
                monster(i,l+1,k)=monster(i,l,k)
             enddo !i
          enddo !l
       enddo !k
       if(nvar.eq.1) then
          do i=1,maxsta
             monster(i,1,1)=ta(i)
          enddo !i
       else
          do i=1,maxsta
             monster(i,1,1)=ta(i)
             monster(i,1,2)=tda(i)
             monster(i,1,3)=ua(i)
             monster(i,1,4)=va(i)  
             monster(i,1,5)=pmsla(i)
             monster(i,1,6)=alta(i)
          enddo !i
       endif
c
       return
       end

