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

      subroutine array_plot(a,b,imax,jmax,NAME,name_array,kmax,cld_hts,s
     1cale)

!     1997 Aug 01 K. Dritz  - Changed NX_L to imax and NY_L to jmax
!     1997 Aug 01 K. Dritz  - Removed include of lapsparms.for

      dimension a(imax,jmax),b(imax,jmax),ia(imax,jmax)
      character*1 c1a_array(imax,jmax),c1b_array(imax,jmax)
      character*1 name_array(imax,jmax)
      real*4 cld_hts(kmax)
      CHARACTER NAME*10
      character*1 c1_cov(0:13)
      data c1_cov
     1   /'.','1','2','3','4','5','6','7','8','9','#','*',')','.'/

c find max and min
      amax=-1.e30
      amin=1.e30
!     sum=0.
      cnt=0.

      do j=1,jmax
      do i=1,imax
          IA(I,J)=0
          if(a(i,j) .gt. 0.05 .or. b(i,j) .gt. 0.05)then
              cnt=cnt+1.
          endif

          if(a(i,j) .ge. 0.00)then
              if (a(i,j).gt.amax) amax=a(i,j)
              if (a(i,j).lt.amin) amin=a(i,j)
!             sum=sum+a(i,j)
          endif

    1     continue
      enddo
      enddo

      IF(CNT.EQ.0) THEN
           write(6,1004) NAME
 1004      FORMAT(1X,'ALL Values in array < 0.05 ' ,A10)
           RETURN
      else
!          write(6,1005) amax,amin,name
 1005      FORMAT(1X,2e12.3,' Max min ',a10)
      ENDIF

!     sum=sum/cnt
!     diffx=amax-sum
!     diffn=sum-amin
!     IF(DIFFN.EQ.0.AND.DIFFx.EQ.0) THEN
!          write(6,1235) CNT,SUM
! 1235 FORMAT(1X,'ALL 'F6.0,' POINTS HAVE EQUAL VALUE OF ',
!     1E12.4)
!         RETURN
!     ENDIF
      fact=1.
      iter=1
      IFLAG=0

      do j = 1,jmax
      do i = 1,imax
        c1a_array(i,j) = c1_cov(int(min(max(A(I,J)*10.*scale,0.),13.)))
        c1b_array(i,j) = c1_cov(int(min(max(B(I,J)*10.*scale,0.),13.)))
        if(NAME(1:4) .eq. 'HORZ')then ! Horizontal Section
            if(name_array(i,j) .ne. ' ')then
                c1a_array(i,j) = name_array(i,j)
                c1b_array(i,j) = name_array(i,j)
            endif
        endif
      enddo
      enddo

!     iplot = min(imax,57)
      iplot = 57

      if(NAME(1:4) .eq. 'VERT')then     ! Vertical Section
          do j = jmax,1,-2
              iarg = nint(cld_hts(j))
              write(6,1001)iarg,
     1   (c1a_array(i,j),i=1,iplot),
     1   (c1b_array(i,j),i=1,iplot)
 1001         FORMAT(1X,i6,2x,57A1,4x,57A1)
          enddo ! j
      elseif(NAME(1:4) .eq. 'HORZ')then ! Horizontal array plot
          do j = jmax,1,-3
              write(6,1002)
     1   (c1a_array(i,j),i=1,iplot),
     1   (c1b_array(i,j),i=1,iplot)
 1002         FORMAT(1X,6x,2x,57A1,4x,57A1)
          enddo ! j
      endif

      RETURN
      END

