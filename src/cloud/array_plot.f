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

      subroutine array_plot(a,b,imax,jmax,NAME,name_array,kmax,cld_hts
     1                     ,scale)

!    ~1990        S. Albers - Ascii plots of cloud fields
!     1997 Aug 01 K. Dritz  - Changed NX_L to imax and NY_L to jmax
!     1997 Aug 01 K. Dritz  - Removed include of lapsparms.for
!     1997 Oct    S. Albers - Make plots work with variable domain sizes.

      integer imax, jmax
      real a(imax,jmax),b(imax,jmax)
      integer ia(imax,jmax)
      character*1 c1a_array(imax,jmax),c1b_array(imax,jmax)
      character*1 name_array(imax,jmax)
      real*4 cld_hts(kmax)
      CHARACTER NAME*10
      character*1 c1_cov(0:13)

      integer max_plot
      parameter(max_plot=65)
      character c_mxp_a*(max_plot)
      character c_mxp_b*(max_plot)

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

      ihigh = imax

      if(imax .gt. 80)then
          iskip = 2
      else
          iskip = 1
      endif

      nplot = (ihigh - 1) / iskip + 1

      if(nplot .gt. max_plot)then ! Prevent lines from getting too long
          nplot = max_plot
          ihigh = 1 + (nplot-1) * iskip
      endif

      nspace = 3

      jskip = 3

      do j = 1,jmax
      do i = 1,ihigh,iskip
        c1a_array(i,j) = c1_cov(int(min(max(A(I,J)*10.*scale,0.),13.)))
        c1b_array(i,j) = c1_cov(int(min(max(B(I,J)*10.*scale,0.),13.)))

        if(NAME(1:4) .eq. 'HORZ')then ! Horizontal Section
            iil = max(i -  iskip    / 2, 1)
            iih = min(i + (iskip-1) / 2, imax)
            jjl = max(j -  jskip    / 2, 1)
            jjh = min(j + (jskip-1) / 2, jmax)

            do ii = iil,iih
            do jj = jjl,jjh
                if(name_array(ii,jj) .ne. ' ')then
                    c1a_array(i,j) = name_array(ii,jj)
                    c1b_array(i,j) = name_array(ii,jj)
                endif

            enddo ! jj
            enddo ! ii

        endif
      enddo
      enddo

      if(NAME(1:4) .eq. 'VERT')then     ! Vertical Section
          do j = jmax,1,-2
              iarg = nint(cld_hts(j))

              do i1 = 1,max_plot
                  i2 = 1 + (i1-1) * iskip
                  if(i2 .le. ihigh)then
                      c_mxp_a(i1:i1) = c1a_array(i2,j)
                      c_mxp_b(i1:i1) = c1b_array(i2,j)
                  endif
              enddo ! i1

              write(6,1001)iarg,c_mxp_a(1:nplot),c_mxp_b(1:nplot)
 1001         FORMAT(1X,i6,2x,a,4x,a)
          enddo ! j

      elseif(NAME(1:4) .eq. 'HORZ')then ! Horizontal array plot
          do j = jmax,1,-jskip

              do i1 = 1,max_plot
                  i2 = 1 + (i1-1) * iskip
                  if(i2 .le. ihigh)then
                      c_mxp_a(i1:i1) = c1a_array(i2,j)
                      c_mxp_b(i1:i1) = c1b_array(i2,j)
                  endif
              enddo ! i1

              write(6,1002)c_mxp_a(1:nplot),c_mxp_b(1:nplot)
 1002         FORMAT(1X,a,3x,a)
          enddo ! j

      endif

      RETURN
      END

