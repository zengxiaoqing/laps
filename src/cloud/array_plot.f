cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
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
      real cld_hts(kmax)
      CHARACTER NAME*10
      character*1 c1_cov(-1:13)

      integer max_plot
      parameter(max_plot=65)
      character c_mxp_a*(max_plot)
      character c_mxp_b*(max_plot)

      data c1_cov
     1   /' ','.','1','2','3','4','5','6','7','8','9','#','*',')','.'/       

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

      iskip = max(imax/48,1)

      nplot = (ihigh - 1) / iskip + 1

      if(nplot .gt. max_plot)then ! Prevent lines from getting too long
          nplot = max_plot
          ihigh = 1 + (nplot-1) * iskip
      endif

      nspace = 3

      jskip = max(jmax/32,2)

      write(6,*)'iskip/jskip/nplot/ihigh',iskip,jskip,nplot,ihigh

      do j = 1,jmax
      do i = 1,ihigh,iskip
        c1a_array(i,j)=c1_cov(int(min(max(A(I,J)*10.*scale,-1.),13.)))
        c1b_array(i,j)=c1_cov(int(min(max(B(I,J)*10.*scale,-1.),13.)))      

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
 1001         FORMAT(i6,2x,a,2x,a)
          enddo ! j

      elseif(NAME(1:4) .eq. 'HORZ')then ! Horizontal array plot
          do j = jmax,1,-jskip

              do i1 = 1,max_plot          ! small grid
                  i2 = 1 + (i1-1) * iskip ! large grid
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

