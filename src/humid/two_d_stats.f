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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine two_d_stats (ii,jj,x,excluded_value)

      implicit none
      integer ii,jj
      real x(ii,jj)
      real excluded_value


c     internal variables

      integer i,j
      real maxa,mina
      integer maxi,maxj, mini,minj
      integer counter
      real x_linear (ii*jj)
      integer istatus
      real ave,adev,sdev,var,skew,curt

c     compute the max and min values of the field

      maxa = -1.e32
      mina = 1.e32

      do i = 1,ii
         do j = 1,jj
            if(x(i,j) .ne. excluded_value) then
               mina = min (mina, x(i,j))
               maxa = max (maxa, x(i,j))
               if(mina.eq.x(i,j)) then
                  mini = i
                  minj = j
               endif
               if(maxa.eq.x(i,j)) then
                  maxi = i
                  maxj = j
               endif
            endif
         enddo
      enddo

      write (6,*)
      write (6,*) 'Computed 2-D statistics'
      write (6,*) 'Max and min values in the field'
      write (6,*) 'Max,i,j', maxa,maxi,maxj
      write (6,*) 'Min,i,j', mina,mini,minj
      write (6,*)

c     begin computation of moment statistics.

      counter = 0

      do i = 1,ii
         do j = 1,jj
            if (x(i,j) .ne. excluded_value) then
               counter = counter + 1
               x_linear(counter) = x(i,j)
            endif
         enddo
      enddo

      call moment_b (x,ii*jj,ave,adev,sdev,var,skew,curt,istatus)

      write (6,*) 'Moment output'

      write (6,*) 'field average ',ave, '+/-',sdev

      return
      end



      
