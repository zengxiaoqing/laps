cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
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



      
