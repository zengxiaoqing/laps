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



      
