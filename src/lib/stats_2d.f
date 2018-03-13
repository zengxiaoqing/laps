cdis
c
c
        subroutine stats_2d(ni,nj,x,y,a_t,b_t,xbar,ybar,r
     1                     ,bias,std,badflag,istatus)
c
c*******************************************************************************
c
c       Routine to calculate regression coefficients from y/x data
c       (formerly t and elev data).

c       Values returned are the coefficients for the regression equations
c       as well as mean of each array, bias, rms   
c
c       Changes:
c               P.A. Stamus     12-01-88        Original (from J. McGinley)
c
c       Inputs/Outputs:
c
c          Variable     Var Type    I/O   Description
c         ----------   ----------  ----- -------------
c          num_sfc         I         I    Number of surface stations.
c          x               RA        I    
c          y               RA        I    
c          b_t             R         O    intercept (y = ax + b)
c          a_t             R         O    slope
c          xbar            R         O    Mean value of the stations.
c          ybar            R         O    Mean value of the gridded field
c
c       User Notes:
c
c       1. Units are not changed in this routine.
c
c*******************************************************************************
c
        real x(ni,nj), y(ni,nj)           

c
c
c.....  Set up storage variables.
c
        cnt = 0.
        sumxy = 0.
        sumx = 0.
        sumy = 0.
        sumx2 = 0.
        sumy2 = 0.
c
c.....  Gather sums and then calculate the 'a' and 'b' for the regression
c.....  equation y = az + b, for both the temperature and dew point.  The
c.....  'b' is the intercept with sea level, and the 'a' is the lapse rate.
c.....  'y' represents the y value and 'z' is x/analyzed
c.....  Also calculate the mean elevation of the stations.
c
        istatus = 0

!       write(6,*)' i  y   x  '
        do i=1,ni
        do j=1,nj
          if(x(i,j).eq.badflag .or. y(i,j).eq.badflag) go to 10
!         write(6,5)i,y(i,j),x(i,j)
5         format(i6,2f8.3)
          sumxy = (x(i,j) * y(i,j)) + sumxy
          sumx = x(i,j) + sumx
          sumx2 = (x(i,j) * x(i,j)) + sumx2
          sumy = y(i,j) + sumy
          cnt = cnt + 1.
          istatus = 1
 10       continue
        enddo
        enddo

        if(cnt .eq. 0.)then
            write(6,*)' no data for stats'
            xbar = 0.
            ybar = 0.
            bias = 0.
            std = 0.
            istatus = 0
            return
        endif

!       Slope    
        denominator = (cnt*sumx2 - sumx*sumx)
        if(denominator .ne. 0.)then
            a_t = (cnt*sumxy - sumx*sumy) / denominator                
        else
            a_t = 0.
            istatus = 0
        endif

!       Intercept
        b_t = (sumy - a_t * sumx) / cnt
c
        xbar = sumx / cnt
        ybar = sumy / cnt

        write(6,*)' cnt/a_t/b_t (slope / intercept) = '
     1           ,int(cnt),a_t,b_t

        write(6,*)' xbar,ybar = ',xbar,ybar

!       Calculate rms (stdev) of the ob-background differences
        cnt = 0
        sumsq = 0.
        do i=1,ni
        do j=1,nj
            if(x(i,j).ne.badflag .and. y(i,j).ne.badflag) then     
                cnt = cnt + 1.
                sumsq = sumsq + (y(i,j)-x(i,j))**2
            endif
        enddo ! i
        enddo ! i

        if(cnt .gt. 0.)then
            std = sqrt(sumsq / cnt)
        else
            std = 0.
        endif

        bias = ybar - xbar

!       if(a_t .lt. 0.1 .OR. a_t .gt. 10.)then
!          write(6,*)' Warning, slope is ill conditioned'
!          istatus = 0
!       endif
c
c.....  End of routine
c
!       Compute correlation coefficient
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        do i=1,ni
        do j=1,nj
            if(x(i,j).ne.badflag .and. y(i,j).ne.badflag) then     
                sum1 = sum1 + ( (x(i,j) - xbar) * (y(i,j) - ybar) )
                sum2 = sum2 + (x(i,j) - xbar)**2
                sum3 = sum3 + (y(i,j) - ybar)**2
            endif
        enddo ! i
        enddo ! i
        denom = sqrt(sum2) * sqrt(sum3)

        if(denom .ne. 0.)then
            r = sum1 / denom
        else
            r = 0.
        endif

        write(6,*)' regression sums = ',sum1,sum2,sum3
        write(6,*)' ratio of points = ',cnt/(float(ni*nj))

        write(6,900)int(cnt),bias,std,r
900     format(2x,' N/bias/rms/r = ',i8,2f9.2,f9.3/)

        return
        end

