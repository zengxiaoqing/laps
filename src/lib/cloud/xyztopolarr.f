
      subroutine xyz_to_polar_r(x,y,z,dec,ra,r)

      IMPLICIT REAL*8(A-Z)

      ATAN3(X,Y)=DMOD((DATAN2(X,Y)+6.2831853071796D0),6.2831853071796D0)

      r=DSQRT(x**2+y**2+z**2)
      dec=ASIN(z/r)
      ra=ATAN3(y,x)

      return
      end
