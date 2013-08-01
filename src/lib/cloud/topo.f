
      SUBROUTINE TOPO(PHI,LON,ut1,TX,TY,TZ)
      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      REAL*8 LST,LON,JD

      include '../../include/astparms.for'

      parameter (R_E_AU = R_E_KM / KM_PER_AU)

!       All arguments are Real * 8
!       Input PHI = latitude in radians
!             LON = longitude in degrees (W is negative)
!             UT1 = Julian Date
!       Output TX,TY,TZ = Equatorial coordinates of point on surface of earth

c     write(6,*)R_E_AU

      FF=(1.D0-1.D0/298.257D0)**2
      CC=1./DSQRT(DCOS(PHI)**2+FF*DSIN(PHI)**2)
      XYG=R_E_AU*CC*DCOS(PHI)

      call sidereal_time(ut1,lon,lst)

      TX=DCOS(LST)*XYG
      TY=DSIN(LST)*XYG
      TZ=R_E_AU*CC*DSIN(PHI)*FF

      RETURN
      END

