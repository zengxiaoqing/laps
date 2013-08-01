        subroutine sidereal_time(ut1,lon,lst)

      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      REAL*8 LST,LON

      include '../../include/astparms.for'

!       UT1 is in days (JD)
!       lon is in degrees (West is negative)
!       lst is in radians

        parameter (rad_per_sec = 2d0 * pi / 86400d0)

      Tu=(ut1-2451545.D0)/36525.d0

        gmst_0ut = 24110.54841d0 +
     1     Tu * (8640184.812866d0  + Tu * (.093104d0 - Tu * 6.2d-6))

        gmst_0ut = gmst_0ut * rad_per_sec

        lst = gmst_0ut + (ut1+.5d0) * 2d0 * pi + lon * rpd

        lst = mod(lst,2d0*pi)

        return
        end
