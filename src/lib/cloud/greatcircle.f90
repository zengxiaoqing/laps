
      subroutine great_circle(alat1,alon1,alat2,alon2,gcdist,gcbearing)


!     http://www.dtcenter.org/met/users/docs/write_ups/gc_simple.pdf

!     solve for bearing of location 2 as seen from location 1
!     north/up is zero bearing with angles counted clockwise
!     angles are in radians

      pi = 3.14159265
      rpd = pi / 180.

!     hav(x) = sin(x/2.)**2
!     hav_inv(x) = ???

      phia = alat1 * rpd
      phib = alat2 * rpd
      deltal = (alon1 - alon2) * rpd

      s = cos(phib) * sin(deltal)
      c = cos(phia) * sin(phib) - sin(phia) * cos(phib) * cos(deltal) 
 
      beta = atan2(s,c)
      gcbearing = beta / rpd

      return
      end
