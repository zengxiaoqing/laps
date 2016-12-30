
      subroutine great_circle(alat1,alon1,alat2,alon2,gcdist,gcbearing)

!     http://www.dtcenter.org/met/users/docs/write_ups/gc_simple.pdf

!     solve for bearing of location 2 as seen from location 1
!     north/up is zero bearing with angles counted clockwise
!     angles are in radians

      hav(x) = sin(x/2.)**2
      hav_inv(x) = 2. * asin(sqrt(x))

      pi = 3.14159265
      rpd = pi / 180.

      phia = alat1 * rpd
      phib = alat2 * rpd
      deltal = (alon1 - alon2) * rpd

      s = cos(phib) * sin(deltal)
      c = cos(phia) * sin(phib) - sin(phia) * cos(phib) * cos(deltal) 
 
      beta = atan2(s,c)
      gcbearing = beta / rpd

      hav_theta = hav(phia-phib) + cos(phia) * cos(phib) * hav(deltal)
      theta = hav_inv(hav_theta)
      gcdist = theta / rpd

      return
      end
