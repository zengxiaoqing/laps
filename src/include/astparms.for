
	real*8 T1950,J2000,OBLIQ1950,pi,c,cinv,k,r_e_km,r_s_km,r_m_km       
        real*8 emrat,km_per_au,dpr,rpd

        parameter (T1950 = 2433282.423357D0,
     1             J2000 = 2451545.0D0,
     1  	   OBLIQ1950 = 23.44578787D0,
     1	           pi = 3.1415926535897932d0,
     1	           c  = 173.1446327d0,
     1	           cinv = 1.d0/173.1446327d0,
     1	           k  = .01720209895d0,
     1             r_e_km = 6378.137d0,
     1             r_s_km = 696000d0,
     1             r_m_km = 1738d0,
     1             emrat = 81.30058827d0,
     1             km_per_au = 149597870.66d0,
     1             dpr = 180d0/pi,
     1             rpd = pi/180d0)

        real*8 zx,zy,zz
      
!       Saturn's Rings
        parameter (zx = -.0912836831D0)
        parameter (zy = -.0724471759D0)
        parameter (zz = -.9931861335D0)
