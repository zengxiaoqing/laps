       function zenith(lat,lon,slat,slon)

c      computes zenith angle give earth lat,lon and satellite subpoint
c      slat, slon.

c      returns zenith angle in degrees for use in SSEC routines

c      input coordinates are all in radians.  E longitude =+

c      authored by Dan Birkenheuer   12/1/95


        real lat,lon,slat,slon
        real zenith
        real cos_of_zen


        cos_of_zen = sin (slat)*sin(lat)+cos(slat)*cos(lat)*
     1      cos(abs(slon-lon))

        zenith = acos(cos_of_zen)

        zenith = zenith + atan2 ( sin(zenith),6.6166-cos_of_zen )

c       convert to degrees

        zenith = zenith * 180./acos(-1.0)

        return

      end




