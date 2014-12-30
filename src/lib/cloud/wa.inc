!       wa.inc

        parameter (nc = 3)

        real wa(nc)                    ! wavelength (um)
        data wa    /.615,.546,.450/

!       Luminance of sunlight from 1AU distance uniformly scattered 
!       throughout a 360 degree sphere. This can be expressed as a
!       function of wavelength by normalizing the radiance at a given
!       wavelength to the radiance of the solar spectrum.
        parameter (day_int0 = 3e9)
