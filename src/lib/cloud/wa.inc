!       wa.inc

        parameter (nc = 3)

        real wa(nc)                    ! wavelength (um)
        data wa    /.615,.546,.450/    ! normal
!       data wa    /.615,.551,.450/    ! dscovr 551nm comparison

!       Luminance of sunlight from 1AU distance uniformly scattered 
!       throughout a 360 degree sphere. This will be defined as the
!       spectral radiance matching that of the solar spectrum.
        parameter (day_int0 = 3e9)
