!       wac.inc

!       parameter (nc = 3)
!       real wa(nc)                    ! wavelength (um)
!       data wa    /.615,.546,.450/

        parameter (nct = 22)

        real wa_tri(nct)               ! wavelength (um)
        data wa_tri /.36,.38,.40,.42,.44,.46,.48,.50,.52,.54,.56,.58,.60&
                    ,.62,.64,.66,.68,.70,.72,.74,.76,.78/

        parameter (gamma = 2.2)

!       http://en.wikipedia.org/wiki/Relative_luminance (luminance)
        rgb_to_y(r,g,b) = 0.2126*R + 0.7152*G + 0.0722*B 

        rgb2luma(R,G,B) = 0.2126*R + 0.7152*G + 0.0722*B ! check luma


