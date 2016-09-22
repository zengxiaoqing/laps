
        subroutine refract(true,apparent,pres)

        include 'trigd.inc'

        real*8 true,apparent,pres

        if(true .gt. -1.0d0)then
            r = 1.02d0/ dtand(true+10.3d0/(true+5.11d0))
            apparent = true + (r/60. * pres/1013.d0)
        else
            apparent = 0.
        endif

        end
