

        function rh_to_cldcv(rh)

!       1999 Steve Albers FSL

!       Convert fractional rh into fractional cloud cover. This is defined
!       as a total cloud cover through a particular cloud layer.
!       It is recommended that the input RH be calculated with subroutine
!       'make_rh' using a 't_ref' between 0C and -10C.

        ramp_thresh = 0.80

        if(rh .lt. ramp_thresh)then
            rh_to_cldcv = 0.0
        elseif(rh .le. 1.00)then
            rh_to_cldcv = (rh - ramp_thresh) / (1.00 - ramp_thresh)  
        else
            rh_to_cldcv = 1.00
        endif

        return
        end
