

        function rh_to_cldcv(rh)

!       1999 Steve Albers FSL

!       Convert fractional rh into fractional cloud cover. This is defined
!       as a total cloud cover through a particular cloud layer.

        ramp_thresh = 0.60

        if(rh .lt. ramp_thresh)then
            rh_to_cldcv = 0.0
        elseif(rh .le. 1.00)then
            rh_to_cldcv = (rh - ramp_thresh) / (1.00 - ramp_thresh)  
        else
            rh_to_cldcv = 1.00
        endif

        return
        end
