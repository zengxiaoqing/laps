
       subroutine get_windob_time_window(c_obstype,i4_window_ob
     1                                  ,istatus)

!      Returns half-width of observation time window in seconds

       character*(*) c_obstype

       logical l_parse

       call get_laps_cycle_time(ilaps_cycle_time,istatus)
       if(istatus .eq. 1)then
           write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
       else
           write(6,*)' Error getting laps_cycle_time'
           return
       endif

       i4_window_ob = ilaps_cycle_time


!      The parse routine should help to handle blanks OK, hopefully it won't
!      get fooled by "look alike" obstypes.

       if(l_parse(c_obstype,'ACARS'))i4_window_ob = ilaps_cycle_time / 2     
       if(l_parse(c_obstype,'VAD'))  i4_window_ob = 1800    

       return
       end


       subroutine get_tempob_time_window(c_obstype,i4_window_ob
     1                                  ,istatus)

!      Returns half-width of observation time window in seconds

       character*(*) c_obstype

       logical l_parse

       call get_laps_cycle_time(ilaps_cycle_time,istatus)
       if(istatus .eq. 1)then
           write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
       else
           write(6,*)' Error getting laps_cycle_time'
           return
       endif

       i4_window_ob = ilaps_cycle_time

!      The parse routine should help to handle blanks OK, hopefully it won't
!      get fooled by "look alike" obstypes.

       if(l_parse(c_obstype,'ACARS'))i4_window_ob = ilaps_cycle_time / 2     

       i4_window_ob = min(i4_window_ob,5400)

       return
       end
