
       subroutine get_windob_time_window(c_obstype,i4_window_ob
     1                                  ,istatus)

!      Returns half-width of observation time window in seconds

       character*(*) c_obstype

       call get_laps_cycle_time(ilaps_cycle_time,istatus)
       if(istatus .eq. 1)then
           write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
       else
           write(6,*)' Error getting laps_cycle_time'
           return
       endif

       i4_window_ob = ilaps_cycle_time

       if(c_obstype .eq. 'ACARS')i4_window_ob = ilaps_cycle_time / 2     
       if(c_obstype .eq. 'VAD')  i4_window_ob = 1800    

       return
       end


       subroutine get_tempob_time_window(c_obstype,i4_window_ob
     1                                  ,istatus)

!      Returns half-width of observation time window in seconds

       character*(*) c_obstype

       call get_laps_cycle_time(ilaps_cycle_time,istatus)
       if(istatus .eq. 1)then
           write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
       else
           write(6,*)' Error getting laps_cycle_time'
           return
       endif

       i4_window_ob = ilaps_cycle_time

       if(c_obstype .eq. 'ACARS')i4_window_ob = ilaps_cycle_time / 2     

       return
       end
