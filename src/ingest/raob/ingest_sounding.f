
!      Driver for sounding ingest

!      Steve Albers      May-1999       Original Version

       character*150 path_to_raw_raob,path_to_local_raob
       character*150 path_to_raw_drpsnd
       character*150 path_to_raw_satsnd
       character*200 path_to_raw_tower
       character*200 path_to_raw_radiometer
       character*8 c8_raob_format, c8_project

       call get_laps_config('nest7grid',istatus)

       call get_snd_parms(path_to_raw_raob,path_to_local_raob
     1                   ,path_to_raw_drpsnd
     1                   ,path_to_raw_satsnd
     1                   ,path_to_raw_tower
     1                   ,path_to_raw_radiometer
     1                   ,istatus)       
       if(istatus .ne. 1)goto999

       call get_c8_project(c8_project,istatus)
       if (istatus .ne. 1) then
          write (6,*) 'Error getting c8_project'
          go to 999
       endif

       call get_grid_dim_xy(ni,nj,istatus)
       if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           stop
       endif

       call get_laps_cycle_time(laps_cycle_time,istatus)
       if(istatus .ne. 1)stop

       lun_out = 21

       if(c8_project .ne. 'RSA')then
           c8_raob_format = c8_project
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_raw_raob,c8_raob_format)

           if(c8_project .eq. 'CWB')then
               c8_raob_format = c8_project
               write(6,*)
               write(6,*)' Call ingest_drpsnd, format=',c8_raob_format       
               call ingest_drpsnd(path_to_raw_drpsnd,c8_raob_format
     1                           ,lun_out)

           else
               c8_raob_format = 'AVAPS'
               write(6,*)
               write(6,*)' Call ingest_drpsnd, format=',c8_raob_format       
               call ingest_drpsnd(path_to_raw_drpsnd,c8_raob_format
     1                           ,lun_out)

           endif

       else ! RSA project
           c8_raob_format = 'WFO'
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_raw_raob,c8_raob_format)

           c8_raob_format = 'RSA'
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_local_raob,c8_raob_format)

           write(6,*)
           write(6,*)' Call tower_driver_sub'
           maxobs = 6000
           itime_before = 900
           itime_after = 900

           call get_max_stations(maxsta, istatus)
           if(istatus .ne. 1)stop

           call tower_driver_sub(ni,nj,lun_out
     1                           ,maxobs,laps_cycle_time
     1                           ,path_to_raw_tower
     1                           ,itime_before,itime_after
     1                           ,maxsta
     1                           ,istatus)

       endif ! RSA project
 

       write(6,*)
       write(6,*)' Call ingest_satsnd'
       call ingest_satsnd(path_to_raw_satsnd)

       write(6,*)
       write(6,*)' Call get_radiometer_data (under construction)'
!      call get_radiometer_data(path_to_raw_radiometer,laps_cycle_time       
!     1                        ,ni,nj,lun_out)

       close(lun_out) ! Output SND file

 999   end


       subroutine get_snd_parms(path_to_raw_raob,path_to_local_raob
     1                         ,path_to_raw_drpsnd
     1                         ,path_to_raw_satsnd
     1                         ,path_to_raw_tower
     1                         ,path_to_raw_radiometer
     1                         ,istatus)

       character*150 path_to_raw_raob,path_to_local_raob
     1              ,path_to_raw_drpsnd
     1              ,path_to_raw_satsnd
     1              ,path_to_raw_radiometer

       character*200 path_to_raw_tower

       namelist /snd_nl/ path_to_raw_raob,path_to_local_raob
     1                  ,path_to_raw_satsnd,path_to_raw_tower
     1                  ,path_to_raw_drpsnd,path_to_raw_radiometer       
 
       character*150 static_dir,filename
 
       call get_directory('nest7grid',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/snd.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,snd_nl,err=901)
       close(1)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading snd_nl in ',filename
       istatus = 0
       return

       end
