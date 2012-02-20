
!      Driver for sounding ingest
       program snd_main

!      Steve Albers      May-1999       Original Version

       use mem_grid, ONLY: topo

!      Access global parameters
       character*8 c8_project
       character*9 a9_time

       call get_grid_dim_xy(ni,nj,istatus)
       if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           stop
       endif

       call get_c8_project(c8_project,istatus)
       if (istatus .ne. 1) then
          write (6,*) 'Error getting c8_project'
          go to 999
       endif

       call get_laps_cycle_time(laps_cycle_time,istatus)
       if(istatus .ne. 1)stop

       call GETENV('LAPS_A9TIME',a9_time)
       call s_len(a9_time,ilen)

       if(ilen .eq. 9)then
           write(6,*)' systime (from env) = ',a9_time
           call i4time_fname_lp(a9_time,i4time_sys,istatus)
       else
           call get_systime(i4time_sys,a9_time,istatus)
           if(istatus .ne. 1)go to 999
           write(6,*)' systime = ',a9_time
       endif

       allocate( topo(ni,nj), STAT=istat_alloc )
       if(istat_alloc .ne. 0)then
           write(6,*)' ERROR: Could not allocate topo'
           stop
       endif

       call read_static_grid(ni,nj,'AVG',topo,istatus)
       if(istatus .ne. 1)then
           stop
       endif

       call snd_sub(ni,nj,c8_project,laps_cycle_time,i4time_sys)          ! I

 999   continue

       write(6,*)' Successful end of ingest_sounding program...'

       deallocate(topo)

       end

       subroutine snd_sub(ni,nj,c8_project,laps_cycle_time,i4time_sys)    ! I

       character*150 path_to_raw_raob,path_to_local_raob
       character*150 path_to_raw_drpsnd
       character*150 path_to_raw_satsnd
       character*150 path_to_raw_poessnd
       character*200 path_to_raw_tower
       character*200 path_to_raw_radiometer
       character*8 c8_raob_format, c8_drpsnd_format, c8_project
       character*13 filename13,cvt_i4time_wfo_fname13

       logical l_fill_ht

       l_fill_ht = .true. ! default value

       call get_snd_parms(path_to_raw_raob,path_to_local_raob
     1                   ,path_to_raw_drpsnd
     1                   ,path_to_raw_satsnd
     1                   ,path_to_raw_poessnd
     1                   ,path_to_raw_tower
     1                   ,path_to_raw_radiometer
     1                   ,l_fill_ht
     1                   ,istatus)       

       if(istatus .ne. 1)goto999

       lun_out = 21

       if(c8_project .ne. 'RSA')then

!          Read RAOB data
           if(c8_project .eq. 'AFTAC')then
             c8_raob_format = 'NIMBUS'
           else
             c8_raob_format = c8_project
           endif
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_raw_raob,c8_raob_format,i4time_sys
     1                     ,l_fill_ht,lun_out)

           call s_len(c8_project,lenc8)

!          Read Dropsonde data
           if(c8_project .eq. 'NIMBUS')then       
               c8_drpsnd_format = c8_project

           elseif(c8_project .eq. 'CWB')then       
               c8_drpsnd_format = c8_project

           elseif(c8_project(1:min(3,lenc8)) .eq. 'AIR')then       
               c8_drpsnd_format = 'SND'

           else
               c8_drpsnd_format = 'AVAPS'

           endif

           write(6,*)
           write(6,*)' Call ingest_drpsnd, format=',c8_drpsnd_format       
           call ingest_drpsnd(path_to_raw_drpsnd,c8_drpsnd_format
     1                       ,l_fill_ht,lun_out)

       else ! RSA project
           c8_raob_format = 'WFO'
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_raw_raob,c8_raob_format,i4time_sys
     1                     ,l_fill_ht,lun_out)

           c8_raob_format = 'RSA'
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_local_raob,c8_raob_format,i4time_sys
     1                     ,l_fill_ht,lun_out)

       endif ! RSA project

       write(6,*)
       write(6,*)' Call tower_driver_sub'
       maxobs = 6000
       itime_before = 900
       itime_after = 900

       call get_max_stations(maxsta, istatus)
       if(istatus .ne. 1)stop

       call tower_driver_sub(ni,nj,lun_out
     1                      ,maxobs,laps_cycle_time
     1                      ,path_to_raw_tower
     1                      ,itime_before,itime_after
     1                      ,maxsta
     1                      ,istatus)

!      Satellite soundings -- added goes (7/26/04 db)

       write(6,*)
       if(c8_project .eq. 'AFWA')then
           write(6,*)' Call ingest_satsnd for AFWA...'
           call ingest_satsnd(path_to_raw_satsnd,c8_project
     1                       ,i4time_sys,lun_out)
       else
           write(6,*)' Call ingest_satsnd for POES...'
           call ingest_satsnd(path_to_raw_poessnd,c8_project
     1                       ,i4time_sys,lun_out)

           write(6,*)
           write(6,*)' Call ingest_goessnd...'
           call ingest_goessnd (path_to_raw_satsnd, i4time_sys, 
     1                          lun_out, istatus)
       endif

!      Radiometer data
       write(6,*)
       i4time_rad = (i4time_sys/3600) * 3600
 900   filename13 = cvt_i4time_wfo_fname13(i4time_rad)
       call s_len(path_to_raw_radiometer,len_path)
       if(len_path .gt. 0)then
           write(6,*)' Call get_radiometer_data, input file at '
     1              ,filename13
           call get_radiometer_data(i4time_sys,laps_cycle_time,ni,nj
     1                             ,i4time_sys - 855/2 ! Radiomtr cyc time / 2
     1                             ,i4time_sys + 855/2 ! Radiomtr cyc time / 2
     1                             ,path_to_raw_radiometer(1:len_path)
     1                                           //'/'//filename13
     1                             ,lun_out,istatus)
           if(istatus .ne. 1 .and. i4time_sys - i4time_rad .lt. 3600
     1                                                           )then      
               write(6,*)' Trying radiometer for an hour earlier'
               i4time_rad = i4time_rad - 3600
               go to 900
           endif
       endif ! len_path > 0

       close(lun_out) ! Output SND file

 999   return
       
       end


       subroutine get_snd_parms(path_to_raw_raob,path_to_local_raob
     1                         ,path_to_raw_drpsnd
     1                         ,path_to_raw_satsnd
     1                         ,path_to_raw_poessnd
     1                         ,path_to_raw_tower
     1                         ,path_to_raw_radiometer
     1                         ,l_fill_ht
     1                         ,istatus)

       character*150 path_to_raw_raob,path_to_local_raob
     1              ,path_to_raw_drpsnd
     1              ,path_to_raw_satsnd
     1              ,path_to_raw_poessnd
     1              ,path_to_raw_radiometer

       character*200 path_to_raw_tower

       logical l_fill_ht

       namelist /snd_nl/ path_to_raw_raob,path_to_local_raob
     1                  ,path_to_raw_satsnd
     1                  ,path_to_raw_poessnd,path_to_raw_tower
     1                  ,path_to_raw_drpsnd,path_to_raw_radiometer 
     1                  ,l_fill_ht      
 
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


      subroutine get_nlevels_snd(pressure_mb,height_m,r_missing_data
     1                          ,nlevels_raw,nlevels_snd) 

      real pressure_mb(nlevels_raw)     
      real height_m(nlevels_raw)     

      do i = 1,nlevels_raw
          if(abs(pressure_mb(i)) .le. 2000. .or.
     1       abs(height_m(i)) .le. 1000000.)then
              nlevels_snd = i
          endif
      enddo ! i

      return
      end

