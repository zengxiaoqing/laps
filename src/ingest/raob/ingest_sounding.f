
!      Driver for sounding ingest

!      Steve Albers      May-1999       Original Version

       character*150 path_to_raw_raob,path_to_local_raob
       character*150 path_to_raw_drpsnd
       character*150 path_to_raw_satsnd
       character*200 path_to_raw_tower
       character*200 path_to_raw_radiometer
       character*8 c8_raob_format, c8_project
       character*9 a9_time
       character*13 filename13,cvt_i4time_wfo_fname13

!      call get_laps_config('nest7grid',istatus)

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

       call get_systime(i4time_sys,a9_time,istatus)
       if(istatus .ne. 1)stop
       write(6,*)' systime = ',a9_time

       call get_snd_parms(path_to_raw_raob,path_to_local_raob
     1                   ,path_to_raw_drpsnd
     1                   ,path_to_raw_satsnd
     1                   ,path_to_raw_tower
     1                   ,path_to_raw_radiometer
     1                   ,istatus)       

       if(istatus .ne. 1)goto999

       lun_out = 21

       if(c8_project .ne. 'RSA')then
           c8_raob_format = c8_project
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_raw_raob,c8_raob_format,lun_out)

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
           call ingest_raob(path_to_raw_raob,c8_raob_format,lun_out)

           c8_raob_format = 'RSA'
           write(6,*)
           write(6,*)' Call ingest_raob, format=',c8_raob_format
           call ingest_raob(path_to_local_raob,c8_raob_format,lun_out)

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
 

!      Satellite soundings
       write(6,*)
       write(6,*)' Call ingest_satsnd...'
       call ingest_satsnd(path_to_raw_satsnd)


!      Radiometer data
       write(6,*)
       write(6,*)' Call get_radiometer_data...'
       i4time_rad = (i4time_sys/3600) * 3600
 900   filename13 = cvt_i4time_wfo_fname13(i4time_rad)
       call s_len(path_to_raw_radiometer,len_path)
       call get_radiometer_data(i4time_rad,laps_cycle_time,ni,nj
     1                         ,i4time_sys - 855/2 ! Radiometer cycle time / 2
     1                         ,i4time_sys + 855/2 ! Radiometer cycle time / 2
     1                         ,path_to_raw_radiometer(1:len_path)
     1                                           //'/'//filename13
     1                         ,lun_out,istatus)
       if(istatus .ne. 1 .and. i4time_sys - i4time_rad .lt. 3600)then
           write(6,*)' Trying radiometer for an hour earlier'
           i4time_rad = i4time_rad - 3600
           go to 900
       endif

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

      subroutine convert_array(array_in,array_out,n,string
     1                        ,r_missing_data,istatus)       

!     QC the observation array and convert units if needed
!     If 'string' is 'none', then do just the QC without conversion

      character*(*) string

      real*4 k_to_c

      real*4 array_in(n),array_out(n)

      do i = 1,n
          if(abs(array_in(n)) .ge. 1e10 .or. 
     1           array_in(n)  .eq. r_missing_data )then
              array_out(n) = r_missing_data
          elseif(string .eq. 'k_to_c')then
              array_out(n) = k_to_c(array_in(n))
          elseif(string .eq. 'pa_to_mb')then
              array_out(n) = array_in(n) / 100.
          elseif(string .eq. 'none')then
              array_out(n) = array_in(n)
          else
              write(6,*)' Unknown operator in convert_array: ',string
              istatus = 0
              return
          endif
      enddo ! i
     
      istatus = 1

      return
      end


      subroutine get_nlevels_snd(pressure_mb,height_m,r_missing_data
     1                          ,nlevels_raw,nlevels_snd)      

      do i = 1,nlevels_raw
          if(abs(pressure_mb(i)) .le. 2000. .or.
     1       abs(height_m(i)) .le. 1000000.)then
              nlevels_snd = i
          endif
      enddo ! i

      return
      end

